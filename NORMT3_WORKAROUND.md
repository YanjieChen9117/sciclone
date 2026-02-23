## 背景：为什么会出现 NORMT3 问题？

在原始的 `sciClone` / `bmm` 生态中：

- 早期版本的 **`bmm` (genome/bmm)** 依赖 **`NORMT3`** 包。
- `sciClone` 的 `DESCRIPTION` 里也把 `NORMT3` 写进了 `Depends`。
- `NORMT3` 已从 CRAN 下架，只能从 Archive 安装：

```bash
wget https://cran.r-project.org/src/contrib/Archive/NORMT3/NORMT3_1.0.4.tar.gz
R CMD INSTALL NORMT3_1.0.4.tar.gz
```

在我们当前的 macOS ARM64 + R 4.5 环境中，上面的安装步骤失败，导致：

- `bmm` 无法从 GitHub 安装（依赖 NORMT3）。
- `sciClone` 无法从 GitHub 安装（依赖 NORMT3 和错误版本的 bmm）。

因此，**NORMT3 是整个安装链路中的“单点失败”**。

---

## 实际报错情况

### 1. 直接安装 NORMT3 失败

命令：

```bash
R CMD INSTALL NORMT3_1.0.4.tar.gz
```

核心报错：

- 缺少 Fortran 编译器：

```text
sh: /opt/gfortran/bin/gfortran: No such file or directory
```

- C 源码不兼容新编译器（旧式声明）：

```text
IPerfcvec.c:32:10: error: too many arguments to function call, expected 0, have 5
IPerfcvec.c:27:6: note: 'IPerfc' declared here: void IPerfc();
...
IPwofzvec.c:30:10: error: too many arguments to function call, expected 0, have 5
```

即使尝试手动修复 C 文件，仍然需要 gfortran 才能编译 Fortran 源码 `toms680-1.f`，在不安装 Fortran toolchain 的前提下，**NORMT3 在该环境基本无法可靠安装**。

### 2. 从 GitHub 安装 genome/bmm 失败

命令：

```r
devtools::install_github("genome/bmm")
```

报错：

```text
ERROR: dependency 'NORMT3' is not available for package 'bmm'
```

原因：`genome/bmm` 的 `DESCRIPTION` 中写有：

```text
Depends: ggplot2, NORMT3
```

### 3. 从 GitHub 安装 genome/sciClone 失败

命令：

```r
devtools::install_github("genome/sciClone")
```

报错：

```text
ERROR: dependency 'NORMT3' is not available for package 'sciClone'
```

原因：`sciClone/DESCRIPTION` 中写有：

```text
Depends: IRanges, bmm, rgl, RColorBrewer, ggplot2, grid, plotrix,
        methods, NORMT3, MKmisc, TeachingDemos, dplyr
```

---

## 我们采用的总体策略

**核心思想：完全避开 NORMT3，而不是强行安装它。**

1. **使用源码形式的旧版 `bmm` (genome/bmm)**，而不是 CRAN 上不兼容的另一个 `bmm` 包。
2. **从 `bmm` 和 `sciClone` 的 `DESCRIPTION` 中移除对 NORMT3 的依赖**。
3. **在 `bmm` 内部自行实现 `erf()` 函数**，替代原本依赖 NORMT3 提供的误差函数。
4. 可选：从 `sciClone` 中移除对 `rgl` 的依赖（仅用于 3D 可视化，在无 X11 的 macOS 上问题多且非核心）。

通过这套方案，我们在 **不安装 NORMT3** 的情况下，成功构建并运行了 `sciClone`，包括：

- 1D 聚类与可视化。
- 使用 `bmm` 的 Beta 混合模型聚类。

---

## 具体修改步骤（可在其他设备复现）

### 步骤 0：前置条件

- 已安装 R（>= 4.x）。
- 已安装开发工具链（Clang，基础 C 编译环境）。
- 已安装 `devtools`、`ggplot2` 等基础 R 包。

```r
install.packages(c("devtools", "ggplot2"), repos = "https://cloud.r-project.org")
```

> 注意：不要求安装 gfortran，也不需要真正安装 NORMT3。

---

### 步骤 1：获取并修正 genome/bmm 源码

在任意工作目录（例如 `/tmp`）：

```bash
cd /tmp
git clone https://github.com/genome/bmm.git
cd bmm
```

编辑 `DESCRIPTION`，将

```text
Depends: ggplot2, NORMT3
```

改为：

```text
Depends: ggplot2
```

然后在 `R/` 目录下添加一个 `erf.R` 文件，实现误差函数：

```r
# 文件: R/erf.R

# Error function implementation
# erf(x) = 2 * pnorm(x * sqrt(2)) - 1

erf <- function(x) {
    2 * stats::pnorm(x * sqrt(2)) - 1
}
```

把 `erf` 导出到 NAMESPACE（在 `NAMESPACE` 文件中添加一行）：

```text
exportPattern("^[^\\.]")
export(erf)
```

构建并安装修正后的 bmm：

```bash
cd /tmp
R CMD build bmm
R CMD INSTALL bmm_0.3.1.tar.gz
```

验证：

```r
library(bmm)
packageVersion("bmm")        # 应为 0.3.1

funcs <- ls("package:bmm")
"init.bmm.hyperparameters" %in% funcs   # TRUE
"bmm.fixed.num.components" %in% funcs   # TRUE
```

至此，我们得到了 **无需 NORMT3 也能正常工作的 bmm 0.3.1**。

---

### 步骤 2：修正本地 sciClone 的 DESCRIPTION（本仓库）

在本仓库根目录（已有 sciClone 源码）：

文件：`DESCRIPTION`

原始内容（关键部分）：

```text
Depends: IRanges, bmm, rgl, RColorBrewer, ggplot2, grid, plotrix,
        methods, NORMT3, MKmisc, TeachingDemos, dplyr
```

修改为（移除 `NORMT3`，可选移除 `rgl`）：

```text
Depends: IRanges, bmm, RColorBrewer, ggplot2, grid, plotrix,
        methods, MKmisc, TeachingDemos, dplyr
```

> 说明：
> - 移除 `NORMT3`：避免强制依赖该包。
> - 移除 `rgl`：在 macOS 上经常需要 XQuartz/GLU，且只影响 3D 图，不影响 1D/2D 主功能。

接着在本仓库根目录构建并安装 sciClone：

```bash
cd /path/to/sciclone           # 本 repo 的根目录
R CMD build sciclone
R CMD INSTALL sciClone_1.1.1.tar.gz
```

验证：

```r
library(sciClone)
library(bmm)
library(VariantAnnotation)

packageVersion("sciClone")   # 1.1.1
packageVersion("bmm")        # 0.3.1
```

---

### 步骤 3：验证 sciClone 聚类功能

在 R 中简单跑一个单样本测试（伪代码示例）：

```r
library(sciClone)

v1 <- read.table("template_vaf_input.txt", header = TRUE, sep = "\t")
cn1 <- read.table("template_copy_number_input.txt", header = TRUE, sep = "\t")

sc <- sciClone(
  vafs           = v1,
  copyNumberCalls= cn1,
  sampleNames    = "Sample1",
  minimumDepth   = 50,
  maximumClusters= 8,
  useSexChrs     = FALSE
)

writeClusterTable(sc, "sciclone_output/test_clusters.txt")
sc.plot1d(sc, "sciclone_output/test_plot_1d.pdf")
```

如果运行成功，说明：

- `bmm` 的 Beta 混合模型部分工作正常。
- `sciClone` 可以完成 1D 聚类与绘图，而无需 NORMT3。

---

## 整体逻辑回顾（便于迁移到其他设备）

1. **问题是如何产生的？**
   - `sciClone` 与旧版 `bmm` 理论上依赖 `NORMT3`。
   - `NORMT3` 已从 CRAN 移除，且在新 macOS + R 版本上编译失败（需要 Fortran 和老式 C 代码修复）。
   - 这使得从 GitHub 安装 `genome/bmm` 和 `genome/sciClone` 均失败。

2. **我们采用的解决思路是什么？**
   - 不强行让 NORMT3 在新环境下编译，而是**从依赖链中去掉它**。
   - 使用 **源码方式** 安装旧版 `bmm`，在其内部自行提供 `erf()` 实现。
   - 修改 `sciClone` 的 `DESCRIPTION`，移除对 NORMT3（以及 rgl）的强依赖。

3. **最终状态如何？**
   - 得到一个：
     - `bmm 0.3.1`（无 NORMT3 依赖，内置 `erf()`）。
     - `sciClone 1.1.1`（依赖本地修正后的 bmm，而非 CRAN bmm，也不再依赖 NORMT3）。
   - 在此基础上，我们成功完成了 TCGA 测试样本的聚类分析与可视化。

---

## 在其他设备上复现的最简步骤（Checklist）

在新机器上，你只需要按顺序做：

1. 安装基础 R 包：

   ```r
   install.packages(c("devtools", "ggplot2"), repos = "https://cloud.r-project.org")
   ```

2. 按本文件 **“步骤 1”** 操作：
   - 从 GitHub 克隆 `genome/bmm`。
   - 修改 `DESCRIPTION` 去掉 `NORMT3`。
   - 添加 `R/erf.R` 和 `NAMESPACE` 导出。
   - `R CMD build` + `R CMD INSTALL`。

3. 按本文件 **“步骤 2”** 操作：
   - 在本地 `sciClone` 源码中修改 `DESCRIPTION`，去掉 `NORMT3`（以及可选的 `rgl`）。
   - `R CMD build` + `R CMD INSTALL`。

4. 用一个简单的单样本数据跑通 `sciClone`（参考 **“步骤 3”**）。

只要上述 4 步完成，你就可以在该设备上运行本仓库中的所有分析脚本，而**不再需要 NORMT3**。


# 🌐 并行程序设计 (Parallel Programming)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-C%2FC%2B%2B-orange.svg)]()
[![CUDA](https://img.shields.io/badge/CUDA-Enabled-green.svg)]()

## 📋 目录 (Table of Contents)

- [项目简介 (Overview)](#-项目简介-overview)
- [项目结构 (Project Structure)](#-项目结构-project-structure)
- [技术栈 (Technology Stack)](#-技术栈-technology-stack)
- [前置要求 (Prerequisites)](#-前置要求-prerequisites)
- [编译与运行 (Compilation & Execution)](#-编译与运行-compilation--execution)
- [实验内容 (Lab Contents)](#-实验内容-lab-contents)
- [硬件环境 (Hardware Specifications)](#-硬件环境-hardware-specifications)
- [性能对比 (Performance Comparison)](#-性能对比-performance-comparison)
- [贡献指南 (Contributing)](#-贡献指南-contributing)
- [许可证 (License)](#-许可证-license)

## 🔍 项目简介 (Overview)

本项目是**并行程序设计课程 (2025 Spring)** 的实验项目集合，基于**高斯消元法 (Gaussian Elimination)** 实现了多种并行计算优化方案。项目展示了从基础串行算法到高级并行技术的完整演进过程，涵盖了现代高性能计算的主要并行编程模型。

### 核心特性

- ✅ **多种并行化策略**：SIMD、多线程、分布式、GPU加速
- ✅ **性能优化实践**：缓存优化、循环展开、向量化、负载均衡
- ✅ **跨平台支持**：x86-64 (SSE/AVX) 和 ARM (NEON) 架构
- ✅ **完整的实验报告**：每个实验都包含详细的PDF报告和性能数据
- ✅ **真实硬件测试**：在华为鲲鹏服务器和NVIDIA Tesla T4 GPU上验证

## 📁 项目结构 (Project Structure)

```
parallel/
├── README.md                          # 项目文档
├── 华为鲲鹏服务器硬件配置.png        # 硬件配置图
│
├── lab0/                              # Lab 0: 并行体系结构调研
│   ├── 2212452-孟启轩-并行体系结构调研.pdf
│   └── TOP500_202411.xlsx             # TOP500超算排名数据
│
├── lab1/                              # Lab 1: 体系结构相关编程
│   ├── 2212452-孟启轩-体系结构相关编程.pdf
│   ├── matrix.cpp/s                   # 矩阵乘法实现
│   ├── sum.cpp/s                      # 求和算法
│   └── data.xlsx                      # 性能测试数据
│
├── lab2/                              # Lab 2: SIMD向量化编程
│   ├── 2212452-孟启轩-SIMD编程.pdf
│   ├── gausscache.cpp                 # 高斯消元缓存优化
│   ├── sse*.cpp                       # SSE指令集实现
│   ├── avx*.cpp                       # AVX/AVX-512指令集实现
│   ├── neon*.cpp                      # ARM NEON指令集实现
│   ├── com*.cpp                       # 各种对比实验
│   ├── data.xlsx                      # 性能测试数据
│   └── 命令2.txt                      # 编译运行命令
│
├── lab3/                              # Lab 3: 多线程并行编程
│   ├── 2212452-孟启轩-多线程编程.pdf
│   ├── col_openmp.cpp                 # OpenMP列划分并行
│   ├── col_barrier.cpp                # 栅障同步实现
│   ├── pthread_neon.cpp               # Pthread + NEON组合
│   ├── avx_omp.cpp                    # AVX + OpenMP组合
│   ├── sse_omp.cpp                    # SSE + OpenMP组合
│   ├── com*.cpp                       # 各种对比实验
│   ├── data.xlsx                      # 性能测试数据
│   └── 命令3.txt                      # 编译运行命令
│
├── lab4/                              # Lab 4: MPI分布式并行编程
│   ├── 2212452-孟启轩-MPI编程.pdf
│   ├── mpi_row.cpp                    # 行划分并行
│   ├── mpi_col*.cpp                   # 列划分并行
│   ├── mpi_2d.cpp                     # 二维划分并行
│   ├── mpi_cache.cpp                  # 缓存优化MPI
│   ├── mpi_omp*.cpp                   # MPI + OpenMP混合
│   ├── mpi_sse.cpp                    # MPI + SSE混合
│   ├── com*.cpp                       # 各种对比实验
│   ├── data.xlsx                      # 性能测试数据
│   └── 命令.txt                       # 编译运行命令
│
├── lab5/                              # Lab 5: GPU并行编程 (CUDA)
│   ├── 2212452-孟启轩-GPU编程.pdf
│   ├── gpu1.cu                        # 基础CUDA实现
│   ├── gpu2.cu                        # 优化CUDA实现
│   ├── gpu3.cu                        # 高级CUDA实现
│   ├── gpu_blocksize.cu               # Block size调优
│   ├── test.cu                        # 测试程序
│   ├── data.xlsx                      # 性能测试数据
│   └── 命令.txt                       # 编译运行命令
│
└── final/                             # 期末研究报告
    ├── 2212452-孟启轩-期末研究报告.pdf
    ├── com.cpp                        # 最终优化版本
    ├── gpu3.cu                        # 最终GPU版本
    └── data.xlsx                      # 综合性能数据
```

## 🛠️ 技术栈 (Technology Stack)

### 编程语言
- C/C++ (C++11标准)
- CUDA C/C++

### 并行编程技术
- **SIMD向量化**
  - Intel SSE (Streaming SIMD Extensions)
  - Intel AVX/AVX2/AVX-512 (Advanced Vector Extensions)
  - ARM NEON
- **多线程并行**
  - POSIX Threads (Pthread)
  - OpenMP (Open Multi-Processing)
- **分布式并行**
  - MPI (Message Passing Interface)
- **GPU加速**
  - NVIDIA CUDA

### 优化技术
- 缓存优化 (Cache Optimization)
- 循环展开 (Loop Unrolling)
- 数据对齐 (Data Alignment)
- 负载均衡 (Load Balancing)
- 内存访问优化 (Memory Access Patterns)

## 📦 前置要求 (Prerequisites)

### 软件要求

#### 基础开发环境
```bash
# GCC/G++ 编译器 (支持C++11)
g++ --version  # 建议 >= 7.0

# 用于SIMD编程
# x86-64: 支持SSE/AVX指令集的处理器
# ARM: 支持NEON指令集的处理器
```

#### Pthread (通常已内置于系统)
```bash
# Linux/Unix系统自带
# 编译时添加 -pthread 标志
```

#### OpenMP
```bash
# Ubuntu/Debian
sudo apt-get install libomp-dev

# CentOS/RHEL
sudo yum install libomp-devel

# macOS (使用 Homebrew)
brew install libomp
```

#### MPI
```bash
# Ubuntu/Debian
sudo apt-get install mpich
# 或
sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# CentOS/RHEL
sudo yum install mpich-devel
# 或
sudo yum install openmpi-devel

# macOS
brew install mpich
# 或
brew install open-mpi

# Windows (推荐使用MS-MPI)
# 从微软官网下载: https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi
```

#### CUDA (用于GPU编程)
```bash
# 需要NVIDIA GPU和CUDA Toolkit
# 从NVIDIA官网下载: https://developer.nvidia.com/cuda-downloads
# 建议版本: CUDA 10.0+ (需根据GPU型号选择合适版本)
# 本项目在CUDA 11.x 和 12.x 上测试通过

# 验证安装
nvcc --version
nvidia-smi
```

### 硬件要求

- **CPU**: x86-64处理器 (支持SSE2及以上) 或 ARM处理器 (支持NEON)
- **内存**: 建议 >= 4GB (大规模矩阵需要更多内存)
- **GPU** (可选): NVIDIA GPU with Compute Capability >= 3.0
- **存储**: >= 1GB 可用空间

## 🚀 编译与运行 (Compilation & Execution)

### Lab 1: 体系结构相关编程

```bash
cd lab1

# 编译矩阵乘法程序
g++ -std=c++11 -O2 matrix.cpp -o matrix
./matrix

# 编译求和程序
g++ -std=c++11 -O2 sum.cpp -o sum
./sum
```

### Lab 2: SIMD编程

```bash
cd lab2

# 基础版本
g++ -std=c++11 -O2 gausscache.cpp -o gausscache
./gausscache

# SSE版本
g++ -std=c++11 -O2 -msse4.2 sse01.cpp -o sse01
./sse01

# AVX版本
g++ -std=c++11 -O2 -mavx2 avx111.cpp -o avx111
./avx111

# AVX-512版本
g++ -std=c++11 -O2 -mavx512f avx512.cpp -o avx512
./avx512

# ARM NEON版本 (需要在ARM平台上编译)
g++ -std=c++11 -O2 -march=armv8-a neon1.cpp -o neon1
./neon1
```

### Lab 3: 多线程编程

```bash
cd lab3

# Pthread版本
g++ -std=c++11 -O2 -pthread com11.cpp -o com11
./com11

# OpenMP版本
g++ -std=c++11 -O2 -fopenmp col_openmp.cpp -o col_openmp
./col_openmp

# OpenMP + AVX混合
g++ -std=c++11 -O2 -fopenmp -mavx2 avx_omp.cpp -o avx_omp
./avx_omp

# OpenMP + SSE混合
g++ -std=c++11 -O2 -fopenmp -msse4.2 sse_omp.cpp -o sse_omp
./sse_omp

# Pthread + NEON混合 (ARM平台)
g++ -std=c++11 -O2 -pthread -march=armv8-a pthread_neon.cpp -o pthread_neon
./pthread_neon
```

### Lab 4: MPI编程

```bash
cd lab4

# Linux/Unix编译示例
# 行划分版本
mpic++ -std=c++11 -O2 mpi_row.cpp -o mpi_row
mpirun -np 4 ./mpi_row

# 列划分版本
mpic++ -std=c++11 -O2 mpi_col.cpp -o mpi_col
mpirun -np 4 ./mpi_col

# 二维划分版本
mpic++ -std=c++11 -O2 mpi_2d.cpp -o mpi_2d
mpirun -np 16 ./mpi_2d

# MPI + OpenMP混合
mpic++ -std=c++11 -O2 -fopenmp mpi_omp.cpp -o mpi_omp
mpirun -np 4 ./mpi_omp

# MPI + SSE混合
mpic++ -std=c++11 -O2 -msse4.2 mpi_sse.cpp -o mpi_sse
mpirun -np 4 ./mpi_sse

# Windows编译示例 (使用MS-MPI)
# 注意: 路径可能因安装位置不同而异，请根据实际安装路径调整
# 64-bit: g++ -O2 mpi_col.cpp -o mpi_col -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include" -L "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -lmsmpi
# 或使用环境变量: g++ -O2 mpi_col.cpp -o mpi_col -I "%MSMPI_INC%" -L "%MSMPI_LIB64%" -lmsmpi
# mpiexec -n 4 mpi_col.exe
```

### Lab 5: GPU编程 (CUDA)

```bash
cd lab5

# 基础版本
nvcc -std=c++11 -O2 gpu1.cu -o gpu1
./gpu1

# 优化版本
nvcc -std=c++11 -O2 gpu2.cu -o gpu2
./gpu2

# 高级优化版本
nvcc -std=c++11 -O2 gpu3.cu -o gpu3
./gpu3

# Block size调优
nvcc -std=c++11 -O2 gpu_blocksize.cu -o gpu_blocksize
./gpu_blocksize

# 查看GPU信息
nvidia-smi
lspci | grep -i nvidia
```

### 期末综合版本

```bash
cd final

# CPU版本
g++ -std=c++11 -O2 -fopenmp -mavx2 com.cpp -o com
./com

# GPU版本
nvcc -std=c++11 -O2 gpu3.cu -o gpu3
./gpu3
```

## 📚 实验内容 (Lab Contents)

### Lab 0: 并行体系结构调研
- **目标**: 调研和分析并行计算体系结构
- **内容**: 
  - TOP500超级计算机排名分析
  - 并行计算体系结构分类
  - 主流并行计算平台对比
- **产出**: 调研报告PDF

### Lab 1: 体系结构相关编程
- **目标**: 理解计算机体系结构对程序性能的影响
- **内容**:
  - 矩阵乘法的缓存优化
  - 向量求和的循环展开
  - 内存访问模式优化
- **关键技术**: Cache优化、循环展开、数据局部性
- **产出**: 实验报告PDF + 源代码 + 性能数据

### Lab 2: SIMD向量化编程
- **目标**: 掌握SIMD指令集进行数据级并行
- **内容**:
  - 高斯消元法的SSE/AVX/AVX-512实现
  - ARM NEON向量化实现
  - 不同SIMD指令集性能对比
- **关键技术**: 
  - Intel SSE (128-bit向量)
  - Intel AVX/AVX2 (256-bit向量)
  - Intel AVX-512 (512-bit向量)
  - ARM NEON (128-bit向量)
- **产出**: 实验报告PDF + 多种SIMD实现 + 性能对比数据

### Lab 3: 多线程并行编程
- **目标**: 掌握共享内存多线程并行编程
- **内容**:
  - Pthread线程池实现
  - OpenMP并行区域和循环并行
  - 行划分、列划分并行策略
  - 栅障同步和临界区
  - SIMD与多线程的混合并行
- **关键技术**:
  - POSIX Threads (Pthread)
  - OpenMP指令和子句
  - 线程同步和负载均衡
- **产出**: 实验报告PDF + 多种多线程实现 + 性能对比数据

### Lab 4: MPI分布式并行编程
- **目标**: 掌握分布式内存并行编程
- **内容**:
  - 行划分、列划分、二维划分并行
  - 点对点通信和集合通信
  - MPI与OpenMP混合编程
  - MPI与SIMD混合编程
  - 负载均衡和通信优化
- **关键技术**:
  - MPI基本通信原语
  - 数据分发和收集
  - 混合并行编程模型
- **产出**: 实验报告PDF + 多种MPI实现 + 性能对比数据

### Lab 5: GPU并行编程 (CUDA)
- **目标**: 掌握GPU大规模并行编程
- **内容**:
  - CUDA编程基础（kernel、thread、block）
  - 全局内存和共享内存优化
  - 线程块大小调优
  - 内存合并访问优化
  - Bank冲突避免
- **关键技术**:
  - CUDA kernel编程
  - GPU内存层次优化
  - 线程组织和同步
- **硬件**: NVIDIA Tesla T4 GPU
- **产出**: 实验报告PDF + 多种CUDA实现 + 性能对比数据

### Final: 期末研究报告
- **目标**: 综合应用所有并行技术，达到最优性能
- **内容**:
  - 所有技术的综合应用
  - 性能调优和瓶颈分析
  - 不同并行方案的深入对比
  - 可扩展性分析
- **产出**: 期末研究报告PDF + 最优实现代码 + 综合性能数据

## 🖥️ 硬件环境 (Hardware Specifications)

本项目在以下硬件平台上进行了测试和性能评估：

### 华为鲲鹏服务器 (ARM架构)
- **处理器**: 华为鲲鹏920 (ARM-based)
- **架构**: ARMv8 64-bit
- **SIMD**: NEON指令集
- **用途**: ARM NEON向量化测试、多线程测试

### x86-64服务器
- **处理器**: Intel Xeon / AMD EPYC (具体型号参见实验报告)
- **SIMD支持**: SSE, AVX, AVX2, AVX-512
- **用途**: x86 SIMD测试、OpenMP/MPI测试

### GPU加速平台
- **GPU型号**: NVIDIA Tesla T4
- **架构**: Turing
- **CUDA Cores**: 2560
- **显存**: 16GB GDDR6
- **Compute Capability**: 7.5
- **用途**: CUDA并行计算测试

详细硬件配置请参见: `华为鲲鹏服务器硬件配置.png`

## 📊 性能对比 (Performance Comparison)

项目中实现的各种并行方案在高斯消元法上的性能对比（具体数据见各实验的data.xlsx）：

### 并行加速比 (Speedup)
1. **串行基准**: 1.0x
2. **缓存优化**: ~2-3x
3. **SIMD (SSE/AVX)**: ~3-8x
4. **多线程 (OpenMP 8线程)**: ~6-10x
5. **SIMD + 多线程**: ~15-30x
6. **MPI (多节点)**: ~10-40x (取决于节点数)
7. **GPU (CUDA)**: ~50-200x (取决于矩阵规模)
8. **混合并行**: ~100-500x (MPI+OpenMP+SIMD)

### 性能优化关键点
- ✅ **缓存优化**: 改善数据局部性，减少Cache Miss
- ✅ **向量化**: 利用SIMD指令实现数据级并行
- ✅ **多线程**: 利用多核CPU实现任务级并行
- ✅ **负载均衡**: 合理划分任务，避免线程空闲
- ✅ **通信优化**: 减少MPI通信开销，重叠计算与通信
- ✅ **GPU优化**: 内存合并访问、共享内存、线程块调优

详细性能数据和分析请参见各实验目录下的PDF报告和Excel数据文件。

## 🤝 贡献指南 (Contributing)

欢迎对本项目进行改进和扩展！

### 如何贡献
1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/YourFeature`)
3. 提交更改 (`git commit -m 'Add some feature'`)
4. 推送到分支 (`git push origin feature/YourFeature`)
5. 创建 Pull Request

### 贡献方向
- 新的并行算法实现
- 性能优化改进
- 文档完善
- Bug修复
- 跨平台支持改进

## 📄 许可证 (License)

本项目采用 MIT 许可证。详见 [LICENSE](LICENSE) 文件。

## 👨‍💻 作者 (Author)

- **姓名**: 孟启轩
- **学号**: 2212452
- **课程**: 并行程序设计 (Parallel Programming)
- **学期**: 2025年春季学期

## 📧 联系方式 (Contact)

如有问题或建议，欢迎通过以下方式联系：
- GitHub Issues: [提交Issue](https://github.com/MengQxuan/parallel/issues)
- Pull Request: 欢迎提交PR

## 🙏 致谢 (Acknowledgments)

- 感谢并行程序设计课程的老师和助教
- 感谢华为提供的鲲鹏服务器测试环境
- 感谢开源社区提供的优秀工具和库

---

**⭐ 如果这个项目对你有帮助，欢迎 Star！**

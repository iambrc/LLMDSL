# LLMDSL

## Getting Started

### 1. Install Gurobi

[Gurobi](https://www.gurobi.com/)

### 2. Clone 

```
git clone https://github.com/iambrc/LLMDSL.git
git submodule update --init --recursive
```

### 3. Build

```
vcpkg install
```


```
cmake -DCMAKE_TOOLCHAIN_FILE="vcpkg/scripts/buildsystems/vcpkg.cmake" -DVCPKG_APPLOCAL_DEPS=ON -DCMAKE_BUILD_TYPE=Release -B build\Release -S .
cmake --build build/Release --config Release
```

### 4. Run
```
build/Release/Release/LLMDSL.exe path\to\yourjsonfile.json 1 1 1 1
```
cmake_minimum_required(VERSION 3.10)

# 定义项目名称 H5，并设置使用 C++ 17 标准
project(H5)

# 设置 C++ 标准为 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置调试和发布模式的编译器标志
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# 设置生成的可执行文件的输出目录到 build/src 目录下
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/src)

# 检查 NphiTotalRange 的源文件，生成可执行文件
add_executable(NphiTotalRange_167
    src/NphiTotalRange_167.cpp  # 主要源文件
)

# 检查 NphiRangeByMass 167 的源文件，生成可执行文件
add_executable(NphiRangeByMass_167
    src/NphiRangeByMass_167.cpp  # 主要源文件
)

# 检查 NphiRangeByMass GM 的源文件，生成可执行文件
add_executable(NphiRangeByMass_GM
    src/NphiRangeByMass_GM.cpp  # 主要源文件
)

# 主程序可执行文件
add_executable(main_program
    src/main.cpp  # 主要源文件
)

# 为 NphiTotalRange 设置头文件目录和编译器标志
target_include_directories(NphiTotalRange_167 PRIVATE ./include)
target_compile_options(NphiTotalRange_167 PRIVATE -Wall -Wextra -Werror)
target_compile_options(NphiRangeByMass_GM PRIVATE -Wall -Wextra -Werror -Wno-unused-parameter)

# 为 NphiRangeByMass 167 设置头文件目录和编译器标志
target_include_directories(NphiRangeByMass_167 PRIVATE ./include)
target_compile_options(NphiRangeByMass_167 PRIVATE -Wall -Wextra -Werror)

# 为 NphiRangeByMass GM 设置头文件目录和编译器标志
target_include_directories(NphiRangeByMass_GM PRIVATE ./include)

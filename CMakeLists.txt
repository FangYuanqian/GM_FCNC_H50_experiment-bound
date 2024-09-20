cmake_minimum_required(VERSION 3.10)

# 定义项目名称 H5，并设置使用 C++ 17 标准
project(H5)

# 设置 C++ 标准为 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置调试和发布模式的编译器标志
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# 检查 NphiTotalRange 的源文件，生成可执行文件
add_executable(NphiTotalRange_167
    src/NphiTotalRange\ 167.cpp  # 主要源文件
)

# 检查 NphiRangeByMass 167 的源文件，生成可执行文件
add_executable(NphiRangeByMass_167
    src/NphiRangeByMass\ 167.cpp  # 主要源文件
)

# 检查 NphiRangeByMass GM 的源文件，生成可执行文件
add_executable(NphiRangeByMass_GM
    src/NphiRangeByMass\ GM.cpp  # 主要源文件
)

# 主程序可执行文件
add_executable(main_program
    src/main.cpp  # 主要源文件
)

# 为 NphiTotalRange 设置头文件目录和编译器标志
target_include_directories(NphiTotalRange_167 PRIVATE /home/fyq/H5/include)
target_compile_options(NphiTotalRange_167 PRIVATE -Wall -Wextra -Werror)

# 为 NphiRangeByMass 167 设置头文件目录和编译器标志
target_include_directories(NphiRangeByMass_167 PRIVATE /home/fyq/H5/include)
target_compile_options(NphiRangeByMass_167 PRIVATE -Wall -Wextra -Werror)

# 为 NphiRangeByMass GM 设置头文件目录和编译器标志
target_include_directories(NphiRangeByMass_GM PRIVATE /home/fyq/H5/include)
target_compile_options(NphiRangeByMass_GM PRIVATE -Wall -Wextra -Werror)

# 为 main_program 设置头文件目录和编译器标志
#target_include_directories(main_program PRIVATE /home/fyq/H5/include)
#target_compile_options(main_program PRIVATE -Wall -Wextra -Werror)

# 如果有外部库，例如数学库或者其他库，添加链接库 (替换为实际的库名称)
# target_link_libraries(NphiTotalRange_167 PRIVATE some_library)
# target_link_libraries(NphiRangeByMass_167 PRIVATE some_library)
# target_link_libraries(NphiRangeByMass_GM PRIVATE some_library)
# target_link_libraries(main_program PRIVATE some_library)

# 为了调试链接器问题，输出更详细的链接器信息（可选）
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -v")

# 添加头文件目录
include_directories("/home/fyq/H5/include/")

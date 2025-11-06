cmake_minimum_required(VERSION 3.22)
project(alglib_gpu LANGUAGES CXX)

option(ALGLIB_GPU_ENABLE_CUDA "Build CUDA backend" ON)
option(ALGLIB_GPU_ENABLE_OPENCL "Build OpenCL fallback backend" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Runtime coerente (evita LNK4098): usar CRT dinâmico e cudart compartilhado
if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
# CUDA 13: '--cudart shared' está obsoleto; usar 'hybrid'.
set(CMAKE_CUDA_RUNTIME_LIBRARY Hybrid)

# Unificar saída de binários na raiz do repo: 'dist-wave' (sem subpastas de configuração)
set(_wave_out "${CMAKE_CURRENT_SOURCE_DIR}/../dist-wave")
file(MAKE_DIRECTORY "${_wave_out}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY         "${_wave_out}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY         "${_wave_out}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY         "${_wave_out}")
foreach(cfg IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${cfg} "${_wave_out}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${cfg} "${_wave_out}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${cfg} "${_wave_out}")
endforeach()

if(ALGLIB_GPU_ENABLE_CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

file(GLOB_RECURSE ALGLIB_SOURCES "src/alglib/*.cpp")

set(CORE_SOURCES
    src/alglib_gpu/GpuBackend.cpp
    src/alglib_gpu/GpuFftDispatcher.cpp
    src/alglib_gpu/alglib_service_runtime.cpp
    src/alglib_gpu/Logging.cpp
)

if(ALGLIB_GPU_ENABLE_CUDA)
  enable_language(CUDA)
  list(APPEND CORE_SOURCES
    src/alglib_gpu/CudaFftExecutor.cu
    src/alglib_gpu/CudaKernels.cu)
endif()

if(ALGLIB_GPU_ENABLE_OPENCL)
  find_package(OpenCL REQUIRED)
  list(APPEND CORE_SOURCES
    src/alglib_gpu/OpenClFftExecutor.cpp
    src/alglib_gpu/OpenClKernels.cpp)
endif()

add_library(alglib_gpu_core STATIC ${CORE_SOURCES})
target_include_directories(alglib_gpu_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib_gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib
)
target_compile_options(alglib_gpu_core PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--cudart=hybrid>)
target_link_options(alglib_gpu_core PRIVATE $<$<LINK_LANGUAGE:CUDA>:--cudart=hybrid>)

target_compile_definitions(alglib_gpu_core
  PRIVATE
    $<$<BOOL:${ALGLIB_GPU_ENABLE_CUDA}>:ALGLIB_GPU_ENABLE_CUDA>
    $<$<BOOL:${ALGLIB_GPU_ENABLE_OPENCL}>:ALGLIB_GPU_ENABLE_OPENCL>
)

# Warnings estritos
target_compile_options(alglib_gpu_core PRIVATE
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/W4 /WX /bigobj>
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-Wall -Wextra -Werror>
)

if(ALGLIB_GPU_ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(_cuda_vendor_runtime "${CMAKE_CURRENT_SOURCE_DIR}/vendor/alglib-gpu/runtime/win64")
  list(APPEND CUDAToolkit_LIBRARY_DIR ${_cuda_vendor_runtime})
  set(_cuda_runtime_targets CUDA::cufft CUDA::cudart CUDA::cuda_driver)
  if(TARGET CUDA::cudadevrt)
    list(APPEND _cuda_runtime_targets CUDA::cudadevrt)
  else()
    set(_cudadevrt_hints
      ${CUDAToolkit_LIBRARY_DIR}
      ${_cuda_vendor_runtime}
      "$ENV{CUDA_PATH}/lib/x64"
      "$ENV{CUDA_PATH}/lib/Win32"
      "$ENV{CUDA_PATH}/targets/x86_64-linux/lib"
      "$ENV{CUDA_PATH}/lib64")
    find_library(CUDA_CUDADEVRT_LIBRARY
      NAMES cudadevrt
      HINTS ${_cudadevrt_hints}
      PATH_SUFFIXES x64 linux)
    if(CUDA_CUDADEVRT_LIBRARY)
      list(APPEND _cuda_runtime_targets ${CUDA_CUDADEVRT_LIBRARY})
    else()
      message(FATAL_ERROR "cudadevrt library not found; ensure CUDA developer runtime is installed")
    endif()
  endif()
  target_link_libraries(alglib_gpu_core PUBLIC ${_cuda_runtime_targets})
  set_target_properties(alglib_gpu_core PROPERTIES
    CUDA_SEPARABLE_COMPILATION OFF
    CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
endif()

if(ALGLIB_GPU_ENABLE_OPENCL)
  target_link_libraries(alglib_gpu_core PUBLIC OpenCL::OpenCL)
endif()

add_library(alglib_gpu SHARED
    src/alglib_gpu/dll_exports.cpp
    ${ALGLIB_SOURCES}
)
target_compile_options(alglib_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--cudart=hybrid>)
target_link_options(alglib_gpu PRIVATE $<$<LINK_LANGUAGE:CUDA>:--cudart=hybrid>)

if(ALGLIB_GPU_ENABLE_CUDA)
  target_sources(alglib_gpu PRIVATE
    src/alglib_gpu/CudaFftExecutor.cu
    src/alglib_gpu/CudaKernels.cu)
endif()

set_target_properties(alglib_gpu PROPERTIES OUTPUT_NAME "alglib")
# Garantir saída em ${_wave_out} sem subpastas
foreach(cfg IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
  set_target_properties(alglib_gpu PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_${cfg} "${_wave_out}"
    LIBRARY_OUTPUT_DIRECTORY_${cfg} "${_wave_out}"
    ARCHIVE_OUTPUT_DIRECTORY_${cfg} "${_wave_out}")
endforeach()

target_include_directories(alglib_gpu
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib_gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib
)

target_compile_definitions(alglib_gpu PRIVATE ALGLIB_GPU_BUILD)

target_link_libraries(alglib_gpu PRIVATE alglib_gpu_core)

if(ALGLIB_GPU_ENABLE_CUDA)
  set_target_properties(alglib_gpu PROPERTIES
    CUDA_SEPARABLE_COMPILATION OFF
    CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
  target_sources(alglib_gpu PRIVATE src/alglib_gpu/CudaFftExecutor.cu)
  target_link_libraries(alglib_gpu PRIVATE ${_cuda_runtime_targets})
endif()

if(ALGLIB_GPU_ENABLE_OPENCL)
  target_link_libraries(alglib_gpu PRIVATE OpenCL::OpenCL)
endif()

if(WIN32)
  target_link_libraries(alglib_gpu PRIVATE ws2_32)
  add_custom_command(TARGET alglib_gpu POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:alglib_gpu>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${_cuda_vendor_runtime}/cudart64_13.dll"
      "${_cuda_vendor_runtime}/cufft64_12.dll"
      "$<TARGET_FILE_DIR:alglib_gpu>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:alglib_gpu>" "${_wave_out}/alglib.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_cuda_vendor_runtime}/cudart64_13.dll" "${_wave_out}/cudart64_13.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_cuda_vendor_runtime}/cufft64_12.dll" "${_wave_out}/cufft64_12.dll"
  )
endif()

install(TARGETS alglib_gpu DESTINATION bin)

add_executable(ALGLIB_DiagService src/alglib_gpu/ALGLIB_DiagService.cpp)
target_link_libraries(ALGLIB_DiagService PRIVATE ws2_32)
set_target_properties(ALGLIB_DiagService PROPERTIES OUTPUT_NAME "alglib_diag_service")
# Garantir saída em ${_wave_out}
foreach(cfg IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
  set_target_properties(ALGLIB_DiagService PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_${cfg} "${_wave_out}")
endforeach()
install(TARGETS ALGLIB_DiagService DESTINATION bin)

# Serviço principal (pipe único) — executável: alglib_service.exe
add_executable(alglib_service src/alglib_gpu/alglib_service.cpp)
target_include_directories(alglib_service PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib_gpu
  ${CMAKE_CURRENT_SOURCE_DIR}/src/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib
)
target_link_libraries(alglib_service PRIVATE alglib_gpu_core ws2_32)
target_compile_options(alglib_service PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX>
  $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall -Wextra -Werror>
)
set_target_properties(alglib_service PROPERTIES OUTPUT_NAME "alglib_service")
# Garantir saída em ${_wave_out}
foreach(cfg IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
  set_target_properties(alglib_service PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_${cfg} "${_wave_out}"
    ARCHIVE_OUTPUT_DIRECTORY_${cfg} "${_wave_out}")
endforeach()
install(TARGETS alglib_service DESTINATION bin)
if(WIN32)
  add_custom_command(TARGET alglib_service POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:alglib_service>" "${_wave_out}/alglib_service.exe"
  )
endif()

 # DLL ponte (transporte + exports tipados): alglib_bridge.dll
add_library(alglib_bridge SHARED
  src/alglib_gpu/alglib_bridge.cpp
  src/alglib_gpu/exports_ops_generated.cpp
)
target_include_directories(alglib_bridge PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib_gpu ${CMAKE_CURRENT_SOURCE_DIR}/src/include)
set_target_properties(alglib_bridge PROPERTIES OUTPUT_NAME "alglib_bridge")
# Garantir saída em ${_wave_out}
foreach(cfg IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
  set_target_properties(alglib_bridge PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_${cfg} "${_wave_out}"
    LIBRARY_OUTPUT_DIRECTORY_${cfg} "${_wave_out}"
    ARCHIVE_OUTPUT_DIRECTORY_${cfg} "${_wave_out}")
endforeach()
if(WIN32)
  target_link_libraries(alglib_bridge PRIVATE ws2_32)
  target_compile_options(alglib_bridge PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/W4 /WX> $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-Wall -Wextra -Werror>)
  add_custom_command(TARGET alglib_bridge POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:alglib_bridge>" "${_wave_out}/alglib_bridge.dll"
  )
endif()
install(TARGETS alglib_bridge DESTINATION bin)

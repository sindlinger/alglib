cmake_minimum_required(VERSION 3.22)
project(alglib_gpu LANGUAGES CXX)

option(ALGLIB_GPU_ENABLE_CUDA "Build CUDA backend" ON)
option(ALGLIB_GPU_ENABLE_OPENCL "Build OpenCL fallback backend" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Unificar saída de binários na raiz do repo: 'dist-wave'
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
    src/dll_core/GpuBackend.cpp
    src/dll_core/GpuFftDispatcher.cpp
    src/service/alglib_service_runtime.cpp
    src/dll_core/Logging.cpp
)

if(ALGLIB_GPU_ENABLE_CUDA)
  enable_language(CUDA)
  list(APPEND CORE_SOURCES
    src/dll_core/CudaFftExecutor.cu
    src/dll_core/CudaKernels.cu)
endif()

if(ALGLIB_GPU_ENABLE_OPENCL)
  find_package(OpenCL REQUIRED)
  list(APPEND CORE_SOURCES
    src/dll_core/OpenClFftExecutor.cpp
    src/dll_core/OpenClKernels.cpp)
endif()

add_library(alglib_core STATIC ${CORE_SOURCES})
target_include_directories(alglib_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dll_core
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dll_bridge
    ${CMAKE_CURRENT_SOURCE_DIR}/src/service
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib_gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib
)

target_compile_definitions(alglib_core
  PRIVATE
    $<$<BOOL:${ALGLIB_GPU_ENABLE_CUDA}>:ALGLIB_GPU_ENABLE_CUDA>
    $<$<BOOL:${ALGLIB_GPU_ENABLE_OPENCL}>:ALGLIB_GPU_ENABLE_OPENCL>
)

target_compile_options(alglib_core PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/bigobj>)

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
  target_link_libraries(alglib_core PUBLIC ${_cuda_runtime_targets})
  set_target_properties(alglib_core PROPERTIES
    CUDA_SEPARABLE_COMPILATION OFF
    CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
endif()

if(ALGLIB_GPU_ENABLE_OPENCL)
  target_link_libraries(alglib_core PUBLIC OpenCL::OpenCL)
endif()

add_library(alglib_bridge SHARED
    src/dll_bridge/dll_exports.cpp
    src/dll_bridge/exports_ops_generated.cpp
    ${ALGLIB_SOURCES}
)

if(ALGLIB_GPU_ENABLE_CUDA)
  target_sources(alglib_bridge PRIVATE
    src/dll_core/CudaFftExecutor.cu
    src/dll_core/CudaKernels.cu)
endif()

# Bridge DLL name
set_target_properties(alglib_bridge PROPERTIES OUTPUT_NAME "alglib_bridge")

target_include_directories(alglib_bridge
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dll_core
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dll_bridge
    ${CMAKE_CURRENT_SOURCE_DIR}/src/service
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib_gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib
)

target_compile_definitions(alglib_bridge PRIVATE ALGLIB_GPU_BUILD)

target_link_libraries(alglib_bridge PRIVATE alglib_core)

if(ALGLIB_GPU_ENABLE_CUDA)
  set_target_properties(alglib_bridge PROPERTIES
    CUDA_SEPARABLE_COMPILATION OFF
    CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
  target_sources(alglib_bridge PRIVATE src/dll_core/CudaFftExecutor.cu)
  target_link_libraries(alglib_bridge PRIVATE ${_cuda_runtime_targets})
endif()

if(ALGLIB_GPU_ENABLE_OPENCL)
  target_link_libraries(alglib_bridge PRIVATE OpenCL::OpenCL)
endif()

if(WIN32)
  target_link_libraries(alglib_bridge PRIVATE ws2_32)
  add_custom_command(TARGET alglib_bridge POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:alglib_bridge>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${_cuda_vendor_runtime}/cudart64_13.dll"
      "${_cuda_vendor_runtime}/cufft64_12.dll"
      "$<TARGET_FILE_DIR:alglib_bridge>")
endif()

add_executable(alglib_service src/service/ALGLIB_Service.cpp)
target_sources(alglib_service PRIVATE ${ALGLIB_SOURCES})
set_target_properties(alglib_service PROPERTIES OUTPUT_NAME "alglib_service")

if(NOT WIN32)
  target_compile_definitions(ALGLIB_Service_app PRIVATE ALGLIB_GPU_SERVICE_STUB)
endif()

target_include_directories(alglib_service
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dll_core
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dll_bridge
    ${CMAKE_CURRENT_SOURCE_DIR}/src/service
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib_gpu
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src/alglib
)

if(WIN32)
  target_link_libraries(alglib_service PRIVATE alglib_core ws2_32)
  add_custom_command(TARGET alglib_service POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:alglib_service>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${_cuda_vendor_runtime}/cudart64_13.dll"
      "${_cuda_vendor_runtime}/cufft64_12.dll"
      "$<TARGET_FILE_DIR:alglib_service>")
else()
  target_link_libraries(alglib_service PRIVATE alglib_core)
endif()

if(ALGLIB_GPU_ENABLE_CUDA)
  target_link_libraries(alglib_service PRIVATE ${_cuda_runtime_targets})
endif()

if(ALGLIB_GPU_ENABLE_OPENCL)
  target_link_libraries(alglib_service PRIVATE OpenCL::OpenCL)
endif()

install(TARGETS alglib_bridge DESTINATION bin)
install(TARGETS alglib_service DESTINATION bin)

add_executable(alglib_diag_service src/service/ALGLIB_DiagService.cpp)
target_link_libraries(alglib_diag_service PRIVATE ws2_32)
set_target_properties(alglib_diag_service PROPERTIES OUTPUT_NAME "alglib_diag_service")
install(TARGETS alglib_diag_service DESTINATION bin)

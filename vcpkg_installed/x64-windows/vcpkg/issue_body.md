Package: eigen3:x64-windows -> 3.4.0#2

**Host Environment**

- Host: x64-windows
- Compiler: MSVC 19.38.33133.0
-    vcpkg-tool version: 2023-09-15-ac02a9f660977426b8ec6392919fbb1d51b10998
    vcpkg-readonly: true
    vcpkg-scripts version: 2c401863dd54a640aeb26ed736c55489c079323b

**To Reproduce**

`vcpkg install `
**Failure logs**

```
CMake Warning at scripts/cmake/vcpkg_buildpath_length_warning.cmake:4 (message):
  eigen3's buildsystem uses very long paths and may fail on your system.

  We recommend moving vcpkg to a short path such as 'C:\src\vcpkg' or using
  the subst command.
Call Stack (most recent call first):
  C:/Users/thoma/AppData/Local/vcpkg/registries/git-trees/250d10d414a5542aaf832350264498ba727c4c03/portfile.cmake:1 (vcpkg_buildpath_length_warning)
  scripts/ports.cmake:147 (include)


-- Using cached libeigen-eigen-3.4.0.tar.gz.
-- Cleaning sources at C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean. Use --editable to skip cleaning for the packages you specify.
-- Extracting source C:/Users/thoma/AppData/Local/vcpkg/downloads/libeigen-eigen-3.4.0.tar.gz
-- Applying patch remove_configure_checks.patch
-- Applying patch fix-vectorized-reductions-half.patch
-- Using source at C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/blds/eigen3/src/3.4.0-74a8d62212.clean
-- Found external ninja('1.11.0').
-- Configuring x64-windows
-- Building x64-windows-dbg
-- Building x64-windows-rel
-- Fixing pkgconfig file: C:/Users/thoma/source/repos/EMinterface/vcpkg_installed/x64-windows/vcpkg/pkgs/eigen3_x64-windows/lib/pkgconfig/eigen3.pc
-- Downloading https://mirror.msys2.org/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst;https://repo.msys2.org/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst;https://mirror.yandex.ru/mirrors/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst;https://mirrors.tuna.tsinghua.edu.cn/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst;https://mirrors.ustc.edu.cn/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst;https://mirror.selfnet.de/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst -> mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst...
[DEBUG] To include the environment variables in debug output, pass --debug-env
[DEBUG] Trying to load bundleconfig from C:\Program Files\Microsoft Visual Studio\2022\Community\VC\vcpkg\vcpkg-bundle.json
[DEBUG] Bundle config: readonly=true, usegitregistry=true, embeddedsha=2c401863dd54a640aeb26ed736c55489c079323b, deployment=VisualStudio, vsversion=17.0
[DEBUG] VS telemetry opted in at SOFTWARE\WOW6432Node\Microsoft\VSCommon\17.0\SQM\\OptIn
[DEBUG] Metrics enabled.
[DEBUG] Feature flag 'binarycaching' unset
[DEBUG] Feature flag 'compilertracking' unset
[DEBUG] Feature flag 'registries' unset
[DEBUG] Feature flag 'versions' unset
[DEBUG] Feature flag 'dependencygraph' unset
Downloading https://mirror.msys2.org/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst
Downloading https://repo.msys2.org/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst
Downloading https://mirror.yandex.ru/mirrors/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst
Downloading https://mirrors.tuna.tsinghua.edu.cn/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst
Downloading https://mirrors.ustc.edu.cn/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst
Downloading https://mirror.selfnet.de/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst
error: Failed to download from mirror set
error: https://mirror.msys2.org/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst: failed: status code 404
error: https://repo.msys2.org/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst: failed: status code 404
error: https://mirror.yandex.ru/mirrors/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst: failed: status code 404
error: https://mirrors.tuna.tsinghua.edu.cn/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst: failed: status code 404
error: https://mirrors.ustc.edu.cn/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst: failed: status code 404
error: https://mirror.selfnet.de/msys2/mingw/mingw32/mingw-w64-i686-pkgconf-1~1.8.0-2-any.pkg.tar.zst: failed: status code 404
[DEBUG] D:\a\_work\1\s\src\vcpkg\base\downloads.cpp(1051): 
[DEBUG] Time in subprocesses: 0us
[DEBUG] Time in parsing JSON: 13us
[DEBUG] Time in JSON reader: 0us
[DEBUG] Time in filesystem: 2681us
[DEBUG] Time in loading ports: 0us
[DEBUG] Exiting after 4.8 s (4748854us)

CMake Error at scripts/cmake/vcpkg_download_distfile.cmake:32 (message):
      
      Failed to download file with error: 1
      If you are using a proxy, please check your proxy setting. Possible causes are:
      
      1. You are actually using an HTTP proxy, but setting HTTPS_PROXY variable
         to `https://address:port`. This is not correct, because `https://` prefix
         claims the proxy is an HTTPS proxy, while your proxy (v2ray, shadowsocksr
         , etc..) is an HTTP proxy. Try setting `http://address:port` to both
         HTTP_PROXY and HTTPS_PROXY instead.
      
      2. If you are using Windows, vcpkg will automatically use your Windows IE Proxy Settings
         set by your proxy software. See https://github.com/microsoft/vcpkg-tool/pull/77
         The value set by your proxy might be wrong, or have same `https://` prefix issue.
      
      3. Your proxy's remote server is out of service.
      
      If you've tried directly download the link, and believe this is not a temporary
      download server failure, please submit an issue at https://github.com/Microsoft/vcpkg/issues
      to report this upstream download server failure.
      

Call Stack (most recent call first):
  scripts/cmake/vcpkg_download_distfile.cmake:270 (z_vcpkg_download_distfile_show_proxy_and_fail)
  scripts/cmake/vcpkg_acquire_msys.cmake:25 (vcpkg_download_distfile)
  scripts/cmake/vcpkg_acquire_msys.cmake:124 (z_vcpkg_acquire_msys_download_package)
  scripts/cmake/vcpkg_acquire_msys.cmake:209 (z_vcpkg_acquire_msys_download_packages)
  scripts/cmake/vcpkg_find_acquire_program(PKGCONFIG).cmake:30 (vcpkg_acquire_msys)
  scripts/cmake/vcpkg_find_acquire_program.cmake:118 (include)
  scripts/cmake/vcpkg_fixup_pkgconfig.cmake:203 (vcpkg_find_acquire_program)
  C:/Users/thoma/AppData/Local/vcpkg/registries/git-trees/250d10d414a5542aaf832350264498ba727c4c03/portfile.cmake:30 (vcpkg_fixup_pkgconfig)
  scripts/ports.cmake:147 (include)



```
**Additional context**

<details><summary>vcpkg.json</summary>

```
{
  "dependencies": [
    "eigen3"
  ]
}

```
</details>

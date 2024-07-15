Example to use rossy_utils in bazel env
```
# ExtTools from rossy_utils
http_archive(
    name = "bazel_rossy",
    sha256 = "4f8fa4efc205b6c1ba2053f7a4ccde6129a47729d246fbba4bd2face053fa6ac",
    strip_prefix = "rossy_utils-bazel-build",
    url = "https://github.com/jeeeunlee/rossy_utils/archive/bazel-build.tar.gz",
)
# eigen
load("@bazel_rossy//thirdparty-bazel/eigen:repository.bzl", "eigen_repository")
eigen_repository(version="3.4.0")
# pinocchio
load("@bazel_rossy//thirdparty-bazel/pinocchio:repository.bzl", "pinocchio_repository")
pinocchio_repository(version="3.1.0")
# Generate local_config_clang repository and configure clang toolchain.
load("@bazel_rossy//thirdparty-bazel/clang:clang_configure.bzl", "clang_autoconf")
clang_autoconf(
    name = "local_config_clang",
    version = "14",
)
```

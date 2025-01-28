# -*- python -*-
#
# Custom build file to use the package with Bazel.

package(default_visibility = ["//visibility:public"])

load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

filegroup(
    name = "srcs",
    srcs = glob(["**"]),
)

cmake(
    name = "osqp",
    lib_source = ":srcs",
    cache_entries = {
        "CMAKE_BUILD_TYPE": "Release",
        "OSQP_BUILD_SHARED_LIB": "ON",
        "OSQP_BUILD_STATIC_LIB": "ON",
        "OSQP_ENABLE_PRINTING": "ON",
        "OSQP_ENABLE_PROFILING": "ON",
        "OSQP_ENABLE_INTERRUPT": "ON",
    },
    # deps = [
    #     "@eigen",
    # ],
    out_include_dir = "include",
    out_static_libs = [
        "libosqpstatic.a",
    ],
    out_shared_libs = [
        "libosqp.so",
    ],
)

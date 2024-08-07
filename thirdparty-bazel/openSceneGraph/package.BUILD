package(
    default_visibility = ["//visibility:public"],
)

load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")


filegroup(
    name = "srcs",
    srcs = glob(["**"]),
)

cmake(
    name = "osg",
    lib_source = ":srcs",
    cache_entries = {
        "BUILD_OSG_EXAMPLES": "OFF",
        "BUILD_OSG_DEPRECATED_SERIALIZERS": "OFF",
        "BUILD_OSG_APPLICATIONS": "OFF",
        "BUILD_OSG_PLUGINS": "OFF",
        "CMAKE_BUILD_TYPE": "Release",
    },
    deps = [
        "@zlib//:zlib",
        # "@eigen",
    ],
    out_include_dir = "include",
    out_shared_libs = [
        "libOpenThreads.so",
        "libosg.so",
        "libosgUtil.so", 
        "libosgDB.so",
        "libosgGA.so",
        "libosgText.so",
        "libosgViewer.so",
        "libosgAnimation.so",
        "libosgFX.so",
        "libosgManipulator.so",
        "libosgParticle.so",
        "libosgUI.so",
        "libosgVolume.so",
        "libosgShadow.so",
        "libosgSim.so",
        "libosgTerrain.so",
        "libosgWidget.so",
        "libosgPresentation.so",
    ],
)

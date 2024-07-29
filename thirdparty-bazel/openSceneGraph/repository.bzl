## repositories.bzl in eigen {:#repositories.bzl-eigen}
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def osg_repository():
    """
    Download repository from GitLab as a ZIP archive, decompress it, and make
    its targets available for binding.

    Args:
        version: version of the library to get.
        sha256: SHA-256 checksum of the downloaded archive.
    """

    http_archive(
        name = "osg",
        urls = [
            "https://github.com/openscenegraph/OpenSceneGraph/archive/refs/tags/OpenSceneGraph-3.6.5.zip",
        ],
        sha256 = "0e9e3e4cc6f463f21a901934a95e9264b231a1d5db90f72dcb4b8cc94b0d1b3b",
        strip_prefix = "OpenSceneGraph-OpenSceneGraph-3.6.5",
        build_file = Label("//rossy_utils/thirdparty-bazel/openSceneGraph:package.BUILD"),
    )

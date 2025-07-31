## repositories.bzl in eigen {:#repositories.bzl-eigen}
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def eigen_repository(
    version = "3.3.4",
    sha256 = ""):
    """
    Download repository from GitLab as a ZIP archive, decompress it, and make
    its targets available for binding.

    Args:
        version: version of the library to get.
        sha256: SHA-256 checksum of the downloaded archive.
    """
    ver2sha = {"3.4.0":"eba3f3d414d2f8cba2919c78ec6daab08fc71ba2ba4ae502b7e5d4d99fc02cda",
                "3.3.4":"b4415a2e16bbaff7568e52853c9e131633248c9bc66d565086e095087e95877d"}
    if(sha256 == ""):
        if version in ver2sha: 
            sha256 = ver2sha[version]
    http_archive(
        name = "eigen",
        urls = [
            "https://gitlab.com/libeigen/eigen/-/archive/{}/eigen-{}.zip".format(version, version),
        ],
        sha256 = sha256,
        strip_prefix = "eigen-{}".format(version),
        build_file = Label("//rossy_utils/thirdparty-bazel/eigen:package.BUILD"),
    )

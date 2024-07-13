## repositories.bzl in pinocchio {:#repositories.bzl-pinocchio}
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def pinocchio_repository(
        version = "2.6.10",
        sha256 = ""):
    """
    Clone repository from GitHub and make its targets available for binding.

    Args:
        version: version of the library to get.
        sha256: SHA-256 checksum of the downloaded archive.
    """
    ver2sha = {"2.6.0":"4ded587192561d434d4158fdae56dd36b0c3a0db65ae4768df94428161bad522" }
    if(sha256 == ""):
        if version in ver2sha: 
            sha256 = ver2sha[version]
    
    http_archive(
        name = "pinocchio",
        urls = [
            "https://github.com/stack-of-tasks/pinocchio/releases/download/v{}/pinocchio-{}.tar.gz".format(version, version),
        ],
        sha256 = sha256,
        strip_prefix = "pinocchio-{}".format(version),
        build_file = Label("//thirdparty-bazel/pinocchio:package.BUILD"),
    )

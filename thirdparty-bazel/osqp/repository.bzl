## repositories.bzl in osqp {:#repositories.bzl-osqp}
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def osqp_repository(
        version = "1.0.0.beta1",
        sha256 = ""):
    """
    Clone repository from GitHub and make its targets available for binding.

    Args:
        version: version of the library to get.
        sha256: SHA-256 checksum of the downloaded archive.
    """
    ver2sha = {"1.0.0.beta1":"648319ecab199e0a2c8d6dcfe92903600381bc5b110c2c3fba57726191eed65f"}
    if(sha256 == ""):
        if version in ver2sha: 
            sha256 = ver2sha[version]
    
    http_archive(
        name = "osqp",
        urls = [
            "https://github.com/osqp/osqp/archive/refs/tags/v{}.tar.gz".format(version),
        ],
        sha256 = sha256,
        strip_prefix = "osqp-{}".format(version),
        build_file = Label("//rossy_utils/thirdparty-bazel/osqp:package.BUILD"),
    )

    

"""Repository rule for clang toolchain configuration"""

def _clang_autoconf_impl(repository_ctx):
    for program in ["clang-format", "clang-tidy"]:
        program_name = "-".join([program, repository_ctx.attr.version])
        program_path = repository_ctx.which(program_name)
        if not program_path:
            fail("{} not found in PATH".format(program_name))
        repository_ctx.symlink(program_path, program)

    repository_ctx.template(
        "BUILD.bazel",
        repository_ctx.attr._build,
        executable = False,
    )

clang_autoconf = repository_rule(
    implementation = _clang_autoconf_impl,
    environ = ["PATH"],
    attrs = {
        "version": attr.string(
            mandatory = True,
            doc = "clang version",
        ),
        "_build": attr.label(
            default = Label("//tools/clang:BUILD.bazel.in"),
            allow_single_file = True,
        ),
    },
    configure = True,
    doc = "Generates repository containing clang toolchain",
)

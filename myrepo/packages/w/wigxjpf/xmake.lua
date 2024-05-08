package("wigxjpf")
    set_homepage("https://fy.chalmers.se/subatom/wigxjpf/")
    set_description("WIGXJPF evaluates Wigner 3j, 6j and 9j symbols accurately using prime factorisation and multi-word integer arithmetic.")

    set_urls("https://fy.chalmers.se/subatom/wigxjpf/wigxjpf-$(version).tar.gz")
    add_versions("1.13", "90ab9bfd495978ad1fdcbb436e274d6f4586184ae290b99920e5c978d64b3e6a")

    on_install(function (package)
        os.vrun("make")
        os.cp("inc/*.h", package:installdir("include"))
        os.cp("lib/*.a", package:installdir("lib"))
        os.cp("bin/*", package:installdir("bin"))
    end)

    on_test(function (package)
        os.vrun("make test")
    end)
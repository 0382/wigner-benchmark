add_rules("mode.debug", "mode.release")

add_requires("gmp")

target("pf_rational")
    set_kind("static")
    set_languages("c11")
    add_files("pf_rational/*.c")
    add_includedirs("pf_rational/")
    add_packages("gmp")
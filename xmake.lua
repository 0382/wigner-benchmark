add_rules("mode.debug", "mode.release")

add_requires("gmp", {system = true})

target("pf_rational")
    set_kind("static")
    set_languages("c11")
    -- add_cflags("-march=native", "-mtune=native", "-fopenmp")
    add_files("pf_rational/*.c")
    add_includedirs("pf_rational/")
    add_packages("gmp")

target("test")
    set_kind("binary")
    set_languages("c11", "cxx17")
    add_files("src/*.cpp")
    add_includedirs("pf_rational/", "WignerSymbol/")
    add_deps("pf_rational")
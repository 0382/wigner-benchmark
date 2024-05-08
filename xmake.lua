add_rules("mode.debug", "mode.release")

add_requires("gnu-gsl")

add_repositories("local-repo myrepo")
add_requires("wigxjpf")

target("error")
    set_kind("binary")
    set_languages("c11", "cxx17")
    add_files("src/error.cpp")
    add_includedirs("WignerSymbol/")
    add_packages("wigxjpf", "gnu-gsl")

target("bench")
    set_kind("binary")
    set_languages("c11", "cxx17")
    add_files("src/bench.cpp")
    add_includedirs("WignerSymbol/")
    add_packages("wigxjpf", "gnu-gsl")

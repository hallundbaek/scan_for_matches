from distutils.core import setup, Extension

extension_mod = Extension(name = "scan_for_matches",
                          sources = ["pymodule.c","../../scan_for_matches.c", "../../parser.c", "../../scanner.c"])

setup(name = "scan_for_matches", version='2.0', ext_modules=[extension_mod])

#!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Ole Schuett

import re
import os
import sys
from os import path
from datetime import datetime

flag_exceptions_re = re.compile(r"__COMPILE_REVISION|__COMPILE_DATE|__COMPILE_ARCH|__COMPILE_HOST|"
                              + r"__INTEL_COMPILER|__cplusplus|\$\{.*\}\$|__.*__")

year = datetime.utcnow().year

BANNER_F = "!--------------------------------------------------------------------------------------------------!\n" \
          +"!   CP2K: A general program to perform molecular dynamics simulations                              !\n" \
          +"!   Copyright (C) 2000 - %d  CP2K developers group                                               !\n"%year \
          +"!--------------------------------------------------------------------------------------------------!\n"

BANNER_Fypp = "#!-------------------------------------------------------------------------------------------------!\n" \
             +"#!   CP2K: A general program to perform molecular dynamics simulations                             !\n" \
             +"#!   Copyright (C) 2000 - %d  CP2K developers group                                              !\n"%year \
             +"#!-------------------------------------------------------------------------------------------------!\n"

BANNER_C = "/*****************************************************************************\n" \
          +" *  CP2K: A general program to perform molecular dynamics simulations        *\n" \
          +" *  Copyright (C) 2000 - %d  CP2K developers group                         *\n"%year \
          +" *****************************************************************************/\n"

#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Usage: analyse_src.py <cp2k-root-dir>")
        print("       This tool checks the source code for violations of the coding conventions.")
        sys.exit(1)

    cp2k_dir = sys.argv[1]

    # check flags and banners
    flags = set()
    for root, _, files in os.walk(path.join(cp2k_dir, "src")):
        #print("Scanning %s ..."%root)
        for fn in files:
            fn_ext = fn.rsplit(".",1)[-1]
            if(fn_ext in ("template", "instantiate")): continue
            content = open(path.join(root, fn)).read()

            # check banner
            if((fn_ext in ("F", ) and not content.startswith(BANNER_F)) or
               (fn_ext in ("fypp", ) and not content.startswith(BANNER_Fypp)) or
               (fn_ext in ("c", "cu", "cpp", "h", "hpp") and not content.startswith(BANNER_C))):
                    print(fn+": Copyright banner malformed")
                    #print('"'+ '"\n"'.join(content.split("\n")[:4]) + '"')

            # find all flags
            for line in content.split("\n"):
                if(len(line) == 0): continue
                if(line[0] != "#"): continue
                if(line.split()[0] not in ("#if","#ifdef","#ifndef","#elif")): continue
                line = line.split("//",1)[0]
                line = re.sub("[|()!&><=*/+-]", " ", line)
                line = line.replace("defined", " ")
                for m in line.split()[1:]:
                    if re.match("[0-9]+[ulUL]*", m): continue # skip numbers
                    if(fn_ext=="h" and fn.upper().replace(".", "_") == m): continue
                    flags.add(m)

    flags = [f for f in flags if not flag_exceptions_re.match(f)]

    #print("Found %d flags."%len(flags))
    #print(flags)

    install_txt = open(path.join(cp2k_dir, "INSTALL.md")).read()
    cp2k_info = open(path.join(cp2k_dir, "src/cp2k_info.F")).read()
    flags_src = re.search(r"FUNCTION cp2k_flags\(\)(.*)END FUNCTION cp2k_flags", cp2k_info, re.DOTALL).group(1)

    for f in sorted(flags):
        if(f not in install_txt):
            print("Flag %s not mentioned in INSTALL.md"%f)
        if(f not in flags_src):
            print("Flag %s not mentioned in cp2k_flags()"%f)

    # check for copies of data files
    data_files = set()
    for _, _, files in os.walk(path.join(cp2k_dir, "data")):
        data_files.update(files)
    data_files.remove("README")
    for root, _, files in os.walk(path.join(cp2k_dir, "tests")):
        d = path.relpath(root, cp2k_dir)
        for c in data_files.intersection(files):
            print("Data file %s copied to %s"%(c, d))

    # check linebreaks and encoding
    for root, _, files in os.walk(cp2k_dir):
        if ("tools/toolchain/build" in root):
            continue
        if any(x in root for x in (".git", "obj", "lib", "exe", "regtesting", "exts")):
            continue
        for fn in files:
            absfn = path.join(root, fn)
            shortfn = absfn[len(cp2k_dir):]
            if(not path.exists(absfn)):
                continue # skip broken symlinks
            content = open(absfn).read()
            if('\0' in content):
                continue # skip binary files
            if("\r\n" in content):
                print("Text file %s contains DOS linebreaks"%shortfn)

            # check for non-ascii chars
            if("# -*- coding: utf-8 -*-" in content):
                continue
            if(not re.search("[\x80-\xFF]", content)):
                continue
            for lineno, line in enumerate(content.split("\n")):
                m = re.search("[\x80-\xFF]", line)
                if(m):
                    print("Found non-ascii char in %s line %d at position %d"%(shortfn, lineno+1, m.start(0)+1))


#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF

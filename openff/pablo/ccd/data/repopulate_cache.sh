#!/usr/bin/env bash

set -euf -o pipefail

cd "$(dirname "$0")"

rm -r ccd_cache || true
mkdir ccd_cache
cd ccd_cache

resnames=(
    ALA
    ARG
    ASN
    ASP
    CYS
    GLN
    GLU
    GLY
    HIS
    ILE
    LEU
    LYS
    MET
    PHE
    PRO
    SER
    THR
    TRP
    TYR
    VAL
    DG
    DA
    DT
    DC
    G
    A
    U
    C
    ACE
    NME
    NA
    CL
    BR
    CS
    IOD
    LI
    RB
    XE
    F
    K
    HOH
)
printf -v joined '%s,' "${resnames[@]}"
curl --remote-name-all --parallel "https://files.rcsb.org/ligands/download/{${joined%,}}.cif"

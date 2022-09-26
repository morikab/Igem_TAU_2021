#!/bin/sh
cd /tamir1/liyamlevi/projects/communique/Igem_TAU_2021/software_analysis/code/constructed_microbiomes/tests/test_local_align_KDVY_first_entry
chmod 777 ./*
for file in ./*; do   ./"$file" ; done

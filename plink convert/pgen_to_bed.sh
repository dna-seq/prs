#!/bin/bash
"""produce a binary fileset"""
plink2 --pfile /path/crossmap_alt --make-bed --out example_out --allow-extra-chr


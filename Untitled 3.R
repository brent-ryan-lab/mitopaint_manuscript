# Title: plate drift correct mitopaint visualisations (mean per well) v1
# R: 4.4.1
# Author: Sarah Franks
# Project: mitopaint manuscript
# Last edit: 11-06-2026

# load packages ####
library(data.table)
library(tidyverse)
library(colorspace)
# set variables ####
file_name_N1
file_name_N2
file_name_N3
dmso_wells
pos_control
pastel_cols <- lighten(c("#440154FF", "#238A8DFF", "#FDE725FF"), amount = 0.3)
# create function to generate plate drift QC plots of dmso ####
# run function to generate plate drift QC plots of dmso ####
# create function to generate plate drift QC plots of dmso + positive control ####
# run function to generate plate drift QC plots of dmso + positive control ####
# create function to calculate proportion of drift corrected features per N ####

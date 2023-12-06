#!/bin/bash

export LC_NUMERIC="en_US.UTF-8"
set -e

# 1) Generate white background movie
for file in _light.*.ppm; 
do 
	echo ${file:0:12}.png
	convert $file -resize 421x221 -transparent white ${file:0:12}.png;
done
img2webp -o ../../../../docs/sphinx/source/tutorials/figures/reactive-silicon-dioxide/hydronium_transfert_light.webp -q 30 -mixed -d 11.11 _light*.png

# 2) Generate black background movie
for file in _dark.*.ppm; 
do 
	echo ${file:0:11}.png
	convert $file -resize 421x221 -transparent black ${file:0:11}.png;
done
img2webp -o ../../../../docs/sphinx/source/tutorials/figures/reactive-silicon-dioxide/hydronium_transfert_dark.webp -q 30 -mixed -d 11.11 _dark*.png

#!/bin/bash

# A script to solve for CCD corrections for multi-spectral images. Must
# have the desired MS image, corresponding PAN image, their cameras.
# Put the ASP bin directory in the path.

# Specify the full path to this tool, as it will look for a python script
# in the same directory.

# This will average the disparity in each column, and save them to
# outputPrefix/run-avg-dx.txt and outputPrefix/run-avg-dy.txt
# Then it will compute the CCD corrections based on these
# using the script ccd_process.py
# and save them in
# outputPrefix/run-corr-dx.txt and outputPrefix/run-corr-dy.txt
#
# If desired to compute these corrections by averaging multiple such
# average disparities, that python script can be invoked directly as:
#
# python ccd_process.py */run-avg-dx.txt --output-prefix combined-corr

# cropWin "0 0 8820 5000"

if [ "$#" -ne 7 ]; then echo Usage: $0 panImage.tif msImage.tif panCamera.xml msCamera.xml msBandNum outputPrefix/run cropWin; exit; fi

panImage=$1
msImage=$2
panCamera=$3
msCamera=$4
msBand=$5
outputPrefix=$6
cropWin=$7

function run_cmd () {
    cmd=$1
    out=$2
    echo ""
    echo $cmd

    if [ "$out" != "" ]; then
        echo Writing the output to $out
        eval $cmd > $out 2>&1
        ans=$?
    else
        eval $cmd
        ans=$?
    fi
    
    if [ "$ans" -ne 0 ]; then
        echo Command failed
        exit 1
    fi
    
    echo ""
}

function gen_dummy_cam {
    name=$1
    echo \
"fu = 1
fv = 1
cu = 1
cv = 1
u_direction = 1 0 0
v_direction = 0 1 0
w_direction = 0 0 1
C = 0 0 0
R = 1 0 0 0 1 0 0 0 1
NULL
" \
> $name

}

if [ "$msBand" -le 1 ] || [ "$msBand" -ge 8 ]; then
    echo "The MS band must be between 1 and 8."
fi

suffix=$(basename $outputPrefix)
outDir=$(dirname $outputPrefix)
if [ "$suffix" = "" ] || [ "$outDir" = "." ]; then
    echo Please specify the output prefix as outputDir/run
    exit 1
fi

cmd="mkdir -p $outDir"
run_cmd "$cmd"

echo Reducing the PAN image resolution to make it comparable to the multispecral image resolution
panImageSub4=$(echo $panImage | perl -p -e "s#\.(TIF|NTF)#_sub4.tif#ig")
panImageSub4=${outputPrefix}-${panImageSub4}
cmd="mkdir -p $(dirname $panImageSub4)"
run_cmd "$cmd"
cmd="gdal_translate -r average -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -outsize 25% 25% $panImage $panImageSub4"
#echo $cmd
run_cmd "$cmd"

echo Separating the desired multispecral band
msImageCurrBand=$(echo $msImage | perl -p -e "s#\.(TIF|NTF)#_b${msBand}.tif#ig")
msImageCurrBand=${outputPrefix}-${msImageCurrBand}
cmd="mkdir -p $(dirname $msImageCurrBand)"
run_cmd "$cmd"
cmd="gdal_translate -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -b $msBand $msImage $msImageCurrBand"
#echo $cmd
run_cmd "$cmd"

leftTsaiCam="$outDir/dummy1.tsai"
rightTsaiCam="$outDir/dummy2.tsai"
echo Creating fake cameras: $leftTsaiCam $rightTsaiCam
gen_dummy_cam $leftTsaiCam
gen_dummy_cam $rightTsaiCam

# Since we do PAN to MS stereo, the correlation and subpixel kernels can be small,
# and in fact, having them small helps in not blurring the disparity average
# too much.
echo Running stereo to find the disparity
cmd="stereo --stop-point 4 --individually-normalize --alignment-method none --left-image-crop-win $cropWin --right-image-crop-win $cropWin --corr-search 5 -5 10 10 --corr-kernel 9 9 --subpixel-kernel 9 9 $panImageSub4 $msImageCurrBand $leftTsaiCam $rightTsaiCam $outputPrefix"
out=${outputPrefix}-stereo-output.txt
#echo $cmd 
run_cmd "$cmd" $out

echo Averaging the disparity
cmd="disp_avg --remove-outliers-params 100 3 --save-no-metadata ${outputPrefix}-RD.tif ${outputPrefix}-avg-dx.txt ${outputPrefix}-avg-dy.txt"
#echo $cmd
run_cmd "$cmd"

echo Process the averaged disparities
ccd_process_path=$(dirname $0)
cmd="$(which python) $ccd_process_path/ccd_process.py ${outputPrefix}-avg-dx.txt --output-prefix ${outputPrefix}-corr --no-plot"
#echo $cmd
run_cmd "$cmd"

echo Separating the channels of the disparity and colorizing them

b=1
cmd="gdal_translate -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -b $b ${outputPrefix}-RD.tif ${outputPrefix}-RD_b${b}.tif"
#echo $cmd
run_cmd "$cmd"

isForward=$(grep -i scan $msCamera | grep -i forward)

if [ "$isForward" != "" ]; then
    cmd="colormap --min 10.5 --max 11.5 ${outputPrefix}-RD_b${b}.tif"
else
    cmd="colormap --min 11.0 --max 12.0 ${outputPrefix}-RD_b${b}.tif"
fi
run_cmd "$cmd"

b=2
cmd="gdal_translate -co TILED=YES -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -b $b ${outputPrefix}-RD.tif ${outputPrefix}-RD_b${b}.tif"
run_cmd "$cmd"

if [ "$isForward" != "" ]; then
    cmd="colormap --min -0.5 --max 0.5 ${outputPrefix}-RD_b${b}.tif"
else
    cmd="colormap --min -1.0 --max 0.0 ${outputPrefix}-RD_b${b}.tif"
fi
run_cmd "$cmd"

echo ""
echo Corrections in x and y:
echo ${outputPrefix}-corr-dx.txt
echo ${outputPrefix}-corr-dy.txt

echo ""
echo Colormaps before the correction:
echo ${outputPrefix}-RD_b1_CMAP.tif
echo ${outputPrefix}-RD_b2_CMAP.tif



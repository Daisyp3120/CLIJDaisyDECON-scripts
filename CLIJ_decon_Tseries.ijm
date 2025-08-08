// To make this script run in Fiji, please activate 
// the clij and clij2 update sites in your Fiji 
// installation. Read more: https://clij.github.io

// Generator version: 2.5.1.6

//Designed to work with OMERO batch process plugin: https://github.com/GReD-Clermont/omero_batch-plugin?tab=readme-ov-file

//Define inputs: path to psf file and number of iterations of Richardson-Lucy Deconvolution
#@ String (visibility=MESSAGE, value="<html>Select PSF file in order for each channel you want to deconvolve.</html>", required=false) msg
#@ File[] (label="Select PSF files", style="files") psf_paths
#@ Integer (label = "Iterations:", min=0, max=5000, value=10) num_iterations
#@ Integer (label = "Save output every x iterations, x:", min=0, max=50, value=1) out_iter
#@ String (visibility=MESSAGE, value="<html>To only output final result set x to 0.</html>", required=false) msg

total_time=getTime();
// Init GPU
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();

Ext.CLIJ2_getGPUProperties(gpu, memory, opencl_version);
print("GPU: " + gpu);
print("Memory in GB: " + (memory / 1024 / 1024 / 1024) );
print("OpenCL version: " + opencl_version);

// Load image
image_stack = getTitle();
Stack.getDisplayMode(mode);
getVoxelSize(dx, dy, dz_stage, unit);
Stack.getDimensions(width, height, image_channels, slices, frames);
//run("Split Channels");

if (lengthOf(psf_paths) != image_channels){
	exit("Number of channels in image and psf are different");
}

//load and push psfs
psfs=newArray(lengthOf(psf_paths));
for (i = 0; i < lengthOf(psf_paths); i++) {
psf_path=psf_paths[i];
open(psf_path);
psfs[i] = getTitle();
print("PSF: "+psfs[i]);
Ext.CLIJ2_push(psfs[i]);
}


//Deconvolve each channel one at a time
//Initialise string for concatenating timepoints
concat_string=" title=Decon_"+image_stack+" open";
for (frame = 1; frame <= frames; frame++) {
	//Initialise string for merging channels
	merge_string="";
for (i = 1; i <= image_channels; i++) {
selectWindow(image_stack);
//Get current LUT of channel
getLut(reds, greens, blues);
Stack.setPosition(i, 1, frame);

print("Number of iterations: "+num_iterations);

time=getTime();
Ext.CLIJ2_pushCurrentZStack(image_stack);
print("Pushing image took " + (getTime() - time) + " msec");

// Richardson Lucy Deconvolution
time = getTime();
if (out_iter == 0) {
	out_iter = num_iterations;
}
num_decon=floor(num_iterations/out_iter);
extra_iterations=num_iterations-(out_iter*num_decon);
image=image_stack;

for (j = 1; j <= num_decon; j++){
	Ext.CLIJx_deconvolveRichardsonLucyFFT(image, psfs[i-1], output_image, out_iter);
	Ext.CLIJ2_release(image);
	//close(image);
	image=Ext.CLIJ2_pull(output_image);
	//run("Duplicate...", "title=CLIJx_Decon_iter"+(out_iter*j)+"_C"+i+"_"+image_stack+" duplicate");
	rename("CLIJx_Decon_iter"+(out_iter*j)+"_C"+i+"_"+image_stack);
	Ext.CLIJ2_push(image);
}
if (extra_iterations > 0) {
	Ext.CLIJx_deconvolveRichardsonLucyFFT(image, psf, output_image, extra_iterations);
}
print("Deconvolution took " + (getTime() - time) + " msec");
Ext.CLIJ2_release(image);
//close(image);

Ext.CLIJ2_pull(output_image);
Ext.CLIJ2_release(output_image);
close(output_image);

rename("CLIJx_Decon_C"+i+"_T"+frame+"-"+image_stack);
//Apply LUT from original image
setLut(reds, greens, blues);
merge_string = merge_string + "c"+i+"=[CLIJx_Decon_C"+i+"_T"+frame+"-"+image_stack+"] ";
}


//Merge deconvolved channels back together
run("Merge Channels...", merge_string+"create");
rename("Decon_T"+frame+"_"+image_stack);
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print("Deconvolved with CLIJx_imageJ2RichardsonLucyDeconvolution on "+year+"-"+month+1+"-"+dayOfMonth);
concat_string=concat_string + " image"+frame+"=[Decon_T"+frame+"_"+image_stack+"]";
}

for (i = 0; i < lengthOf(psfs); i++) {
Ext.CLIJ2_release(psfs[i]);
close(psfs[i]);
}

if (frames > 1) {
	concat_string=concat_string + " image"+frames+1+"=[-- None --]";
	run("Concatenate...", concat_string);
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+image_channels+" slices="+slices+" frames="+frames+" display=Color");
	Stack.setDisplayMode(mode);
}
setVoxelSize(dx, dy, dz_stage, unit);
print("Decon total time " + (getTime() - total_time) + " msec");

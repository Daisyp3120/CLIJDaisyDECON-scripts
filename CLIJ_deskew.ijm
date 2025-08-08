//Haase, R., Royer, L.A., Steinbach, P. et al. CLIJ: GPU-accelerated image processing for everyone. Nat Methods 17, 5–6 (2020). https://doi.org/10.1038/s41592-019-0650-
#@ String (choices={"3i","Zeiss"}, style="listBox") scope
#@ Boolean (label="Rotate/coverslip correct") rotate

total_time=getTime();
getVoxelSize(dx, dy, dz_stage, unit);
Stack.getDimensions(width, height, channels, slices, frames);
Stack.getDisplayMode(mode);
orig_title=getTitle();

run("Duplicate...", "duplicate");
title=getTitle();
image_stack=title;
if (scope == "3i") {
	angle=32.8;
	dz=sin(angle*PI/180)*dz_stage;
	deskewfactor = (cos(angle*PI/180)*dz_stage)/dx;
	padding=width+slices*deskewfactor;
	run("Canvas Size...", "width="+padding+" height="+height+" position=Center-Left zero");
} else if (scope == "Zeiss") {
	angle=30;
	dz=sin(angle*PI/180)*dz_stage;
	deskewfactor = (cos(angle*PI/180)*dz_stage)/dx;
	padding=height+slices*deskewfactor;
	run("Canvas Size...", "width="+width+" height="+padding+" position=Top-Center zero");
} else {
	exit("Scope not recognised");
}

// affine transform
run("CLIJ2 Macro Extensions","cl_device=");
concat_string=" title=Deskew_"+image_stack+" open";
for (frame = 1; frame <= frames; frame++) {
	//Initialise string for merging channels
	merge_string="";
	for (i = 1; i <= channels; i++) {
		selectWindow(image_stack);
		//Get current LUT of channel
		getLut(reds, greens, blues);
		Stack.setPosition(i, 1, frame);

		image1 = title;
		Ext.CLIJ2_pushCurrentZStack(image1);
		if (scope == "3i") {
			image2 = title+"_shearXZ";
			transform = "shearXZ="+-deskewfactor;
			if (rotate){
				transform = transform + " rotateY="+angle;
			}
		}
		if (scope == "Zeiss") {
			image2 = title+"_shearYZ";
			transform = "shearYZ="+-deskewfactor;
		}
		Ext.CLIJ2_affineTransform3D(image1, image2, transform);
		Ext.CLIJ2_pull(image2);
		rename("CLIJ2_Deskew_C"+i+"_T"+frame+"_"+image_stack);
		//close(image1);
		merge_string = merge_string + "c"+i+"=[CLIJ2_Deskew_C"+i+"_T"+frame+"_"+image_stack+"] ";
	}
	run("Merge Channels...", merge_string+"create");
	rename("Deskew_T"+frame+"_"+image_stack);
	concat_string=concat_string + " image"+frame+"=Deskew_T"+frame+"_"+image_stack;
}

if (frames > 1) {
	concat_string=concat_string + " image"+frames+1+"=[-- None --]";
	run("Concatenate...", concat_string);
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+slices+" frames="+frames+" display=Color");
	Stack.setDisplayMode(mode);
}
rename("Deskew_"+orig_title);
setVoxelSize(dx, dy, dz, unit);
close(title);
print("Deskew total time " + (getTime() - total_time) + " msec");

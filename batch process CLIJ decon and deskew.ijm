// CLIJ2 Deconvolution + Deskew Macro

#@ File (label = "Input folder", style = "directory") input
#@ File (label = "Output folder", style = "directory") output
#@ File[] (label="PSF files (one per channel)", style="files") psf_paths
#@ Integer (label = "Iterations", value = 10) num_iterations
#@ String (choices={"3i", "Zeiss"}, label = "Microscope Type") scope
#@ String (label = "File suffix", value = ".tif") suffix

print("Starting batch process...");
print("Input:  " + input);
print("Output: " + output);

function processFolder(dir) {
    files = getFileList(dir);
    for (i = 0; i < files.length; i++) {
        item = files[i];
        full = dir + File.separator + item;
        if (File.isDirectory(full)) {
            processFolder(full);
        } else if (endsWith(item, suffix)) {
            processFile(dir, output, item);
        }
    }
}

function processFile(inDir, outDir, fn) {
    print("\n Processing: " + fn);
    op = inDir + File.separator + fn;
    open(op);
    imageTitle = getTitle();

    getVoxelSize(dx, dy, dzs, unit);
    Stack.getDimensions(w, h, ch, sl, fr);
    print("  → Channels: " + ch + ", Slices: " + sl + ", Frames: " + fr);

    // Validate PSFs
    if (lengthOf(psf_paths) < ch) {
        print("ERROR: Not enough PSFs for all channels.");
        close(); return;
    }

    // Init CLIJ
    run("CLIJ2 Macro Extensions", "cl_device=");
    Ext.CLIJ2_clear();

    // Deconvolve per channel
    deconTitles = "";
        // Deskew
    angle = 32.8;
    dz = sin(angle * PI / 180) * dzs;
    df = (cos(angle * PI / 180) * dzs) / dx;
    padding = w + sl * df;
    for (c = 1; c <= ch; c++) {
    	selectWindow(imageTitle);
        Stack.setPosition(c, 1, 1);
        waitForUser;
        Ext.CLIJ2_pushCurrentZStack(imageTitle);

        open(psf_paths[c - 1]);
        psfTitle = getTitle();
        Ext.CLIJ2_push(psfTitle);

        dn = "decon_c" + c;
        Ext.CLIJx_deconvolveRichardsonLucyFFT(imageTitle, psfTitle, dn, num_iterations);
        //Find equivalent CLIJ command
        //run("Canvas Size...", "width=" + padding + " height=" + h + " position=Center-Left zero");
        Ext.CLIJ2_affineTransform3D(dn, deskewed, "shearXZ=" + (-df));
    Ext.CLIJ2_pull(deskewed);
    deskewedTitle = "deskewed_" + dn +"_"+ fn;
    rename(deskewedTitle);
        deconTitles += "c" + c + "=" + deskewedTitle + " ";
        print(deconTitles);
        close(psfTitle);
    }

    // Merge all deconvolved channels
    run("Merge Channels...", deconTitles + "create");
    merged = "merged_" + fn;
    rename(merged);
    print("   Decon done → " + merged);



    

    // Save
    outPath = outDir + File.separator + deskewedTitle;
    print("  Saving to: " + outPath);
    saveAs("Tiff", outPath);

    close("*");
    Ext.CLIJ2_clear();
    print( "Done: " + fn);
}

// Start
processFolder(input);
print("All files processed.");


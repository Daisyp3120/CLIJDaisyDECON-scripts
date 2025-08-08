#@ File (label = "Input folder", style = "directory") input
#@ File (label = "Output folder", style = "directory") output
#@ File[] (label="PSF files (one per channel)", style="files") psf_paths
#@ Integer (label = "Iterations", value = 10) num_iterations
#@ String (choices={"3i", "Zeiss"}, label = "Microscope Type") scope
#@ String (label = "File suffix", value = ".tif") suffix

print("Starting batch process...");
print("Input:  " + input);
print("Output: " + output);

all_pre_times = newArray();
all_mid_times = newArray();
all_post_times = newArray();
all_total_times = newArray();
all_filenames = newArray();

start_batch = getTime();

// Run processing
processFolder(input);

// Summary
print("\nBatch Summary");
total_batch_time = (getTime() - start_batch) / 1000;
print("Total batch time: " + total_batch_time + " s");


function processFolder(dir) {
    files = getFileList(dir);
    for (i = 0; i < files.length; i++) {
        item = files[i];
        full = dir + File.separator + item;
        if (File.isDirectory(full)) {
            processFolder(full);
        } else {
            itemLower = toLowerCase(item);
            suffixLower = toLowerCase(suffix);
            if (endsWith(itemLower, suffixLower)) {
                print("  Processing file: " + item);
                times = processFile(dir, output, item);

                print("Returned times array received.");

                // no typeof check here because ImageJ macro doesn't support typeof()

                all_pre_times = Array.concat(all_pre_times, times[0]);
                all_mid_times = Array.concat(all_mid_times, times[1]);
                all_post_times = Array.concat(all_post_times, times[2]);
                all_total_times = Array.concat(all_total_times, times[3]);
                all_filenames = Array.concat(all_filenames, times[4]);
            } else {
                print("Skipping: " + item);
            }
        }
    }

    Array.getStatistics(all_pre_times, min_pre, max_pre, mean_pre, std_pre);
    Array.getStatistics(all_mid_times, min_mid, max_mid, mean_mid,std_mid);
    Array.getStatistics(all_post_times, min_post, max_post, mean_post, std_post);
    Array.getStatistics(all_total_times, min_total, max_total, mean_total, std_total);

    print("Pre Time mean: " + mean_pre);
    print("Pre Time stdDev: " + std_pre);
    print("Mid Time mean: " + mean_mid);
    print("Mid Time stdDev: " + std_mid);
    print("Post Time mean: " + mean_post);
    print("Post Time stdDev: " + std_post);
    print("Total Time mean: " + mean_total);
    print("Total Time stdDev: " + std_total);

    print("All files processed.");
}

function processFile(inDir, outDir, fn) {
    file_start = getTime();
    op = inDir + File.separator + fn;

    // Pre-processing
    pre_start = getTime();
    run("Bio-Formats Importer", "open=["+op+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    imageTitle = getTitle();
    getVoxelSize(dx, dy, dzs, unit);
    Stack.getDimensions(w, h, ch, sl, fr);
    print("  → Channels: " + ch + ", Slices: " + sl + ", Frames: " + fr);
    print("  → PSF files provided: " + lengthOf(psf_paths));

    if (lengthOf(psf_paths) < ch) {
        print("  ERROR: Not enough PSFs for all channels.");
        close();
        return newArray(0,0,0,0,fn);
    }

    run("CLIJ2 Macro Extensions", "cl_device=");
    Ext.CLIJ2_clear();
    pre_end = getTime();
    pre_time = (pre_end - pre_start) / 1000;

    // Mid-processing
    mid_start = getTime();
    deconTitles = newArray();
    angle = 32.8;
    dz = sin(angle * PI / 180) * dzs;
    df = (cos(angle * PI / 180) * dzs) / dx;
    padding = w + sl * df;

    // Array for each timepoint to hold channel images for merging later
    merged_per_time = newArray();

    for (t = 1; t <= fr; t++) {
        per_ch_titles = newArray();

        for (c = 1; c <= ch; c++) {
            selectWindow(imageTitle);
            Stack.setPosition(c, 1, t);
            Ext.CLIJ2_pushCurrentZStack(imageTitle);

            open(psf_paths[c - 1]);
            psfTitle = getTitle();
            Ext.CLIJ2_push(psfTitle);

            decon = "decon_c" + c + "_t" + t;
            Ext.CLIJx_deconvolveRichardsonLucyFFT(imageTitle, psfTitle, decon, num_iterations);
            Ext.CLIJ2_pull(decon);
            rename(decon);
            per_ch_titles = Array.concat(per_ch_titles, decon);
            close(psfTitle);
        }

        // Merge channels for this timepoint if more than 1 channel
        if (lengthOf(per_ch_titles) > 1) {
            merge_cmd = "";
            for (i = 0; i < lengthOf(per_ch_titles); i++) {
                merge_cmd += " c" + (i+1) + "=[" + per_ch_titles[i] + "]";
            }
            run("Merge Channels...", merge_cmd + " create keep");
            mergedTitle = "merged_T" + t;
            rename(mergedTitle);

            // Close individual channel windows
            for (i = 0; i < lengthOf(per_ch_titles); i++) {
                if (isOpen(per_ch_titles[i])) close(per_ch_titles[i]);
            }

            merged_per_time = Array.concat(merged_per_time, mergedTitle);
        } else if (lengthOf(per_ch_titles) == 1) {
            // Only one channel, no merge needed
            merged_per_time = Array.concat(merged_per_time, per_ch_titles[0]);
        } else {
            print("  WARNING: No channels found for timepoint " + t);
        }
    }
    mid_end = getTime();
    mid_time = (mid_end - mid_start) / 1000;

    // Concatenate merged images across time if more than one
    if (lengthOf(merged_per_time) > 1) {
        concat_cmd = "title=FinalStack";
        for (i = 0; i < lengthOf(merged_per_time); i++) {
            concat_cmd += " image" + (i + 1) + "=[" + merged_per_time[i] + "]";
        }
        run("Concatenate...", concat_cmd);
        selectWindow("FinalStack");
    } else if (lengthOf(merged_per_time) == 1) {
        selectWindow(merged_per_time[0]);
        rename("FinalStack");
    } else {
        print("ERROR: No images to concatenate.");
        close("*");
        Ext.CLIJ2_clear();
        post_time = 0;
        total_time = (getTime() - file_start) / 1000;
        return newArray(pre_time, mid_time, post_time, total_time, fn);
    }

    // Fix hyperstack dimensions
    run("Stack to Hyperstack...",
        "order=xyczt(default) channels=1 slices=" + sl + " frames=" + fr + " display=Color");

    // Save output
    outPath = outDir + File.separator + "output_decon_" + fn;
    print("Saving output to: " + outPath);
    saveAs("Tiff", outPath);
    close("*");
    Ext.CLIJ2_clear();

    post_end = getTime();
    post_time = (post_end - mid_end) / 1000;
    total_time = (post_end - file_start) / 1000;

    print("Pre:   " + pre_time + " s");
    print("Mid:   " + mid_time + " s");
    print("Post:  " + post_time + " s");
    print("Total: " + total_time + " s");
    print("Done:  " + fn);

    return newArray(pre_time, mid_time, post_time, total_time, fn);
}

// Function to compute standard deviation (not currently used)
function stddev(arr, mean_val) {
    n = lengthOf(arr);
    if (n == 0) return NaN;
    sum_sq = 0;
    for (i = 0; i < n; i++) {
        sum_sq += pow(arr[i] - mean_val, 2);
    }
    return sqrt(sum_sq / n);
}

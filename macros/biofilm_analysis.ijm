print("Starting BEL/FEUP macro for biofilm analysis");

print("Please select a working directory");

input = getDirectory("Choose a Directory")
list = getFileList(input);


print("Please select a C-scan image (3D)");

Dialog.create("Analyze biofilm");
Dialog.addChoice("Open File:", list);
Dialog.show();
filename = Dialog.getChoice();

print("Opening file " + filename);
open(filename);


print("Please ajust the threshold for the binarization process");

Dialog.create("Adjust threshold");
Dialog.addNumber("Min Intensity:", 200);
Dialog.addNumber("Max Intensity:", 255);
Dialog.show();

mincolor = Dialog.getNumber();
maxcolor = Dialog.getNumber();

print("Starting binarization process");
setThreshold(mincolor, maxcolor);
call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
run("Convert to Mask", "method=Moments background=Dark calculate black");
print("Binarization process concluded");

print("Starting reslicing process");
getPixelSize(unit, pw, ph, pd);
run("Reslice [/]...", "output="+ph+" start=Top avoid");
print("Reslicing process concluded");

print("Please define the parameters");
suportslice = nSlices;
bulkslice = 1;
Dialog.create("Adjust the parameters");
Dialog.addNumber("Porosity Threshold:", 0.975);
Dialog.addNumber("Number of samples:", 10);
Dialog.addNumber("Suport slice:", suportslice);
Dialog.addNumber("Last Bulk slice:", bulkslice);
Dialog.show();
porosity_threshold = Dialog.getNumber();
n = Dialog.getNumber();
suportslice = Dialog.getNumber();
bulkslice = Dialog.getNumber();



print("Starting analysis process");


title1 = "Biofilm Analysis per Sample"; 
title2 = "["+title1+"]"; 
f=title2; 
run("New... ", "name="+title2+" type=Table"); 
print(f,"\\Headings:Sample\tX\tY\tporosity\tthickness"); 

thickness_array = newArray(0);
porosity_array = newArray(0);
v0_array = newArray(0);
v1_array = newArray(0);

width = getWidth;
height = getHeight;
getPixelSize(unit, pw, ph, pd);

searchValue_0 = 0;
searchValue_1 = 255;

countValue_all_0 = 0;
countValue_all_1 = 0;

j = 0;
for( ix = 1; ix <= n; ix++ ) {

	for( iy = 1; iy <= n; iy++ ) {

		j += 1;

		nx = round(getWidth()/n);
		ix_ini = (ix-1)*nx;
		ix_fini = (ix)*nx;

		ny = round(getHeight()/n);
		iy_ini = (iy-1)*ny;
		iy_fini = (iy)*ny;

		makeRectangle(ix_ini, iy_ini, ix_fini-ix_ini, iy_fini-iy_ini);

		thickness = 0;

        countValue_all_0 = 0;
        countValue_all_1 = 0;
		for( i = suportslice; i >= bulkslice; i-- ) {

		    setSlice(i);

	        countValue_0 = 0;
	        countValue_1 = 0;
			for( x = ix_ini; x <= ix_fini; x++ ) {
				for( y = iy_ini; y <= iy_fini; y++ ) {
					    if( getPixel(x,y) == searchValue_0 ) {
                                countValue_0 += 1;
                        }
                        if( getPixel(x,y) == searchValue_1 ) {
                                countValue_1 += 1;
                        }
                }

			}

	        slice_porosity = countValue_0/(countValue_0+countValue_1);

			if( slice_porosity > porosity_threshold ){
				break;
			}

		    thickness += pd;
	        countValue_all_0 += countValue_0;
			countValue_all_1 += countValue_1;

		}
		local_porosity = countValue_all_0/(countValue_all_0+countValue_all_1);

		thickness_array = Array.concat(thickness_array,thickness);
		porosity_array = Array.concat(porosity_array,local_porosity);
		v0_array = Array.concat(v0_array,countValue_all_0);
		v1_array = Array.concat(v1_array,countValue_all_1);

		print(f,j+"\t"+ix+"\t"+iy+"\t"+local_porosity+"\t"+thickness); 

	}
}

Array.getStatistics(thickness_array, thickness_min, thickness_max, thickness_mean, thickness_stdDev)
Array.getStatistics(porosity_array, porosity_min, porosity_max, porosity_mean, porosity_stdDev)

roughness = 0
for(i=0;i<thickness_array.length;i++){
    roughness += abs(thickness_array[i] - thickness_mean);
}
roughness = roughness/thickness_array.length


title1 = "Biofilm Statistics"; 
title2 = "["+title1+"]"; 
f=title2;
run("New... ", "name="+title2+" type=Table");
print(f,"\\Headings:Parameter or Variable\tValue\tUnit"); 
print(f,"Porosity Threshold"+"\t"+porosity_threshold+"\t"+"-"); 
print(f,"Samples"+"\t"+n*n+"\t"+"-");
print(f,"Width"+"\t"+width*pw+"\t"+unit);
print(f,"Length"+"\t"+height*ph+"\t"+unit);
print(f,"Height"+"\t"+pd*(suportslice-bulkslice+1)+"\t"+unit);
print(f,"Biofilm Porosity"+"\t"+porosity_mean+"+/-"+porosity_stdDev+"\t"+"-"); 
print(f,"Biofilm Thickness"+"\t"+thickness_mean+"+/-"+thickness_stdDev+"\t"+unit); 
print(f,"Biofilm Roughness"+"\t"+roughness+"+/-"+thickness_stdDev+"\t"+unit); 

print("Analysis process concluded");


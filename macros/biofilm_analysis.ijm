print("Starting BEL/FEUP macro for biofilm analysis");

print("Please select a working directory");

input = getDirectory("Choose a Directory")
function filter(i, name) {
    if (endsWith(name,"/")) return false;
    if (!endsWith(name,".tif")) return false;
    return true;
}
list = getFileList(input);
list2 = newArray(0);
for (i=0; i<list.length; i++) {
    showProgress(i, list.length);
    if (filter(i, list[i])) {
    	list2 = Array.concat(list2,list[i]);
    }
}

print("Please select a C-scan image (3D)");

Dialog.create("Analyze biofilm");
Dialog.addChoice("Open File:", list2);
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
Dialog.addNumber("Porosity Threshold:", 0.99);
Dialog.addNumber("Number of samples:", 10);
Dialog.addNumber("Suport slice:", suportslice);
Dialog.addNumber("Last Bulk slice:", bulkslice);
Dialog.show();
porosity_threshold = Dialog.getNumber();
n = Dialog.getNumber();
suportslice = Dialog.getNumber();
bulkslice = Dialog.getNumber();


k = round(random*999)

print("Starting analysis process (#" + k + ")");


title1 = "Biofilm Analysis per Sample (#" + k + ")"; 
title2 = "["+title1+"]"; 
f=title2; 
run("New... ", "name="+title2+" type=Table"); 
print(f,"\\Headings:Sample\tX\tY\tZI\tZF\tporosity\tthickness\ttransitions"); 

thickness_array = newArray(0);
porosity_array = newArray(0);
v0_array = newArray(0);
v1_array = newArray(0);
x_array = newArray(0);
y_array = newArray(0);

width = getWidth;
height = getHeight;
getPixelSize(unit, pw, ph, pd);

searchValue_0 = 0;
searchValue_1 = 255;

countValue_all_0 = 0;
countValue_all_1 = 0;


// Defining position sets
setSlice(suportslice);
x_set = newArray(0);
for( ix = 1; ix <= n; ix++ ) {
	nx = round(getWidth()/n);
	ix_ini = (ix-1)*nx;
	ix_fini = (ix)*nx;
	x_set = Array.concat(x_set,pw*(ix_ini+ix_fini)*0.5);
}
y_set = newArray(0);
for( iy = 1; iy <= n; iy++ ) {
	ny = round(getHeight()/n);
	iy_ini = (iy-1)*ny;
	iy_fini = (iy)*ny;
	y_set = Array.concat(y_set,ph*(iy_ini+iy_fini)*0.5);
}

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

		z_ini = 0.0;
		z_fini = 0.0;
		z = 0.0;
		transitions = 0;

		// position: 0 for under biofilm, 1 for in biofilm and 2 for over biofilm
		position = 0;
		

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

			if( ( position == 0 ) & ( slice_porosity < porosity_threshold ) ){
				position = 1;
				z_ini = z;
				transitions += 1;
			}else {
				if( ( position == 1 ) & ( slice_porosity >= porosity_threshold ) ){
					position = 2;
					z_fini = z;
			        countValue_buffer_0 = 0;
			        countValue_buffer_1 = 0;
					transitions += 1;
				}else {
					if( ( position == 2 ) & ( slice_porosity < porosity_threshold ) ){
						position = 1;
						countValue_all_0 += countValue_buffer_0;
						countValue_all_1 += countValue_buffer_1;
						transitions += 1;
					}
				}
			}

		    z += pd;

		    if ( position == 0 ){
	        	countValue_all_0 += 0;
				countValue_all_1 += 0;
		    }
		    if ( position == 1 ){
	        	countValue_all_0 += countValue_0;
				countValue_all_1 += countValue_1;
		    }
			if( position == 2 ){
		        countValue_buffer_0 += countValue_0;
		        countValue_buffer_1 += countValue_1;
				// break;
			}

			print(f,j+"\t"+ix+"\t"+iy+"\t"+z+"\t"+z+"\t"+slice_porosity+"\t-\t"+transitions); 


		}

		if ( (countValue_all_0+countValue_all_1)>0 ){
			
			local_porosity = countValue_all_0/(countValue_all_0+countValue_all_1);
	
	        thickness = z_fini - z_ini;
	
			thickness_array = Array.concat(thickness_array,thickness);
			porosity_array = Array.concat(porosity_array,local_porosity);
			v0_array = Array.concat(v0_array,countValue_all_0);
			v1_array = Array.concat(v1_array,countValue_all_1);
			transitions_array = Array.concat(transitions_array,transitions);
	
			x_array = Array.concat(x_array,pw*(ix_ini+ix_fini)*0.5);
			y_array = Array.concat(y_array,ph*(iy_ini+iy_fini)*0.5);
	
			print(f,j+"\t"+ix+"\t"+iy+"\t"+z_ini+"\t"+z_fini+"\t"+local_porosity+"\t"+thickness+"\t"+transitions); 
		}
	}
}
saveAs("Text", "biofilm_data_" + k + ".csv");

Array.getStatistics(x_array, x_min, x_max, x_mean, x_stdDev);
Array.getStatistics(y_array, y_min, y_max, y_mean, y_stdDev);
Array.getStatistics(thickness_array, thickness_min, thickness_max, thickness_mean, thickness_stdDev);
Array.getStatistics(porosity_array, porosity_min, porosity_max, porosity_mean, porosity_stdDev);

total_number_points = j

roughness = 0
for(i=0;i<thickness_array.length;i++){
    roughness += abs(thickness_array[i] - thickness_mean);
}
roughness = roughness/thickness_array.length;
roughness_stdDev = thickness_stdDev * sqrt(thickness_array.length) / thickness_array.length;


title1 = "Biofilm Statistics (#" + k + ")"; 
title2 = "["+title1+"]"; 
f=title2;
run("New... ", "name="+title2+" type=Table");
print(f,"\\Headings:Parameter or Variable\tValue\tUnit"); 
print(f,"Min Intensity (0-255)"+"\t"+mincolor+"\t"+"-"); 
print(f,"Max Intensity (0-255)"+"\t"+maxcolor+"\t"+"-"); 
print(f,"Porosity Threshold"+"\t"+porosity_threshold+"\t"+"-"); 
print(f,"Samples"+"\t"+n*n+"\t"+"-");
print(f,"Width"+"\t"+width*pw+"\t"+unit);
print(f,"Length"+"\t"+height*ph+"\t"+unit);
print(f,"Height"+"\t"+pd*(suportslice-bulkslice+1)+"\t"+unit);
print(f,"Biofilm Porosity"+"\t"+porosity_mean+"+/-"+porosity_stdDev+"\t"+"-"); 
print(f,"Biofilm Thickness"+"\t"+thickness_mean+"+/-"+thickness_stdDev+"\t"+unit); 
print(f,"Biofilm Maximum Thickness"+"\t"+thickness_max+"+/-"+thickness_stdDev+"\t"+unit); 
print(f,"Biofilm Minimum Thickness"+"\t"+thickness_min+"+/-"+thickness_stdDev+"\t"+unit); 
print(f,"Biofilm Roughness"+"\t"+roughness+"+/-"+roughness_stdDev+"\t"+unit); 
saveAs("Text", "biofilm_statistics_" + k + ".csv");


requires("1.52f");
values = newArray(2,3.01,3,3,4,4);
Plot.create("Biofilm Thickness (#" + k + ")", "Thickness ("+unit+")", "Frequency");
Plot.setColor("blue", "#ccccff");
binWidth = (thickness_max - thickness_min)/(0.5*n);//use 0 for auto-binning
binCenter = 0;
Plot.addHistogram(thickness_array, binWidth, binCenter);
Plot.setLimits(NaN, NaN, 0, NaN);
Plot.show;

requires("1.52f");
values = newArray(2,3.01,3,3,4,4);
Plot.create("Biofilm Porosity (#" + k + ")", "Porosity", "Frequency");
Plot.setColor("red", "#FF9100");
binWidth = (porosity_max - porosity_min)/(0.5*n);//use 0 for auto-binning
binCenter = 0;
Plot.addHistogram(porosity_array, binWidth, binCenter);
Plot.setLimits(NaN, NaN, 0, NaN);
Plot.show;

requires("1.52f");
Plot.create("Porosity and Thickness (#" + k + ")", "Porosity", "Thickness ("+unit+")");
Plot.setColor("black");
Plot.add("circle", porosity_array, thickness_array);
Plot.setColor("blue");
x_aux = newArray(1);
y_aux = newArray(1);
x_aux[0] = porosity_mean;
y_aux[0] = thickness_mean;
Plot.add("circle", x_aux, y_aux);
x_aux = newArray(2);
y_aux = newArray(2);
x_aux[0] = porosity_mean-porosity_stdDev;
x_aux[1] = porosity_mean+porosity_stdDev;
y_aux[0] = thickness_mean;
y_aux[1] = thickness_mean;
Plot.add("line", x_aux, y_aux);
x_aux = newArray(2);
y_aux = newArray(2);
x_aux[0] = porosity_mean;
x_aux[1] = porosity_mean;
y_aux[0] = thickness_mean-thickness_stdDev;
y_aux[1] = thickness_mean+thickness_stdDev;
Plot.add("line", x_aux, y_aux);
Plot.setLimits(porosity_min, porosity_max, thickness_min, thickness_max);
Plot.show();


requires("1.52f");
name = "Porosity profiles 1 (#" + k + ")";
Plot.create(name, "Width ("+unit+")", "Porosity (-)");
Plot.setColor("black");
numSlices = 0;
for(i=0;i<y_set.length;i++){
	numSlices += 1;
	x_aux = newArray(0);
	y_aux = newArray(0);
	for(j=0;j<porosity_array.length;j++){
	    if ( y_array[j] == y_set[i] ){
	    	x_aux = Array.concat(x_aux,x_array[j]);
	    	y_aux = Array.concat(y_aux,porosity_array[j]);
	    }
	}
	Plot.setColor("red", "#FF9100");
	Plot.add("filled", x_aux, y_aux);
	Plot.setLimits(x_min, x_max, 0, 1);
	Plot.setXYLabels("Width ("+unit+")", "Porosity (-)"+" for length " + d2s(y_set[i], 5) + " " + unit);
	Plot.appendToStack;
}
Plot.show();
rename(name);
run("Stack to Hyperstack...", "slices=&numSlices frames=1")	


requires("1.52f");
name = "Porosity profiles 2 (#" + k + ")";
Plot.create(name, "Length ("+unit+")", "Porosity (-)");
Plot.setColor("black");
numSlices = 0;
for(i=0;i<x_set.length;i++){
	numSlices += 1;
	x_aux = newArray(0);
	y_aux = newArray(0);
	for(j=0;j<porosity_array.length;j++){
	    if ( x_array[j] == x_set[i] ){
	    	x_aux = Array.concat(x_aux,y_array[j]);
	    	y_aux = Array.concat(y_aux,porosity_array[j]);
	    }
	}
	Plot.setColor("red", "#FF9100");
	Plot.add("filled", x_aux, y_aux);
	Plot.setLimits(y_min, y_max, 0, 1);
	Plot.setXYLabels("Length ("+unit+")", "Porosity (-)"+" for width " + d2s(x_set[i], 5) + " " + unit);
	Plot.appendToStack;
}
Plot.show();
rename(name);
run("Stack to Hyperstack...", "slices=&numSlices frames=1")	


requires("1.52f");
name = "Thickness profiles 3 (#" + k + ")";
Plot.create(name, "Width ("+unit+")", "Thickness ("+unit+")");
Plot.setColor("black");
numSlices = 0;
for(i=0;i<y_set.length;i++){
	numSlices += 1;
	x_aux = newArray(0);
	y_aux = newArray(0);
	for(j=0;j<thickness_array.length;j++){
	    if ( y_array[j] == y_set[i] ){
	    	x_aux = Array.concat(x_aux,x_array[j]);
	    	y_aux = Array.concat(y_aux,thickness_array[j]);
	    }
	}
	Plot.setColor("blue", "#ccccff");
	Plot.add("filled", x_aux, y_aux);
	Plot.setLimits(x_min, x_max, 0, thickness_max+0.05*thickness_mean);
	Plot.setXYLabels("Width ("+unit+")", "Thickness ("+unit+")"+" for length " + d2s(y_set[i], 5) + " " + unit);
	Plot.appendToStack;
}
Plot.show();
rename(name);
run("Stack to Hyperstack...", "slices=&numSlices frames=1")	


requires("1.52f");
name = "Thickness profiles 4 (#" + k + ")
Plot.create(name, "Length ("+unit+")", "Thickness ("+unit+")");
Plot.setColor("black");
numSlices = 0;
for(i=0;i<x_set.length;i++){
	numSlices += 1;
	x_aux = newArray(0);
	y_aux = newArray(0);
	for(j=0;j<thickness_array.length;j++){
	    if ( x_array[j] == x_set[i] ){
	    	x_aux = Array.concat(x_aux,y_array[j]);
	    	y_aux = Array.concat(y_aux,thickness_array[j]);
	    }
	}
	Plot.setColor("blue", "#ccccff");
	Plot.add("filled", x_aux, y_aux);
	Plot.setLimits(y_min, y_max, 0, thickness_max+0.05*thickness_mean);
	Plot.setXYLabels("Length ("+unit+")", "Thickness ("+unit+")"+" for width " + d2s(x_set[i], 5) + " " + unit);
	Plot.appendToStack;
}
Plot.show();
rename(name);
run("Stack to Hyperstack...", "slices=&numSlices frames=1")	



print("Analysis process concluded");

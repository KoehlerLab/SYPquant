### Import required packages
import sys
import glob
import os
import ast
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import tifffile as TIFFfile
from skimage.io import imread
from skimage.measure import label
from skimage.measure import regionprops
from skimage.measure import marching_cubes
from skimage.measure import mesh_surface_area

from scipy import ndimage as ndi
exec("from skimage.filters import threshold_" + sys.argv[5])

### Parameters
filename='/home/xxx_518.tif' ## full path to one-channel image of axis staining
AXIS_channel='518' #Channel you want to segment (for Olympus channel number for airyscan excitation wavelength)
SC_channel='603' #Channel you want to quantify intensity (for Olympus channel number for airyscan excitation wavelength)
DAPI_channel='455'
threshold_method='mean' 
stitched='True' #or False
Stitched_MIP_folder='/home/MaxTiffs' # directory to *TileConfiguration.registered.txt file for the image if using stitched files
### Calculate thresholding mask per nucleus (requires cellpose nuclei segmentation)

print('Analysing file:'+ filename)

im=imread(filename)

path=os.path.dirname(filename)
name=filename.replace(path + "/", "")
name=name.replace("_" + AXIS_channel + ".tif", "")

VSI_Directory=os.path.dirname(os.path.dirname(os.path.dirname(path)))

Position=os.path.dirname(path).replace(os.path.dirname(os.path.dirname(path)) + "/", "")

im_mask=imread(path + "/" + name + "_cp_masks.tif")

raw_mask=np.load(path + "/" + name + "_segm_raw.npy")

htp3_mask=np.zeros(im.shape)
for i, j in enumerate(np.unique(im_mask[im_mask>0])):
    t=eval("threshold_"+threshold_method)(im[im_mask==j])
    nucleus_mask=np.where(np.logical_and(im_mask==j, im>t),1,0)
    htp3_mask[nucleus_mask>0]=1



TIFFfile.imwrite(path + "/" + name + "_SC_mask.tif", htp3_mask.astype('uint16'))

print('Finished generating SC mask for ' + filename)

### Calculate SYP intensity on and off the HTP3 mask per nucleus (normalise to SYP background)
#Load images
SYP_image=(path + "/" + name + "_" + SC_channel + ".tif")

im_cp_mask=im_mask
im_SC_mask=htp3_mask
im_SYP=imread(SYP_image)

## Calculate background in SYP channel by applying a dilation of 50 on cellpose segmented binary image
#and substracting a dilation of 5 of cellpose binary image
#to remove pixels that are still true signal near the objects

cp_mask_binary=np.copy(raw_mask)
cp_mask_binary[cp_mask_binary>0]=1

cp_mask_binary_dilation1=ndi.binary_dilation(cp_mask_binary, iterations=10, border_value=0)
cp_mask_binary_dilation2=ndi.binary_dilation(cp_mask_binary, iterations=50, border_value=0)

gonad_Bkg=np.zeros(im_cp_mask.shape)
gonad_Bkg=np.where(np.logical_and(cp_mask_binary_dilation2==True, cp_mask_binary_dilation1==False),1,0)

TIFFfile.imwrite(path + "/" + name + "_cp_mask_dilation1.tif", cp_mask_binary_dilation1.astype('uint16'))
TIFFfile.imwrite(path + "/" + name + "_cp_mask_dilation2.tif", cp_mask_binary_dilation2.astype('uint16'))


print('Finished SC mask dilations: dilation 1: 10 iterations, dilation 2: 50 iterations')

MaxBackground=np.max(im_SYP[gonad_Bkg==1])
MinBackground=np.min(im_SYP[gonad_Bkg==1])
MeanBackground=np.mean(im_SYP[gonad_Bkg==1])
MedianBackground=np.median(im_SYP[gonad_Bkg==1])
    
print('Mean SYP Background:', MeanBackground)
print('Median SYP Background:', MedianBackground)
    
### Calculate DAPI intensity to correct for z plane
#Load images
DAPI_image=(path + "/" + name + "_" + DAPI_channel + ".tif")

im_DAPI=imread(DAPI_image)
MaxDAPIBkg=np.max(im_DAPI[gonad_Bkg==1])
MinDAPIBkg=np.min(im_DAPI[gonad_Bkg==1])
MeanDAPIBkg=np.mean(im_DAPI[gonad_Bkg==1])
MedianDAPIBkg=np.median(im_DAPI[gonad_Bkg==1])

print('Mean DAPI Background:', MeanDAPIBkg)
print('Median DAPI Background:', MedianDAPIBkg)


#Read Metadata and get pixel size
with TIFFfile.TiffFile(SYP_image) as tif:
    metadata=tif.imagej_metadata
info=metadata['Info']

start_i=info.find('PhysicalSizeX')+16
pixelSizeX=float(info[start_i:(info[start_i:].find(',')+start_i)])
start_i=info.find('PhysicalSizeY')+16
pixelSizeY=float(info[start_i:(info[start_i:].find(',')+start_i)])
start_i=info.find('PhysicalSizeZ')+16
pixelSizeZ=float(info[start_i:(info[start_i:].find(',')+start_i)])

VoxelV=pixelSizeX*pixelSizeY*pixelSizeZ
Voxel_array=[pixelSizeZ, pixelSizeY, pixelSizeX]

print('Voxel size (z,y,x):', Voxel_array)

#Generate csv file containing all the information regarding each segmented nucleus (nucleus id, centroid position, nuclear volume, etc) 

results=pd.DataFrame()
regions=regionprops(im_mask)
h, w, c = im_cp_mask.shape

for i,j in enumerate(np.unique(im_cp_mask[im_cp_mask>0])):
    tempM=np.copy(im_cp_mask)
    tempM[im_cp_mask!=j]=0

    vts, fs, ns, cs = marching_cubes(tempM, level=None,gradient_direction='ascent',spacing=Voxel_array,method='lewiner')
    volume=regions[i].area*VoxelV
    surface_a=mesh_surface_area(vts,fs)

    
    filter_SYP=im_SYP[im_cp_mask==j] #sc pixels in cellpose mask

    filter_SC=np.zeros(im_SYP.shape)
    filter_SC=np.where(np.logical_and(im_cp_mask==j, im_SC_mask==1), 1, 0)
    filter_SYP_SC=im_SYP[filter_SC==1]  #sc pixels within the nucleus (cellpose mask) AND on the SC (htp3 mask)

    filter_Nuc=np.zeros(im_SYP.shape)
    filter_Nuc=np.where(np.logical_and(im_cp_mask==j, im_SC_mask==0), 1, 0)
    filter_SYP_Nuc=im_SYP[filter_Nuc==1] #sc pixels within the nucleus (cellpose mask)  AND outside the SC (htp3_mask)

    filter_DAPI=im_DAPI[im_cp_mask==j] #dapi pixels in cellpose mask
    filter_DAPI_SC=im_DAPI[filter_SC==1]  #sc pixels within the nucleus (cellpose mask) AND on the SC (htp3 mask)
    

    zmin,ymin,xmin,zmax,ymax,xmax=regions[i].bbox
    ymin_tile=0
    xmin_tile=0
    ymax_tile=w
    xmax_tile=c
    if (ymin in range(ymin_tile,ymin_tile+2)
        or xmin in range(xmin_tile,xmin_tile+2)
        or ymax in range(ymax_tile-2,ymax_tile)
        or xmax in range(xmax_tile-2,xmax_tile)):
    
        filter_bbox='Yes'
    else:
        filter_bbox='No'

    

    temp=pd.DataFrame(
        {'Folder': [VSI_Directory],
         'Position': [Position],
         'Tile': [(name + '.vsi')],
         'TileShape': [im_cp_mask.shape],
         'Nucleus_ID': [regions[i].label],
         'Nucleus_BBox': [regions[i].bbox],
         'Nucleus_BBox_TileEdge': [filter_bbox], #if Yes it means it is a cropped nucleus and should be removed during filtering
         'Cz': [regions[i].centroid[0]],
         'Cy': [regions[i].centroid[1]],
         'Cx': [regions[i].centroid[2]],
         'VolumeNucleus(um3)': [regions[i].area*VoxelV], #a good cutoff could be >15 <35
         'Sphericity': [(((math.pi)**(1/3))*((6*volume)**(2/3)))/(surface_a)], #Ivana uses 0.75 but I may have to go down
         'RadiusNucleus(um)':[np.power(3*(regions[i].area*VoxelV)/(4*np.pi),1/3)],  #a good cutoff could be >1.5 
         'VolumeSC(um3)': [len(filter_SYP_SC)*VoxelV],
         'MaxI_SC_SYP': [np.max(filter_SYP_SC)],
         'MinI_SC_SYP': [np.min(filter_SYP_SC)],
         'MeanI_SC_SYP': [np.mean(filter_SYP_SC)],
         'MedianI_SC_SYP': [np.median(filter_SYP_SC)],
         'MaxI_Nuc_SYP': [np.max(filter_SYP_Nuc)],
         'MinI_Nuc_SYP': [np.min(filter_SYP_Nuc)],
         'MeanI_Nuc_SYP': [np.mean(filter_SYP_Nuc)],
         'MedianI_Nuc_SYP': [np.median(filter_SYP_Nuc)],
         'MaxI_TOTAL_SYP': [np.max(filter_SYP)],
         'MinI_TOTAL_SYP': [np.min(filter_SYP)],
         'MeanI_TOTAL_SYP': [np.mean(filter_SYP)],
         'MedianI_TOTAL_SYP': [np.median(filter_SYP)],
         'I_Total_SYP': [np.sum(filter_SYP)],
         'I_Total_SYPSC': [np.sum(filter_SYP_SC)],
         'I_Total_SYPNuc': [np.sum(filter_SYP_Nuc)],
         '#CP_Pixels': [len(filter_SYP)],
         '#SC_Pixels': [len(filter_SYP_SC)],
         '#Nuc_Pixels': [len(filter_SYP_Nuc)],
         'MaxI_TOTAL_DAPI': [np.max(filter_DAPI)],
         'MinI_TOTAL_DAPI': [np.min(filter_DAPI)],
         'MeanI_TOTAL_DAPI': [np.mean(filter_DAPI)],
         'MedianI_TOTAL_DAPI': [np.median(filter_DAPI)],
         'I_Total_DAPI': [np.sum(filter_DAPI)],
         'I_Total_DAPISC': [np.sum(filter_DAPI_SC)]
     })
    
    results = pd.concat([results, temp])


results['MaxI_SYP_Background']=MaxBackground
results['MinI_SYP_Background']=MinBackground
results['MeanI_SYP_Background']=MeanBackground
results['MedianI_SYP_Background']=MedianBackground
results['MaxI_SC_SYP-Bkg']=results['MaxI_SC_SYP']-results['MaxI_SYP_Background']
results['MinI_SC_SYP-Bkg']=results['MinI_SC_SYP']-results['MinI_SYP_Background']
results['MeanI_SC_SYP-Bkg']=results['MeanI_SC_SYP']-results['MeanI_SYP_Background']
results['MedianI_SC_SYP-Bkg']=results['MedianI_SC_SYP']-results['MedianI_SYP_Background']
results['MaxI_Nuc_SYP-Bkg']=results['MaxI_Nuc_SYP']-results['MaxI_SYP_Background']
results['MinI_Nuc_SYP-Bkg']=results['MinI_Nuc_SYP']-results['MinI_SYP_Background']
results['MeanI_Nuc_SYP-Bkg']=results['MeanI_Nuc_SYP']-results['MeanI_SYP_Background']
results['MedianI_Nuc_SYP-Bkg']=results['MedianI_Nuc_SYP']-results['MedianI_SYP_Background']
results['MaxI_TOTAL_SYP-Bkg']=results['MaxI_SC_SYP-Bkg']+results['MaxI_Nuc_SYP-Bkg']
results['MinI_TOTAL_SYP-Bkg']=results['MinI_SC_SYP-Bkg']+results['MinI_Nuc_SYP-Bkg']
results['MeanI_TOTAL_SYP-Bkg']=results['MeanI_SC_SYP-Bkg']+results['MeanI_Nuc_SYP-Bkg']
results['MedianI_TOTAL_SYP-Bkg']=results['MedianI_SC_SYP-Bkg']+results['MedianI_Nuc_SYP-Bkg']
results['I_Total_SYP-Bkg']=results['I_Total_SYP']-(results['MeanI_SYP_Background']*results['#CP_Pixels'])
results['I_Total_SYPSC-Bkg']=results['I_Total_SYPSC']-(results['MeanI_SYP_Background']*results['#SC_Pixels'])
results['I_Total_SYPNuc-Bkg']=results['I_Total_SYPNuc']-(results['MeanI_SYP_Background']*results['#Nuc_Pixels'])
results['Mean_Ratio_SCvsNuc']=results['MeanI_SC_SYP-Bkg']/results['MeanI_Nuc_SYP-Bkg'] ##Concentration of protein on axis versus nucleoplasm
results['Mean_Ratio_SCvsTotal']=results['MeanI_SC_SYP-Bkg']/results['MeanI_TOTAL_SYP-Bkg'] ##Concentration of protein on axis versus total nucleus
results['Sum_Ratio_SCvsNuc']=results['I_Total_SYPSC-Bkg']/results['I_Total_SYPNuc-Bkg'] ##Fraction of protein on axis versus nucleoplasm
results['Sum_Ratio_SCvsTotal']=results['I_Total_SYPSC-Bkg']/results['I_Total_SYP-Bkg'] ##Fraction of protein on axis versus total nucleus
results['MaxI_DAPI_Bkg']=MaxDAPIBkg
results['MinI_DAPI_Bkg']=MinDAPIBkg
results['MeanI_DAPI_Bkg']=MeanDAPIBkg
results['MedianI_DAPI_Bkg']=MedianDAPIBkg
results['MaxI_TOTAL_DAPI-Bkg']=results['MaxI_TOTAL_DAPI']-results['MaxI_DAPI_Bkg']
results['MinI_TOTAL_DAPI-Bkg']=results['MinI_TOTAL_DAPI']-results['MinI_DAPI_Bkg']
results['MeanI_TOTAL_DAPI-Bkg']=results['MeanI_TOTAL_DAPI']-results['MeanI_DAPI_Bkg']
results['MedianI_TOTAL_DAPI-Bkg']=results['MedianI_TOTAL_DAPI']-results['MedianI_DAPI_Bkg']
results['I_Total_DAPI-Bkg']=results['I_Total_DAPI']-(results['MeanI_DAPI_Bkg']*results['#CP_Pixels'])
results['I_Total_DAPISC-Bkg']=results['I_Total_DAPISC']-(results['MeanI_DAPI_Bkg']*results['#SC_Pixels'])


##Re-calculate centroid coordinates in stitched gonad to remove duplicated nuclei

print(stitched)
if stitched=="True":
    if len(Stitched_MIP_folder)>1:
        folder=Stitched_MIP_folder
    else:
        folder=VSI_Directory + "/MaxTiffs/"
        
    print('Found ', folder)

    config_file=glob.glob(folder+"*TileConfiguration.registered.txt")
    Registration=open(config_file[0], "r")
    Registration_lines=Registration.readlines()
    for line in Registration_lines:
        print('Registration_line:',line)
        if name in line:
            start_x=line.find(name + ".tif") + len(name+".tif") + 5
            print('start_x:',start_x)
            x_image=float(line[start_x:(line[start_x:].find(',')+start_x)])
            print('x_image:',x_image)
            start_y=line[start_x:].find(',')+start_x + 1
            print('start_y:',start_y)
            y_image=float(line[start_y:(line[start_y:].find(')')+start_y)])
            print('y_image:',y_image)

    results.insert(4,'TileYCorr',y_image)
    results.insert(5,'TileXCorr',x_image)
    results.insert(12,'CorrCy',results['Cy']+y_image)
    results.insert(13,'CorrCx',results['Cx']+x_image)
    results.insert(14,'PixelSizeZ', pixelSizeZ)
    results.insert(15,'PixelSizeY', pixelSizeY)
    results.insert(16,'PixelSizeX', pixelSizeX)

##Save resuls table as csv in tif image folder
results.to_csv(path + "/" + name +'_ouput.csv') 

print('Finished calculating object properties for ' + filename)

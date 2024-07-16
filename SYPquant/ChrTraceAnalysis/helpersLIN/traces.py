# Ignore warnings in notebook
import warnings
warnings.filterwarnings('ignore')

import sys
import pandas as pd 
import numpy as np
import math
import scipy
import scipy.misc
import re
import tifffile as tifff

from skimage import measure, data, filters
from skimage.filters import *
from skimage.measure import label, regionprops, regionprops_table 
import skimage as skimage
import cv2
import operator
from scipy import spatial
import math
from datetime import date as dateFun

from scipy.interpolate import splprep, splev
from scipy.spatial.distance import cdist

import re, glob, csv, ast, os

#My functions

def sphericity_IC(masks,label,marchcubesspacing=(1,1,1),level=None,gradient_direction='ascent'):
    #heavily based on: https://github.com/scikit-image/scikit-image/issues/3797
    # marchcubesspacing use voxel size in correct order of the axis in masks
    tempM=np.copy(masks)
    tempM[masks!=label] = 0
    vts, fs, ns, cs =  skimage.measure.marching_cubes_lewiner(tempM, level=level,gradient_direction=gradient_direction,spacing=marchcubesspacing)
    ##STILL TO ADD  allow_degenerate=False, and option to use expanded mask as mask option ()
    ##maybe   skimage.measure.find_contours
    volume=np.sum(masks==label)*math.prod(marchcubesspacing)
    surface_a=skimage.measure.mesh_surface_area(vts, fs)
    sphericity=(((math.pi)**(1/3))*((6*volume)**(2/3)))/(surface_a)
    return sphericity, surface_a, volume


def get_SubfolderFiles(mainfolder,subfoler_pattern,file_ending=""):
    subfolder_files=glob.glob(mainfolder+subfoler_pattern+"*/*"+file_ending)
    return subfolder_files

def get_INTsubfolder_order(subfolder_files,numberingPattern="Position_[0-9]+"):
    positions=set(re.findall(numberingPattern, ''.join(subfolder_files)))
    nums=re.findall("[0-9]+", ''.join(positions))
    numsI=np.array(list(map(int, nums)))
    out=numsI[np.argsort(numsI)]
    return out

def is_in(elem,vector):
    out=[elem==vector_elem for vector_elem in vector]
    occurance=sum(out)
    return occurance>0

def is_in_fulloutput(elem,vector):
    out=[elem==vector_elem for vector_elem in vector]
    occurance=sum(out)
    indices= [i for i, x in enumerate(out) if x]
    return out, indices, occurance, occurance>0
    
def adjust_shape_ZYX(obj,newshape):
    out=np.zeros(newshape,dtype='uint16')
    current=obj.shape
    out[0:current[0],0:current[1],0:current[2]]=obj
    return out

def extract_array_field_allROIs(array_field,roiset,keepROIOrder=False):
    initiate_loc=True
    for roi in roiset:
        if keepROIOrder:
            mask_value=int(roi["Order"][0][-1])
        else:
            mask_value=int(roi["Name"][0][-1])
        if initiate_loc:
            extracted_field=np.c_[ roi[array_field], np.ones(np.shape(roi[array_field])[0])*mask_value, list(range(np.shape(roi[array_field])[0])) ]
            initiate_loc=False
        else:
            temp=np.c_[ roi[array_field], np.ones(np.shape(roi[array_field])[0])*mask_value, list(range(np.shape(roi[array_field])[0])) ]
            extracted_field=np.r_[extracted_field,temp]

    return extracted_field

def combine_field_from_allROIs(roiset,internalpath='["Spline"]',fullname=False,columns = ['Z','Y','X']):
    
    initiate_loc=True
    for roi in roiset:
        if fullname:
            mask_value=roi["Name"][0]
        else:
            mask_value=int(roi["Name"][0][-1])
                
        if initiate_loc:
            extracted_field=eval("roi"+internalpath)
            
            if columns is None:
                out=pd.DataFrame(extracted_field)
                # trace_name=pd.DataFrame([mask_value] * len(out.index), columns="trace_name")
                out["trace_name"]=mask_value
                # out = pd.concat([out, trace_name.reindex(out.index)], axis=1)

            else:
                out=pd.DataFrame(extracted_field, columns = columns)
                # trace_name=pd.DataFrame([mask_value] * len(out.index), columns="trace_name")
                out["trace_name"]=mask_value
                # out = pd.concat([out, trace_name.reindex(out.index)], axis=1)
            initiate_loc=False
        else:
            if columns is None:
                tmp=pd.DataFrame(extracted_field)
                # trace_name=pd.DataFrame([mask_value] * len(tmp.index), columns="trace_name")
                # tmp = pd.concat([tmp, trace_name.reindex(tmp.index)], axis=1)
                tmp["trace_name"]=mask_value
                out = pd.concat([out, tmp], ignore_index=True, sort=False)
                

            else:
                tmp=pd.DataFrame(extracted_field, columns = columns)
                # trace_name=pd.DataFrame([mask_value] * len(tmp.index), columns="trace_name")
                # tmp = pd.concat([tmp, trace_name.reindex(tmp.index)], axis=1)
                tmp["trace_name"]=mask_value
                out = pd.concat([out, tmp], ignore_index=True, sort=False)

    return out



def get_shift(i1,i2,axis=None,corr_mode="same",corr_method="auto"):
    corr=scipy.signal.correlate(i1,i2,mode=corr_mode,method=corr_method)
    out=abs(np.array(np.array(i1.shape).astype("double")*(1/2)).astype("int")-np.unravel_index(np.argmax(corr), np.shape(corr), order='C'))
    if axis is not None:
        out=abs(np.array(np.array(i1.shape).astype("double")*(1/2)).astype("int")-np.unravel_index(np.argmax(corr), np.shape(corr), order='C'))[axis]
    return out

def combine_pandasfield_from_allROIs(roiset,internalpath='["Spline"]',fullname=False,keepROIOrder=True):
    
    initiate=True
    for roi in roiset:
        if fullname:
            mask_value=roi["Name"][0]
        else:
            mask_value=int(roi["Name"][0][-1])

        
        if initiate:
                        
            out=pd.DataFrame.from_dict(eval("roi"+internalpath))
            # trace_name=pd.DataFrame([mask_value] * len(out.index), columns=["trace_name"])
            out["trace_name"]=mask_value
            # out = pd.concat([out, trace_name.reindex(out.index)], axis=1)
            initiate=False
        else:
            tmp=pd.DataFrame.from_dict(eval("roi"+internalpath))
            # trace_name=pd.DataFrame([mask_value] * len(tmp.index), columns=["trace_name"])
            # tmp = pd.concat([tmp, trace_name.reindex(tmp.index)], axis=1)
            tmp["trace_name"]=mask_value
            out = pd.concat([out, tmp], ignore_index=True, sort=False)

    return out
def distance_to_centre(centre,current_coord,zscale):
    out=abs(current_coord-centre)
    out=math.sqrt(out[0]*out[0]*zscale+out[1]*out[1]+out[2]*out[2])
    return out

def arclength_simple(coord_3D):
    out=0
    for i in range(np.shape(coord_3D)[0]-1):
        d=np.linalg.norm(coord_3D[i,:] - coord_3D[i+1,:])
        out +=d
    return out

def TNB_frame(R):
    #TNB (Frenet–Serret) frame implementation as in MATLAB
    #https://tutorial.math.lamar.edu/Classes/CalcIII/VectorFunctions.aspx
    #https://ch.mathworks.com/matlabcentral/fileexchange/80981-tnb
    # R is matrix of 3D coordinates of ordered points
    #add option diff or gradient
    
    dR=np.gradient(R,axis=0)
    ds=np.sqrt(np.sum(np.square(dR),axis=1),dtype="longdouble")
    unit_Tangent=dR/ds[:,None] #unit tangent vector
    dT=np.gradient(unit_Tangent,axis=0)
    dTds=dT/ds[:,None]
    kappa=np.sqrt(np.sum(np.square(dTds),axis=1),dtype="longdouble") 
    unit_Normal =dTds/kappa[:,None]
    unit_Binormal=np.cross(unit_Tangent,unit_Normal)
    dB=np.gradient(unit_Binormal,axis=0)
    dBdsneg=-dB/ds[:,None]
    tau=[np.dot(dBdsneg[index,:],unit_Normal[index,:]) for index in range(np.shape(unit_Normal)[0])]
    out={"T": unit_Tangent, "N": unit_Normal, "B": unit_Binormal, "curvature": kappa, "torsion": tau, "arclength": sum(ds)}
    return out
    
def interparc(nt,coords,linear=False, method="CubicSpline", axis=0, bc_type='not-a-knot',diff="diff",wrap=False,reltol=1e-9):
    # Defininig interparc function - translation from Matlab code by John D'Errico
    #https://ch.mathworks.com/matlabcentral/fileexchange/34874-interparc
    #exporting function handle, second derivative, and nt=3 checkpoints are still missing. Still to test are also ppval and np.polyval possible differences
    # also wrapping into a closed curve (where the first end is added to the end) is optional  as here csape is stil not optional

    #to fix: returns nt+1 points, check for identical spots in the output, increase precision
    #LINEAR=TRUE - linear inerpolation only returns the first derivative as in original code
    #method options:CubicSpline, PchipInterpolator, Akima1DInterpolator. I might add CubicHermiteSpline in future.
       #For CubicSpline bc_type argument can be addded: 'not-a-knot', 'periodic, 'natural', 'clamped'
    ##details: https://docs.scipy.org/doc/scipy/reference/interpolate.html
    import decimal

    t=np.linspace(0,1,nt,dtype="longdouble")
    to=np.copy(t)
    n=np.shape(coords)[0]
    ndim=np.shape(coords)[1]
    coords=coords.astype("longdouble")
    p1 = coords[0,:]
    pend = coords[-1,:]
    
    #Wrap the vector function
    if np.linalg.norm(p1-pend)>10*np.spacing(np.linalg.norm(np.amax(np.absolute(coords),axis=0))) and wrap:
        coordsW=np.r_[coords, p1[None,:]]
        n=n+1
        #nt=nt+1 #once csape is available uncomment and comment above one 
    else:
        coordsW=np.copy(coords)

    dR=np.diff(coordsW,axis=0)
    ds=np.sqrt(np.sum(np.square(dR,dtype="longdouble"),axis=1),dtype="longdouble")
    chordlen=ds/sum(ds)
    cumarc=np.cumsum(chordlen,dtype="longdouble")
    cumarc=np.insert(cumarc,0,0).astype("longdouble")
    if linear:
        tbins=np.digitize(t,cumarc)
        tbins[(tbins<=0) | (t<=0)]=1
        tbins[(tbins >= n) | (t >= 1)]=n
        tbins=tbins-1 #-1 for 0-based

        s=(t-cumarc[tbins])/chordlen[tbins-1]
        pt = coordsW[tbins,:] + (coordsW[tbins+1,:] - coordsW[tbins,:])*np.tile(s[:,None],(1,ndim))
        #derive
        dudt = (coordsW[tbins+1,:] - coordsW[tbins,:])/np.tile(chordlen[tbins,None],(1,ndim))
        #get a function to create more points - missing
        out={"pt": pt, "dudt": dudt}
    else:
        
        spl=[]
        spld=[]
        spld_i=[]
        spldd=[]
        spldd_i=[]
        splddd=[]
        splddd_i=[]
        
        # diffarray=np.array([[3, 0,0],[0,2,0],[0,0,1],[0,0,0]],dtype="longdouble")
        # ddiffarray=np.array([[2, 0],[0,1],[0,0]],dtype="longdouble")
        # dddiffarray=np.array([[1],[0]],dtype="longdouble")
        #just c spline, later I can ass bsplines and scipy
        for i in range(ndim):
            if method=="CubicSpline":
                pp=eval("scipy.interpolate."+method+"(cumarc,coordsW[:,i],axis=axis,bc_type=bc_type)")
            else:
                pp=eval("scipy.interpolate."+method+"(cumarc,coordsW[:,i],axis=axis)")
            spl.append(pp)
            #check if they all have 4 parameters
            # dd={"c": [], "order": []}
            # dd["c"].append(np.matmul(np.transpose(pp.c),diffarray))
            # dd["order"].append(3)
            
            spld_i.append(pp.derivative(nu=1))#append(scipy.interpolate.PPoly(np.transpose(spld[i]["c"][0].astype("longdouble")),cumarc))
            spld.append(np.transpose(spld_i[i].c))#.append(dd)
            spldd_i.append(pp.derivative(nu=2))#append(scipy.interpolate.PPoly(np.transpose(spldd[i].astype("longdouble")),cumarc))
            spldd.append(np.transpose(spldd_i[i].c))#append(np.matmul(spld[i]["c"][0].astype("longdouble"),ddiffarray.astype("longdouble"),dtype="longdouble"))
            splddd_i.append(pp.derivative(nu=3))#append(scipy.interpolate.PPoly(np.transpose(splddd[i].astype("longdouble")),cumarc))
            splddd.append(np.transpose(splddd_i[i].c))#.append(np.matmul(spldd[i].astype("longdouble"),dddiffarray.astype("longdouble"),dtype="longdouble"))

        # % catch the case where there were exactly three points
        # % in the curve, and spline was used to generate the
        # % interpolant. In this case, spline creates a curve with
        # % only one piece, not two
        if len(cumarc) == 3:
            cumarc = spl[1].x
            n = len(cumarc)
            chordlen = sum(chordlen)


        # Generate the total arclength along the curve
        # by integrating each segment and summing the
        # results. The integration scheme does its job
        # using an ode solver.
        # polyarray here contains the derivative polynomials
        # for each spline in a given segment
        polyarray = np.zeros((ndim,3),dtype="longdouble")
        seglen = np.zeros((n-1,1),dtype="longdouble")
        def segkernel(t,y,ndim,polyarray):
            #sqrt((dx/dt)^2 + (dy/dt)^2 + ...)
            val = np.zeros(np.shape(t))
            for k in range(ndim):
                val = val + np.polyval(polyarray[k,:],t)**2

            val = np.sqrt(val,dtype="longdouble")
            return [val]

        for i in range(n-2):############-1 for number of segments and -1 for  0-based
            for j in range(ndim):
                polyarray[j,:]=spld[j][i,:].astype("longdouble")

        #  integrate the arclength for the i'th segment
        #  using ode45 for the integral. I could have
        #  done this part with quad too, but then it
        #  would not have been perfectly (numerically)
        #  consistent with the next operation in this tool.
            vdp1 = lambda t,y: segkernel(t,y,ndim,polyarray)
            sol = scipy.integrate.solve_ivp(vdp1, [0,chordlen[i].astype("longdouble")], [0],rtol=reltol,atol=1e-10,method="BDF")

            tout = sol.t
            yout = sol.y
            seglen[i] = yout[0][-1]

        #and normalize the segments to have unit total length
        totalsplinelength = sum(seglen.astype("longdouble"))
        cumseglen=np.insert(np.cumsum(seglen,dtype="longdouble"),0,0)

        #which interval did each point fall into, in
        #terms of t, but relative to the cumulative
        #arc lengths along the parametric spline?
        tbins=np.digitize(to*totalsplinelength,cumseglen)
        #catch any problems at the ends
        tbins[(tbins<=0) | (to<=0)]=1
        tbins[(tbins >= n) | (to >= 1)]=n
        tbins=tbins-1 #-1 for 0-based

        #Do the fractional integration within each segment
        #for the interpolated points. t is the parameter
        #used to define the splines. It is defined in terms
        #of a linear chordal arclength. This works nicely when
        #a linear piecewise interpolant was used. However,
        #what is asked for is an arclength interpolation
        #in terms of arclength of the spline itself. Call s
        #the arclength traveled along the spline.
        s = totalsplinelength*to

        def ode_events(t,y):
            #ode event trap, looking for zero crossings of y.
            terminal = np.ones(np.shape(y), dtype=bool)
            direction= np.ones(np.shape(y))
            return y[0]#
        # ode_events.terminal=True
        # ode_events.direction=1

        ti = np.copy(to)
        for i in range(nt-1):
            # si is the piece of arc length that we will look
            # for in this spline segment.
            si = s[i].astype("longdouble") - cumseglen[tbins[i]].astype("longdouble")
        
            # extract polynomials for the derivatives
            # in the interval the point lies in
            for z in range(ndim):
                polyarray[z,:] = spld[z][tbins[i],:].astype("longdouble")

            # we need to integrate in t, until the integral
            # crosses the specified value of si. Because we
            # have defined totalsplinelength, the lengths will
            # be normalized at this point to a unit length.
            #
            # Start the ode solver at -si, so we will just
            # look for an event where y crosses zero.
            #[tout,yout,te,ye] 
            vdp2 = lambda t,y: segkernel(t,y,ndim,polyarray)
            
            sol_2 = scipy.integrate.solve_ivp(vdp2, [0,chordlen[tbins[i]].astype("longdouble")], [-si],events=[ode_events],method="BDF",rtol=reltol,atol=1e-10)#
            #sol_2= ode45(@(t,y) segkernel(t,y),[0,chordlen(tbins(i))],-si,opts)
            tout = sol_2.t
            yout = sol_2.y
            te=sol_2.t_events
            ye=sol_2.y_events
            
            # we only need that point where a zero crossing occurred
            # if no crossing was found, then we can look at each end.
            if not ((te is None) or np.size(te, axis=None)==0):
                ti[i] = te[0].astype("longdouble") + cumarc[tbins[i]].astype("longdouble")
            else:
                # a crossing must have happened at the very
                # beginning or the end, and the ode solver
                # missed it, not trapping that event.
                if abs(yout[0][0]) < abs(yout[0][-1]):
                # the event must have been at the start.
                    ti[i] = tout[0].astype("longdouble") + cumarc[tbins[i]].astype("longdouble")
                else:
                # the event must have been at the end.
                    ti[i] = tout[-1].astype("longdouble") + cumarc[tbins[i]].astype("longdouble")

        # Interpolate the parametric splines at ti to get
        # our interpolated value.
        pt=np.empty((nt,ndim),dtype="longdouble")
        dudt = np.zeros((nt,ndim),dtype="longdouble")
        dduddt=np.zeros((nt,ndim),dtype="longdouble")
        d3ud3t=np.zeros((nt,ndim),dtype="longdouble")
        coefs=[]

        for ll in range(ndim):
            pt[:,ll] = spl[ll](ti).astype("longdouble")
            dudt[:,ll] = spld_i[ll](ti).astype("longdouble")
            dduddt[:,ll] = spldd_i[ll](ti).astype("longdouble")
            d3ud3t[:,ll] = splddd_i[ll](ti).astype("longdouble")
            coefs.append(spl[ll].c)
        #create a function handle for evaluation, passing in the splines - missing
        indexList = [np.any(i) for i in np.isnan(pt)]
        # delete all the rows with any NaN value usually last point if points were generated with Akima1D interpoplator
        pt = np.delete(pt, indexList, axis=0)
        dudt = np.delete(dudt, indexList, axis=0)
        dduddt = np.delete(dduddt, indexList, axis=0)
        d3ud3t = np.delete( d3ud3t, indexList, axis=0)
        ti=ti[[not elem for elem in indexList]]

        kappa=np.linalg.norm(np.cross(dudt,dduddt, axis=1),axis=1)/(np.linalg.norm(dudt,axis=1)**3).astype("longdouble")
        scndXthrd=np.cross(dduddt,d3ud3t, axis=1).astype("longdouble")
        tau=[np.dot(dudt[index,:],scndXthrd[index,:]).astype("longdouble") for index in range(np.shape(scndXthrd)[0])]/(np.linalg.norm(np.cross(dudt,dduddt, axis=1),axis=1)**2).astype("longdouble")
        
        def foft(t,spl):
            return Warning("function handle export is not working yet")
        out={"pt": pt, "dudt": dudt, "d2ud2t":dduddt,"d3ud3t": d3ud3t, "fofthandle": [foft(t,spl)], "curvature": kappa, "torsion":tau, "t": ti, "arclength": totalsplinelength,"coefs": coefs}
    return out


def profile_along_3D_polyline(input_image,polyline, polyline_transform=None,expanded_mask=None,cval=None, expansion_radius=1,I_profiles_function="max",I_profiles_function_option="axis=0",weight_dist_power=2,sample_more_points=False,factor=2,return_point_index=False):
    if sample_more_points:
        polyline=interparc(factor*np.shape(polyline)[0],polyline)["pt"]
    #polyline cols: coordinates, rows: points
    #polyline_transform needed if polyline coordinates do not fit pixels of the image (e.g is polyline is in microns or the image has been scaled in some way)
    if polyline_transform is not None:
        polyline=np.matmul(polyline,polyline_transform).astype('int')
    if (I_profiles_function_option=="") or (I_profiles_function_option is None):
            I_profiles_function_option="axis=0"
    if cval is None:
        cval=1

    polyline_mask=np.zeros(np.shape(input_image))
    polyline_mask[(tuple(polyline.T))]=cval

    if expanded_mask is None:
        expanded_mask=skimage.segmentation.expand_labels(polyline_mask, distance=expansion_radius)

    sampled_pixels=np.array(np.where(expanded_mask==cval)).T
    sampled_pixels_intensities=input_image[tuple(sampled_pixels.T)]
    sampled_pixels_distances=cdist(sampled_pixels,polyline)
    sampled_pixels_order=np.argmin(sampled_pixels_distances,axis=1)
    sampled_pixels_distances=np.amin(sampled_pixels_distances,axis=1)
    sampled_pixels_distances[sampled_pixels_distances==0]=1e-8
    sampled_pixels_weigths=1/np.power(sampled_pixels_distances,weight_dist_power)
    out=np.zeros([len(np.unique(sampled_pixels_order)),1],dtype="float")
    for loc_idx, polyline_point in enumerate(np.unique(sampled_pixels_order)):
        if I_profiles_function=="w_average":
            out[loc_idx,...]=eval("np.average"+"(sampled_pixels_intensities[np.where(sampled_pixels_order==polyline_point)], weights=sampled_pixels_weigths[np.where(sampled_pixels_order==polyline_point)],"+re.sub('axis=[0-9]', 'axis=0', I_profiles_function_option)+")")
        else:
            out[loc_idx,...]=eval("np."+I_profiles_function+"(sampled_pixels_intensities[np.where(sampled_pixels_order==polyline_point)],"+re.sub('axis=[0-9]', 'axis=0', I_profiles_function_option)+")")
    if return_point_index:
        return out.T, np.unique(sampled_pixels_order)
    else:
        return out.T

def find_minimal_distances(point_set1,point_set2,return_pairing=False):
    distances=cdist(point_set1,point_set2)
    order=np.argmin(distances,axis=1)
    out=np.amin(distances,axis=1)
    # out[out==0]=1e-8
    if return_pairing:
        return out, order
    else:
        return out

def straightener_3D(points, polyline_transform, interparc_output,Voxel,selection_radius,interparc_spacing, axis=[0,0,1],use_ssc=False,overhang_factor=0,return_polyline=False):
    polyline=np.matmul(interparc_output["pt"],polyline_transform).astype('int')
    # nt=polyline.shape[0]
    # ndim=polyline.shape[1]
    t=np.copy(interparc_output["t"])
    cumarc=t*interparc_output["arclength"]/Voxel[0]
    selR = max(selection_radius*1.5,1.5*interparc_spacing/Voxel[0])#preselection to speed up a bit
    newpoints=np.copy(points)
    finin=np.zeros([np.shape(points)[0],1],dtype=bool)
    # curr_coefs=np.copy(interparc_output['coefs'])
    dist2O=np.ones((np.shape(points)[0],1))#*selR
    for p in np.array(range(1,polyline.shape[0],2)):
        # if p==0:
        #     continue
        Or=polyline[p,:]
        #np.array([scipy.interpolate.PPoly(c, t, extrapolate=None, axis=0)(input[1]) for c in curr_coefs])/Voxel[0]
        start=polyline[p-1,:]
        if p==polyline.shape[0]-1:
            stop=Or
        elif p==polyline.shape[0]-2:
            stop=polyline[-1,:]
        else:
            stop=polyline[p+1,:]
          
        Tan=stop-start
        n=interparc_output['d2ud2t'][p]/Voxel[0]
        if p==polyline.shape[0]-1:
            # Tan=interparc_output["dudt"][p]/Voxel[0]
            n=interparc_output['d2ud2t'][p-1]/Voxel[0]#-1
        r=Tan/np.linalg.norm(Tan)
        n=n/np.linalg.norm(n)
        bn=np.cross(r,n)
        #n=np.cross(bn, n)/np.linalg.norm(np.cross(bn, n))
        br=np.cross(bn, r)/np.linalg.norm(np.cross(bn, r))
        overhang=overhang_factor*interparc_spacing/Voxel[0]
        startD=min(np.dot(start-Or,r),np.dot(stop-Or,r))-overhang #will be negative, 30 is margin to avoid loosing points
        stopD=max(np.dot(start-Or,r),np.dot(stop-Or,r))+overhang #will be positive
        if p==1:
            if np.dot(start-Or,r)<0:
                startD=-(2*selection_radius+interparc_spacing/Voxel[0])
            else:
                stopD=(2*selection_radius+interparc_spacing/Voxel[0])
            selR=selR+2*interparc_spacing/Voxel[0]
        elif (p>=polyline.shape[0]-2):
            if np.dot(stop-Or,r)>0: ##### CHECK!!!!
                stopD=(polyline.shape[0]-p)*selection_radius+interparc_spacing/Voxel[0]
            else:
                startD=-(polyline.shape[0]-p)*selection_radius+interparc_spacing/Voxel[0]
            selR=selR+(polyline.shape[0]-p)*interparc_spacing/Voxel[0]#abs(polyline.shape[0]-p)
        
        indexH=insphere(points,Or,selR*1.5)[0] #preselection
        index=np.zeros((np.shape(points)[0],1),dtype=bool)
        for h in range(len(indexH)):
            if indexH[h]:
                D=np.dot(points[h,:]-Or,r)
                Do=abs(np.linalg.norm(points[h,:]-Or))
                #Dt=np.dot(points[h,:]-Or,br)
                if D>=startD and D<=stopD:# and (np.linalg.norm(points[h,:]-Or)<=selR or np.linalg.norm(points[h,:]-start)<=selR or np.linalg.norm(points[h,:]-start)<=selR): #:# and abs(Dt)<=selR 
                    if finin[h]:
                        if dist2O[h]>Do:
                            index[h]=True
                            dist2O[h]=Do
                        else:
                            index[h]=False
                    else:
                        index[h]=True
                        dist2O[h]=Do
        angle= -math.acos(np.dot(r,axis))#- means rotate clockwise!! %atan2d(norm(cross(r,axis)),dot(r,axis)) %LINK missing to matlab and wiki page explaining this
        u=np.cross(r,axis)
        
        if abs(np.linalg.norm(u))==0:
            u=u
            rotMat=np.eye(3)
            if np.dot(axis,r)<0:
                rotMat[np.where(np.array(axis)==1)[0],np.where(np.array(axis)==1)[0]]=-1
        #elif dot(r,axis)==0:
        else:
            u=u/np.linalg.norm(u)
            if use_ssc:
                ssc = [[0, -u[2], u[1]], [u[2],0,-u[0]], [-u[1],u[0],0]] #040423 https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
                rotMat=np.eye(3) + ssc + np.matmul(ssc,ssc)*(1+np.dot(axis,r))/(np.linalg.norm(u)**2) / 1
            else:
                rotMat=np.eye(3)*math.cos(angle)+math.sin(angle)*np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])+(1-math.cos(angle))*np.array([[u[0]**2, u[0]*u[1], u[0]*u[2]], [u[0]*u[1], u[1]**2, u[1]*u[2]], [u[0]*u[2], u[1]*u[2], u[2]**2] ])
        dist=cumarc[p]
        # print(u)
                
        
        newpoints[np.where(index)[0],:]=np.transpose(np.matmul(rotMat,np.transpose(points[np.where(index)[0],:]-Or)))
        newpoints[np.where(index)[0],np.where(np.array(axis)==1)[0]]=newpoints[np.where(index)[0],np.where(np.array(axis)==1)[0]]+dist
        finin[np.where(index)[0]]=True
    if return_polyline:
        return newpoints, finin, polyline
    else:
        return newpoints, finin

def TNB_frame_splines(pt, makeEq=False, method="CubicSpline", axis=0, bc_type='not-a-knot',diff="diff"):
    #do not change axis and diff parameter yet
    #pt is 3d coordinate matrix with ZYX coordinates where each row is a point
    ndim=np.shape(pt)[1]
    n=np.shape(pt)[0]
    t=np.linspace(0,1,n)
    dR=np.diff(pt,axis=0)
    ds=np.sqrt(np.sum(np.square(dR,dtype="longdouble"),axis=1),dtype="longdouble")
    chordlen=ds/sum(ds)
    cumarc=np.cumsum(chordlen,dtype="longdouble")
    cumarc=np.insert(cumarc,0,0).astype("longdouble")
    nt=len(cumarc)
    if makeEq:
        cumarc=np.copy(t)

    spl=[]
    spld_i=[]
    spldd_i=[]
    splddd_i=[]

    dudt = np.zeros((nt,ndim),dtype="longdouble")
    dduddt=np.zeros((nt,ndim),dtype="longdouble")
    d3ud3t=np.zeros((nt,ndim),dtype="longdouble")
    
    for i in range(ndim):
        if method=="CubicSpline":
            pp=eval("scipy.interpolate."+method+"(cumarc,pt[:,i],axis=axis,bc_type=bc_type)")
        else:
            pp=eval("scipy.interpolate."+method+"(cumarc,pt[:,i],axis=axis)")
        
        spl.append(pp)
        #check if they all have 4 parameters
        spld_i.append(pp.derivative(nu=1))
        spldd_i.append(pp.derivative(nu=2))
        splddd_i.append(pp.derivative(nu=3))

        dudt[:,i] = spld_i[i](cumarc).astype("longdouble")
        dduddt[:,i] = spldd_i[i](cumarc).astype("longdouble")
        d3ud3t[:,i] = splddd_i[i](cumarc).astype("longdouble")
            #create a function handle for evaluation, passing in the splines - missing
            
    kappa=np.linalg.norm(np.cross(dudt,dduddt, axis=1),axis=1)/(np.linalg.norm(dudt,axis=1)**3).astype("longdouble")
    scndXthrd=np.cross(dduddt,d3ud3t, axis=1).astype("longdouble")
    tau=[np.dot(dudt[index,:],scndXthrd[index,:]).astype("longdouble") for index in range(np.shape(scndXthrd)[0])]/(np.linalg.norm(np.cross(dudt,dduddt, axis=1),axis=1)**2).astype("longdouble")
    ds=np.sqrt(np.sum(np.square(dudt),axis=1))
    unit_Tangent=dudt/ds[:,None] #unit tangent vector
    dT=np.gradient(unit_Tangent,axis=0)
    dTds=dT/ds[:,None]
    kappa_TNB=np.sqrt(np.sum(np.square(dTds),axis=1)) 
    unit_Normal =dTds/kappa[:,None]
    unit_Binormal=np.cross(unit_Tangent,unit_Normal)
    dB=np.gradient(unit_Binormal,axis=0)
    dBdsneg=-dB/ds[:,None]
    tau_TNB=[np.dot(dBdsneg[index,:],unit_Normal[index,:]) for index in range(np.shape(unit_Normal)[0])]
    out={"T": unit_Tangent, "N": unit_Normal, "B": unit_Binormal, "curvature": kappa, "torsion": tau, "arclength": sum(ds), "kappa_TNB":kappa_TNB, "tau_TNB": tau_TNB}
    return out

#Oane's functions edited (_IC)
def get_rois_IC(fname):
    # print(fname)
    try:
        ifs = open(fname, mode='r',encoding="utf-8")
        ifs.readlines()
        ifs.close()
        ifs = open(fname, mode='r',encoding="utf-8")
    except:
        ifs = open(fname, mode='r',encoding="cp1252")

    # print(ifs.readlines())
    rois = [] # list of all csv rows per ROI
    genericinfo = [] # should contain all data that is outside of specific ROIs, might not fully work with groups yet
    generic = True
    csvreader = csv.reader(ifs)
    for row in csvreader:
        # print(row)

        if 'BT_Roi' in row:
            rois.append({'Order': row[-1],'Type': [],'Name': [],'GroupInd': [], 'csvlist': []})
            generic = False  
        if generic: 
            genericinfo.append(row)
        else:
            if row[0] in ['Type','Name','GroupInd']: 
                rois[-1][row[0]].append(row[-1].lower())##########################
            else:
                rois[-1]['csvlist'].append(row)

    ifs.close()
    for roi in rois:

        if roi["Type"][0]=="linetrace":
            index=0
            while index<len(roi["csvlist"]):
                if roi["csvlist"][index][0] == "Vertices":
                    if roi["csvlist"][index][0] not in roi:
                        roi[roi["csvlist"][index][0]] = roi["csvlist"][index][1]
                    try:
                        roi["User_Points"]=np.array([list(map(float,x)) for x in roi["csvlist"][index+1:index+1+int(roi["csvlist"][index][1])]])
                    except:
                        temp=np.array(roi["csvlist"][index+1:index+1+int(roi["csvlist"][index][1])])
                        temp[temp==""]='0'
                        roi["User_Points"]=temp.astype("float")
                    index = index + int(roi["csvlist"][index][1])+1
                    
                elif roi["csvlist"][index][0] == "Segment":
                    try:
                        points=np.array([list(map(float,x)) for x in roi["csvlist"][index+1:index+1+int(roi["csvlist"][index][3])]]) # backward compatible
                    except:
                        temp=np.array(roi["csvlist"][index+1:index+1+int(roi["csvlist"][index][3])])
                        temp[temp==""]='0'
                        points=temp.astype("float")
                    points[:,3]=int(roi["csvlist"][index][1])
                    if "Points" not in roi:
                        roi["Points"] = points
                    else:
                        roi["Points"]=np.concatenate((roi["Points"], points))
                        
                    index =index + int(roi["csvlist"][index][3]) +1
                else: 
                    if roi["csvlist"][index][0] not in roi:
                        roi[roi["csvlist"][index][0]] = roi["csvlist"][index][-1]
                    index +=1 
            ###Remove redundant points
            selectRows=[tuple(roi["Points"][Row,0:2])!=tuple(roi["Points"][Row+1,0:2]) for Row in range(len(roi["Points"])-1)]
            selectRows.append(True)
            roi["Points"]=roi["Points"][selectRows,:]

        elif roi["Type"][0]=="polyline":
            for index, element in enumerate(roi["csvlist"]):
                if element[0] == "Vertices":
                    if element[0] not in roi:
                        roi[element[0]] = element[1]
                    roi["User_Points"]=np.array([list(map(float,x)) for x in roi["csvlist"][index+1:index+1+int(element[1])]])
                    ### Iinterpolate more points as -- not implemented yet, not needed
                    roi["Points"]=np.copy(roi["User_Points"])
                    break
                else: 
                    if element[0] not in roi:
                        roi[element[0]] = element[1]

        else:
            print("Only ROIs annotated as Polyline or LineTrace can be used.")
        
        roi.pop('csvlist')
    Voxel=[]
    for ind, field in enumerate(genericinfo):
        if field[0] in ["ImageVoxelWidth", "ImageVoxelHeight", "ImageVoxelDepth"]:
            Voxel.append(float(field[-1]))
    return rois, genericinfo, Voxel

def DoG_sharpening(imageI,initial_sigma=1,sigma=None,gaussian_options="",only_blur_inital=False):
    if sigma is None:
        sigma=2*initial_sigma
        print("Warning in traces.DoG_sharpening: sigma parameter is not defined and will be set to double of the initial_sigma.")
    if isinstance(sigma,list) or isinstance(sigma,np.ndarray):
        if len(sigma)!=len(np.shape(imageI)):
            print("Warning in traces.DoG_sharpening: shape of the sigma parameter does not correspond to the image shape. First value will be used")
    if isinstance(initial_sigma,list) or isinstance(initial_sigma,np.ndarray):
        if len(initial_sigma)!=len(np.shape(imageI)):
            print("Warning in traces.DoG_sharpening: shape of the sigma parameter does not correspond to the image shape. First value will be used")
    if np.any(np.array(sigma)<=np.array(initial_sigma)): ##Needs modifying in case you have list and not initigers
        sys.exit("Error in traces.DoG_sharpening: sigma parameter must be larger than the initial sigma in all corresponding dimensions.")
    
    #If there is an empty slice, sharpening will create weird artefacts
    z=np.where(np.sum(imageI,axis=(1,2))!=0)[0]
    x=np.where(np.sum(imageI,axis=(0,1))!=0)[0]
    y=np.where(np.sum(imageI,axis=(0,2))!=0)[0]
    image_c=np.copy(imageI[np.ix_(z,y,x)])

    if only_blur_inital:
        if gaussian_options=="":
            out =skimage.filters.gaussian(image_c, sigma=initial_sigma)
        else:
            out=eval("skimage.filters.gaussian(image_c, sigma=initial_sigma,"+gaussian_options+")") 

    else:
        if gaussian_options=="":
            blurred_sigma = skimage.filters.gaussian(image_c, sigma=initial_sigma)        
            blurred_res = skimage.filters.gaussian(blurred_sigma, sigma=sigma)#zyx_resolution_pxl
        else:
            blurred_sigma = eval("skimage.filters.gaussian(image_c, sigma=initial_sigma,"+gaussian_options+")")     
            blurred_res = eval("skimage.filters.gaussian(blurred_sigma, sigma=sigma,"+gaussian_options+")")    
        out = blurred_sigma - blurred_res
    del image_c
    # print(np.max(out))
    image_o=np.zeros_like(imageI).astype("float")
    image_o[np.ix_(z,y,x)]=out.astype("float")
    return image_o

def DoG_sharpening_std(imageI,initial_sigma=1,sigma=None,gaussian_options="",only_blur_inital=False):
    if sigma is None:
        sigma=2*initial_sigma
        print("Warning in traces.DoG_sharpening: sigma parameter is not defined and will be set to double of the initial_sigma.")
    if isinstance(sigma,list) or isinstance(sigma,np.ndarray):
        if len(sigma)!=len(np.shape(imageI)):
            print("Warning in traces.DoG_sharpening: shape of the sigma parameter does not correspond to the image shape. First value will be used")
    if isinstance(initial_sigma,list) or isinstance(initial_sigma,np.ndarray):
        if len(initial_sigma)!=len(np.shape(imageI)):
            print("Warning in traces.DoG_sharpening: shape of the sigma parameter does not correspond to the image shape. First value will be used")
    if np.any(np.array(sigma)<=np.array(initial_sigma)): ##Needs modifying in case you have list and not initigers
        sys.exit("Error in traces.DoG_sharpening: sigma parameter must be larger than the initial sigma in all corresponding dimensions.")
    #If there is an empty slice, sharpening will create weird artefacts
    z=np.where(np.sum(imageI,axis=(1,2))!=0)[0]
    x=np.where(np.sum(imageI,axis=(0,1))!=0)[0]
    y=np.where(np.sum(imageI,axis=(0,2))!=0)[0]
    image_c=np.copy(imageI[np.ix_(z,y,x)])

    if only_blur_inital:
        if gaussian_options=="":
            out =skimage.filters.gaussian(image_c, sigma=initial_sigma)
        else:
            out=eval("skimage.filters.gaussian(image_c, sigma=initial_sigma,"+gaussian_options+")") 
    else:
        if gaussian_options=="":
            blurred_sigma = skimage.filters.gaussian(image_c, sigma=initial_sigma)        
            blurred_res = skimage.filters.gaussian(image_c, sigma=sigma)#zyx_resolution_pxl
        else:
            blurred_sigma = eval("skimage.filters.gaussian(image_c, sigma=initial_sigma,"+gaussian_options+")")     
            blurred_res = eval("skimage.filters.gaussian(image_c, sigma=sigma,"+gaussian_options+")")    
        out = blurred_sigma - blurred_res
    del image_c
    
    image_o=np.zeros_like(imageI,   dtype="float")
    image_o[np.ix_(z,y,x)]=out

    # print(np.max(image_o))
    return image_o

def sfloat_to_uint(image,mask=None,bit_depth=16,ql=0,qu=1,final_multiplier=None,remove_outside=True):
    if mask is None:
        mask=np.ones(np.shape(image))
    out=np.copy(image)
    if remove_outside:
        out[mask==0]=np.min(out)
    out=(out-np.quantile(out[mask!=0],q=ql))/(np.quantile(out[mask!=0],q=qu)-np.quantile(out[mask!=0],q=ql))
    out[out<0]=0
    out[out>1]=1
    if final_multiplier is None:
        final_multiplier=2**int(bit_depth)
    out=out*final_multiplier
    return out

def threshold_mine(image,mask=None,id=0,threshold_mth="otsu",binary_output=True,preThresh=None,return_normalized=False,upper_quantile=1,sharpen_prior=False,sharpen_sigma=0.75,sigma_res=[1.875, 3.132, 3.132],stdDog=False):
    #Add threshold by value
    image_loc=np.copy(image)

    if sharpen_prior:
        if stdDog:
            image_loc=DoG_sharpening_std(image_loc,sharpen_sigma,sigma_res)
        else:
            image_loc=DoG_sharpening(image_loc,sharpen_sigma,sigma_res)#np.copy(DoG)#*image_loc
    if preThresh is not None:
        image_loc[image<preThresh]=0
        image_loc[image==0]=0

    
    if mask is not None:
        if id==0:
            image_loc[mask==id]=0
        else:
            image_loc[mask!=id]=0

    

    out=image_loc.astype("float")
    out[image==0]=0
    
    if threshold_mth=="simple":
        trHold=preThresh
        out[out<=trHold]=0
    else:
        cntrl=np.sum(out!=0)
        if cntrl==0:
            trHold=1
        else:
            trHold=eval("skimage.filters.threshold_"+threshold_mth+"(out[out!=0])")
        out[out<trHold]=0

    if binary_output:
        out[out>=trHold]=1
        out[image==0]=0
    elif return_normalized and cntrl!=0:
        out=out.astype("float")/np.quantile(out[out>0].astype("float"),upper_quantile)
        out[out>1]=1
    return out

def threshold_mine_v1(image,mask=None,id=0,threshold_mth="otsu",binary_output=True):
    if mask is not None:
        image[mask!=id]=0

    trHold=eval("skimage.filters.threshold_"+threshold_mth+"(image[image!=0])")
    out=np.copy(image)
    out[out<trHold]=0

    if binary_output:
        out[out>=trHold]=1
        
    return out

def get_shift_v2(i1,i2,axis=None,corr_mode="same",corr_method="auto"):
    corr=scipy.signal.correlate(i1,i2,mode=corr_mode,method=corr_method)
    out=np.array(np.array(i1.shape).astype("double")*(1/2)).astype("int")-np.unravel_index(np.argmax(corr), np.shape(corr), order='C')
    if axis is not None:
        out=out[axis]
    return out

def metadata_dict__spotMAXtiff(path):
    #Reads metadata from Tiffs that are saved for cellpose and spotmax and converts it to a dictionary
    with tifff.TiffFile(path) as tiff:
        imagej_metadata = tiff.imagej_metadata
    out=dict((key,val) for key,val in imagej_metadata.items() if key.lower()!='info')
    info_dict=dict((VKpair[0].strip(), VKpair[1].strip()) for VKpair in (re.split('=|:',line.replace(",",'')) for line in imagej_metadata["Info"].replace(" ",'').split(",") if '=' in line or ':' in line))
    for key,val in info_dict.items():
        if key in out.keys():
            out[key+'_info']=val
        else:
            out[key]=val
        if "array" in val:
            out[key]=np.array(re.split(",",re.findall("\[.*\]",re.split('=|:',re.findall(key+"[=|:]array\(\[.*\]\)",imagej_metadata["Info"].replace(" ",''))[0])[1])[0][1:][:-1])).astype("float")
    return out

def resolution_info_spotMAXtiff(path):
    meta_dct=metadata_dict__spotMAXtiff(path)
    try:
        out={"NA": float([val for key,val in meta_dct.items() if "LensNA" in key][0]),
            "emission": float([val for key,val in meta_dct.items() if "emissionwavelength" in key.lower()][0]),
            "voxel": [float(meta_dct["PhysicalSizeZ"]),float(meta_dct["PhysicalSizeY"]),float(meta_dct["PhysicalSizeX"])], "axis": 'ZYX',
            "unit":[val for key,val in meta_dct.items() if "unit" == key.lower()][0]}
    except:
        sys.exit("Error: Information about resolution can not be find in the metadata. Function: resolution_info_spotMAXtiff()")
    
    out["rayleigh_nm"]=1.22*out["emission"]/(2*out["NA"])
    out["abbe_nm"]=out["emission"]/(2*out["NA"])
    #Rayleigh (used): 1.22*lambda/(2*NA), Abbe:0.5*lambda/NA, Sparrow:0.47*lambda/NA float([val for key,val in meta_dct.items() if "emissionwavelength" in key.lower()][0])
    return out

def get_sharpening_input(path,Z_limit_um=None,resolution_multiplier=1,emission=None,output_res_info=False):
    image=tifff.imread(glob.glob(path))
    res_info=resolution_info_spotMAXtiff(path)

    if emission is not None:
        res_info["emission"]=emission
        res_info["rayleigh_nm"]=res_info["rayleigh_nm"]=1.22*emission/(2*res_info["NA"])
    if res_info["unit"].lower()!="um":
        print("Warning: formula in sharpening part needs to be adjusted. Currently it is assumed that the unit of tiff is um, but something else is save in the metadata. Everything will still run through")
    if Z_limit_um is None:
        res_ZYX=[1,0.001*resolution_multiplier*res_info["rayleigh_nm"]/res_info["voxel"][1],0.001*resolution_multiplier*res_info["rayleigh_nm"]/res_info["voxel"][2]]
    else:
        res_ZYX=[Z_limit_um/res_info["voxel"][0],0.001*resolution_multiplier*res_info["rayleigh_nm"]/res_info["voxel"][1],0.001*resolution_multiplier*res_info["rayleigh_nm"]/res_info["voxel"][2]]
    if output_res_info:
        return image, res_ZYX, res_info
    else:
        return image, res_ZYX
    
def insphere(data,position, dist):
    index=[np.square(data[:,0]-position[0])+np.square(data[:,1]-position[1])+np.square(data[:,2]-position[2])<=dist**2]
    return index
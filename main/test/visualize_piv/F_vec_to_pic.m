function [data]=F_vec_to_pic(X,loc_normal,x,y)
    datavec=NaN(y*x,1);
    datavec(loc_normal)=X;
    data=reshape(datavec,[y,x]);
end
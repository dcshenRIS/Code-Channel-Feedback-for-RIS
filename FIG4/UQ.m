function [h_UQ,err_quant] = UQ(h,B,hmin,hmax)%将角度的三角函数量化
%% descriptions of this function
% This function is to quan. and norm. the spacial angle with B bits to the
% range in [hmin,hmax]
% ---------------- input descriptions -------------------------------------
%   "h" is the spacial angles. 
%   "B" is the quan. bits.
% ---------------- output descriptions ------------------------------------
%   "h_UQ" is the quantized angles.
m=length(h);
h_UQ=zeros(size(h));
for i=1:m
if h(i)>hmax
    h_UQ(i)=hmax;
else if h(i)<hmin
        h_UQ(i)=hmin;
    end
delta=(hmax-hmin)/2^B;
h_UQ(i) = hmin + delta*round((h(i)-hmin)/delta);
end
end
err_quant=norm(h_UQ-h)^2/norm(h)^2;
end
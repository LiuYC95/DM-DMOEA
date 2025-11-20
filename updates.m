function [chrom,obj]=updates(weight,chrom,obj,newpoint,zz,z)
weight(weight==0)=0.00001; 
% 采取切比雪夫聚合方法，先计算个体的各个子目标函数值和参考点值的差值的绝对值，
% 然后分别乘以对应的权重，选取其中最大的值作为聚合后的目标函数值
part=abs(zz-z);  % 计算得到的新解在当前权重下的聚合后函数值
newobj = max(weight.*part); 
part=abs(obj-z); % 计算当前个体的聚合后函数值
oldobj=max(weight.*part); 
if newobj<oldobj % 若新解更优，替换
    chrom=newpoint;
    obj=zz;
end
end


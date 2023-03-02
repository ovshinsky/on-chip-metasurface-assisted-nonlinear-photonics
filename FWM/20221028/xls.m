function I = xls( path )
%EXCELREAD Summary of this function goes here
%   Detailed explanation goes here
            I=xlsread(path,2);
            I(:,1)=I(:,1)*1e9;
            for k=2:size(I,2)
                I(:,k)=-1*I(:,k);
            end
end


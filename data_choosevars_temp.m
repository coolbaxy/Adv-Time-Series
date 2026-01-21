

function [Y, X, vars, lev]=data_choosevars_temp(begin_date,end_date,esty1,estq1,esty2,estq2,nlags,choose_vars,samp,data,macrovars)

macro = data;   
macro_select=choose_vars;

%%
daty1=begin_date;       % First Year of Data Set
datq1=1;                % First Quarter of Data Set
n1 = (esty1-daty1)*4 + ( estq1 - datq1 + 1);          
n2 = (esty2-daty1)*4 + ( estq2 - datq1 + 1);        
data = macro(n1:n2,macro_select);
vars = macrovars(macro_select,:);   

%set prior for first lag of each variable (if estimation with Minnesota
%prior)
[rowy,coly] = size(data);
lev = ones(coly,1);            
%compute nlags lags
[T,nvars]=size(data);
for p=1:nlags;
    X(:,1+(p-1)*nvars:p*nvars)=data((nlags+1-p):(T-p),:);                                           
end;

%rescaling variables since we loose nlags observations through the lagging 
Y=data((nlags+1):T,:);
%%
% label for industries:
%1.Food And Kindred Products	 
 %2.Tobacco Products	 
 %3.Textile Mill Products	 
 %4.Apparel And Other Finished Products Made From Fabrics And Similar Materials	 
 %5.Lumber And Wood Products, Except Furniture	 
 %6.Furniture And Fixtures	 
 %7.Paper And Allied Products	 
 %8.Printing, Publishing, And Allied Industries	 
 %9.Chemicals And Allied Products	 
 %10.Petroleum Refining And Related Industries	 
 %11Rubber And Miscellaneous Plastics Products	 
 %12.Leather And Leather Products	 
 %13.Stone, Clay, Glass, And Concrete Products	
 %14.Primary Metal Industries	 
 %15.Fabricated Metal Products, Except Machinery And Transportation Equipment	
 %16.Industrial And Commercial Machinery And Computer Equipment	 
 %17.Electronic And Other Electrical Equipment And Components, Except Computer Equipment	 
 %18.Transportation Equipment	 
 %19.Measuring, Analyzing, And Controlling Instruments; Photographic, Medical And Optical Goods; Watches And Clocks	
 %20.Miscellaneous Manufacturing Industries

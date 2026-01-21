function c = correl(x,y) ;

c = corrcoef(x,y);
c = c(1,2);
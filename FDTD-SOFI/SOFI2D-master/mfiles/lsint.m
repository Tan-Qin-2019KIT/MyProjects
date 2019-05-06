function [xold,fold,para,how]=lsint(xnew,fnew,para)
%LSINT	Function to initialize LEASTSQ routine.

%	Copyright (c) 1990 by the MathWorks, Inc.
%	Andy Grace 7-9-90.
xold=xnew;
fold=fnew;
para=foptions(para);
if para(14)==0, para(14)=length(xnew)*100;end 
if para(1)>0,
	disp('')
	if para(5)>0
	disp('f-COUNT      RESID    STEP-SIZE      GRAD/SD  LINE-SEARCH')
else
	disp('f-COUNT      RESID    STEP-SIZE      GRAD/SD        LAMBDA LINE-SEARCH')
	end
end

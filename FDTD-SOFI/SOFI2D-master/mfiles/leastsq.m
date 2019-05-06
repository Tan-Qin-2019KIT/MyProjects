function [x,OPTIONS,f,JOCOB] = leastsq(FUN,x,OPTIONS,GRADFUN,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10)
%LEASTSQ Solves non-linear least squares problems.
% 	LEASTSQ solves problems of the form:
%	min  sum {FUN(X).^2}    where FUN and X may be vectors of matrices.   
%             x
%
%	X=LEASTSQ('FUN',X0) starts at the matrix X0 and finds a minimum to the
%	sum of squares of the functions described in FUN. FUN is usually
%	an M-file which returns a matrix of objective functions: F=FUN(X).
%
%	X=LEASTSQ('FUN',X0,OPTIONS) allows a vector of optional parameters to
%	be defined. OPTIONS(2) is a measure of the precision required for the 
%	values of X at the solution. OPTIONS(3) is a measure of the precision
%	required of the objective function at the solution. See HELP FOPTIONS. 
%
%	X=LEASTSQ('FUN',X0,OPTIONS,'GRADFUN') enables a function'GRADFUN'
%	to be entered which returns the partial derivatives of the functions,
%	dF/dX, (stored in columns) at the point X: gf = GRADFUN(X).
%	

%	[X,OPTIONS,F,J]=LEASTSQ('FUN',X0,...) returns, F, the residuals and
%	J the Jacobian of the function FUN at the solution. 
%	F and J can be used with to calculate the variance and
%	confidence intervals using CONFINT(X,S,J).

%	Copyright (c) 1990 by the MathWorks, Inc.
%	Andy Grace 7-9-90.

%	The default algorithm is the Levenberg-Marquardt method with a 
%	mixed quadratic and cubic line search procedure.  A Gauss-Newton
%	method is selected by setting  OPTIONS(5)=1. 
%
%	X=LEASTSQ('FUN',X,OPTIONS,GRADFUN,P1,P2,..) allows
%	coefficients, P1, P2, ... to be passed directly to FUN:
%	F=FUN(X,P1,P2,...). Empty arguments ([]) are ignored.

% ------------Initialization----------------
XOUT=x(:);
[nvars]=length(XOUT);

evalstr = [FUN];
if ~any(FUN<48)
	evalstr=[evalstr, '(x'];
	for i=1:nargin - 4
		evalstr = [evalstr,',P',int2str(i)];
	end
	evalstr = [evalstr, ')'];
end


if nargin < 3, OPTIONS=[]; end
if nargin < 4, GRADFUN=[]; end

if length(GRADFUN)
	evalstr2 = [GRADFUN];
	if ~any(GRADFUN<48)
		evalstr2 = [evalstr2, '(x'];
		for i=1:nargin - 4
		    evalstr2 = [evalstr2,',P',int2str(i)];
		end
		evalstr2 = [evalstr2, ')'];
	end
end

f = eval(evalstr);
f=f(:);
nfun=length(f);
GRAD=zeros(length(XOUT),nfun);
OLDX=XOUT;
MATX=zeros(3,1);
MATL=[f'*f;0;0];
OLDF=f'*f;
FIRSTF=f'*f;
[OLDX,OLDF,OPTIONS]=lsint(XOUT,f,OPTIONS);
EstSum=0.5;
GradFactor=0; 
CHG = 1e-7*abs(XOUT)+1e-7*ones(nvars,1);

OPTIONS(10)=1;
status=-1;

while status~=1

% Work Out Gradients
	if ~length(GRADFUN) | OPTIONS(9)
		OLDF=f;
		CHG = sign(CHG+eps).*min(max(abs(CHG),OPTIONS(16)),OPTIONS(17));
		for gcnt=1:nvars
		    temp = XOUT(gcnt);
		    XOUT(gcnt) = temp +CHG(gcnt);
		    x(:) = XOUT;
		    f(:) = eval(evalstr);
		    GRAD(gcnt,:)=(f-OLDF)'/(CHG(gcnt));
		    XOUT(gcnt) = temp;
		end
		f = OLDF;
		OPTIONS(10)=OPTIONS(10)+nvars;
% Gradient check
		if OPTIONS(9) == 1
		    GRADFD = GRAD;
		    x(:)=XOUT; GRAD = eval(evalstr2);
		    graderr(GRADFD, GRAD, evalstr2);
		    OPTIONS(9) = 0;
		end
	else
		x(:) = XOUT;
		OPTIONS(11)=OPTIONS(11)+1;
		GRAD = eval(evalstr2);
	end
	% Try to set difference to 1e-8 for next iteration
	if nfun==1
		CHG = nfun*1e-8./GRAD;
	else
		CHG = nfun*1e-8./sum(abs(GRAD)')';
	end

	gradf = 2*GRAD*f;
	fnew = f'*f;
%---------------Initialization of Search Direction------------------
	if status==-1
		if cond(GRAD)>1e8
		    SD=-(GRAD*GRAD'+(norm(GRAD)+1)*(eye(nvars,nvars)))\(GRAD*f);
		    if OPTIONS(5)==0, GradFactor=norm(GRAD)+1; end
		    how='COND';
		else
%		SD=GRAD'\(GRAD'*X-F)-X;
		    SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*f);
		end
		FIRSTF=fnew;
		OLDG=GRAD;
		GDOLD=gradf'*SD;
		OPTIONS(18)=1;
		if OPTIONS(1)>0
		    disp([sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),fnew,OPTIONS(18)),sprintf('%12.3g  ',GDOLD)]);
		end
		XOUT=XOUT+OPTIONS(18)*SD;
		if OPTIONS(5)==0
		    newf=GRAD'*SD+f;
		    GradFactor=newf'*newf;
		    SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*f); 
		end
		newf=GRAD'*SD+f;
		XOUT=XOUT+OPTIONS(18)*SD;
		EstSum=newf'*newf;
		status=0;
		if OPTIONS(7)==0; PCNT=1; end
		
	else
%-------------Direction Update------------------
		gdnew=gradf'*SD;
		if OPTIONS(1)>0, 
		    num=[sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),fnew,OPTIONS(18)),sprintf('%12.3g  ',gdnew)];
		end
		if gdnew>0 & fnew>FIRSTF

% Case 1: New function is bigger than last and gradient w.r.t. SD -ve
% ... interpolate. 
		    how='inter';
		    [stepsize]=cubici1(fnew,FIRSTF,gdnew,GDOLD,OPTIONS(18));
		    OPTIONS(18)=0.9*stepsize;
		elseif fnew<FIRSTF

%  New function less than old fun. and OK for updating 
%         .... update and calculate new direction. 
		    [newstep,fbest] =cubici3(fnew,FIRSTF,gdnew,GDOLD,OPTIONS(18));
		    if fbest>fnew,fbest=0.9*fnew; end 
		    if gdnew<0
		        how='incstep';
		        if newstep<OPTIONS(18),  newstep=(2*OPTIONS(18)+1e-4); how=[how,'IF']; end
		        OPTIONS(18)=abs(newstep);
		    else 
		        if OPTIONS(18)>0.9
		            how='int_step';
		            OPTIONS(18)=min([1,abs(newstep)]);
		        end
		    end
% SET DIRECTION.
% Gauss-Newton Method    
		    temp=1;
		    if OPTIONS(5)==1 
		        if OPTIONS(18)>1e-8 & cond(GRAD)<1e8
		            SD=GRAD'\(GRAD'*XOUT-f)-XOUT;
		            if SD'*gradf>eps,how='ERROR- GN not descent direction',  end
		            temp=0;
		        else
		            if OPTIONS(1) > 0
		                disp('Conditioning of Gradient Poor - Switching To LM method')
		            end
		            how='CHG2LM';
		            OPTIONS(5)=0;
		            OPTIONS(18)=abs(OPTIONS(18));				
		        end
		    end
		    
		    if (temp)      
% Levenberg_marquardt Method N.B. EstSum is the estimated sum of squares.
%                                 GradFactor is the value of lambda.
% Estimated Residual:
		        if EstSum>fbest
		            GradFactor=GradFactor/(1+OPTIONS(18));
		        else
		            GradFactor=GradFactor+(fbest-EstSum)/(OPTIONS(18)+eps);
		        end
		        SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*f); 
		        OPTIONS(18)=1; 
		        estf=GRAD'*SD+f;
		        EstSum=estf'*estf;
		        if OPTIONS(1)>0, num=[num,sprintf('%12.6g ',GradFactor)]; end
		    end
		    gdnew=gradf'*SD;

		    OLDX=XOUT;
% Save Variables
		    FIRSTF=fnew;
		    OLDG=gradf;
		    GDOLD=gdnew;	

		    % If quadratic interpolation set PCNT
		    if OPTIONS(7)==0, PCNT=1; MATX=zeros(3,1);  MATL(1)=fnew; end
		else 
% Halve Step-length
		    how='Red_Step';
		    if fnew==FIRSTF,
		        if OPTIONS(1)>0,
		            disp('No improvement in search direction: Terminating')
		        end
		        status=1;
		    else
		        OPTIONS(18)=OPTIONS(18)/8;
		        if OPTIONS(18)<1e-8
		            OPTIONS(18)=-OPTIONS(18);
		        end
		    end
		end
		XOUT=OLDX+OPTIONS(18)*SD;
		if OPTIONS(1)>0,disp([num,how]),end

	end %----------End of Direction Update-------------------
	if OPTIONS(7)==0, PCNT=1; MATX=zeros(3,1);  MATL(1)=fnew; end
% Check Termination 
	if max(abs(SD))< OPTIONS(2) & (gradf'*SD) < OPTIONS(3) & max(abs(gradf)) < 10*(OPTIONS(3)+OPTIONS(2))
		if OPTIONS(1) > 0
		    disp('Optimization Terminated Successfully')  
		end
		status=1; 
	elseif OPTIONS(10)>OPTIONS(14)
		disp('maximum number of iterations has been exceeded');
		if OPTIONS(1)>0
		    disp('Increase OPTIONS(14)')
		end
		status=1;
	else

% Line search using mixed polynomial interpolation and extrapolation.
		if PCNT~=0
		    while PCNT > 0
		        x(:) = XOUT; f(:)  = eval(evalstr); OPTIONS(10)=OPTIONS(10)+1;
		        fnew = f'*f;
		        if fnew<OLDF'*OLDF, OX = XOUT; OLDF=f; end
		        [PCNT,MATL,MATX,steplen,fnew,how]=searchq(PCNT,fnew,OLDX,MATL,MATX,SD,GDOLD,OPTIONS(18),how);
		        OPTIONS(18)=steplen;
		        XOUT=OLDX+steplen*SD;
		        if fnew==FIRSTF,  PCNT=0; end
		    end
		    XOUT = OX;
		    f=OLDF;
		else
		    x(:)=XOUT; f(:) = eval(evalstr); OPTIONS(10)=OPTIONS(10)+1;
		end
	end
end
OPTIONS(8) = fnew;
XOUT=OLDX;
x(:)=XOUT;


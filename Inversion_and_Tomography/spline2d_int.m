function [z, l] = spline2d_int (x_out, y_out, arg_3, arg_4, arg_5, arg_6 )
%
% SPLINE2D	Gridding using Green's function for splines in tension
%
%	SPLINE2D will find a spline-based surface using continuous curvature splines
%	in tension (if set).  The algorithm uses the Green's function for the spline.
%	You can supply data constrains, slope constrains, or a mix of both.
%	Solution can be evaluated on a grid or at arbitrary locations
%
%	Use one of the following call formats:
%
%	z = spline2d_int (x_out, y_out, x_data, y_data, z_data)
%	z = spline2d_int (x_out, y_out, x_data, y_data, z_data, t)
%
% The input parameters are:
%
%	x_out	- Desired output x positions (can be vector or matrix)
%	y_out	- Desired output y positions (can be vector or matrix)
%
%	x_data	- x-coordinates of points with data constraints
%	y_data	- y-coordinates of points with data constraints
%	z_data	- data constraints at the above points
%
%	t	- tension to use, 0 <= t < 1
%		  if t is a vector of length 2 the second value is taken as the lengthscale
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%       Math. Geol., 30, 77-93.


IM = sqrt (-1);
length_scale = abs (max(x_out(:)) - min(x_out(:)) + IM * (max(y_out(:)) - min(y_out(:)))) / 50;

if (nargin == 5 | nargin == 10 | nargin == 7)	% No tension selected, set default
	t = 0;
	n_args = nargin;
else
	n_args = nargin - 1;
	t = eval (['arg_' int2str(nargin)]);
	if (length(t) == 2)	% User gave both tension and lengthscale
		length_scale = t(2);
		t = t(1);
	end
	if (t < 0.0 | t >= 1.0)
		error ('spline2d: tension must be 0 <= t < 1 !')
	end
end

% Now figure out what was passed and assign the data accordingly

if (n_args == 5 | n_args == 10)	% z_data supplied
	x_data = arg_3;
	y_data = arg_4;
	z_data = arg_5;
end
if (n_args == 7)	% only z_slope supplied
	x_slope = arg_3;
	y_slope = arg_4;
	i_slope = arg_5;
	j_slope = arg_6;
	z_slope = arg_7;
elseif (n_args == 10)	% z_slope supplied
	x_slope = arg_6;
	y_slope = arg_7;
	i_slope = arg_8;
	j_slope = arg_9;
	z_slope = arg_10;
end

% Misc initializations

p = sqrt (t / (1 - t));
p = p / length_scale;
n0 = 0;
n1 = 0;

% First we must enforce the use of column vectors for the data constrainst

if (n_args == 5 | n_args == 10)	% z_data supplied; check if we must transpose
	[m,n] = size (x_data); if (m < n), x_data = x_data'; end
	[m,n] = size (y_data); if (m < n), y_data = y_data'; end
	[m,n] = size (z_data); if (m < n), z_data = z_data'; end
	n0 = length (z_data);
end
if (n_args == 7 | n_args == 10)	% z_slope supplied; check if we must transpose
	[m,n] = size (x_slope); if (m < n), x_slope = x_slope'; end
	[m,n] = size (y_slope); if (m < n), y_slope = y_slope'; end
	[m,n] = size (i_slope); if (m < n), i_slope = i_slope'; end
	[m,n] = size (j_slope); if (m < n), j_slope = j_slope'; end
	[m,n] = size (z_slope); if (m < n), z_slope = z_slope'; end
	n1 = length (z_slope);
end

% Assembly final xp, yp, and zp vectors (possibly combination of data and slopes)

if (n_args == 10)	% z_data and z_slope supplied; put xyz in general point vector
	xp = [x_data; x_slope];
	yp = [y_data; y_slope];
	zp = [z_data; z_slope];
elseif (n_args == 7)	% z_slope supplied; put xyz in general point vector
	xp = x_slope;
	yp = y_slope;
	zp = z_slope;
else 			% z_data supplied; put xyz in general point vector
	xp = x_data;
	yp = y_data;
	zp = z_data;
end

% Now build the square n x n linear system that must be solved for the alpha's

disp ('build matrix')
tic
n = n0 + n1;
A = zeros (n, n);
if (n_args == 5 | n_args == 10)	% z_data supplied; build data matrix rows
	for i = 1:n0
		r = (abs ((xp(i) - xp) + IM * (yp(i) - yp)));
		A(i,:) = (spline2d_green (r, p))';
	end
end
if (n_args == 7 | n_args == 10)	% z_slope supplied; build slope matrix rows
	for i = 1:n1
		j = i + n0;
		dx = xp(j) - xp;
		dy = yp(j) - yp;
		A(j,:) = (spline2d_grad (dx, dy, i_slope(i), j_slope(i), p))';
	end
end
toc
disp ('solve matrix')
tic
% Done building square linear system, now solve it

alpha = A \ zp;
toc
if (nargout == 2)	% Return eigenvalues
	disp ('find eigenvalues')
	tic
	l = svd (A);
	toc
end

disp ('evaluate')
tic
% Now evaluate final solution at output locations

z = zeros (size(x_out));
for i = 1:length(alpha)
	r = abs ((x_out - xp(i)) + IM * (y_out - yp(i)));
	z = z + alpha(i) * spline2d_green (r, p);
end

toc
disp ('done')


function G = spline2d_green (x, c)
%       $Id: spline2dgreen.m,v 1.2 2006/03/31 23:11:02 pwessel Exp $
% SPLINE2DT_GREEN       Green function for 2-D spline in tension
%
%       g = spline2dt_green (x, c)
%
%       x       - abscissa values
%       c       - tension parameter (c = sqrt (t/(1-t))
%
%  See Wessel, P, D. Bercovici, 1998, Gridding with Splines in Tension : A Green function Approach,
%       Math. Geol., 30, 77-93.

% Translated to Matlab from C.
% spline2d_green computes the Green function for a 2-d spline possibly
% in tension, G(u) = Ko(u) + log(u), where u = c * x and c = sqrt (t/(1-t)).
% The modified Bessel function Ko of order zero is based on Num. Rec.
% All x must be >= 0.  When c = 0 it degenerates to x^2 * log(x)
% 

ic = 1/c;
g0 = 0.115931515658412420677337;        %/* log(2) - 0.5772156... */
mask = find(x == 0);
if length(mask)>0, x(mask) = ones(length(mask),1); end        
id = find(x<=2*ic);
cx= c*x(id);
t = cx.*cx;
y = 0.25*t;
z = t/14.0625;
G = zeros(size(x));
G(id) = (-log(0.5*cx) .* (z .* (3.5156229 + z .* (3.0899424 + z .* (1.2067492 + z .* (0.2659732 + z .* (0.360768e-1 + z .* 0.45813e-2))))))) + (y .* (0.42278420 + y .* (0.23069756 + y .* (0.3488590e-1 + y .* (0.262698e-2 + y .* (0.10750e-3 + y .* 0.74e-5))))));
id = find(x>2*ic);    
y = 2*ic./x(id);
cx = c*x(id);
G(id) = (exp (-cx) ./ sqrt (cx)) .* (1.25331414 + y .* (-0.7832358e-1 + y .* (0.2189568e-1 + y .* (-0.1062446e-1 + y .* (0.587872e-2 + y .* (-0.251540e-2 + y .* 0.53208e-3)))))) + log (cx) - g0;
if length(mask)>0, G(mask) = zeros(length(mask),1); end        



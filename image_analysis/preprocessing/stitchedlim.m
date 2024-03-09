function lowhigh = stitchedlim(varargin)
%Find limits to contrast stretch an image, while ignoring zero-padding
%Adapted from the STRETCHLIM function
%   LOW_HIGH = stitchedlim(I,TOL) returns a pair of gray values that can be
%   used by IMADJUST to increase the contrast of an image. TOL specifies
%   the percentage of values to saturate at high and low levels
%
%   LOW_HIGH = stitchedlim(I,TOL,BV) allows the zero (background) value to
%   be specified manually as BV
%
%   LOW_HIGH = stitchedlim(I,opts) may also be passed, where opts is a
%   struct that can have fields 'Tol' and 'blackvalue' corresponding to TOL
%   and BV
%
%   TOL = [LOW_FRACT HIGH_FRACT] specifies the fraction of the image to
%   saturate at low and high pixel values. The fractions are taken over
%   only nonzero values in the image
%
%   If TOL is a scalar, TOL = LOW_FRACT, and HIGH_FRACT = 1 - LOW_FRACT,
%   which saturates equal fractions at low and high pixel values.
%
%   If you omit the argument, TOL defaults to [0.01 0.99], saturating 2%,
%   and BV defaults to 0.
%
%   If TOL = 0, LOW_HIGH = [min(I(:)); max(I(:))].
%
%   To specify BV but use default TOL, use LOW_HIGH = stitchedlim(I,[],BV)
%
%   LOW_HIGH = STRETCHLIM(RGB,TOL) returns a 2-by-3 matrix of pixel value
%   pairs to saturate each plane of the RGB image. TOL specifies the same
%   fractions of saturation for each plane.
%
%   Class Support
%   -------------
%   The input image can be uint8, uint16, int16, double, or single, and must
%   be real and nonsparse. The output limits are double and have values
%   between 0 and 1.
%
%   Note
%   ----
%   If TOL is too big, such that no pixels would be left after saturating
%   low and high pixel values, then STRETCHLIM returns [0; 1].
%
%   Example
%   -------
%       I = imread('pout.tif');
%       J = imadjust(I,stretchlim(I),[]);
%       figure, imshow(I), figure, imshow(J)
%
%   See also BRIGHTEN, DECORRSTRETCH, HISTEQ, IMADJUST.

%   Copyright 1999-2014 The MathWorks, Inc.

[img,tol,bv] = ParseInputs(varargin{:});

if isa(img,'uint8')
    nbins = 256;
else
    nbins = 65536;
end

tol_low = tol(1);
tol_high = tol(2);
 
p = size(img,3);

if tol_low < tol_high
    ilowhigh = zeros(2,p);
    for i = 1:p                          % Find limits, one plane at a time
        im = img(:,:,i);
        N = imhist(im(im > bv),nbins);
        if sum(N) > 0
            cdf = cumsum(N)/sum(N); %cumulative distribution function
            ilow = find(cdf > tol_low, 1, 'first');
            ihigh = find(cdf >= tol_high, 1, 'first');
        else %if image is all zero (or <= bv)
            ilow = 0; ihigh = 0;
        end
        if ilow == ihigh   % this could happen if img is flat
            ilowhigh(:,i) = [1;nbins];
        else
            ilowhigh(:,i) = [ilow;ihigh];
        end
    end
    lowhigh = (ilowhigh - 1)/(nbins-1);  % convert to range [0 1]

else
    %   tol_low >= tol_high, this tolerance does not make sense. For example, if
    %   the tolerance is .5 then no pixels would be left after saturating
    %   low and high pixel values. In all of these cases, STRETCHLIM
    %   returns [0; 1]. See gecks 278249 and 235648.
    lowhigh = repmat([0;1],1,p);
end


%-----------------------------------------------------------------------------
function [img,tol,bv] = ParseInputs(varargin)

narginchk(1, 3);

img = varargin{1};
validateattributes(img, {'uint8', 'uint16', 'double', 'int16', 'single'}, {'real', ...
    'nonsparse','nonempty'}, mfilename, 'I or RGB', 1);
if (ndims(img) > 3)
    error(message('images:stretchlim:dimTooHigh'))
end

%default input values
defaulttol = [0.01 0.99];
defaultbv = 0;

if nargin == 1
    tol = defaulttol; %default
    bv = defaultbv; %default
elseif nargin == 2
    if isstruct(varargin{2})
        opts = varargin{2};
        if isfield(opts,'Tol')
            tol = opts.Tol;
        else
            tol = defaulttol;
        end
        if isfield(opts,'blackvalue')
            bv = opts.blackvalue;
        else
            bv = defaultbv;
        end
    elseif isnumeric(varargin{2})
        tol = varargin{2};
        bv = defaultbv;
    else
        error('second argument is of the wrong class (struct or numeric)')
    end
elseif nargin == 3
    tol = varargin{2};
    bv = varargin{3};
end
switch numel(tol)
    case 1
        tol(2) = 1 - tol;

    case 2
        if (tol(1) >= tol(2))
            error(message('images:stretchlim:invalidTolOrder'))
        end
    otherwise
        error(message('images:stretchlim:invalidTolSize'))
end

if ( any(tol < 0) || any(tol > 1) || any(isnan(tol)) )
    error(message('images:stretchlim:tolOutOfRange'))
end





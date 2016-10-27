function [Erz1, Erz2, Erz3] = rsBackPropMud( img1, img2, Hq1, Hq2, Hq3 )
% Implemented according to:
% [1]O. Mudanyali, D. Tseng, C. Oh, S. O. Isikman, I. Sencan, W. Bishara, C. Oztoprak, S. Seo, B. Khademhosseini, and A. Ozcan, “Compact, Light-weight and Cost-effective Microscope based on Lensless Incoherent Holography for Telemedicine Applications,” Lab Chip, vol. 10, no. 11, pp. 1417–1428, Jun. 2010.
% [2]S.-H. Lee and D. G. Grier, “Holographic microscopy of holographically trapped three-dimensional structures,” Optics express, vol. 15, no. 4, pp. 1505–1512, 2007.
% [3]F. C. Cheong, B. J. Krishnatreya, and D. G. Grier, “Strategies for three-dimensional particle tracking with holographic video microscopy,” Opt. Express, vol. 18, no. 13, pp. 13563–13573, 2010.
% Note: In [1] They use spatial frequency and in [2,3] spatial angular
% frequency
%
% The derivation of the formulae for the propagator (without refractive
% index n) can be found in (p. 60):
% [4]J. W. Goodman, Introduction to Fourier optics. McGraw-Hill, 1996.

img1 = imfilter(img1, fspecial('gaussian', 4, 0.5));
img2 = imfilter(img2, fspecial('gaussian', 4, 0.5));

Bq1 = fft2(img1);
Bq2 = fft2(img2);

Erz1 = Hq1.*Bq1;
Erz2 = Hq2.*Bq1;
Erz3 = Hq3.*Bq2;

Erz1 = ifft2(Erz1);
Erz2 = ifft2(Erz2);
Erz3 = ifft2(Erz3);

Erz1 = abs(Erz1);
Erz2 = abs(Erz2);
Erz3 = abs(Erz3);

end


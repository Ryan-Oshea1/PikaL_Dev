%convert png to actual value image
function output_image = png_value_convertor(image_path)

hyperspectralImageAE = imread(image_path);
hyperspectralImageAE_OG = hyperspectralImageAE;
hyperspectralImageAEHigh=bitshift(hyperspectralImageAE, 16-12);
hyperspectralImageAElow = bitshift(hyperspectralImageAE, 16-28);
hyperspectralImageAE =hyperspectralImageAEHigh+ hyperspectralImageAElow;
hyperspectralImageAE = hyperspectralImageAE;
output_image = hyperspectralImageAE;

end
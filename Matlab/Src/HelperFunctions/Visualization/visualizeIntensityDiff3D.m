function visualizeIntensityDiff3D(Image)
    im = Image;
    figure;
    mesh(im);
    grid on;
    xlabel("x"); ylabel("y"); zlabel("Intensities");
    title("Image Intensity Difference Profile");
end
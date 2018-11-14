function visualize(force1, force2, X, Y, Diff)
% Visualize Force Vector Field At each Location for the Constituent
    % image
    quiver(X, Y, force1, force2,"r"); axis image
    hax = gca; %get the axis handle
    imagesc(hax.XLim,hax.YLim,Diff); %plot the image within the axis limits
    hold on; %enable plotting overwrite
    quiver(X, Y, force1, force2,"--g",'LineWidth',2 ); %plot the quiver on top of the image (same axis limits)
    %imagesc(Diff);
end


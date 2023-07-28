function [] = Image_realScale(img1, tit1, scan0, varargin)

    vrange = [-60 0];

    if nargin==5
        scan = varargin{1};
        DataType = varargin{2};
        
        switch DataType
        case 'iq'
            env = abs(reshape(img1,[scan0.Nz, scan0.Nx]));
        case 'rf'
            reshaped_bf_data = reshape(img1,[scan0.Nz, scan0.Nx]);
            env = interp1(scan0.z_axis, tools.envelope(double(reshaped_bf_data)), scan.z_axis,'linear',0);
        end
        
        im = 20*log10(env./max(env(:)));
        x_lim = [min(scan.x_matrix(:)) max(scan.x_matrix(:))]*1e3; 
        z_lim = [min(scan.z_matrix(:)) max(scan.z_matrix(:))]*1e3;

        imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,im); 
        shading flat; colormap gray; caxis(vrange); hold on; axis equal manual;
        set(gca,'YDir','reverse');
        set(gca,'fontsize',10); axis([x_lim z_lim]); title([tit1,' env+log']); 

%         imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,im); 
%         shading flat; colormap gray; caxis(vrange); colorbar; hold on; axis equal manual;
%         xlabel('x [mm]'); ylabel('z [mm]'); set(gca,'YDir','reverse');
%         set(gca,'fontsize',10); axis([x_lim z_lim]); title([tit1,' env+log']);         


    else 
        im = reshape(abs(img1),[scan0.Nz,scan0.Nx]);
        im = 20*log10(im./max(im(:)));
        %-- setting axis limits (mm)
        x_lim = [min(scan0.x_matrix(:)) max(scan0.x_matrix(:))]*1e3; 
        z_lim = [min(scan0.z_matrix(:)) max(scan0.z_matrix(:))]*1e3;
        
%         imagesc((scan0.x_axis)*1e3,(scan0.z_axis)*1e3, im); 
%         shading flat; colormap gray; caxis(vrange); hold on; axis equal manual;  
%         xlabel('x [mm]'); ylabel('z [mm]'); set(gca,'YDir','reverse');
%         set(gca,'fontsize',10); axis([x_lim z_lim]); title(tit1);

        imagesc((scan0.x_axis)*1e3,(scan0.z_axis)*1e3, im); 
        shading flat; colormap gray; caxis(vrange);  hold on; axis equal manual;  
        set(gca,'YDir','reverse'); colorbar;
        set(gca,'fontsize',10);  axis([x_lim z_lim]); title(tit1); 

    end

end


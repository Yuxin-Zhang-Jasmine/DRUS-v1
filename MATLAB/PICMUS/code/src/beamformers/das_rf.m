function image = das_rf(scan,dataset,pw_indices)

    %-- Function which implements the conventional Delay And Sum (DAS) beamform technique with apodization in reception
    %-- The corresponding code is dedicated to the reconstrucion of dataset (rawdata) saved in RF format

    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    assert(isempty(dataset.modulation_frequency)||dataset.modulation_frequency==0,'The supplied dataset is not RF');

    %-- select the plane waves that will be used in each frame
    if nargin < 3
        pw_indices{1} = 1:dataset.firings;
    end

    %-- define scan based on time axis
    time = (0:(size(dataset.data,1)-1)).'/dataset.sampling_frequency+dataset.initial_time;
    z_axis= time*dataset.c0/2;
    rf_scan = linear_scan(scan.x_axis,z_axis);

    %-- receive apodization
    %-- dynamically expanding receive aperture with hanning apodization
    rx_f_number = 1.75;
    rx_aperture = rf_scan.z/rx_f_number;
    rx_aperture_distance = abs(rf_scan.x*ones(1,dataset.channels)-ones(rf_scan.pixels,1)*dataset.probe_geometry(:,1).');
    receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey25');

    %-- angular apodization -> no apodization
    angular_apodization = ones(rf_scan.pixels,dataset.firings); 

    %-- beamforming loop
    beamformed_data = zeros(rf_scan.pixels,length(pw_indices));
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    wb=waitbar(0,'DAS beamforming');

    for f=1:length(pw_indices)
        waitbar(f/length(pw_indices),wb,sprintf('DAS-RF beamforming %0.0f%%',f/length(pw_indices)*100));
        for pw=pw_indices{f}
            %-- transmit delay
            transmit_delay = rf_scan.z*cos(dataset.angles(pw))+rf_scan.x*sin(dataset.angles(pw));
            for nrx=1:dataset.channels
                %-- receive delay
                receive_delay = sqrt((dataset.probe_geometry(nrx,1)-rf_scan.x).^2+(dataset.probe_geometry(nrx,3)-rf_scan.z).^2);
                %-- total delay
                delay = (transmit_delay+receive_delay)/dataset.c0;
                %-- beamformed data
                beamformed_data(:,f) = beamformed_data(:,f)+angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
            end
            clc;
            disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])           
        end
    end
    close(wb);

    beamformed_data(isnan(beamformed_data))=0;

    %-- reshape
    reshaped_beamformed_data = reshape(beamformed_data,[numel(rf_scan.z_axis) numel(rf_scan.x_axis)  length(pw_indices)]);

    %-- compute envelope
    envelope_beamformed_data = tools.envelope(reshaped_beamformed_data);

    %-- interpolate the requested grid
    resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));
    for f=1:length(pw_indices)
        resampled_envelope_beamformed_data(:,:,f) = interp1(rf_scan.z_axis,envelope_beamformed_data(:,:,f),scan.z_axis,'linear',0);
    end

    %-- declare an us_image object to store the beamformed data
    image = us_image('DAS-RF beamforming');
    image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
    image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
    image.algorithm = 'Delay-and-Sum (RF version)';
    image.scan = scan;
    image.number_plane_waves = cellfun('length',pw_indices);
    image.data = resampled_envelope_beamformed_data;
    image.transmit_f_number = 0;
    image.receive_f_number = rx_f_number;
    image.transmit_apodization_window = 'none';
    image.receive_apodization_window = 'Tukey 25%';

end

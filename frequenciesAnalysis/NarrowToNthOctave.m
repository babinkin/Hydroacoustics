function bands = NarrowToNthOctave(p,freq,fexact_l_c_u)


% Compute the acoustic variebl per N octave band
len_band = size(fexact_l_c_u,1);
bands = zeros(len_band,1);
for band_i = 1:len_band
    bands(band_i) = 0;
    idx = find( freq >= fexact_l_c_u(band_i,1) & freq < fexact_l_c_u(band_i,3) );
    % If we have no 'measurement' point in this band:
    if ( isempty(idx) )
      fprintf('Warning: no point found in band centered at %4.3f\n',fexact_l_c_u(band_i,2));
    % If we have only 1 'measurement' point in this band:
    elseif ( length(idx) == 1 )
      fprintf('Warning: only one point found in band centered at %4.3f\n',fexact_l_c_u(band_i,2));
      bands(band_i) = p(idx);
    % If we have more than 1 'measurement' point in this band averge integrate:
    elseif ( length(idx) > 1 )
      for i = 1:(length(idx)-1)
        bands(band_i) = bands(band_i) + ...
		  ( freq(idx(1) + i) - freq(idx(1) + i - 1) )* ...
		  abs( p(idx(1) + i) + p(idx(1) + i - 1) )/2;
      end
      bands(band_i) = bands(band_i) / ( freq(idx(end)) - freq(idx(1)) );
    end
end


end
function[out_path] = elobes_writewav_wrapper(z,fs,out_path)

% forces writewav to scale the output and save as 32 bit floating point 
% values. Use readwave with 'g' option to recover the absolute values
mode_str = 'gesvL'; 

% scale the values so that this value maps to +/-1. Since floating point is 
% used there will be no clipping if this value is exceeded but it makes
% auditioning the files in players which don't understand voicebox's 'g'
% chunk easier
ref_amplitude = 10;

% ensure the folder structure is in place
[outdir,~] = fileparts(out_path);
check_output_dir_exists(outdir)

try
    v_writewav(z,fs,out_path,mode_str,[],[],ref_amplitude);
    %fprintf('\n\nProcessed audio saved to %s\n',out_path);
catch err
    fprintf('Error saving file to %s',out_path);
    rethrow(err)
end



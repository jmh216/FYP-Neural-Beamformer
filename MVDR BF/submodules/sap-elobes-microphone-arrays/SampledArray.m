classdef SampledArray < ElobesMicArray
    % SampledArray provides methods for interpolation for arrays whose 
    % manifold has been sampled
    % 
    properties (SetAccess = protected)
        maxSensorRadius; % distance of the microphone furthest from origin
        H = [];          % [fix(nfft/2)+1 nSensors nSrcPos] matrix containing rfft of sampled data
        srcUnitVec;      % [nSrc 3] cartesian direction of arrival of measurement samples
        interp_fcn;      % function handle to for obtaining response at interpolated direction of arrival
    end
    methods (Abstract)
        [ir,src_pos,fs] = loadSampledData(obj)
        % returns the raw impulse responses and their positions
        %       ir: impulse responses in [nSamples,nMic,nPos] matrix
        %  src_pos: source positions as [nPos x 3] matrix in cartesian
        %           co-ordinates
        %       fs: sample rate of the data (in Hz)
    end

    methods
        function[obj] = SampledArray()
            % call superclass's constructor
            obj = obj@ElobesMicArray();
            
        end
        % preparing data depends on selected interpolation method
        % - none/nearest neighbour - nothing to be done to the data
        % - others to be added
        function[] = prepareData(obj,req_fs,varargin)
            % acquire and resample as necessary
            obj.fs = req_fs;
            [ir,srcPos,in_fs] = obj.loadSampledData();
            [nSamples,nMic,nPos] = size(ir);
            
            % resampling must be done as 2D matrix
            ir = resample(reshape(ir,nSamples,nMic*nPos),req_fs,in_fs);
            nSamples = size(ir,1);
            ir = reshape(ir,nSamples,nMic,nPos);
            obj.nfft = size(ir,1);
            obj.H = rfft(ir,[],1);
            obj.f = (0:((obj.nfft+2)/2 - 1)).' * obj.fs/obj.nfft; %[nf 1] where nf is number of bins up to nyquist
            
            % can't deal with different distances yet
            srcRadius = sqrt(sum(srcPos.^2,2));           
            if max(abs(srcRadius-mean(srcRadius))) > 20*eps
                error('Require fixed measurement radius, for now')
            end
            obj.srcUnitVec = bsxfun(@rdivide, srcPos, srcRadius);
            
            switch obj.interpolationMethod
                case 'none'
                	obj.interp_fcn = @obj.interp_exact_match; 
                case 'nearest_neighbour'
                    obj.interp_fcn = @obj.interp_nearest_neighbour;
                otherwise
                    error('Unknown interpolation method')
            end
            
            
%             
%             % define the response at the origin (if it hasn't been
%             % specified) which will form the basis of the other
%             % microphone responses
%             
%             % f0 of sinc determines the bandwidth
%             % len_filt determines the resolution/roll off
%             f0=0.95*obj.fs/2;
%             len_filt = 99;
%             half_len = (len_filt-1)/2;
%             % the pulse at the origin
%             h = hamming(len_filt) .* (sinc(2*f0*(-half_len:half_len).'/obj.fs));
%             
%             % need to zero pad both sides to allow for relative propagation
%             % delay to furtherst sensor
%             pad_len = ceil(obj.maxSensorRadius / obj.c * obj.fs);
%             h = [zeros(pad_len,1); h; zeros(pad_len,1)];
%             obj.nfft = size(h,1);
%             obj.t0 = pad_len+half_len+1; % sample offset of main peak
%             obj.H0 = rfft(h);
            
            
        end
        function[H,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc(obj,src_az,src_inc,varargin)
             % account for the rotation (if any)
             [rel_src_az,rel_src_inc] = obj.planeWaveDoaWrtArray(src_az,src_inc);
             H = obj.interp_fcn(rel_src_az,rel_src_inc);
        end
                

        function[solidAngles] = targetAngleWrtSamples(obj,targetAz,targetInc)
            
            targetUnitVec = mysph2cart(targetAz,targetInc,ones(size(targetAz)));
            solidAngles = acos(obj.srcUnitVec*targetUnitVec.');
        end  
            
        function[H] = interp_nearest_neighbour(obj,src_az,src_inc)
            solidAngles = obj.targetAngleWrtSamples(src_az,src_inc);
            [~,imin] = min(solidAngles,[],1);
            H = obj.H(:,:,imin);
        end    
        function[H] = interp_exact_match(obj,src_az,src_inc)
            threshold = 2*eps;
            solidAngles = obj.targetAngleWrtSamples(src_az,src_inc);
            [angle,imin] = min(solidAngles,1);
            if any(angle>threshold)
                error('No exact match for requested source direction')
            end
            H = obj.H(:,:,imin);
        end    
    end
end
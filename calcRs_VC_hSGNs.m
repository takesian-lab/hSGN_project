function [Rin, Rs, Cm, error, tau]=calcRs_VC_hSGNs(data, samplingRate, pulseAmp, pulseStart, pulseLength, baselineStart, baselineEnd)
		Rin=0;
		Rs=0;
		Cm=0;
		error=0; %will turn into 1 if Rin, rs and cm can't be calculated

    %e    waveo('calcRsFit', nan(1,10) , 'xscale', [pulseStart samplingRate]);

		dataPulseStart = 1+round(pulseStart/samplingRate);
		dataPulseEnd = 1+round((pulseStart+pulseLength)/samplingRate);
		
%         if dataPulseStart>size(data, 2)
% 			disp('calcRs: start of pulse is beyond end of data')
% 			error=1;
% 			return
% 		end
% 		
% 	   if dataPulseEnd>size(data, 2)
% 			disp('calcRs: end of pulse is beyond end of data.  Fitting to available data')
% 			dataPulseEnd=size(data, 2);
% 	   end
		
		if pulseAmp>0
			[peak, peakloc]=max(data(dataPulseStart:dataPulseEnd));
		else
			[peak, peakloc]=min(data(dataPulseStart:dataPulseEnd));
		end
		peakloc=peakloc+round(pulseStart/samplingRate)-1;
		
		delta=round(0.1*pulseLength/samplingRate);
		endline=mean(data(round((pulseStart+pulseLength*0.7)/samplingRate)+1:round((pulseStart+pulseLength)/samplingRate)));
		baseline=mean(data(round(baselineStart/samplingRate)+1:round(baselineEnd/samplingRate)-1));
		
		if pulseAmp>0
			above=find(data(dataPulseStart:round((pulseStart+pulseLength)/samplingRate))>(peak-endline)/5+endline);
			if isempty(above)
				error=1;
				return;
			end
			last=above(end)+round(pulseStart/samplingRate);
		else
			above=find(data(dataPulseStart:round((pulseStart+pulseLength)/samplingRate))<(peak-endline)/5+endline);
			if isempty(above)
				error=1;
				return;
			end
			last=above(end)+round(pulseStart/samplingRate);
		end

		delta=round((last-peakloc)/2);
		clear above last
		
		peak1=peak-endline;
		peak2=data(delta+peakloc)-endline;
		peak3=data(2*delta+peakloc)-endline;
		peakloc=peakloc-round(pulseStart/samplingRate)-1;
		
        if (peak1-peak2)*(peak2-peak3)<0
 			disp('calcRs: 3 points are non-monotonic.  No fit possible')
            error=1;
            return
        end
        
		tau=delta*(1/log(peak1/peak2)+1/log(peak2/peak3)+2/log(peak1/peak3))/3;
		amp=(peak1/exp(-peakloc/tau)+peak2/exp(-(peakloc+delta)/tau)+peak3/exp(-(peakloc+2*delta)/tau))/3;
		
		Rs=1000*pulseAmp/amp;
		Rin=1000*pulseAmp/(endline-baseline)-Rs;
		Cm=1000*tau*samplingRate/Rs;

        if any(~isreal([tau amp Rs Rin Cm]))
  			disp('calcRs: complex number returned.  No fit possible')
            Rin=0;
            Rs=0;
            Cm=0;
            error=1;
            return
        end
        
	%	waveo('calcRsFit', endline+amp*exp(-[0:round(pulseLength/samplingRate)]/tau), 'xscale', [pulseStart samplingRate]);
		%clear tau amp peak peak1 peak2 peak3 peakloc endlkine baseline
        clear amp peak peak1 peak2 peak3 peakloc endlkine baseline    

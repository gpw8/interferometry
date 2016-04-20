function stfun=stfunc(nsamp,delta,tp,ts,wavetype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is based on a subroutine in Ohminato and Chouet, BSSA, 1998
% nsamp - number of samples
% delta - time between evenly-spaced samples (1/sample rate)
% tp - time constant to define width of wavelet
% ts -time shift
% wavetype - type of wavelet
%  wavetype = 1: Gaussian               2: Ricker
%             3: 1-cos(t) (One cycle)   4: 1-cos(t) (Half cycle)
%             5: Smoothed step          6: Ramp
%             7: Triangle               8: Others
%
%
%
stfun=zeros(1,nsamp);
nt=1:nsamp;
t0all=delta*(nt-1);
for nt=1:nsamp
    t0=t0all(nt);
    if(wavetype == 1)
        %       tp=0.4
        %       ts=2.0
        t=t0-ts;
        stfun(nt)=exp(-t^2/tp^2);
        %
        %                         ricker wavelet
        %
    elseif(wavetype == 2)
        %       tp=1.83
        %       ts=2.0
        wc=2*pi/tp;
        t=t0-ts;
        dum=(wc*t/2)^2;
        stfun(nt)=(2*dum-1)*exp(-dum);
        %
        %			  1-cos(t)(One cycle)
        %
    elseif(wavetype == 3)
        %       tp=1
        t=t0-(ts-tp/2);
        if((t<0) || (t>tp))
            stfun(nt)=0;
        else
            stfun(nt)=(1-cos(2*pi*t/tp))/2;
        end
        %
        %			  1-cos(t)(Half cycle)
        %
    elseif(wavetype == 4)
        %       tp=3.d0
        t=t0-(ts-tp/2);
        if(t<0)
            stfun(nt)=0;
        elseif(t<tp)
            stfun(nt)=(1.-cos(pi*t/tp))/2;
        else
            stfun(nt)=1.0;
        end
        %
        %                         smoothed step function
        %
    elseif(wavetype == 5)
        %       tp=0.25
        %       ts=2
        t=t0-ts;
        stfun(nt)=(1.+tanh(t/tp))/2;
        %
        %                         ramp function
        %
    elseif(wavetype == 6)
        %       tp=0.5d0
        %       ts=2.0
        t=t0-(ts-tp/2);
        if(t<0)
            stfun(nt)=0;
        elseif(t<tp)
            stfun(nt)=t/tp;
        else
            stfun(nt)=1;
        end
        %
        %                         Triangle
        %
    elseif(wavetype == 7)
        %       tp=0.5d0
        t=t0-(ts-tp/2);
        if(t<0)
            stfun(nt)=0;
        elseif(t<tp/2)
            stfun(nt)=2*t/tp;
        elseif(t<tp)
            stfun(nt)=2*(1-t/tp);
        else
            stfun(nt)=0;
        end
        %
        %                         Others
        %
    elseif(wavetype == 8)
        %       tp=1.0d0
        %       ts=2.0
        t=t0-(ts-tp/2);
        stfun(nt)=-2.*t/(tp^2)*exp(-t^2/tp^2);
        %       if(t<0)
        %         stfun(nt)=0
        %       elseif(t<tp)
        %         stfun(nt)=1.0/tp
        %       else
        %         stfun(nt)=0
        %       end
    end
end
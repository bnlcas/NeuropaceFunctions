function ecog=applyLineNoiseNotch_60HzHarmonics(ecog,sampFreq)

% attempted to use GPU
%{
keyboard
%decide if you can use GPU
if parallel.gpu.GPUDevice.isAvailable
    props = gpuDevice();
    data = ecog.data;
    p = whos('data');
    useGPU = props.FreeMemory > p.bytes;
else
    useGPU = 0;
end


tic;

if useGPU
    ecog.data = gpuArray(ecog.data);
end
%}

notchFreq=60;
while notchFreq<sampFreq/2
    fprintf(['Notch Filters: ' int2str(notchFreq) 'Hz.. \n'])
    [b,a]=fir2(1000,[0 notchFreq-1 notchFreq-.5 notchFreq+.5 notchFreq+1 sampFreq/2]/(sampFreq/2),[1 1 0 0 1 1 ]);
    ecog.data=filtfilt(b,a,ecog.data')';
    notchFreq=notchFreq+60;
end
fprintf('\nNotch Filters Applied\n')
%{
if useGPU
    ecog.data = gather(ecog.data);
end
toc;
%}
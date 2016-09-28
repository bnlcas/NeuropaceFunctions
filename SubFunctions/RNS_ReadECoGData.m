function [ECoG_hdr, ECoG_data] = RNS_ReadECoGData(DATFile, LAYFile)

%Example Usage
%
%DATFile = '100315/127392201131430000.DAT';
%LAYFile = '100315/127392201131430000.lay';
%[ECoG_hdr, ECoG_data] = ReadECoGData(DATFile, LAYFile)
%

lay = textread(LAYFile,'%s','delimiter','\n','whitespace','');

ECoG_hdr.NPConfigString = ReadHeader('NPConfigStr', lay, 's');
ECoG_hdr.PatientInitials = ReadHeader('PatientInitials', lay, 's');
ECoG_hdr.DeviceID = ReadHeader('DeviceSerialNumber', lay, 'n');
ECoG_hdr.TimeStampPatientLocalString = ReadHeader('ECoGTimeStampAsLocalTime', lay, 's');
ECoG_hdr.SamplingRate = ReadHeader('SamplingRate', lay, 'n');
ECoG_hdr.WaveformCount = ReadHeader('WaveformCount', lay, 'n');
ECoG_hdr.ChannelMap = ReadChannelMap(lay);
ECoG_hdr.TriggerReason = ReadHeader('TriggerReason', lay, 's');
ECoG_hdr.Annotations = ReadHeader('Annotations', lay, 's');
if strcmp(ReadHeader('AmplifierChannel1ECoGStorageEnabled',lay, 's'), 'On')
    ECoG_hdr.EnabledChannels(1) = 1; else ECoG_hdr.EnabledChannels(1) = 0; 
end
if strcmp(ReadHeader('AmplifierChannel2ECoGStorageEnabled',lay, 's'), 'On')
    ECoG_hdr.EnabledChannels(2) = 1; else ECoG_hdr.EnabledChannels(2) = 0; 
end
if strcmp(ReadHeader('AmplifierChannel3ECoGStorageEnabled',lay, 's'), 'On')
    ECoG_hdr.EnabledChannels(3) = 1; else ECoG_hdr.EnabledChannels(3) = 0; 
end
if strcmp(ReadHeader('AmplifierChannel4ECoGStorageEnabled',lay, 's'), 'On')
    ECoG_hdr.EnabledChannels(4) = 1; else ECoG_hdr.EnabledChannels(4) =0; 
end
if findstr(ECoG_hdr.TriggerReason, 'USER_SAVED'); 
    ECoG_hdr.EnabledChannels = [1 1 1 1];
end



fid = fopen(DATFile);

dat = fread(fid,'int16');
fclose(fid);   
                         
% populate channels
ChannelNum = 0;
for ChannelIndex = 1:4
     if ECoG_hdr.EnabledChannels(ChannelIndex)
         ChannelNum = ChannelNum + 1;
%      	 ECoG_data{ChannelIndex} = dat(ChannelNum:ECoG_hdr.WaveformCount:end)'-512;
        ECoG_data(ChannelIndex,:) = dat(ChannelNum:ECoG_hdr.WaveformCount:end)'-512;
     else
%          ECoG_data{ChannelIndex} = [];
        ECoG_data(ChannelIndex,:) = [];
     end
end

%--------------------------------------------------------------------------
% function [value] = ReadHeader(VariableName, lay, type)
% function to read variables from .lay header
%--------------------------------------------------------------------------


function [value] = ReadHeader(VariableName, lay, type)
switch type
    case 's'
        myvalue{1} = '';
    case 'n'
        myvalue{1} = 0;
end
VariableNum = 0;
for i = 1:length(lay)
    if ~isempty(strfind(lay{i}, [VariableName '=']))
        VariableNum = VariableNum+1;
        switch type
            case 'n'
                myvalue{VariableNum} = str2num(lay{i}(length(VariableName) + 2:end));
            case 's'
                myvalue{VariableNum} = lay{i}(length(VariableName) + 2:end);
        end
        
    end
end
if VariableNum <= 1
    value = myvalue{1};
else
    value = myvalue;
end

%--------------------------------------------------------------------------
% function [ChannelMap] = ReadChannelMap(lay)
% function to read ChannelMap from .lay header
%--------------------------------------------------------------------------

function [ChannelMap] = ReadChannelMap(lay)
ChannelMap{1} = ''; ChannelMap{2} = ''; ChannelMap{3} = ''; ChannelMap{4} = '';
for i = 1:length(lay)
    if ~isempty(findstr(lay{i},'[ChannelMap]'))
        for j = 1:4
            if ~isempty(findstr(lay{i+j},'='))
              ChannelNumber = str2num(lay{i+j}(findstr(lay{i+j},'=')+1:end));
              ChannelMap{ChannelNumber} = lay{i+j}(1:findstr(lay{i+j},'=')-1);
            end
        end        
    end
end

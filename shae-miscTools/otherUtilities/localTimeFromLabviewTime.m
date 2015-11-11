function localTime = localTimeFromLabviewTime(inputTime);
% there could be issues with daylight savings time

% define constants
LV_EPOCH_START = datenum('12:00 am Jan 1 1904');
JV_EPOCH_START = datenum('12:00 am Jan 1 1970');

% current time UTC 
utc_time_ms = java.lang.System.currentTimeMillis;
% current time in matlab
currentMatlabTimePartialDays = now();
currentMatlabTimeHours = currentMatlabTimePartialDays*24;
% find the difference (in hours)
utc_time_partialDays = utc_time_ms./(1000*60*60*24);
utc_time_partialDaysSinceMatlabStart = utc_time_partialDays+JV_EPOCH_START;
utc_time_HoursSinceMatlabStart = utc_time_partialDaysSinceMatlabStart*24;

utcOffsetHours = round(utc_time_HoursSinceMatlabStart-currentMatlabTimeHours);
utcOffsetDays = utcOffsetHours/24;


inputTimePartialDays = inputTime./(60*60*24);
inputTimePartialDaysMatlabStart = inputTimePartialDays+LV_EPOCH_START;

localTime = inputTimePartialDaysMatlabStart-utcOffsetDays;




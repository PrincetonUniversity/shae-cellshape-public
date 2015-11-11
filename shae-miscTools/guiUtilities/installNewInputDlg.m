function installNewInputDlg
% this function will retrieve and install newid, the new input dialog
% where pressing enter is the same as clicking the Okay button

inputdlgpaths = which('inputdlg.m','-all');
if length(inputdlgpaths)~=1
    error('installNewInputDlg:unknownDestination','Incorrect number of destinations, please install manually');
end
%%
urlwrite('http://www.mathworks.com/matlabcentral/answers/uploaded_files/1727/newid.m',...
    [tempdir,filesep,'newid.m'])
[inputdlgpathParent] = fileparts(inputdlgpaths{1});
%%
movefile([tempdir,filesep,'newid.m'],inputdlgpathParent);
rehash toolboxcache;
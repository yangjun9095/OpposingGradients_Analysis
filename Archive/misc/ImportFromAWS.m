function ImportFromAWS(Prefix)

%Get folders
    [~,UserProcPath,UserDynResPath,~,UserPreProcPath] = ...
        DetermineLocalFolders(Prefix);
    UserPreProcPath_Prefix = [UserPreProcPath,filesep,Prefix];
    UserProcPath_Prefix = [UserProcPath,filesep,Prefix,'_'];
    UserDynResPath_Prefix =[UserDynResPath,filesep,Prefix];
    [~, username] = system('echo %USERNAME%');
    username = strrep(username, sprintf('\n'),''); %removes new line
    HGlabDataFolder = ['E:\HGlab\Dropbox\',username,'\LivemRNA\Data'];
    HGlabPreProcPath_Prefix = [HGlabDataFolder,filesep,'PreProcessedData\',Prefix];
    HGlabProcPath_Prefix = [HGlabDataFolder,filesep,'ProcessedData\',Prefix,'_'];
    HGlabDynResPath_Prefix = [HGlabDataFolder,filesep,'DynamicsResults\',Prefix];
    
%Move files from HGlab Dropbox to users Data folder
%Need to find a way that will transfer the files in these folders to the folders they were copied from...
    movefile(HGlabPreProcPath_Prefix, UserPreProcPath_Prefix) % try UserProcPath_AWS
    movefile(HGlabProcPath_Prefix, UserProcPath_Prefix) % try UserPreProcPath_AWS
    movefile(HGlabDynResPath_Prefix, UserDynResPath_Prefix) % try UserDynResPath_AWS

end
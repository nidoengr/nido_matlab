myUrl = ...'http://google.com';
'https://www.google.com/search?q=speed+test&rlz=1C1CHBF_enAE925US932&oq=speed&aqs=chrome.0.35i39j69i57j35i39j0i67j46i67i175i199j46i67j0i67j0i67i131i433j0i67l2.1144j0j7&sourceid=chrome&ie=UTF-8';
% web(myUrl,'-browser');

% [stat,h,url] = web(myUrl,'-browser') %<--This didnt work
% % Warning: [STAT,H,URL] = WEB(___) does not return a handle
% % or URL for pages that open in the system browser. Use STAT
% % = WEB(___) instead. 
% % > In web>displayWarningMessage (line 438)
% % In web (line 99)
% % In running_average_internet_speed_test (line 5) 

% % This just opens the url above - no straightforward way to
% interact
% stat = web(myUrl,'-browser'); 

[system_status, cmd_out] = system('speed-test','-echo')

% strfind(cmd_out,'Download')
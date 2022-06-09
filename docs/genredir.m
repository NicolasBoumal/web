pattern = '<html>\r\n<head>\r\n<title>Redirection to %s.pdf</title>\r\n<meta http-equiv="refresh" content="0; URL=%s.pdf">\r\n\r\n<script type="text/javascript">\r\nvar _gaq = _gaq || [];\r\n_gaq.push([''_setAccount'', ''UA-19509417-1'']);\r\n_gaq.push([''_trackPageview'']);\r\n(function() {\r\nvar ga = document.createElement(''script''); ga.type = ''text/javascript''; ga.async = true;\r\nga.src = (''https:'' == document.location.protocol ? ''https://ssl'' : ''http://www'') + ''.google-analytics.com/ga.js'';\r\nvar s = document.getElementsByTagName(''script'')[0]; s.parentNode.insertBefore(ga, s);\r\n})();\r\n</script></head>\r\n<body>\r\nIf you are not redirected, please click on <a href="%s.pdf">%s.pdf</a>\r\n\r\n</body>\r\n</html>\r\n';

lst = cellstr(ls('*.pdf'));
N = length(lst);
for i = 1 : N
    file = lst{i};
    name = file(1:end-4);
    fid = fopen([name '.htm'], 'w+');
    fprintf(fid, pattern, name, name, name, name);
    fclose(fid);
end

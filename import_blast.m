function out = import_blast(process, coln, ifname, ofname, del)
global fatal_error
global blast_length

unix_cmd = ['awk -F"\t" ''{ print $' num2str(coln) ' }'' ' ifname ' > ' ...
    ofname];
[st, w] =  unix(unix_cmd);
if st ~= 0
    fatal_error = 1;
    fatal_msg(process, ...
        {['Failed reading bit scores of alignments from "' ifname '"'], ...
        w});
    return
end
[fid, w] = fopen(ofname, 'r');
if fid == -1
    fatal_error = 1;
    fatal_msg(process, {'Failed loading tmp file', w});
    return
end
out = textscan(fid, del, 'delimiter', '\n', 'bufsize', blast_length);
out = out{1};
fclose(fid);
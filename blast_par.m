function blast_par(process)
global last_process
global folder
global blast_dir
%global blast_db
global blast_length
global blast_cutoff
global blast_cpus
global organism
global fatal_error
global OS

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'BLAST parameters loaded');
    write_log(process);
    return
end

%% Import target sequences from fasta file
%blast_dir = read_config('blast_dir', process);
if blast_dir(end) ~= '/'
    blast_dir = [blast_dir '/'];
end
%blast_db = read_config('blast_db', process);
%fasta_db = read_config('fasta_db', process);
%blast_cpus = uint32(str2double(read_config('blast_cpus', process)));
blast_length = uint32(str2double(read_config('blast_length', process)));
blast_cutoff = uint32(str2double(read_config('blast_cutoff', process)));
blast_cutoff = check_value(process, blast_cutoff, ...
    'blast_cutoff', [200, 100, 1000], [100, 2000]);
organism = read_config('organism', process);

%% Get OS type
OS = get_os_type(process);

if OS == 1
	nr_mac_cpus(process)
elseif OS == 2
    nr_linux_cpus(process)
end

%% Get number of computer cores
if OS > 0
    fid = fopen([folder 'nr_cpus.out'], 'r');
        if fid ~= -1
            st = 0;
        else
            st = 1;
        end
        error_msg(process, st, {'Failed getting the system type'});
    	cpus = textscan(fid, '%d');
        cpus = cpus{1};
    fclose(fid);
    delete([folder 'nr_cpus.out']);
else
    cpus = 0;
end

if cpus < 1
    cpus = 1;
    error_msg(process, 1, ...
        {'Could not acquire the number of available processor cores', ...
        'Assuming there is one CPU in your computer'});
end

%% Check if any of the parameters are left empty
if isempty(blast_cpus) || ~(blast_cpus > 0) || blast_cpus > cpus ...
        || length(blast_cpus) > 1
    blast_cpus = cpus;
    error_msg(process, st, ...
        {'You have not provided valid number of CPUs', ...
        ['Will use ' num2str(cpus) ' CPUs']});
end

%% Check maximum blast length allowance
if isempty(blast_cpus) || blast_length < 1000 || length(blast_length) > 1
    blast_length = 1000;
    error_msg(process, st, ...
        {'Your maximum target length is too short', ...
        ['Will use maximum target length of ' ...
        num2str(blast_length) ' bp']});
end

%% Check maximum blast length allowance
if isempty(blast_cpus) || blast_length > 5000000 || ...
        length(blast_length) > 1
    blast_length = 5000000;
    error_msg(process, st, ...
        {'Your maximum target length is too long', ...
        ['Will use maximum target length of ' ...
        num2str(blast_length) ' bp']});
end

%% Check is organism is set
if isempty(organism)
    fatal_error = 1;
    fatal_msg(process, {'You have not chosen valid organism'});
    return
end

save([folder 'mats/fp' num2str(process) '.mat'], 'blast_dir', ...
    'organism', 'blast_length', 'blast_cpus', 'OS', 'blast_cutoff');
write_log(process, 'Imported BLAST parameters');

%% Process finished
write_log(process);
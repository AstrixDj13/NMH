filename = 'userInput.txt';

opts = detectImportOptions(filename);

opts = setvartype(opts,{'u0','gamma','points','scheme'},...
                        {'double','double','int16','string'});

userInput = readtable(filename,opts);

velocity = table2array(userInput(:,'u0'));
gamma = table2array(userInput(:,'gamma'));
points = table2array(userInput(:,'points'));
scheme = string(table2array(userInput(:,'scheme')));


%Copyright (C) 2022 by Frida Viset

if ispc
    run('tools\load_tools.m');
    run('ProcessData\RunFilters.m');
else
    run('tools/load_tools.m');
    run('ProcessData/RunFilters.m');
end

function [wgts, spheres] = load_wgts_spheres(EEG, icatype, ica_path)
% load weights & spheres

init_dir = pwd;

switch icatype
    case 'amica'
        EEG = pop_loadmodout10(EEG,ica_path);
        wgts = EEG.icaweights;
        spheres = EEG.icasphere;
    case 'cudaica'
        cd(ica_path)
        weightsfile = dir('*.wts');
        spherefile = dir('*.sph');
        nchans = EEG.nbchan;
        ncomps = nchans;
        wgts = floatread(weightsfile.name,[ncomps Inf],[],0);
        spheres = floatread(spherefile.name,[nchans Inf],[],0);
end

cd(init_dir)
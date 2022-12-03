package org.exbio.tfprio.configs;

import org.exbio.pipejar.configs.ConfigModuleCollection;
import org.exbio.tfprio.configs.modules.DeSeq2;
import org.exbio.tfprio.configs.modules.InputConfigs;
import org.exbio.tfprio.configs.modules.MixOptions;

public class Configs extends ConfigModuleCollection {
    public static final InputConfigs inputConfigs = new InputConfigs();
    public static final MixOptions mixOptions = new MixOptions();
    public static final DeSeq2 deSeq2 = new DeSeq2();
}
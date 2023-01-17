package org.exbio.tfprio.steps.report;

import org.exbio.pipejar.configs.ConfigTypes.FileTypes.InputFile;
import org.exbio.pipejar.configs.ConfigTypes.FileTypes.OutputFile;
import org.exbio.pipejar.pipeline.ExecutableStep;

import java.util.Collection;
import java.util.HashSet;
import java.util.concurrent.Callable;

import static org.exbio.pipejar.util.ScriptExecution.executeAndWait;

public class ReportCreation extends ExecutableStep {
    public final OutputFile outputFile = addOutput("report");
    private final InputFile reportDirectory;

    public ReportCreation(OutputFile reportDirectory) {
        super(false, reportDirectory);

        this.reportDirectory = addInput(reportDirectory);
    }

    @Override
    protected Collection<Callable<Boolean>> getCallables() {
        return new HashSet<>() {{
            add(() -> {
                String command = "ng build --prod --output-path " + reportDirectory.getAbsolutePath();

                executeAndWait(command, true);

                return true;
            });
        }};
    }
}
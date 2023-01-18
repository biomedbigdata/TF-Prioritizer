package org.exbio.tfprio.steps.report;

import org.apache.commons.lang3.tuple.Pair;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

import static org.exbio.pipejar.util.FileManagement.extend;
import static org.exbio.pipejar.util.FileManagement.softLink;

public class TfGroup {
    private final String symbol;
    private final Collection<TranscriptionFactor> transcriptionFactors;
    private final Map<String, Double> hmDcg = new HashMap<>();
    private final File outDir;
    private final Path base;
    private final Map<String, Map<String, String>> hmPairingHeatmap = new HashMap<>();
    private String biophysicalLogo;
    private Map<String, String> tfSequence;

    public TfGroup(String symbol, Collection<TranscriptionFactor> transcriptionFactors, File srcDir) {
        this.symbol = symbol;
        this.transcriptionFactors = transcriptionFactors;
        this.outDir = extend(srcDir, "assets", "tfs", symbol);
        this.base = srcDir.toPath();
    }

    public void setBiophysicalLogo(File biophysicalLogo) throws IOException {
        File out = extend(outDir, "biophysicalLogo.png");
        softLink(out, biophysicalLogo);
        this.biophysicalLogo = base.relativize(out.toPath()).toString();
    }

    public void setTfSequence(Collection<File> tfSequence) {
        File outDir = extend(this.outDir, "tfSequence");
        this.tfSequence = tfSequence.stream().map(f -> {
            String name = f.getName().substring(0, f.getName().lastIndexOf('.'));
            File out = extend(outDir, f.getName());
            try {
                softLink(out, f);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
            return Pair.of(name, base.relativize(out.toPath()).toString());
        }).collect(Collectors.toMap(Pair::getKey, Pair::getValue));
    }

    public void setHeatmap(String pairing, String hm, File file) {
        File out = extend(this.outDir, "heatmaps", hm, pairing, "heatmap.png");
        try {
            softLink(out, file);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        hmPairingHeatmap.computeIfAbsent(hm, k -> new HashMap<>()).put(pairing,
                base.relativize(out.toPath()).toString());
    }


    public void setRank(String hm, Double rank) {
        hmDcg.put(hm, rank);
    }

    public Collection<TranscriptionFactor> getTranscriptionFactors() {
        return transcriptionFactors;
    }

    public String getSymbol() {
        return symbol;
    }

    public JSONObject toJSON() {
        return new JSONObject() {{
            put("symbol", symbol);
            put("hmDcg", hmDcg);
            put("transcriptionFactors", transcriptionFactors.stream().map(TranscriptionFactor::toJSON).toList());
            put("biophysicalLogo", biophysicalLogo);
            put("tfSequence", tfSequence);
            put("heatmaps", hmPairingHeatmap);
        }};
    }
}

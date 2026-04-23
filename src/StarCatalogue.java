import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Stream;

/**
 * Temporary catalogue loader that simply forwards rows from {@code data.csv}.
 */
public final class StarCatalogue {
    private final Path csvPath;

    public StarCatalogue() {
        this(Path.of("data.csv"));
    }

    public StarCatalogue(Path csvPath) {
        this.csvPath = csvPath;
    }

    /**
     * Returns every row from {@code data.csv} without applying transformation.
     */
    public List<String> pumpData() throws IOException {
        return Files.readAllLines(csvPath);
    }

    /**
     * Streams rows from {@code data.csv} without applying transformation.
     */
    public Stream<String> pumpDataStream() throws IOException {
        return Files.lines(csvPath);
    }
}

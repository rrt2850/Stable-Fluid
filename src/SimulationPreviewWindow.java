import javax.swing.*;
import java.awt.image.BufferedImage;

/**
 * Small Swing window for displaying the latest simulation frame while the solver runs.
 */
public final class SimulationPreviewWindow {
    private final JFrame frame;
    private final JLabel imageLabel;

    public SimulationPreviewWindow(int width, int height) {
        try {
            final JFrame[] holder = new JFrame[1];
            final JLabel[] labelHolder = new JLabel[1];
            SwingUtilities.invokeAndWait(() -> {
                JFrame previewFrame = new JFrame("Stable Fluid - Live Preview");
                JLabel previewLabel = new JLabel();
                previewLabel.setPreferredSize(new java.awt.Dimension(width, height));
                previewFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                previewFrame.setContentPane(previewLabel);
                previewFrame.pack();
                previewFrame.setLocationByPlatform(true);
                previewFrame.setVisible(true);
                holder[0] = previewFrame;
                labelHolder[0] = previewLabel;
            });
            this.frame = holder[0];
            this.imageLabel = labelHolder[0];
        } catch (Exception exception) {
            throw new RuntimeException("Failed to initialize preview window.", exception);
        }
    }

    public void update(BufferedImage image, int step) {
        if (image == null) {
            return;
        }

        SwingUtilities.invokeLater(() -> {
            imageLabel.setIcon(new ImageIcon(image));
            frame.setTitle("Stable Fluid - Live Preview (Step " + step + ")");
        });
    }

    public void close() {
        try {
            SwingUtilities.invokeAndWait(frame::dispose);
        } catch (Exception exception) {
            throw new RuntimeException("Failed to close preview window cleanly.", exception);
        }
    }
}
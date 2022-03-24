/**
 *
 */
package org.theseed.shared.meta;

/**
 * This is a simple interface for reporting progress in a metabolic search.  Note that in a GUI environment
 * we expect these methods to be called from a background thread, so they will have to use a runnable to
 * update a display.  For example, the following is a typical override for showStatus:
 *
    @Override
    public void showStatus(String msg) {
        Platform.runLater(new Runnable() {
            @Override
            public void run() {
                txtMessageBuffer.setText(msg);
            }
        });
    }
 *
 * @author Bruce Parrello
 *
 */
public interface IProgressReporter {

    /**
     * Denote the current progress.
     *
     * @param p		fraction between 0 (no progress) and 1 (completed)
     */
    public void showProgress(double p);

    /**
     * Display a progress message.
     *
     * @param msg		message to display
     */
    public void showStatus(String msg);

    /**
     * Denote we are done.
     */
    public void showCompleted();

}

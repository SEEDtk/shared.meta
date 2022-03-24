/**
 *
 */
package org.theseed.shared.meta;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This is a default progress reporter that logs the messages and does nothing with the progress value.
 *
 * @author Bruce Parrello
 *
 */
public class LogProgressReporter implements IProgressReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(LogProgressReporter.class);


    @Override
    public void showProgress(double p) {
    }

    @Override
    public void showStatus(String msg) {
        log.info(msg);
    }

    @Override
    public void showCompleted() {
    }

}

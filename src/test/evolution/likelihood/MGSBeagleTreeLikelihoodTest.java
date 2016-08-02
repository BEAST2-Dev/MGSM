package test.evolution.likelihood;
import java.io.File;


import org.junit.Test;

import beast.core.Logger;
import beast.util.LogAnalyser;
import beast.util.XMLParser;
import junit.framework.TestCase;

public class MGSBeagleTreeLikelihoodTest extends TestCase{

	@Test
	public void testBeagleTreeLikelihood() throws Exception {
        XMLParser parser = new XMLParser();
        beast.core.Runnable runable = parser.parseFile(new File("src/test/0rg.xml"));
        Logger.FILE_MODE = Logger.LogFileMode.overwrite;
        runable.run();

        String logFile = "0rg.log";
        System.out.println("\nAnalysing log " + logFile);
        LogAnalyser logAnalyser = new LogAnalyser(logFile, 0);
        Double [] beagleLikelihood = logAnalyser.getTrace("treeLikelihood.0");
        Double [] javaLikelihood = logAnalyser.getTrace("treeLikelihood.1");
        for (int i = 0; i < beagleLikelihood.length; i++) {
        	assertEquals(beagleLikelihood[i], javaLikelihood[i], 1e-8);
        }

	}
}

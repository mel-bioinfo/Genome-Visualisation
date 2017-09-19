package gv15;

import analysis.PackCreator;
import data.Assembly;
import data.Consensus;
import data.Contig;
import data.IReadManager;
import data.auxiliary.Feature;
import io.AssemblyFile;
import io.AssemblyFileHandler;
import io.TabletFile;
import io.TabletFileHandler;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 *
 * @author ranasi01
 */
public class TabletDataHandler {
    
    private String cachePath;
    private ArrayList<String> loadedReference;
    private IReadManager reads;
    private HashMap<String,ArrayList<Feature>> loadedFeatures;

    public TabletDataHandler(String cachePath) {
        this.cachePath = cachePath;
        loadedFeatures = new HashMap();
    }
    
    //@E read in data form the correct contig within the coordinates
    public void ExtractDataAtCoordinates(String[] fileNames,int startCoordinate,
            int endCoordinate,int contigNumber,String sampleName) throws Exception{
        
              
        File testBam = new File(fileNames[0]);
        File testReference = new File(fileNames[1]);
        
        if (!testBam.exists()) {
            System.out.println("ERROR: The BAM file '" + fileNames[0] + "' was not found. Check file location and address in prefs.txt file");
            System.exit(0);
        }
        
        if (!testReference.exists()) {
            System.out.println("ERROR: The human reference genome file '" + fileNames[1] + "' was not found. Check file location and address in prefs.txt file");
            System.exit(0);
        }
        
        
        TabletFile tabletFile;
        tabletFile = TabletFileHandler.createFromFileList(fileNames); //@E address to the bam and fasta file
        
        AssemblyFile[] files = tabletFile.getFileList(); //@E address to the bam and fasta file
        File cacheDir = new File(cachePath); //@E cache path does not exist
        
        AssemblyFileHandler assemblyFileHandler = new AssemblyFileHandler(files, cacheDir); //@E stores the bam and fasta addresses
        assemblyFileHandler.runJob(0);
        
        TabletFileHandler.addAsMostRecent(tabletFile); //@E address to the bam and fasta file
        Assembly assembly = assemblyFileHandler.getAssembly(); //@E bam file location and cacheID? and bai file

        assembly.getAssemblyStatistics().setAssembly(tabletFile.assembly); //@E bam file location
        assembly.getAssemblyStatistics().setReference(tabletFile.assembly);
        
        //Loading data from contigs
        Contig selectedCotig = assembly.getContig(contigNumber); //@E contig, length. 

        //Set Location
        assembly.getBamBam().setSize(endCoordinate-startCoordinate);
        assembly.getBamBam().setBlockStart(selectedCotig, startCoordinate);

        selectedCotig.clearCigarFeatures();
        assembly.getBamBam().loadDataBlock(selectedCotig);
        assembly.getBamBam().indexNames();

        //Extracting Reference Data
        Consensus consensus = selectedCotig.getConsensus();            
        byte[] referenceData = consensus.getRange(startCoordinate-1, endCoordinate-1);
        loadedReference = new ArrayList();
        
        int counts =0;
        
        for(int i = 0;i<referenceData.length;i++){  // read in bases // time taken between these two breaks
            loadedReference.add(UtilityFunctions.getInstance().GetBaseFromVal(referenceData[i]));
            counts = i+1;
        }
                
        //Sorting Reads
        selectedCotig.getReads().trimToSize();
        Collections.sort(selectedCotig.getReads());
        selectedCotig.calculateOffsets(assembly);
        
        //Packing Reads
        PackCreator packCreator = new PackCreator(selectedCotig, false);
        packCreator.runJob();      
        
        //Packing Reads                    
        reads = selectedCotig.getPackManager();
        
        loadedFeatures.put(sampleName, selectedCotig.getFeatures());
    }
    
    public IReadManager getReads(){
        return reads;
    }
    
    public ArrayList<String> getLoadedReference(){
        return loadedReference;
    }
    
    public HashMap<String,ArrayList<Feature>> getLoadedFeatures(){
        return loadedFeatures;
    }
}

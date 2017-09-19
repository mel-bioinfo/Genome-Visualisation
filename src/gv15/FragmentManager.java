package gv15;

import data.Read;
import data.auxiliary.CigarEvent;
import data.auxiliary.CigarFeature;
import data.auxiliary.CigarInsertEvent;
import data.auxiliary.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;


/**
 *
 * @author ranasi01
 */
public class FragmentManager {
    //DebugInfo red = new DebugInfo();
    
    private int maxReadCount;
    private String dataPath;
    private String cachePath;
    private ReadManager readManager;
    
    public FragmentManager(String dataPath,String cachePath){
        this.dataPath = dataPath;
        this.cachePath = cachePath;
        this.maxReadCount = 0;
    }
    
    public void ProcessFragments(int j, HashMap<String,ArrayList<Phenotype>> phenotypes,ReferenceManager referenceManager,
            PanelManager panelManager,int flank,VariantContext currentVariant) throws Exception{
        //@E j added sent from i
        int startCoord = currentVariant.getStart() - flank;
        int endCoord  = currentVariant.getStart() + flank;
        
        int loopCount =0;
        //DebugInfo loopC = new DebugInfo(); // new object
        
        //Load the Read data for all Samples
        readManager = new ReadManager(referenceManager.getReferencePath(),
                dataPath,cachePath);
        readManager.LoadDataFromSamples(phenotypes, startCoord, endCoord,referenceManager,currentVariant); //@E Load in bam files
        readManager.CreateInsertionArrays(phenotypes,startCoord); // Create the insertion arrays
        
        //Adjust reference sequence for insertions. Add the insertions to this
        referenceManager.AdjustReferences(j, readManager.InsertionArrays); //@E this method adds the ins. to the reference J added
        
        //Create the Panel with Fragments using the extracted Data
        for(String type:phenotypes.keySet()){ // pheotype loop
            //if(type.equals("Neg_Control")){
                 
                int maxReadCountForPhenotype = -1;
                
                int panelSize = referenceManager.AdjustedReferenceData.get(type).size(); //@E the number of insertions at a position
                panelManager.GetPanelFromPhenotype(type).Fragments = new HashMap[panelSize]; 
                
                Map<String,FragmentNode>[] tempFrags = new HashMap[panelSize]; //@E becomes the fragment, before the connections are made
                //Create each panel column
                int totalInsertColumns = 0;
                int panelColumnNo = 0;
                for(int columnNo = 0;columnNo<readManager.InsertionArrays.get(type).length;columnNo++){

                    int insertCount = readManager.InsertionArrays.get(type)[columnNo];
                    int readCountForColumn = 0;
                    //Loop through all the samples for the phenotype
                    for(int sampleNo = 0;sampleNo<phenotypes.get(type).size();sampleNo++){

                        //Extract Reads for sample
                        int sampleReadCount = readManager.GetReadsForSample(
                                phenotypes.get(type).get(sampleNo).FileName).size();
                        readCountForColumn+=sampleReadCount;
                        //Loop through all the reads for the sample
                        for(int readNo = 0;readNo<sampleReadCount;readNo++){ // each base from bam file??
                            gv15.Read currentRead = readManager.GetReadsForSample(
                                phenotypes.get(type).get(sampleNo).FileName).get(readNo);
                            
                            //DebugInfo testBase = new DebugInfo();
                            

                            String[] readBases = currentRead.BaseValues;
                            
                            //Ensure that the read is within the target region
                            int baseIndex = (startCoord - (currentRead.StartPosition+1)) + columnNo; 

                            if(baseIndex >= 0 && baseIndex < currentRead.Length){
                                String baseVal = readBases[baseIndex];
                                //@E
                                //testBase.addText(baseVal);
                                //testBase.createFile(columnNo +""+ readNo +""+sampleNo);
                                
                                
                                FragmentNode tempFragNode = new FragmentNode();  
                                    
                                if(tempFrags[panelColumnNo] == null)
                                    tempFrags[panelColumnNo] = new HashMap();
                                    
                                if(!tempFrags[panelColumnNo].containsKey(baseVal))
                                    tempFrags[panelColumnNo].put(baseVal, tempFragNode);
                                    
                                //Increment the Read Count
                                tempFrags[panelColumnNo].get(baseVal).ReadCount++;

                                //Get the Insert features of the current Read
                                ArrayList<InsertFeature> insertFeatures = readManager.GetInsertsForReadAtPosition(currentRead,
                                    phenotypes.get(type).get(sampleNo).FileName,columnNo,startCoord);
                                
                                //The last fragment does not have any connected fragments
                                if(panelColumnNo < tempFrags.length-1){
                                if(tempFrags[panelColumnNo].get(baseVal).ConnectedFragments == null)
                                    tempFrags[panelColumnNo].get(baseVal).ConnectedFragments = new HashMap();
                               
                                
                                //IF No inserts for this read therefore connect directly with the next Base
                                if(insertFeatures.isEmpty()){
                                    //Add connected Fragments
                                    if( (baseIndex+1) < currentRead.Length){
                                        String nextBaseVal = readBases[baseIndex+1];

                                        if(!tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.containsKey(nextBaseVal))
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.put(nextBaseVal, new HashSet());

                                        int connectionColumn = panelColumnNo+readManager.InsertionArrays.get(type)[columnNo]+1;
                                        if(connectionColumn < tempFrags.length)
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.get(nextBaseVal).
                                                add(connectionColumn);

                                    }                                        
                                }else{
                                    for(InsertFeature insFeature:insertFeatures){
                                        //Add the inserted Fragments
                                        String finalConnectedBase = null;
                                        if( (baseIndex+1) < currentRead.Length)
                                            finalConnectedBase = readBases[baseIndex+1];

                                        if(insFeature.InsertedBases.size() == 6)
                                            System.err.println("");
                                        
                                        AddInsertedBases(tempFrags, panelColumnNo+1, insFeature.InsertedBases, //@E call method at the bottom
                                                finalConnectedBase,panelColumnNo+readManager.InsertionArrays.get(type)[columnNo]+1);
                                        
                                        //Connect the current Fragment to the inserted Fragment
                                        String nextBaseVal = insFeature.InsertedBases.get(0);
                                       
                                        if(!tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.containsKey(nextBaseVal))
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.put(nextBaseVal, new HashSet());

                                        if( (panelColumnNo+1) < tempFrags.length)
                                            tempFrags[panelColumnNo].get(baseVal).ConnectedFragments.get(nextBaseVal).
                                                    add(panelColumnNo+1);
                                    }
                                }
                            }
                            }
                        }//End Read loop
                    }//End sample Loop
                    
                    //if(insertCount > 0)
                        //System.err.println("");
                    panelColumnNo+=(insertCount+1);
                    totalInsertColumns+=insertCount;
                    
                    if(readCountForColumn > maxReadCountForPhenotype)
                        maxReadCountForPhenotype = readCountForColumn;
                    
                }//End Column Loop    
                
                if(maxReadCountForPhenotype > maxReadCount)
                    maxReadCount = maxReadCountForPhenotype;
                
                //Add the fragments to the panel fragments //@E type is the control, CIN3 etc. 
                panelManager.GetPanelFromPhenotype(type).Fragments = tempFrags;//@E fragments passed to the pannel
                panelManager.MaxReadCount = maxReadCount; 


                int adjustedPos = 0; //@E start from the beginning 
                int addedVal = 0;
                for(int i = 0;i<flank+1;i++){ //@E 15ish.  +1 for the target base conversion identification
                    while(referenceManager.AdjustedReferenceData.get(type).get(adjustedPos).equals("INS")){
                        addedVal++; //@E number increments while the adjusted data equals INS
                        adjustedPos++;   //@E experiment with this value to see if the postion of the red line changes
                    }
                    adjustedPos++; //@E advance the position place for every base pased not just the ins
                    //@E the while is never reached again if the the adjusted Pos is not incrmented
                    
                    loopCount = loopCount + 1;
                }
                referenceManager.ShiftVals.put(type, addedVal); //@E stores the RED line value
                
            //}//End Type Check
        }//End Phenotype Loop
        //red.createFile("RED"); // print the file
    }
    
    //@E called from process fragments
    public void AddInsertedBases(Map<String,FragmentNode>[] fragments, int index,
            ArrayList<String> insertedBases,String finalConnectedBase,int finalConnectedColumn){
        
        for(int insertIndex = 0;insertIndex<insertedBases.size();insertIndex++){
        
            //Add the inserted Base
            if(fragments[index+insertIndex] == null)
                fragments[index+insertIndex] = new HashMap();

            FragmentNode tempFragNode = new FragmentNode();  
            String insertedBase = insertedBases.get(insertIndex);
            
            if(!fragments[index+insertIndex].containsKey(insertedBase))
                fragments[index+insertIndex].put(insertedBase, tempFragNode);
                                    
            //Increment the Read Count
            fragments[index+insertIndex].get(insertedBase).ReadCount++;
            
            if(fragments[index+insertIndex].get(insertedBase).ConnectedFragments == null)
                fragments[index+insertIndex].get(insertedBase).ConnectedFragments = new HashMap();

            //Add the next Base as the connected Fragment
            if(insertIndex < insertedBases.size()-1){
                String nextInsertedBase = insertedBases.get(insertIndex+1);
                                                                       
                if(!fragments[index+insertIndex].get(insertedBase).ConnectedFragments.containsKey(nextInsertedBase))
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.put(nextInsertedBase, new HashSet());
                
                if( (index+insertIndex+1) < fragments.length )
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.get(nextInsertedBase).
                        add(index+insertIndex+1);
                
            }else if (insertIndex == insertedBases.size()-1){
                //Add the next non-insert base as the connected fragment of the last Insert
                
                if(!fragments[index+insertIndex].get(insertedBase).ConnectedFragments.containsKey(finalConnectedBase))
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.put(finalConnectedBase, new HashSet());

                if( finalConnectedColumn < fragments.length )
                    fragments[index+insertIndex].get(insertedBase).ConnectedFragments.get(finalConnectedBase).
                        add(finalConnectedColumn);                
            }

        }
                  
    }    

    public void FragmentPrinter(Map<String,FragmentNode>[] fragments){
        System.out.println("Printing Fragments\n");
        
        for(int index = 0;index<fragments.length;index++){
            if(fragments[index]!=null){
                for(String baseVal:fragments[index].keySet()){
                    System.out.print(baseVal + "(" + fragments[index].get(baseVal).ReadCount + ") ");
                    
                    for(String connectedFrag:fragments[index].
                            get(baseVal).ConnectedFragments.keySet()){
                        System.out.print(connectedFrag);
                        for(Object connectedIndex:fragments[index].get(baseVal).ConnectedFragments.get(connectedFrag)){
                            System.out.print("["+connectedIndex+"]");
                        }
                    }
                    
                    System.out.println("");
                }
            }
            System.out.println("");
        }
    }
}

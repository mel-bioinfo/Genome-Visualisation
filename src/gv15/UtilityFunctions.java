package gv15;

import io.CigarParser;
import javafx.scene.paint.Color;

/**
 *
 * @author ranasi01
 */
public class UtilityFunctions {
    
    public int ReadCountRenderThreshold;
    public boolean InsertionsOnlyAtVariant = true;
    public String ReadColour_Unvaried;
    public String ReadColour_Varied;
    public String ReadColour_Insertion;
    
    public CigarParser CigarParser = null;
    
    private static UtilityFunctions instance = new UtilityFunctions();
    
    private UtilityFunctions(){}
        
    public static UtilityFunctions getInstance(){
        return instance;
    }
 
    public String RowNumberToBaseType(int num){  //@E 1 to 5 and return the base on the the graph
        switch(num){
            case 0: return "A";
            case 1: return "C";
            case 2: return "G";
            case 3: return "T";
            case 4: return "N";
        }        
        return "N";
    }
    
    public String GetBaseFromVal(int val){
        switch (val){
            case 6: return "A";
            case 7: return "A";
            case 8: return "C";
            case 9: return "C";
            case 10: return "G";
            case 11: return "G";
            case 12: return "T";
            case 13: return "T";
            
            default: break;
        }
        
        return "N";
    }  
}

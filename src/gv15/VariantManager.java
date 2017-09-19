package gv15;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.util.ArrayList;
/**
 *
 * @author ranasi01
 */
public class VariantManager {
    private VariantContext selectedVariant;
    private ArrayList<VariantContext> variants;
    
    public int TotalVariantCount = 0;
    
    public VariantManager(File vcfFile){
        variants = new ArrayList();
        VCFFileReader reader = new VCFFileReader(vcfFile, false);
        for(VariantContext context:reader.iterator().toList()){ 
            variants.add(context);
            TotalVariantCount++;
        }
    }
    
    public VariantContext getSelectedVariant(){
        return selectedVariant;
    }
    
    public void setVariant(int index){
        selectedVariant = variants.get(index);
    }
}

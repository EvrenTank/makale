package makale;

import java.awt.*;
import java.awt.font.FontRenderContext;

public class deneme2 {



    public static void main(String[] args) {
 liquids liquids = new liquids();
String name[] = {"CH2F2_difluoromethane","C2HF5_pentafluoroethane"};
double w_array[][] = {{0.50,0.50}};


        for(int i = 0;i< w_array.length;i++){
            liquids.weight_fraction_to_mole_fraction(name,w_array[i]);
        }


 //System.out.println( liquids.vis_P("C3H8_propane",273,5.5,20));
        //System.out.println(liquids.k_mix_Flippov(0.210,0.152,0.4,0.6));
       /* System.out.println(liquids.Parachor("C10H22_decane"));
        String names [] = {"C3H8O_propylalcohol","C4H11N_diethylamine"};
        double x1[] = {0.0,1.0};
        double x2[] = {0.1,1-0.1};
        double x3[] = {0.3,1-0.3};
        double x4[] = {0.5,1-0.5};
        double x5[] = {0.7,1-0.7};
        double x6[] = {0.9,1-0.9};
        double x7[] = {1.0,0.0};
        double T = 100; // farketmez*/

  /*      liquids.surten_mix_ZuoandStendby_Kays(names,x1,T);
        liquids.surten_mix_ZuoandStendby(names,x1,T);
        liquids.surten_mix_ZuoandStendby_Kays(names,x2,T);
        liquids.surten_mix_ZuoandStendby(names,x2,T);
        liquids.surten_mix_ZuoandStendby_Kays(names,x3,T);
        liquids.surten_mix_ZuoandStendby(names,x3,T);
        liquids.surten_mix_ZuoandStendby_Kays(names,x4,T);
        liquids.surten_mix_ZuoandStendby(names,x4,T);
        liquids.surten_mix_ZuoandStendby_Kays(names,x5,T);
        liquids.surten_mix_ZuoandStendby(names,x5,T);
        liquids.surten_mix_ZuoandStendby_Kays(names,x6,T);
        liquids.surten_mix_ZuoandStendby(names,x6,T);
        liquids.surten_mix_ZuoandStendby_Kays(names,x7,T);
        liquids.surten_mix_ZuoandStendby(names,x7,T); */






   }
}

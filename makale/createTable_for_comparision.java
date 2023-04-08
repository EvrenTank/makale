package makale;

import javax.swing.*;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;

//  0. T(Celsius)
//  1.T(Kelvin)                      //  2.Pressure(MPa)                 //  3.Density liquid (kg/m^3)      //  4.Specific volume vapor(m^3/kg)      //  5.Liquid enthalpy( kJ/kg)
//  6.Vapor enthalpy(kJ/kg)          //  7.Liquid entropy(kJ/(kgK))      //  8.Vapor entropy(kJ/(kgK))      //  9.Liquid specific heat (kJ/(kgK))    //  10.Vapor specific heat (kJ/(kgK))
// 11.cp/cv   (Birimi yok.)          //  12.Liquid vel. of sound (m/s)	 //  13.Vapor vel. of sound (m/s)   //  14.Liquid  viscosity ( mikroPa.s)    //  15.Vapor viscosity ( mikroPa.s)
// 16.Liquid therm. cond. ( mW/mK))  //17.Vapor therm. cond. ( mW/mK))	 //18.Surface  Tension ( mN/m)      //  19. T( C )

// Yukaridaki birimler tablolarda olan birimleri gosteriyor. Ama metotlard acogunluk olarak kullanilan birimler neyse bu birimler de onlara donusturulecek ve
// hata analizi icin olusturulan tablolarda bu birimler kullanilacak.

public class createTable_for_comparision extends JPanel {

    /* Bu sınıfta özelliklerin hesaplanması için kullanılan farklı yöntemlerden elde edilen
    değerler tablo halinde verilecek ve kıyaslamaları yapılacaktır.
     */
    JTable table;
    JLabel label = new JLabel("");
    /*String liquids[] = { "CH4_methane","C3H8_propane","C3H6_propylene","CF4_carbontetrafluoride","CCl4_carbontetrachloride","C3H6O_acetone","CH4O_methanol","C4H10_butane","C4H10_isobutane","C7H16_heptane","C3H8O3_glycerol",
                         "CHCl3_chloroform","CH3Cl_methylchloride","C6H6_benzene","C2H6_ethane","C2H6O_ethylalcohol","CO2_carbondioxide","C7H8_toluene","C8H18_octane","C9H20_nonane","C10H22_decane","Ar_argon","Br2_bromine",
                         "N2_nitrogen","NH3_ammonia","O2_oxygen","He_helium4",
                         "Hg_mercury","H2O2_hydrogenperoxide","Bi_bismuth","Pb_lead","Na_sodium","K_potassium" };*/
  /*  String properties[] = {"Density (kg/m^3)","Surface tension (N/m)","Thermal conductivity (W/(mK))","Viscosity (Pa.s)",
            "Specific heat (kJ/(kmolK))","h_evaporation (kJ/(kmol))","deltah (kJ/(kmol))","deltas (kJ/(kmolK))","Pvapor (kPa)"};*/
    String properties[] = {"Density (kg/m^3)","Surf. tension (N/m)","Thermal cond. (W/(mK))","Viscosity (Pa.s)",
            "Specific heat (kJ/(kmolK))","Evaporation enthalpy (kJ/(kmol))","deltah (kJ/(kmol))","deltas (kJ/(kmolK))","Pvapor (kPa)"};
    String liquids[] = {"Ar_argon","CH4_methane","C2H2F4_1112tetrafluoroethane","C2H3F3_111trifluoroethane","C2H4F2_11difluoroethane",
            "C2HF5_pentafluoroethane","C2HClF4_2chloro1112tetrafluoroethane","C2HCl2F3_22dichloro111trifluoroethane",
            "C3H8_propane","C3H6_propylene","C4H10_butane","C4H10_isobutane","CCl2F2_dichlorodifluoromethane","CH2F2_difluoromethane",
            "CHClF2_chlorodifluoromethane","CHF3_fluoroform","C2H6_ethane","CO2_carbondioxide","N2_nitrogen","NH3_ammonia",
            "O2_oxygen","He_helium4","H2O_water"};
    JComboBox <String> liquid_list=new JComboBox<String>(liquids);
    JComboBox <String> property_list=new JComboBox<String>(properties);
    String liquid = "Ar_argon";
    String property = "Density (kg/m^3)";
    String column[] ={"T(Kelvin)","Referans","Katsayilar"," % Error","Rackett","% Error","Yamada ve Gunn","% Error","HBT","% Error"};;
    Object row[][];
    JScrollPane scrollPane;

    public createTable_for_comparision(){
        this.setLayout(null);
        liquid_list.setBounds(20,50,300,50);
        property_list.setBounds(350,50,300,50);
        this.add(liquid_list);
        this.add(property_list);
         double x_values[],y_values[];
        row = density_values_for_Table(liquid);
        label.setBounds(35,130,500,50);
        label.setText("Liquid: "+liquid + "  Property: "+property);
        this.add(label);
        table=new JTable(row,column);
        resizeTableColumnWidth();
        table.setPreferredScrollableViewportSize(new Dimension(700,600));
        scrollPane = new JScrollPane(table);
        scrollPane.setBounds(20,200,650,600);
        this.add(scrollPane);
        ActionListener actionListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                createTable_for_comparision.this.remove(label);
                createTable_for_comparision.this.remove(scrollPane);
                liquid = (String) liquid_list.getSelectedItem();
                property = (String) property_list.getSelectedItem();
                label.setText("Liquid: "+liquid + "  Property: "+property);
                label.setBounds(35,130,500,50);
                createTable_for_comparision.this.add(label);
                row = calculate(liquid,property);
                table = new JTable(row,column);
                resizeTableColumnWidth();
                table.setPreferredScrollableViewportSize(new Dimension(700,600));
                scrollPane = new JScrollPane(table);
                scrollPane.setBounds(20,200,650,600);
                createTable_for_comparision.this.add(scrollPane);
                createTable_for_comparision.this.revalidate();
            }
        };
        liquid_list.addActionListener(actionListener);
        property_list.addActionListener(actionListener);
    }
    public Object[][] calculate(String liquid,String property){
        Object row[][] ;

        if(property == "Density (kg/m^3)"){
            row = density_values_for_Table(liquid);
        }
        else if(property == "Surf. tension (N/m)"){
            row = surten_values_for_Table(liquid);
        }
        else if(property == "Thermal cond. (W/(mK))"){
            row = thermalconductivity_values_for_Table(liquid);
        }
        else if(property == "Viscosity (Pa.s)"){
            row = viscosity_values_for_Table(liquid);
        }
        else if(property == "Specific heat (kJ/(kmolK))"){
            row = cp_values_for_Table(liquid);
        }
        else if(property == "Evaporation enthalpy (kJ/(kmol))"){
            row = hbuharlasma_values_for_Table(liquid);
        }
        else if(property == "deltah (kJ/(kmol))"){
            row = deltah_values_for_Table(liquid);
        }
        else if(property == "deltas (kJ/(kmolK))"){
            row = deltas_values_for_Table(liquid);
        }
        else if(property == "Pvapor (kPa)"){
            row = Pvapor_values_for_Table(liquid);
        }
        else{
            row = density_values_for_Table(liquid);
        }

        return  row;
    }
    public void resizeTableColumnWidth(){
     // Loop through all columns
        for (int i = 0; i < table.getColumnCount(); i++) {
            TableColumn column = table.getColumnModel().getColumn(i);

            // Calculate the maximum width of the content in the column
            int maxColumnWidth = 0;
            for (int j = 0; j < table.getRowCount(); j++) {
                TableCellRenderer cellRenderer = table.getCellRenderer(j, i);
                Object value = table.getValueAt(j, i);
                Component component = cellRenderer.getTableCellRendererComponent(table, value, false, false, j, i);
                maxColumnWidth = Math.max(maxColumnWidth, component.getPreferredSize().width);
            }

            // Calculate the maximum width of the header text for the column
            TableCellRenderer headerRenderer = column.getHeaderRenderer();
            if (headerRenderer == null) {
                headerRenderer = table.getTableHeader().getDefaultRenderer();
            }
            Object headerValue = column.getHeaderValue();
            Component headerComponent = headerRenderer.getTableCellRendererComponent(table, headerValue, false, false, 0, i);
            int maxHeaderWidth = headerComponent.getPreferredSize().width;

            // Set the preferred width of the column to the maximum of the content and header widths
            int preferredWidth = Math.max(maxColumnWidth, maxHeaderWidth);
            column.setPreferredWidth(preferredWidth);
        }


    }


    public Object[][] surten_values_for_Table (String name){
        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        double surten_c[]=values.getsurtension(name);
        String metot_names[]= {"T(Kelvin)","Referans","Katsayilar","% Error","Brock ve Bird","% Error","Pitzer","% Error","Zuo ve Stendby",
                               "% Error","Sastri ve Rao","% Error","Macleod ve Sugden","% Error"};
        column = metot_names;
        Object row[][]=new Object[table_Values.length][metot_names.length]; // Tablolara eklenecek olan satırlar
        double sigma; // Yüzey gerilimi Birimi: N/m
        double sigma_referans;
        double percent_error;
        double T;
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        NumberFormat formatter3 = new DecimalFormat("#0.0000000",symbol);
        for(int i=0;i<table_Values.length;i++){
            //row[i][1]=String.format("%,.5f", sigma); Boyle de formatlanabilir. Bircok yontem var.
            T = table_Values[i][1];
            row[i][0]=T;
            sigma = table_Values[i][18]; // Birimi mN/m. Bunu N/m yapmak icin 1000'e bolecegim.
            sigma = sigma/1000; // Birimini degistirdim.
            sigma_referans = sigma;
            row[i][1]=formatter3.format(sigma);
            sigma = liquids.sur_tension(name,T);
            percent_error = (sigma - sigma_referans) / sigma_referans*100;
            row[i][2]=formatter.format(sigma) ;
            row[i][3]=formatter2.format(percent_error);
            sigma = liquids.surten_BrockandBird(name,T);
            percent_error = (sigma-sigma_referans)/sigma_referans*100;
            row[i][4]=formatter.format(sigma) ;
            row[i][5]=formatter2.format(percent_error);
            sigma = liquids.surten_Pitzer(name,T);
            percent_error = (sigma-sigma_referans)/sigma_referans*100;
            row[i][6]=formatter.format(sigma);
            row[i][7]=formatter2.format(percent_error);
            sigma = liquids.surten_ZuoandStendby(name,T);
            percent_error = (sigma-sigma_referans)/sigma_referans*100;
            row[i][8]=formatter.format(sigma) ;
            row[i][9]=formatter2.format(percent_error);
            sigma = liquids.surten_SastriandRao(name,T);
            percent_error = (sigma-sigma_referans)/sigma_referans*100;
            row[i][10]=formatter.format(sigma) ;
            row[i][11]=formatter2.format(percent_error);
            sigma = liquids.surten_MacleodandSugden(name,T,"double");
            percent_error = (sigma-sigma_referans)/sigma_referans*100;
            System.out.println("Macleod and Sugden="+sigma);
            row[i][12]=formatter.format(sigma) ;
            row[i][13]=formatter2.format(percent_error);
        }
        return row;
    }
    public Object[][] viscosity_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        double vis_c[]=values.getvis(name);
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar","% Error","Przezdziecki ve Sridhar","% Error"};
        column = metot_names;
        Object row[][]=new Object[table_Values.length][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        double T;
        double vis; // Pa.s
        double vis_referans;
        double percent_error;
        for(int i=0;i<table_Values.length;i++){
            T = table_Values[i][1];
            row[i][0] = T;
            vis = table_Values[i][14];//Birimi mikroPa.s . Bunu Pa.s yapmak icin 1000000'e bolecegim.
            vis = vis / 1000000; // Birimini cevirdim.
            vis_referans = vis;
            row[i][1]=formatter.format(vis);
            vis = liquids.vis(name,T);
            percent_error = (vis-vis_referans)/vis_referans*100;
            row[i][2]=formatter.format(vis) ;
            row[i][3]=formatter2.format(percent_error);
            vis = liquids.vis_Przezdziecki_and_Sridhar(name,T,"double");
            percent_error = (vis-vis_referans)/vis_referans*100;
            row[i][4]=formatter.format(vis);
            row[i][5]=formatter2.format(percent_error);


        }
        return row;
    }
    public Object[][] cp_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        double cp_c[]=values.getcp(name);
        double critical_values[] = values.get_critical(name);
        double M = critical_values[0];
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar"," % Error","CSP","% Error"};
        column = metot_names;
        Object row[][]=new Object[table_Values.length][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.0000",symbol);
        double T;
        double cp;
        double cp_referans;
        double percent_error;
        for(int i=0;i<table_Values.length;i++){
            T = table_Values[i][1];
            row[i][0] = T;
            cp = table_Values[i][9]; // Birimini degistirecegim. Su an kJ/(kgK) ben bunu kJ/(kmolK) yapacagim.
            cp = cp * M ;
            cp_referans = cp;
            row[i][1]=formatter.format(cp);
            cp = liquids.cp(name,T);
            percent_error = (cp-cp_referans)/cp_referans*100;
            row[i][2]=formatter.format(cp);
            row[i][3]=formatter2.format(percent_error);
            cp = liquids.cp_CSP(name,T,"double");
            percent_error = (cp-cp_referans)/cp_referans*100;
            row[i][4]=formatter.format(cp);
            row[i][5]=formatter2.format(percent_error);
        }
        return row;
    }
    public Object[][] density_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        double cp_c[]=values.getro(name);
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar"," % Error","Rackett","% Error","Yamada ve Gunn","% Error","HBT","% Error"}; // HBT: Hankinson and Thomson
        column = metot_names;
        Object row[][]=new Object[table_Values.length][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        double T;
        double ro; // Pa.s
        double ro_referans;
        double percent_error;
        for(int i=0;i<table_Values.length;i++){
            T = table_Values[i][1];
            row[i][0] = T;
            ro = table_Values[i][3];
            ro_referans = ro;
            row[i][1]=formatter.format(ro);
            ro = liquids.ro(name,T);
            percent_error = (ro-ro_referans)/ro_referans*100;
            // Simdilik katsayi ile hesaplanann degerler referans olarak kullanilacak.
            row[i][2]=formatter.format(ro) ;
            row[i][3]=formatter2.format(percent_error);
            ro = liquids.ro_Rackett(name,T);
            percent_error = (ro-ro_referans)/ro_referans*100;
            row[i][4]=formatter.format(ro);
            row[i][5]=formatter2.format(percent_error);
            ro = liquids.ro_Yamada_Gunn(name,T);
            percent_error = (ro-ro_referans)/ro_referans*100;
            row[i][6]=formatter.format(ro);
            row[i][7]=formatter2.format(percent_error);
            ro = liquids.ro_HBT(name,T);
            percent_error = (ro-ro_referans)/ro_referans*100;
            row[i][8]=formatter.format(ro);
            row[i][9]=formatter2.format(percent_error);
        }
        return row;
    }
    public Object[][] thermalconductivity_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        double cp_c[]=values.getk(name);
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar","% Error","Latini","% Error","Sastri","% Error"};
        column = metot_names;
        Object row[][]=new Object[table_Values.length][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        double T;
        double k;
        double k_referans;
        double percent_error;
        for(int i=0;i<table_Values.length;i++){
            T = table_Values[i][1];
            row[i][0] = T;
            k = table_Values[i][16]; // Birimi mW/(mK). Bunu W/(mK) yapmak icin 1000'e bolecegim.
            k = k/1000;
            k_referans = k;
            row[i][1]=formatter.format(k);
            k = liquids.k(name,T);
            percent_error = (k-k_referans)/k_referans*100;
            row[i][2]=formatter.format(k);
            row[i][3]=formatter2.format(percent_error);
            k = liquids.k_Latini(name,T);
            percent_error = (k-k_referans)/k_referans*100;
            row[i][4]=formatter.format(k);
            row[i][5]=formatter2.format(percent_error);
            k = liquids.k_Sastri(name,T);
            percent_error = (k-k_referans)/k_referans*100;
            row[i][6]=formatter.format(k);
            row[i][7]=formatter2.format(percent_error);
        }
        return row;
    }

    public Object[][] hbuharlasma_values_for_Table(String name) {
        // Degerlerin okundugu tablodaki birim kJ/kg
        // Hata analizi Tablosundaki degerler kJ/kmol biriminden olacak. Onun icin M'yi kullanmam gerekli.
        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar","% Error"};
        column = metot_names;
        double critical_values[] = values.get_critical(name);
        double M = critical_values[0];
        Object row[][]=new Object[table_Values.length][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        double T;
        double hf,hg;
        double hvap; //
        double hvap_referans;
        double percent_error;
        for(int i=0;i<table_Values.length;i++){
            T = table_Values[i][1];
            row[i][0] = T;
            hf = table_Values[i][5];
            hg =  table_Values[i][6];
            hvap = (hg - hf); // Birimi kJ/kg. Ben bunu kJ/kmol yapacagim.
            hvap = hvap * M;
            hvap_referans = hvap;
            row[i][1]=formatter.format(hvap);
            hvap = liquids.hvap(name,T);
            percent_error = (hvap-hvap_referans)/hvap_referans*100;
            row[i][2]=formatter.format(hvap);
            row[i][3]=formatter2.format(percent_error);
        }
        return row;
    }
    public Object[][] deltah_values_for_Table(String name) {
        // Degerlerin okundugu tablodaki birim kJ/kg
        // Hata analizi Tablosundaki degerler kJ/kmol biriminden olacak. Onun icin M'yi kullanmam gerekli.
        //Referans h degerlerinin farkli olmasi h degerlerinin de farkli olmasina neden olacagi icin ben hata analizi direkt olarak deltah
        // yani (h(T2) - h(T1)) uzerinden yapacagim ki href degerlerinin onemi kalmasin. Entropide de ayni sekilde yapacagim.

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar","% Error"};
        column = metot_names;
        double critical_values[] = values.get_critical(name);
        double M = critical_values[0];
        Object row[][]=new Object[table_Values.length-1][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        double T1,T2;
        double h1,h2;
        double deltah; //
        double deltah_referans;
        double percent_error;
        for(int i=0;i<table_Values.length-1;i++){
            T1 = table_Values[i][1];
            T2 = table_Values[i+1][1];
            row[i][0] = (T1 + " - "+ T2);
            h1 = table_Values[i][5];
            h2 =  table_Values[i+1][5];
            deltah = (h2 - h1); // Birimi kJ/kg. Ben bunu kJ/kmol yapacagim.
            deltah = deltah * M;
            deltah_referans = deltah;
            row[i][1]=formatter.format(deltah);
            h1 = liquids.h(name,T1);
            h2 = liquids.h(name,T2);
            deltah = h2-h1;// Bunun birimi zaten kJ/kmol. O yuzden birimine dokunmuyorum.
            percent_error = (deltah-deltah_referans)/deltah_referans*100;
            row[i][2]=formatter.format(deltah);
            row[i][3]=formatter2.format(percent_error);
        }
        return row;
    }
    public Object[][] deltas_values_for_Table(String name) {
        // Degerlerin okundugu tablodaki birim kJ/kg
        // Hata analizi Tablosundaki degerler kJ/kmol biriminden olacak. Onun icin M'yi kullanmam gerekli.
        //Referans s degerlerinin farkli olmasi s degerlerinin de farkli olmasina neden olacagi icin ben hata analizi direkt olarak deltah
        // yani (s(T2) - s(T1)) uzerinden yapacagim ki sref degerlerinin onemi kalmasin. Entalpide de ayni sekilde yapacagim.

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar","% Error"};
        column = metot_names;
        double critical_values[] = values.get_critical(name);
        double M = critical_values[0];
        Object row[][]=new Object[table_Values.length-1][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        double T1,T2;
        double s1,s2;
        double deltas; //
        double deltas_referans;
        double percent_error;
        for(int i=0;i<table_Values.length-1;i++){
            T1 = table_Values[i][1];
            T2 = table_Values[i+1][1];
            row[i][0] = (T1 + " - "+ T2);
            s1 = table_Values[i][7];
            s2 =  table_Values[i+1][7];
            deltas = (s2 - s1); // Birimi kJ/kg. Ben bunu kJ/kmol yapacagim.
            deltas = deltas * M;
            deltas_referans = deltas;
            row[i][1]=formatter.format(deltas);
            s1 = liquids.s(name,T1);
            s2 = liquids.s(name,T2);
            deltas = s2-s1;// Bunun birimi  kJ/(kgK).  Ben bunu kJ/kmol yapacagim.
            deltas = deltas * M;
            percent_error = (deltas-deltas_referans)/deltas_referans*100;
            row[i][2]=formatter.format(deltas);
            row[i][3]=formatter2.format(percent_error);
        }
        return row;
    }

    public Object[][] Pvapor_values_for_Table(String name) { // Pbuhar = Pdoyma
        // Degerlerin okundugu tablodaki birim MPa
        // Hata analizi Tablosundaki degerler kPa biriminden olacak.

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues2(name);
        liquids liquids = new liquids();
        String metot_names[] = {"T(Kelvin)","Referans","Katsayilar","% Error"};
        column = metot_names;
        Object row[][]=new Object[table_Values.length][metot_names.length]; // Tablolara eklenecek olan satırlar
        DecimalFormatSymbols symbol= new DecimalFormatSymbols();
        symbol.setDecimalSeparator('.');
        NumberFormat formatter = new DecimalFormat("#0.0000000",symbol);
        NumberFormat formatter2 = new DecimalFormat("#0.00",symbol);
        double Pvapor;
        double Pvapor_referans;
        double percent_error;
        double T;
        for(int i=0;i<table_Values.length;i++){
            T = table_Values[i][1];
            row[i][0] = (T);
            Pvapor = table_Values[i][2]; // Birimi MPa. 1000 ile çarparak kPa yapacagim.
            Pvapor = Pvapor * 1000;
            Pvapor_referans = Pvapor;
            row[i][1]=formatter.format(Pvapor);
            Pvapor = liquids.Pvapor(name,T,"double");
            percent_error = (Pvapor-Pvapor_referans)/Pvapor_referans*100;
            row[i][2]=formatter.format(Pvapor);
            row[i][3]=formatter2.format(percent_error);
        }
        return row;
    }

    public static void main(String[] args) {
        JFrame frame = new JFrame("Farklı Metotlar ile Elde Edilen Değerlerin Tablo Değerleri ile Karşılaştırılması");
        frame.setSize(800,800);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        createTable_for_comparision panel = new createTable_for_comparision();
        frame.add(panel);
        frame.setVisible(true);
    }
}

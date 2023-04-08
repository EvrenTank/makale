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

//  0.P(MPa)
//  1.T(Kelvin) Bubble               //  2. T(K) Dew                         //  3.Density liquid (kg/m^3)        //  4.Specific volume vapor(m^3/kg)
//  5.Liquid enthalpy( kJ/kmol)      //  6.Vapor enthalpy(kJ/kmol)           //  7.Liquid entropy(kJ/(kmolK))     //  8.Vapor entropy(kJ/(molK))
//  9.Liq. spec. heat (kJ/(kmolK))   //  10.Vapor specific heat (kJ/(kmolK)) //  11.cp/cv   (Birimi yok.)         //  12.Liquid vel. of sound (m/s)
//  13.Vapor vel. of sound (m/s)     //  14.Liquid  viscosity ( mikroPa.s)   //  15.Vapor viscosity ( mikroPa.s)  //  16.Liquid therm. cond. ( mW/mK))
//  17.Vapor therm. cond. ( mW/mK))	 //  18.Surface  Tension ( mN/m)         //  19. P(MPa)

// Yukaridaki birimler tablolarda olan birimleri gosteriyor. Ama metotlard acogunluk olarak kullanilan birimler neyse bu birimler de onlara donusturulecek ve
// hata analizi icin olusturulan tablolarda bu birimler kullanilacak.

public class comparisonTablesforRefMixtures extends JPanel {

    /* Bu sınıfta özelliklerin hesaplanması için kullanılan farklı yöntemlerden elde edilen
    değerler tablo halinde verilecek ve kıyaslamaları yapılacaktır.
     */
    JTable table;
    JLabel label = new JLabel("");
    String properties[] = {"Density (kg/m^3)","Surf. tension (N/m)","Thermal cond. (W/(mK))","Viscosity (Pa.s)","Specific heat (kJ/(kmolK))","Pvapor (kPa)"};
    String liquids[] = {"R404A","R407C","R410A","R507A"};
    JComboBox <String> liquid_list=new JComboBox<String>(liquids);
    JComboBox <String> property_list=new JComboBox<String>(properties);
    String liquid = "R404A";
    String property = "Density (kg/m^3)";
    String column[] ={"T(Kelvin)","Referans","Total mass/total volume"," % Error","Spencer ve Danner","% Error"};
    Object row[][];
    JScrollPane scrollPane;
    String liquids_in_mixture[];
    double mole_numbers[];

    public comparisonTablesforRefMixtures(){
        this.setLayout(null);
        liquid_list.setBounds(20,50,300,50);
        property_list.setBounds(350,50,300,50);
        this.add(liquid_list);
        this.add(property_list);
        double x_values[],y_values[];
        liquids_in_mixture = new String[]{"C2HF5_pentafluoroethane","C2H3F3_111trifluoroethane"};
        mole_numbers = new double[]{0.411839,0.588161};
        row = density_values_for_Table(liquid);
        label.setBounds(35,130,500,50);
        label.setText("Liquid: "+liquid + "  Property: "+property);
        this.add(label);
        table=new JTable(row,column);
        resizeTableColumnWidth();
        table.setPreferredScrollableViewportSize(new Dimension(700,600));
        scrollPane = new JScrollPane(table);
        scrollPane.setBounds(20,200,700,600);
        this.add(scrollPane);
        ActionListener actionListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                comparisonTablesforRefMixtures.this.remove(label);
                comparisonTablesforRefMixtures.this.remove(scrollPane);
                liquid = (String) liquid_list.getSelectedItem();
                property = (String) property_list.getSelectedItem();
                if(liquid.equals("R404A")){
                     liquids_in_mixture = new String[]{"C2HF5_pentafluoroethane","C2H3F3_111trifluoroethane","C2H2F4_1112tetrafluoroethane"};
                     mole_numbers =new double[] {0.358,0.604,0.038};
                }
                else if( liquid.equals("R407C") ){
                   liquids_in_mixture = new String[]{"CH2F2_difluoromethane","C2HF5_pentafluoroethane","C2H2F4_1112tetrafluoroethane"};
                    mole_numbers = new double[]{0.381110,0.179558,0.439332};
                }
                else if( liquid.equals("R410A") ){
                    liquids_in_mixture = new String[]{"CH2F2_difluoromethane","C2HF5_pentafluoroethane"};
                    mole_numbers = new double[]{ 0.697616,0.302384};
                }
                else if(liquid.equals("R507A")) { // 507A icin burasi
                    liquids_in_mixture = new String[]{"C2HF5_pentafluoroethane","C2H3F3_111trifluoroethane"};
                    mole_numbers = new double[]{0.411839,0.588161};
                }
                else { // 507A icin burasi
                    liquids_in_mixture = new String[]{"C2HF5_pentafluoroethane","C2H3F3_111trifluoroethane"};
                    mole_numbers = new double[]{0.411839,0.588161};
                }
                label.setText("Liquid: "+liquid + "  Property: "+property);
                label.setBounds(35,130,500,50);
                comparisonTablesforRefMixtures.this.add(label);
                row = calculate(liquid,property);
                table = new JTable(row,column);
                resizeTableColumnWidth();
                table.setPreferredScrollableViewportSize(new Dimension(700,600));
                scrollPane = new JScrollPane(table);
                scrollPane.setBounds(20,200,700,600);
                comparisonTablesforRefMixtures.this.add(scrollPane);
                comparisonTablesforRefMixtures.this.revalidate();
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
        double table_Values[][] = values.getTableValues(name);
        liquids liquids = new liquids();
        double surten_c[]=values.getsurtension(name);
        String metot_names[]= {"T(Kelvin)","Referans","Hadden","% Error","WK & MS","% Error","WK & HW",
                "% Error","ZS & Kays","% Error","ZS & RiceTeja","% Error"};
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

            try{
                sigma = Double.parseDouble(liquids.surten_mix_Hadden(liquids_in_mixture,mole_numbers,T,-1.0));
                percent_error = (sigma - sigma_referans) / sigma_referans*100;
                row[i][2]=formatter.format(sigma) ;
                row[i][3]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][2]=liquids.surten_mix_Hadden(liquids_in_mixture,mole_numbers,T,-1.0);
                row[i][3]="Hesaplanamiyor";
            }
            try{
                sigma =  Double.parseDouble(liquids.surten_mix_WeinaugKatz_MacleodSugden(liquids_in_mixture,mole_numbers,T));
                percent_error = (sigma-sigma_referans)/sigma_referans*100;
                row[i][4]=formatter.format(sigma) ;
                row[i][5]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                row[i][4]=liquids.surten_mix_WeinaugKatz_MacleodSugden(liquids_in_mixture,mole_numbers,T);
                row[i][5]="Hesaplanamadi";
            }
            try{
                sigma =  Double.parseDouble(liquids.surten_mix_WeinaugKatz_HugillandWelsenes(liquids_in_mixture,mole_numbers,T));
                percent_error = (sigma-sigma_referans)/sigma_referans*100;
                row[i][6]=formatter.format(sigma) ;
                row[i][7]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                row[i][6]=liquids.surten_mix_WeinaugKatz_HugillandWelsenes(liquids_in_mixture,mole_numbers,T);
                row[i][7]="Hesaplanamadi";
            }
            try{
                sigma =  Double.parseDouble(liquids.surten_mix_ZuoandStendby_Kays(liquids_in_mixture,mole_numbers,T));
                percent_error = (sigma-sigma_referans)/sigma_referans*100;
                row[i][8]=formatter.format(sigma) ;
                row[i][9]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                row[i][8]=liquids.surten_mix_ZuoandStendby_Kays(liquids_in_mixture,mole_numbers,T);
                row[i][9]="Hesaplanamadi";
            }
            try{
                sigma =  Double.parseDouble(liquids.surten_mix_ZuoandStendby_RiceTeja(liquids_in_mixture,mole_numbers,T));
                percent_error = (sigma-sigma_referans)/sigma_referans*100;
                row[i][10]=formatter.format(sigma) ;
                row[i][11]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                row[i][10]=liquids.surten_mix_ZuoandStendby_RiceTeja(liquids_in_mixture,mole_numbers,T);
                row[i][11]="Hesaplanamadi";
            }
        }
        return row;
    }
    public Object[][] viscosity_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues(name);
        liquids liquids = new liquids();
        double vis_c[]=values.getvis(name);
        String metot_names[] = {"T(Kelvin)","Referans","Teja ve Rice","% Error"};
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

            try{
                vis = Double.parseDouble(liquids.vis_mix_Teja_and_Rice(liquids_in_mixture,mole_numbers,T,100.0));
                percent_error = (vis-vis_referans)/vis_referans*100;
                row[i][2]=formatter.format(vis) ;
                row[i][3]=formatter2.format(percent_error);

            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][2]=liquids.vis_mix_Teja_and_Rice(liquids_in_mixture,mole_numbers,T,100.0);
                row[i][3]="Hesaplanamadi";
            }
        }
        return row;
    }
    public Object[][] cp_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues(name);
        liquids liquids = new liquids();
        double cp_c[]=values.getcp(name);
        double critical_values[] = values.get_critical(name);
        double M = critical_values[0];
        String metot_names[] = {"T(Kelvin)","Referans","Molar oran"," % Error","Teja","% Error"};
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
            cp_referans = cp;
            row[i][1]=formatter.format(cp);

            try{
                cp = Double.parseDouble(liquids.cp_mix_molarfraction(liquids_in_mixture,mole_numbers,T));
                percent_error = (cp-cp_referans)/cp_referans*100;
                row[i][2]=formatter.format(cp);
                row[i][3]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][2]=liquids.cp_mix_molarfraction(liquids_in_mixture,mole_numbers,T) ;
                row[i][3]="Hesaplanamadi";
            }
            try{
                cp = Double.parseDouble(liquids.cp_mix_Teja(liquids_in_mixture,mole_numbers,T));
                percent_error = (cp-cp_referans)/cp_referans*100;
                row[i][4]=formatter.format(cp);
                row[i][5]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][4]=liquids.cp_mix_Teja(liquids_in_mixture,mole_numbers,T) ;
                row[i][5]="Hesaplanamadi";
            }

        }
        return row;
    }
    public Object[][] density_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues(name);
        liquids liquids = new liquids();
        double cp_c[]=values.getro(name);
        String metot_names[] = {"T(Kelvin)","Referans","Total mass/total volume"," % Error","Spencer ve Danner","% Error"};
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

            try{
                ro = Double.parseDouble(liquids.ro_mix_molar(liquids_in_mixture,mole_numbers,T));
                percent_error = (ro-ro_referans)/ro_referans*100;
                row[i][2]=formatter.format(ro) ;
                row[i][3]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][2]=liquids.ro_mix_molar(liquids_in_mixture,mole_numbers,T) ;
                row[i][3]="Hesaplanamadi";
            }
            try{
                ro = Double.parseDouble(liquids.ro_mix_Spencer_and_Danner(liquids_in_mixture,mole_numbers,T));
                percent_error = (ro-ro_referans)/ro_referans*100;
                row[i][4]=formatter.format(ro) ;
                row[i][5]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][4]=liquids.ro_mix_Spencer_and_Danner(liquids_in_mixture,mole_numbers,T) ;
                row[i][5]="Hesaplanamadi";
            }

        }
        return row;
    }
    public Object[][] thermalconductivity_values_for_Table(String name) {

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues(name);
        liquids liquids = new liquids();
        double cp_c[]=values.getk(name);
        String metot_names[] = {"T(Kelvin)","Referans","Filippov","% Error","Baroncini","% Error","Li","% Error","Power Law","% Error"};
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
            try{
                k = Double.parseDouble(liquids.k_mix_Filippov_x(T,liquids_in_mixture,mole_numbers));
                percent_error = (k-k_referans)/k_referans*100;
                row[i][2]=formatter.format(k);
                row[i][3]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][2]=liquids.k_mix_Filippov_x(T,liquids_in_mixture,mole_numbers) ;
                row[i][3]="Hesaplanamadi";
            }
            try{
                k = Double.parseDouble(liquids.k_mix_Baroncini(liquids_in_mixture,mole_numbers,T));
                percent_error = (k-k_referans)/k_referans*100;
                row[i][4]=formatter.format(k);
                row[i][5]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][4]=liquids.k_mix_Baroncini(liquids_in_mixture,mole_numbers,T) ;
                row[i][5]="Hesaplanamadi";
            }
            try{
                k = Double.parseDouble(liquids.k_mix_Li(liquids_in_mixture,mole_numbers,T));
                percent_error = (k-k_referans)/k_referans*100;
                row[i][6]=formatter.format(k);
                row[i][7]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][6]=liquids.k_mix_Li(liquids_in_mixture,mole_numbers,T) ;
                row[i][7]="Hesaplanamadi";
            }
            try{
                k = Double.parseDouble(liquids.k_mix_PowerLaw(liquids_in_mixture,mole_numbers,T));
                percent_error = (k-k_referans)/k_referans*100;
                row[i][8]=formatter.format(k);
                row[i][9]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][8]=liquids.k_mix_PowerLaw(liquids_in_mixture,mole_numbers,T) ;
                row[i][9]="Hesaplanamadi";
            }
        }
        return row;
    }


    public Object[][] Pvapor_values_for_Table(String name) { // Pbuhar = Pdoyma
        // Degerlerin okundugu tablodaki birim MPa
        // Hata analizi Tablosundaki degerler kPa biriminden olacak.

        liquid_values values = new liquid_values();
        double table_Values[][] = values.getTableValues(name);
        liquids liquids = new liquids();
        String metot_names[] = {"T(Kelvin)","Referans","Aalto","% Error"};
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
            Pvapor = table_Values[i][0]; // Birimi MPa. 1000 ile çarparak kPa yapacagim.
            Pvapor = Pvapor * 1000;
            Pvapor_referans = Pvapor;
            row[i][1]=formatter.format(Pvapor);
            try{
                Pvapor = Double.parseDouble(liquids.Pvapor_mix_Aalto(T,liquids_in_mixture,mole_numbers));
                percent_error = (Pvapor-Pvapor_referans)/Pvapor_referans*100;
                row[i][2]=formatter.format(Pvapor);
                row[i][3]=formatter2.format(percent_error);
            }
            catch (NumberFormatException e){
                e.printStackTrace();
                row[i][2]=liquids.Pvapor_mix_Aalto(T,liquids_in_mixture,mole_numbers);
                row[i][3]="Hesaplanamadi";
            }
        }
        return row;
    }

    public static void main(String[] args) {
        JFrame frame = new JFrame("Farklı Metotlar ile Elde Edilen Değerlerin Tablo Değerleri ile Karşılaştırılması");
        frame.setSize(800,800);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        comparisonTablesforRefMixtures panel = new comparisonTablesforRefMixtures();
        frame.add(panel);
        frame.setVisible(true);
    }
}

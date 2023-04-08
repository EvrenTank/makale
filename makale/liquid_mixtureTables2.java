package makale;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

public class liquid_mixtureTables2 extends JPanel {
    liquids sample=new liquids();
    String liquid_name1="CO";
    String liquid_name2="CO";
    JLabel label_l1,label_l2,label_T,label_P,label_m1,label_m2;// label_l1,label_l2 mean liquid1 and liquid2
    JTable table;
    JLabel l,l2; JComboBox<String> c; JTextField t;
    Object row[][]={{"","","",""}};
    JScrollPane sp;
    JTextField field_T,field_P,field_m1,field_m2;
    String T="300.0";
    String P="100.0";//kPa
    double m1,m2;
    ArrayList<JLabel> labels=new ArrayList<JLabel>();
    ArrayList<JLabel> labels2=new ArrayList<JLabel>();
    ArrayList<JComboBox> comboboxes=new ArrayList<JComboBox>();
    ArrayList<JTextField> textfields=new ArrayList<JTextField>();
    public liquid_mixtureTables2()
    {
        label_T=new JLabel("T(K):");
        label_T.setBounds(20,100,30,30);
        this.add(label_T);
        field_T=new JTextField(""+T,10);
        field_T.setBounds(60,100,50,30);
        Border blackline=BorderFactory.createLineBorder(Color.BLACK);
        field_T.setBorder(blackline);
        this.add(field_T);
        liquid_names l_names= new liquid_names();
        String isimler [] = l_names.get_names();
        this.setLayout(null);
        JPanel panel_compounds=new JPanel();// Bu paneli malzeme sayisina gore olusturulacak olan alanlar icin kullanacagim. GridLayout yaparak
        // direkt olarak elemanlari buna ekleyip, bu paneli de ana panele ekleyecegim.
        String column[]= {"Property"," Value"," Unit","Available temp. range"};
        //String column[]= {"ozellik"," Deger","Birim","Gecerli sicaklik araligi"};
        table=new JTable(row , column);
        resizeTableColumnWidth();
        table.setPreferredScrollableViewportSize(new Dimension(600,600));
        JLabel label=new JLabel("Enter the comp. numbers in the mixture:");
        //JLabel label=new JLabel("Karisimdaki sivi sayisini giriniz:");
        JLabel label_P=new JLabel("P(kPa):");
        field_P = new JTextField(""+P,10);
        label_P.setBounds(120,100,50,30);
        field_P.setBounds(180,100,50,30);
        this.add(label_P);
        this.add(field_P);
        label.setBorder(blackline);
        label.setBounds(240,100,230,30 );
        this.add(label);
        JTextField field=new JTextField(10);
        field.setBounds(480,100,30,30);
        this.add(field);
        sp=new JScrollPane(table);
        sp.setBounds(600,20,600,600);
        liquid_mixtureTables2.this.add(sp);
        JButton hesapla_butonu=new JButton("Calculate");
        hesapla_butonu.setBounds(260,150,140,20);
        this.add(hesapla_butonu);

        ActionListener get_values=new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                T=field_T.getText();
                P=field_P.getText();
                System.out.println(P);
                String liquid_names[] = new String[comboboxes.size()];
                double mole[] = new double[comboboxes.size()];
                for(int i=0;i<comboboxes.size();i++){
                    liquid_names[i] = labels.get(i).getText();
                    mole[i] = Double.parseDouble(textfields.get(i).getText());
                }
                row=sample.calculate_values_for_mixtures(liquid_names,mole,Double.parseDouble(T),Double.parseDouble(P));
                liquid_mixtureTables2.this.remove(sp);
                table=new JTable(row , column);
                resizeTableColumnWidth();
                table.setPreferredScrollableViewportSize(new Dimension(600,600));
                sp=new JScrollPane(table);
                sp.setBounds(600,20,600,600);
                liquid_mixtureTables2.this.add(sp);
                //liquid_mixtureTables2.this.repaint();
                liquid_mixtureTables2.this.revalidate();



            }
        };

        hesapla_butonu.addActionListener(get_values);

        field.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int compound_number=Integer.parseInt(field.getText());
                Component[] components=liquid_mixtureTables2.this.getComponents();

                for(int j=8;j<components.length;j++){
                    liquid_mixtureTables2.this.remove(components[j]);
                }

                labels.clear();
                comboboxes.clear();
                textfields.clear();
                labels.clear();
                int x=20;int y=180;
                for(int i=0;i<compound_number;i++){

                    c=new JComboBox<String>(isimler);
                    c.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            String liquid_name=(String)((JComboBox)e.getSource()).getSelectedItem();
                            for(int i=0; i<comboboxes.size();i++ ){
                                if( e.getSource() == comboboxes.get(i)){
                                    labels.get(i).setText(liquid_name);
                                }
                            }
                        }
                    });
                    c.setBounds(x+105,y+30*i,300,20);
                    comboboxes.add(c);
                    liquid_mixtureTables2.this.add(c);

                    l=new JLabel((String)c.getSelectedItem());
                    l.setBounds(x,y+30*i,100,20);
                    labels.add(l);
                    liquid_mixtureTables2.this.add(l);

                    l2=new JLabel("Mol sayisi:");
                    l2.setBounds(x+410,y+30*i,70,20);
                    labels2.add(l2);
                    liquid_mixtureTables2.this.add(l2);

                    t=new JTextField();
                    t.setText("1");
                    t.setBounds(x+470,y+30*i,45,20);
                    textfields.add(t);
                    liquid_mixtureTables2.this.add(t);

                }

                liquid_mixtureTables2.this.repaint();

                liquid_mixtureTables2.this.revalidate();
            }
        });

        label_P=new JLabel("P:");
        label_m1=new JLabel("m1");
        label_m2=new JLabel("m2");

        field_m1=new JTextField("1",5); // mol sayilarini girmek icin olan alan
        field_m2=new JTextField("1",5); // mol sayilarini girmek icin olan alan

        JComboBox <String> isim_listesi=new JComboBox<String>(isimler);
        JComboBox <String> isim_listesi2=new JComboBox<String>(isimler);
        label_l1=new JLabel(" Label L1 ");
        label_l2=new JLabel("Label L2");
        isim_listesi.setBounds(200, 100, 120, 150);
        isim_listesi2.setBounds(200, 100, 120, 150);

        liquids liquid=new liquids();

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


    public static void main(String[] args) {

        liquid_mixtureTables2 t1=new liquid_mixtureTables2();
        JFrame frame=new JFrame("Mixture icin olan table");
        frame.setSize(1200,800);
        frame.setResizable(false);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(t1);
        frame.setVisible(true);

    }



}

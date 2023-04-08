package makale;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.*;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

public class liquids_Table extends JPanel{
    String liquid_name="CO";
    liquid_names l_names= new liquid_names();
    JLabel label;
    JTable table;
    Object row[][];
    JScrollPane sp;
    JTextField field_T,field_P;
    String T="300.0";
    String P="100.0";//kPa
    public liquids_Table()
    {

        JLabel label_T=new JLabel("T(K):");
        JLabel label_P=new JLabel("P(kPa):");
        field_T=new JTextField(""+T,10);
        field_T.setBackground(this.getBackground());
        field_P=new JTextField(""+P,10);
        field_P.setBackground(this.getBackground());
        liquids liquid=new liquids();

        row=liquid.calculate_values_for_pure("H2O_water",Double.parseDouble(T),Double.valueOf(P));

        String isimler[] = l_names.get_names();
        JComboBox <String> isim_listesi=new JComboBox<String>(isimler);
        isim_listesi.setSelectedItem("H2O_water");
        label=new JLabel(" H2O_water ");
        this.add(label);

        isim_listesi.setBounds(200, 100, 100, 50);
        this.add(isim_listesi);
        this.add(label_T);
        this.add(field_T);
        this.add(label_P);
        this.add(field_P);

        String column[]= {"Property"," Value"," Unit"," Available temp. range"};

        table=new JTable(row,column);
        resizeTableColumnWidth();
        table.setPreferredScrollableViewportSize(new Dimension(600,600));
        sp=new JScrollPane(table);
        this.add(sp);

        ActionListener ac_lis=new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // TODO Auto-generated method stub
                liquids_Table.this.remove(sp);
                T=(field_T.getText());
                P=(field_P.getText());
                liquid_name=(String) isim_listesi.getSelectedItem();
                //label.setText(this.getLayout().getClass().getName());
                label.setText(liquid_name);
                row=liquid.calculate_values_for_pure(liquid_name,Double.parseDouble(T),Double.valueOf(P));
                table=new JTable(row, column);
                resizeTableColumnWidth();
                table.setPreferredScrollableViewportSize(new Dimension(600,600));
                sp=new JScrollPane(table);
                liquids_Table.this.add(sp);
                liquids_Table.this.revalidate();
//		        /*Object comp[]=this.getComponents();
//		        for(int i=0;i<comp.length;i++) {
//		        	System.out.println(comp[i]);
//		        }*/

            }
        };

        field_P.addActionListener(ac_lis);
        field_T.addActionListener(ac_lis);
        isim_listesi.addActionListener(ac_lis);

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

    public static void main(String args[]) {
        JFrame frame=new JFrame(" Properties of Liquids:");
        JTabbedPane tabbedPane=new JTabbedPane(JTabbedPane.TOP);
        liquids_Table panel=new liquids_Table();
        liquid_mixtureTables2 panel2=new liquid_mixtureTables2();
        grafikPanel_last panel3 = new grafikPanel_last();
        createTable_for_comparision panel4 = new createTable_for_comparision();
        tabbedPane.addTab("Pure liquids",panel);
        tabbedPane.addTab("Mixtures",panel2);
        tabbedPane.addTab("Property graphics",panel3);
        tabbedPane.addTab("Error calculations",panel4);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setExtendedState(JFrame.MAXIMIZED_BOTH);
        frame.add(tabbedPane);
        frame.setVisible(true);

    }

}

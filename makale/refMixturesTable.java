package makale;

import javax.swing.*;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class refMixturesTable extends JPanel {
    String liquid_name="CO";
    liquid_names l_names= new liquid_names();
    JLabel label;
    JTable table;
    Object row[][];
    JScrollPane sp;
    JTextField field_T,field_P;
    String T="300.0";
    String P="100.0";//kPa
    public refMixturesTable()
    {

        JLabel label_T=new JLabel("T(K):");
        JLabel label_P=new JLabel("P(kPa):");
        field_T=new JTextField(""+T,10);
        field_T.setBackground(this.getBackground());
        field_P=new JTextField(""+P,10);
        field_P.setBackground(this.getBackground());
        liquids liquid=new liquids();

        row=liquid.calculate_values_for_refMixtures("R404A",Double.parseDouble(T),Double.valueOf(P));

        String isimler[] = {"R404A","R407C","R410A","R507A"};
        JComboBox <String> isim_listesi=new JComboBox<String>(isimler);
        isim_listesi.setSelectedItem("R404A");
        label=new JLabel("R404A");
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
                refMixturesTable.this.remove(sp);
                T=(field_T.getText());
                P=(field_P.getText());
                liquid_name=(String) isim_listesi.getSelectedItem();
                //label.setText(this.getLayout().getClass().getName());
                label.setText(liquid_name);
                row=liquid.calculate_values_for_refMixtures(liquid_name,Double.parseDouble(T),Double.valueOf(P));
                table=new JTable(row, column);
                resizeTableColumnWidth();
                table.setPreferredScrollableViewportSize(new Dimension(600,600));
                sp=new JScrollPane(table);
                refMixturesTable.this.add(sp);
                refMixturesTable.this.revalidate();
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
        refMixturesTable panel1 = new refMixturesTable();
        comparisonTablesforRefMixtures panel2 = new comparisonTablesforRefMixtures();
        tabbedPane.addTab("Refrigerant Mixtures",panel1);
        tabbedPane.addTab("Refrigerant Mixtures Error Analyse",panel2);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setExtendedState(JFrame.MAXIMIZED_BOTH);
        frame.add(tabbedPane);
        frame.setVisible(true);

    }

}

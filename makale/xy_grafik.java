package makale;


import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.AffineTransform;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

    public class xy_grafik extends JPanel implements ActionListener {

        public ArrayList<double[]> xlist=new ArrayList<double[]>();
        public ArrayList<double[]> ylist=new ArrayList<double[]>();
        public String [] curve_names;
        private boolean tagsVisible=false;
        private JButton buton;
        private String xlabel="deneme";
        private String ylabel="deneme";

        double[] xseries,xseries2;
        double[] yseries,yseries2;
        //private NumberFormat numberformat=NumberFormat.getInstance(); // Virgülden sonra kaç rakam
        // gösterileceğini ayarlamam için gerekli
        public xy_grafik(String x_label,String y_label) {
            setLayout(null);
            xlabel = x_label;
            ylabel = y_label;
            this.setBounds(0,0,300,300);
            //this.setSize(300,300);
            buton=new JButton("Veri etiketlerini gizle");
            buton.setBounds(700, 300, 170, 30);
            buton.addActionListener(this);
            buton.setFocusPainted(false);
            add(buton);
            setBackground(Color.gray.brighter());
/*		JFrame frame=new JFrame("XY Grafiği");
		frame.setLayout(null);
		frame.setBounds(0, 0, 900, 700);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.add(this);
		frame.setVisible(true);*/

        }


        public void setValues(double x[],double y[]) {

        xlist.clear();
        ylist.clear();
        xlist.add(x);
        ylist.add(y);
        }
        public void setCurves(String [] names){
            curve_names = names;
        }

        public void setValues(ArrayList<double[]> xlist,ArrayList<double[]> ylist) {
            this.xlist.clear();
            this.ylist.clear();
            this.xlist= xlist;
            this.ylist= ylist;
        }

        public void createGraphic(String name,String graphic_name,double P){
            // name : liquid name'i ifade eder.
            Object o[];
            if(graphic_name.equalsIgnoreCase("viscosity")){
                o = viscosity_values(name,P);
            }
            else if(graphic_name.equalsIgnoreCase("surface tension")){
                o = surten_values(name,P);
            }
            else if(graphic_name.equalsIgnoreCase("thermal conductivity")){
                o= k_values(name,P);
            }
            else if(graphic_name.equalsIgnoreCase("density")){
                o= ro_values(name,P);
            }
            else {
                o = cp_values(name,P);
            }
            ArrayList<double[]> x ;
            ArrayList<double[]> y ;
            x = (ArrayList<double[]>) o[0];
            y = (ArrayList<double[]>) o[1];
            setValues(x,y);
            repaint();
     /*       for(int i=0;i<x.size();i++){
                for(int k=0;k<x.get(i).length;k++){

                }
            }*/
        }
        public Object[] cp_values(String name,double P) {
            liquid_values values = new liquid_values();
            liquids liquids = new liquids();
            double cp_c[]=values.getcp(name);
            ArrayList<double[]> x = new ArrayList<double[]>();
            ArrayList<double[]> y = new ArrayList<double[]>();
            double x_ekseni[]  = new double[20];
            double x_ekseni2[]  = new double[20];
            double y_ekseni[]  = new double[20];
            double y_ekseni2[]  = new double[20];
            double Tmin= cp_c[4];
            double Tmax= cp_c[5];
            double T;
            double cp; // Özgül ısı, Birimi kJ/(kmolK)
            for(int i=0;i<20;i++){
                T = Tmin+(Tmax-Tmin)/19*i;
                x_ekseni[i] = T;
                x_ekseni2[i] = T;
                cp = liquids.cp(name,T);
                y_ekseni[i] = cp;
                try {
                    cp = Double.parseDouble(liquids.cp_CSP(name,T));
                }
                catch (NumberFormatException e){
                    e.printStackTrace();
                    cp = 0.0;
                }
                y_ekseni2[i] = cp;
            }
            x.add(x_ekseni);
            y.add(y_ekseni);
            x.add(x_ekseni2);
            y.add(y_ekseni2);
            String curve_names[] = {"Katsayilar","CSP metodu"};
            setCurves(curve_names);
            Object [] object =new Object[2];
            object[0] = x;
            object[1] = y;
            return object;
        }
        public Object[] ro_values(String name,double P) {
            liquid_values values = new liquid_values();
            liquids liquids = new liquids();
            double ro_c[]=values.getro(name);
            ArrayList<double[]> x = new ArrayList<double[]>();
            ArrayList<double[]> y = new ArrayList<double[]>();
            double x_ekseni[]  = new double[20];
            double y_ekseni[]  = new double[20];
            double x_ekseni2[]  = new double[20];
            double y_ekseni2[]  = new double[20];
            double x_ekseni3[]  = new double[20];
            double y_ekseni3[]  = new double[20];
            double x_ekseni4[]  = new double[20];
            double y_ekseni4[]  = new double[20];
            double x_ekseni5[]  = new double[20];
            double y_ekseni5[]  = new double[20];
            double x_ekseni6[]  = new double[20];
            double y_ekseni6[]  = new double[20];
            double x_ekseni7[]  = new double[20];
            double y_ekseni7[]  = new double[20];
            double Tmin= ro_c[4];
            double Tmax= ro_c[5];
            double T;
            double ro; // Yoğunluk, Birimi kg/m^3
            for(int i=0;i<20;i++){
                T = Tmin+(Tmax-Tmin)/19*i;
                x_ekseni[i] = T;
                x_ekseni2[i] = T;
                x_ekseni3[i] = T;
                x_ekseni4[i] = T;
                x_ekseni5[i] = T;
                x_ekseni6[i] = T;
                ro = liquids.ro(name,T);
                y_ekseni[i] = ro;
                ro = liquids.ro_Rackett(name,T);
                y_ekseni2[i] = ro;
                ro = liquids.ro_Yamada_Gunn(name,T);
                y_ekseni3[i] = ro;
                ro = liquids.ro_HBT(name,T);
                y_ekseni4[i] = ro;
                ro= liquids.ro_Tait(name,T,P);
                y_ekseni5[i] = ro;
                ro = liquids.ro_Chang_and_Zhao(name,T,P);
                y_ekseni6[i] = ro;
            }
            x.add(x_ekseni);
            y.add(y_ekseni);
            x.add(x_ekseni2);
            y.add(y_ekseni2);
            x.add(x_ekseni3);
            y.add(y_ekseni3);
            x.add(x_ekseni4);
            y.add(y_ekseni4);
            x.add(x_ekseni5);
            y.add(y_ekseni5);
            x.add(x_ekseni6);
            y.add(y_ekseni6);
            String curve_names[] = {"Katsayilar","Rackett ","Yamada-Gunn ","HBT ","Tait","Chang ve Zhao"};
            setCurves(curve_names);
            Object [] object =new Object[2];
            object[0] = x;
            object[1] = y;
            return object;
        }
        public Object[] k_values(String name, double P) {
            liquid_values values = new liquid_values();
            liquids liquids = new liquids();
            double k_c[]=values.getk(name);
            ArrayList<double[]> x = new ArrayList<double[]>();
            ArrayList<double[]> y = new ArrayList<double[]>();
            ArrayList<double[]> x_Latini = new ArrayList<double[]>();
            ArrayList<double[]> y_Latini = new ArrayList<double[]>();// Latini yöntemi ile hesaplanan değerler.
            double x_ekseni[]  = new double[20];
            double y_ekseni[]  = new double[20];
            double x_ekseni2[]  = new double[20];
            double y_ekseni2[]  = new double[20];
            double x_ekseni3[]  = new double[20];
            double y_ekseni3[]  = new double[20];
            double x_ekseni4[]  = new double[20];
            double y_ekseni4[]  = new double[20];
            double x_ekseni5[]  = new double[20];
            double y_ekseni5[]  = new double[20];
            double Tmin= k_c[3];
            double Tmax= k_c[4];
            double T;
            double k; // Isıl iletkenlik katsayısı: N/m
            for(int i=0;i<20;i++){
                T = Tmin+(Tmax-Tmin)/19*i;
                x_ekseni[i] = T;
                x_ekseni2[i] = T;
                x_ekseni3[i] = T;
                x_ekseni4[i] = T;
                x_ekseni5[i] = T;
                k = liquids.k(name,T);
                y_ekseni[i] = k;
                k = liquids.k_Latini(name,T);
                y_ekseni2[i] = k;
                k=liquids.k_Sastri(name,T);
                y_ekseni3[i] = k;
                k=liquids.k_Missenard(name,T,P,"double");
                y_ekseni4[i] = k;
                k = liquids.k_Latini_and_Baroncini(name,T,P,"double");
                y_ekseni5[i] = k;
            }
            x.add(x_ekseni);
            x.add(x_ekseni2);
            x.add(x_ekseni3);
            x.add(x_ekseni4);
            x.add(x_ekseni5);
            y.add(y_ekseni);
            y.add(y_ekseni2);
            y.add(y_ekseni3);
            y.add(y_ekseni4);
            y.add(y_ekseni5);

            String curve_names[] = {"Katsayilar","Latini vd.","Sastri","Missenard","Lat. ve Baron."};
            setCurves(curve_names);
/*
            for(int i=0;i<20;i++){
                T = Tmin+(Tmax-Tmin)/19*i;
                x_ekseni[i] = T;
                k = liquids.k_Latini(name,T);
                y_ekseni[i] = k;
            }
            x.add(x_ekseni);
            y.add(y_ekseni);
*/
            Object [] object =new Object[2];
            object[0] = x;
            object[1] = y;

            return object;
        }
        public Object[] surten_values(String name, double P) {
            liquid_values values = new liquid_values();
            liquids liquids = new liquids();
            double surten_c[]=values.getsurtension(name);
            ArrayList<double[]> x = new ArrayList<double[]>();
            ArrayList<double[]> y = new ArrayList<double[]>();
            double x_ekseni[]  = new double[20];
            double x_ekseni2[]  = new double[20];
            double x_ekseni3[]  = new double[20];
            double x_ekseni4[]  = new double[20];
            double x_ekseni5[]  = new double[20];
            double x_ekseni6[]  = new double[20];
            double y_ekseni[]  = new double[20];
            double y_ekseni2[]  = new double[20];
            double y_ekseni3[]  = new double[20];
            double y_ekseni4[]  = new double[20];
            double y_ekseni5[]  = new double[20];
            double y_ekseni6[]  = new double[20];
            double Tmin= surten_c[3];
            double Tmax= surten_c[4];
            double T;
            double sigma; // Yüzey gerilimi Birimi: N/m
            for(int i=0;i<20;i++){
                T = Tmin+(Tmax-Tmin)/19*i;
                x_ekseni[i] = T;
                x_ekseni2[i] = T;
                x_ekseni3[i] = T;
                x_ekseni4[i] = T;
                x_ekseni5[i] = T;
                x_ekseni6[i] = T;
                sigma = liquids.sur_tension(name,T);
                y_ekseni[i] = sigma;
                sigma = liquids.surten_BrockandBird(name,T);
                y_ekseni2[i] = sigma;
                sigma = liquids.surten_Pitzer(name,T);
                y_ekseni3[i] = sigma;
                sigma = liquids.surten_ZuoandStendby(name,T);
                y_ekseni4[i] = sigma;
                sigma = liquids.surten_SastriandRao(name,T);
                y_ekseni5[i] = sigma;
                sigma = liquids.surten_MacleodandSugden(name,T,"double");
                y_ekseni6[i] = sigma;
            }
            x.add(x_ekseni);
            y.add(y_ekseni);
            x.add(x_ekseni2);
            y.add(y_ekseni2);
            x.add(x_ekseni3);
            y.add(y_ekseni3);
            x.add(x_ekseni4);
            y.add(y_ekseni4);
            x.add(x_ekseni5);
            y.add(y_ekseni5);
            x.add(x_ekseni6);
            y.add(y_ekseni6);
            boolean isNan = false;
            /* Bazı sıvılar için sur_tension2 metodunda kullanılan kritik değerler bilinmediği için bu değerler 0 olarak kabul ediliyor.
            Bu da sigma değerinin NaN sonucu vermesine neden oluyor. Bu durumda grafikte sorun çıkacağı için eğer değerler arasında NaN
            varsa onun grafiğini çizdirmeyecek şekilde ayarladım.
            */
/*            for(int j=0;j<y_ekseni2.length;j++){
                if(Double.isNaN(y_ekseni2[j]) ){
                    isNan = true;
                    break;
                }
            }
            if(isNan == false){
                x.add(x_ekseni2);
                y.add(y_ekseni2);
            }

 */
            String curve_names[] = {"Katsayilar","Brock-Bird","Pitzer","Zuo-Stendby","Sastri-Rao","Macl.-Sugden"};
           // String curve_names[] = {"Katsayılar","Brock-Bird","Pitzer","Zuo-Stendby","Sastri-Rao","Macleod and Sugden"};
            setCurves(curve_names);
            Object [] object =new Object[2];
            object[0] = x;
            object[1] = y;
            return object;
        }

        public Object[] viscosity_values(String name,double P) {
            liquid_values values = new liquid_values();
            liquids liquids = new liquids();
            double vis_c[]=values.getvis(name);
            ArrayList<double[]> x = new ArrayList<double[]>();
            ArrayList<double[]> y = new ArrayList<double[]>();
            double x_ekseni[]  = new double[20];
            double y_ekseni[]  = new double[20];
            double x_ekseni2[]  = new double[20];
            double y_ekseni2[]  = new double[20];
            double x_ekseni3[]  = new double[20];
            double y_ekseni3[]  = new double[20];

            double Tmin= vis_c[4];
            double Tmax= vis_c[5]-1.0;
            double T;
            double vis; // Pa.s

            for(int i=0;i<20;i++){
                T = Tmin+(Tmax-Tmin)/19*i;
                x_ekseni[i] = T;
                x_ekseni2[i] = T;
                x_ekseni3[i] = T;
               // x_ekseni2[i] = T;
                vis = liquids.vis(name,T);
                y_ekseni[i] = vis;
                vis = liquids.vis_Przezdziecki_and_Sridhar(name,T,"double");
                y_ekseni2[i] = vis;
                vis = liquids.vis_Lucas(name,T,P,"double");
                y_ekseni3[i] = vis;

            }
            x.add(x_ekseni);
            y.add(y_ekseni);
            x.add(x_ekseni2);
            y.add(y_ekseni2);
            x.add(x_ekseni3);
            y.add(y_ekseni3);
            String curve_names[] = {"Katsayilar","Przez. ve Sridhar","Lucas"};
            setCurves(curve_names);
Object [] object =new Object[2];
object[0] = x;
object[1] = y;
return object;
        }

        public void myPaint(Graphics g){

            //boolean tagsVisible=true;
            JButton buton;

            int alanx=60;
            int alany=40;
            int alanwidth=240;
            int alanheight=200;
            double[] xseries,xseries2;
            double[] yseries,yseries2;
            // y ekseni için kullanıyorum
            NumberFormat numberformat1=NumberFormat.getInstance(); // Virgülden sonra kaç rakam
            numberformat1.setMaximumFractionDigits(5);
            numberformat1.setMinimumFractionDigits(5);

            // x ekseni için kullanıyorum.
            NumberFormat numberformat2=NumberFormat.getInstance(); // Virgülden sonra kaç rakam
            numberformat2.setMaximumFractionDigits(2);
            numberformat2.setMinimumFractionDigits(2);

            Graphics2D g2=(Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setColor(Color.white); // Grafik arka planı rengi
            g2.fillRect(alanx, alany, alanwidth, alanheight);

            //double x_min = 0.0;  // Bunları 0 yaparsam grafik tam bir grafik gibi olur. Ama daha anlaşılır bir grafik için böyle yapılabilir.
            double x_min = 10000;  // Bunları 0 yaparsam grafik tam bir grafik gibi olur. Ama daha anlaşılır bir grafik için böyle yapılabilir.
            double x_max = 0.0;
           //double y_min = 0.0;  // Bunları 0 yaparsam grafik tam bir grafik gibi olur. Ama daha anlaşılır bir grafik için böyle yapılabilir.
            double y_min = 10000;  // Bunları 0 yaparsam grafik tam bir grafik gibi olur. Ama daha anlaşılır bir grafik için böyle yapılabilir.
            double y_max = 0.0;
            double x[];
            double y[];
            Color colors[] = {Color.black,Color.red,Color.blue.brighter(),Color.green.darker(),Color.magenta.darker(),Color.yellow,Color.green};

            for ( int j1=0;j1<xlist.size();j1++){
                double xd[] = xlist.get(j1);
                double yd[] = ylist.get(j1);
                for(int k=0;k<xd.length;k++){
                    if(xd[k] < x_min){ x_min = xd[k];}
                    if(xd[k] > x_max){ x_max = xd[k];}
                    if(yd[k] < y_min){ y_min = yd[k];}
                    if(yd[k] > y_max){ y_max = yd[k];}
                }

            }
            //g2.drawString("8",-(alany+alanheight)/2-140,alanx-50);

            for ( int j=0;j<xlist.size();j++)
            {
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setStroke(new BasicStroke(1.0f));
                 x = xlist.get(j);
                 y = ylist.get(j);



                int boyut=Math.min(x.length,y.length);
                xseries=new double[boyut];
                yseries=new double[boyut];
                xseries2=new double[boyut];
                yseries2=new double[boyut];
                for(int i=0;i<boyut;i++) {
                    xseries[i]=x[i];
                    yseries[i]=y[i];
                    xseries2[i]=x[i];
                    yseries2[i]=y[i];
                }

                g2.setColor(Color.black); // Grafik eğrisinin üzerindeki veri noktalarının rengi
                Arrays.sort(xseries2); // Minimum değeri bulmak içinn sıraladım
                Arrays.sort(yseries2);

                double xrange=x_max-x_min; // En büyük değer - En küçük değer
                double yrange=y_max-y_min; // En büyük değer - En küçük değer


                double katsayix;
                double katsayiy;
                if(xrange == 0.0){
                    katsayix = 0;
                    katsayiy = 0;
                }
                else {
                    katsayix=alanwidth/xrange;
                    katsayiy=alanheight/yrange;
                }

                for(int i=0;i<xseries.length;i++) {
                    g2.fillOval((int)((xseries[i]-x_min)*katsayix+alanx-1),(int)(alany+alanheight-(yseries[i]-y_min)*katsayiy-1), 2,2);


                }

                g2.setColor(Color.red);
                g2.setFont(new Font("arial", Font.PLAIN, 8));
                if(tagsVisible) {

                    for(int i=0;i<xseries.length;i++) {
                        g2.drawString("("+numberformat2.format(xseries[i])+";"+numberformat2.format((yseries[i]))+")"
                                ,(int)((xseries[i]-x_min)*katsayix+alanx+2),(int)(alany+alanheight-(yseries[i]-y_min)*katsayiy));

                    }

                }

                g2.setColor(colors[j]);

                // Veri noktalarını bağlayan çizgiler ( Eğrinin kendisi yani )
                for(int i=0; i<xseries.length-1;i++) {

                    g2.drawLine((int)((xseries[i]-x_min)*katsayix+alanx),(int)(alany+alanheight-(yseries[i]-y_min)*katsayiy),
                            (int)((xseries[i+1]-x_min)*katsayix+alanx),(int)(alany+alanheight-(yseries[i+1]-y_min)*katsayiy));
                }
                g2.setColor(Color.black);
                g2.setStroke(new BasicStroke(1.5f));
                // x çizgisini belirginleştirmek için kullanıldı.
                //g2.drawLine(alanx, (int)(alany+alanheight-(-y_min)*katsayiy),alanx+alanwidth, (int)(alany+alanheight-(-y_min)*katsayiy));
                // y çizgisini belirginleştirmek için kullanıldı.
                //g2.drawLine((int)(-x_min*katsayix+alanx),alany,(int)(-x_min*katsayix+alanx), (alany+alanheight));
                g2.setColor(colors[j]);
                g2.setFont(new Font("arial", Font.PLAIN, 8));

                g2.drawLine(alanx+alanwidth+10,alany+alanheight/2-40+j*20 ,alanx+alanwidth+30,alany+alanheight/2-40+j*20 );
                g2.drawString(curve_names[j],alanx+alanwidth+35,alany+alanheight/2-40+j*20 );

                g2.setColor(Color.black);
                g2.setFont(new Font("arial", Font.PLAIN, 8));
                if( j==0){ // Bu kısmın bir defa çizdirilmesi yeterli olduğu için bu şekilde yaptım.Sadece bir defa çizdirsin diye.
                    for(int i=0;i<=5;i++) {
                        double xdegeri=(double)(xrange/5*i+x_min);
                        String xi= ""+(numberformat2.format(xdegeri)); // x ekseninde yazılacak değer
                        int xi_length = g2.getFontMetrics().stringWidth(xi);
                        g2.drawString(xi, alanx+(i*(alanwidth/5))-xi_length/2, alany+alanheight+15);
                        g2.drawRect(alanx+(i*(alanwidth/5)), alany+alanheight-1, 1, 2);
                        String yi = ""+numberformat1.format((double)(yrange/5*i+y_min)); // y ekseninde yazılacak değer
                        int yi_length = g2.getFontMetrics().stringWidth(yi);
                        g2.drawString(yi, alanx-yi_length-5, alany+alanheight-(i*(alanheight/5))+3);
                        //g2.drawString(yi.replace(".",""), alanx-yi_length-5, alany+alanheight-(i*(alanheight/5))+3); // Dilenirse bu şekilde de kullanılabilir.
                        g2.drawRect(alanx-1,alany+alanheight-(i*(alanheight/5)), 2, 1);
                    }
                }

            }


            g2.setFont(new Font("arial", Font.BOLD, 10));
            g2.drawString(xlabel, (alanx+alanwidth)/2+10, alany+alanheight+30);

            /*AffineTransform at=new AffineTransform();
            at.rotate(-Math.PI/2);
            g2.setTransform(at);
            g2.setFont(new Font("arial", Font.BOLD, 15));
            g2.drawString(ylabel,-(alany+alanheight)/2-140,alanx-50);*/
            int ylabel_length = g2.getFontMetrics().stringWidth(ylabel); // bir stringin pixel cinsinden uzunluğunu verir.
            g2.drawString(ylabel,alanx-ylabel_length/2,alany-15);


        }

        @Override
        public void paint(Graphics g) {
            // TODO Auto-generated method stub
            super.paint(g);

            myPaint(g);

        }

        @Override
        public void actionPerformed(ActionEvent e) {
            // TODO Auto-generated method stub

            if(e.getSource()==buton) {
                tagsVisible = !tagsVisible;
                if(buton.getText()=="Veri etiketlerini gizle") {
                    buton.setText("Veri etiketlerini göster");

                }

                else {
                    buton.setText("Veri etiketlerini gizle");
                }

                repaint();
            }


        }


        public static void main(String[] args) {
            ArrayList<double[]>x=new ArrayList<double[]>();
            ArrayList<double[]>y=new ArrayList<double[]>();

            JFrame frame = new JFrame();

            xy_grafik grafik=new xy_grafik("x label","Viskozite");
            xy_grafik grafik2=new xy_grafik("x label","Yoğunluk");
            grafik.setBackground(Color.green);
            grafik.setBounds(20,20,400,300);
            grafik2.setBackground(Color.pink);
            grafik2.setBounds(450,20,400,300);
            //grafik.createGraphic("CO");
            double x1[]=new double[11];
            double y1[]=new double[11];
            double x2[]=new double[11];
            double y2[]=new double[11];

            for(int i=0;i<11;i++){
                int a = i-5;
                x1[i] = a;
                x2[i] = a;
                y1[i] = a*a*a;
                y2[i] = a*a;
            }
            x.add(x1);
            x.add(x2);
            y.add(y1);
            y.add(y2);
            grafik2.repaint();
            grafik.repaint();

            grafik.setValues(x,y);
            grafik2.setValues(x,y);
            frame.add(grafik);
            frame.add(grafik2);
            frame.setSize(700,700);
            frame.setLayout(null);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        }




    }






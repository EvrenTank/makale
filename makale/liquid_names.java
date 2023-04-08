package makale;

import java.io.File;  // Import the File class
import java.io.FileNotFoundException;  // Import this class to handle errors
import java.util.ArrayList;
import java.util.Scanner; // Import the Scanner class to read text files

public class liquid_names {

    public liquid_names(){

    }

    public  String[] get_names() {
        ArrayList <String> isimler = new ArrayList<String>();
        try {
            File myObj = new File("D:\\Kullanicilar-Lenovo-silme\\eclipse-workspace\\Bitirme Tezi\\src\\bitirme_tezi\\isimler_liste.txt");
            Scanner myReader = new Scanner(myObj);
            while (myReader.hasNextLine()) {
                String data = myReader.nextLine();
                isimler.add(data);
                //System.out.println(data);

            }
            myReader.close();
        } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

        String isim_dizisi [] =new String[isimler.size()] ;
        for (int i=0;i<isimler.size();i++){
            isim_dizisi[i] = isimler.get(i);

        }
       /* for ( int k=0;k<isim_dizisi.length;k++){
            System.out.println(isim_dizisi[k]);
        }*/

        return  isim_dizisi;
    }

    public static void main(String[] args) {
        ArrayList <String> isimler = new ArrayList<String>();
        try {
            File myObj = new File("D:\\Kullanicilar-Lenovo-silme\\eclipse-workspace\\Bitirme Tezi\\src\\bitirme_tezi\\isimler_liste.txt");
            Scanner myReader = new Scanner(myObj);
            while (myReader.hasNextLine()) {
                String data = myReader.nextLine();
                isimler.add(data);
                //System.out.println(data);

            }
            myReader.close();
        } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

        String isim_dizisi [] =new String[isimler.size()] ;
        for (int i=0;i<isimler.size();i++){
            isim_dizisi[i] = isimler.get(i);

        }

        for ( int k=0;k<isim_dizisi.length;k++){
            //System.out.println(isim_dizisi[k]);
        }
        //System.out.println(isimler.size());
        //System.out.println(isim_dizisi.length);

    }
}



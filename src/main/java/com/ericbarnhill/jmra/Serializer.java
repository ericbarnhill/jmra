package com.ericbarnhill.jmra;

import java.io.*;

public class Serializer {

    public static double[][][] loadData(File file) {
        double[][][] data = new double[0][][];
        try {
            FileInputStream fileStream = new FileInputStream(file);
            ObjectInputStream objectStream = new ObjectInputStream(fileStream);
            data = (double[][][])objectStream.readObject();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return data;
    }

    public static void writeData(double[][][] data, File file) {
        try {
            FileOutputStream fileStream = new FileOutputStream(file);
            ObjectOutputStream objectStream = new ObjectOutputStream(fileStream);
            objectStream.writeObject(data);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}

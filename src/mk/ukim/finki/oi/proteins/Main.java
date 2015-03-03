package mk.ukim.finki.oi.proteins;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;

public class Main {

    public static final String DATA_DIR = "data/";
    public static final String FILE_PROTEIN_LINKS = "9606.protein.links.detailed.v9.1.txt";
    public static final String FILE_MAPPING_ENTREZ = "entrez_gene_id.vs.string.v9.05.28122012.txt";
    public static final String FILE_HI2 = "HI-II-14.tsv";
    public static final String FILE_LIT = "Lit-BM-13.tsv";
    public static final String FILE_VENKA = "Venkatesan-09.tsv";
    public static final String FILE_YU = "Yu-11";

    public static void main(String[] args) {
        try {
            long t1 = System.currentTimeMillis();

            // #Entrez_Gene_ID	STRING_Locus_ID
            HashMap<String, Integer> namingMap = new HashMap<String, Integer>();
            HashMap<Integer, String> namingMapInv = new HashMap<Integer, String>();
            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_MAPPING_ENTREZ)));

            br.readLine(); //The column header line
            String line = null;
            while((line = br.readLine()) != null) {
                StringTokenizer tok = new StringTokenizer(line);
                int entrezId = Integer.parseInt(tok.nextToken());
                String locusId = tok.nextToken();

                namingMap.put(locusId, entrezId);
                namingMapInv.put(entrezId, locusId);
            }
            long t2 = System.currentTimeMillis();
            System.out.println(String.format("Map built in %d ms: %d", t2-t1, namingMap.size()));

            processInteractions();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private static void processInteractions() throws IOException {
        long t1 = System.currentTimeMillis();
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_PROTEIN_LINKS)));
        String line = null;
        long idx = 0;
        br.readLine(); //The column header line

        // protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
        HashSet<String> allProteins = new HashSet<String>();

        while((line = br.readLine()) != null) {
            StringTokenizer tok = new StringTokenizer(line);

            String protein1 = tok.nextToken();
            String protein2 = tok.nextToken();
            int neighborhood = Integer.parseInt(tok.nextToken());
            int fusion = Integer.parseInt(tok.nextToken());
            int cooccurence = Integer.parseInt(tok.nextToken());
            int coexpression = Integer.parseInt(tok.nextToken());
            int experimental = Integer.parseInt(tok.nextToken());
            int database = Integer.parseInt(tok.nextToken());
            int textmining = Integer.parseInt(tok.nextToken());
            int combinedScore = Integer.parseInt(tok.nextToken());

            allProteins.add(protein1);
            allProteins.add(protein2);
            idx++;
        }
        long t2 = System.currentTimeMillis();
        System.out.println(String.format("File read in %d ms: %d", t2-t1, idx));
        System.out.println("Number of distinct proteins: " + allProteins.size());
    }
}

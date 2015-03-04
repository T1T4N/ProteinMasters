package mk.ukim.finki.oi.proteins;

import java.io.*;
import java.util.*;
import java.util.concurrent.RunnableFuture;

public class Main {

    public static final String DATA_DIR = "data/";
    public static final String FILE_PROTEIN_LINKS = "9606.protein.links.detailed.v9.1.txt";
    public static final String FILE_MAPPING_ENTREZ = "entrez_gene_id.vs.string.v9.05.28122012.txt";
    public static final String FILE_HI2 = "HI-II-14.tsv";
    public static final String FILE_LIT = "Lit-BM-13.tsv";
    public static final String FILE_VENKA = "Venkatesan-09.tsv";
    public static final String FILE_YU = "Yu-11";

    static HashMap<String, Integer> entrezMap;
    static HashMap<Integer, String> entrezMapInv;
    private static HashMap<String, Integer> proteinIndexMap;
    private static String[] sortedProteins;
    private static HashMap<StringPair, ProteinInteraction> proteinInteractions;

    public static void main(String[] args) {
        try {

            long t1 = System.currentTimeMillis();
            buildEntrezMap();
            long t2 = System.currentTimeMillis();
            System.out.println(String.format("Map built in %d ms: %d", t2 - t1, entrezMap.size()));


            t1 = System.currentTimeMillis();
            processInteractions();
            t2 = System.currentTimeMillis();
            System.out.println(String.format("Protein interactions processed in %d ms: %d", t2-t1, proteinInteractions.size()));
            System.out.println("Number of distinct proteins: " + sortedProteins.length);


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void buildEntrezMap() throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_MAPPING_ENTREZ)));

        // #Entrez_Gene_ID	STRING_Locus_ID
        entrezMap = new HashMap<String, Integer>();
        entrezMapInv = new HashMap<Integer, String>();

        br.readLine(); // The column header line
        String line;
        while((line = br.readLine()) != null) {
            StringTokenizer tok = new StringTokenizer(line);
            int entrezId = Integer.parseInt(tok.nextToken());
            String locusId = tok.nextToken();

            entrezMap.put(locusId, entrezId);
            entrezMapInv.put(entrezId, locusId);
        }
        br.close();
    }
    private static void processInteractions() throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_PROTEIN_LINKS)));
        br.readLine(); // The column header line

        String line;
        // protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
        HashSet<String> allProteins = new HashSet<String>();
        proteinInteractions = new HashMap<StringPair, ProteinInteraction>();

        // int notFound = 0;
        while((line = br.readLine()) != null) {
            StringTokenizer tok = new StringTokenizer(line);

            String protein1 = tok.nextToken();
            String protein2 = tok.nextToken();
            /*
            // Adds + 3 sec because of branch prediction
            if(entrezMap.get(protein1) == null || entrezMap.get(protein2) == null) {
                notFound = notFound + (entrezMap.get(protein1) == null ? 1 : 0);
                notFound = notFound + (entrezMap.get(protein2) == null ? 1 : 0);
                // continue;
            }
            */

            int neighborhood = Integer.parseInt(tok.nextToken());
            int fusion = Integer.parseInt(tok.nextToken());
            int cooccurence = Integer.parseInt(tok.nextToken());
            int coexpression = Integer.parseInt(tok.nextToken());
            int experimental = Integer.parseInt(tok.nextToken());
            int database = Integer.parseInt(tok.nextToken());
            int textmining = Integer.parseInt(tok.nextToken());
            int combinedScore = Integer.parseInt(tok.nextToken());

            ProteinInteraction proteinInteraction = new ProteinInteraction(protein1, protein2, neighborhood, fusion, cooccurence, coexpression, experimental, database, textmining, combinedScore);
            StringPair stringKey = new StringPair(protein1, protein2);

            // Very slow
            proteinInteractions.put(stringKey, proteinInteraction);

            allProteins.add(protein1);
            allProteins.add(protein2);
        }
        br.close();
        proteinIndexMap = new HashMap<String, Integer>(allProteins.size());
        sortedProteins = allProteins.toArray(new String[0]);
        Arrays.sort(sortedProteins);
        for (int i = 0; i < sortedProteins.length; i++) {
            proteinIndexMap.put(sortedProteins[i], i);
        }
    }
}

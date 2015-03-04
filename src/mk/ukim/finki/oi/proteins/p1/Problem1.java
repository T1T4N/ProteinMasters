package mk.ukim.finki.oi.proteins.p1;

import java.io.*;
import java.util.*;

public class Problem1 {

    public static final String DATA_DIR = "data/";
    public static final String FILE_PROTEIN_LINKS = "9606.protein.links.detailed.v9.1.txt";
    public static final String FILE_MAPPING_ENTREZ = "entrez_gene_id.vs.string.v9.05.28122012.txt";
    public static final String FILE_HI2 = "HI-II-14.tsv";
    public static final String FILE_LIT = "Lit-BM-13.tsv";
    public static final String FILE_VENKA = "Venkatesan-09.tsv";
    public static final String FILE_YU = "Yu-11.tsv";

    private HashMap<String, Integer> entrezMap;
    private HashMap<Integer, String> entrezMapInv;
    private HashMap<String, Integer> proteinIndexMap;
    private String[] sortedProteins;
    private HashMap<StringPair, ProteinInteraction> proteinInteractionsMap;
    private HashMap<String, HashSet<String>> proteinAssociations;

    public void execute() {
        try {
            long t1 = System.currentTimeMillis();
            buildEntrezMap();
            long t2 = System.currentTimeMillis();
            System.out.println(String.format("Map built in %d ms: %d", t2 - t1, entrezMap.size()));


            t1 = System.currentTimeMillis();
            processInteractions();
            t2 = System.currentTimeMillis();
            System.out.println(String.format("Protein interactions processed in %d ms: %d", t2 - t1, proteinInteractionsMap.size()));
            System.out.println(String.format("Number of distinct proteins: %d", sortedProteins.length));

            t1 = System.currentTimeMillis();
            processDataSet(1);  // 1 = FILE_HI2
            processDataSet(2);  // 2 = FILE_LIT
            processDataSet(3);  // 3 = FILE_VENKA
            processDataSet(4);  // 4 = FILE_YU
            t2 = System.currentTimeMillis();
            System.out.println(String.format("Data sets processed in %d ms", t2-t1));

            // No need for this loop since there is no combinedScore change if
            // No change was made to the proteinInteraction
            /*
            for (ProteinInteraction pi : proteinInteractionsMap.values()){
                int oldScore = pi.combinedScore;
                int newScore = pi.calculateCombinedScore();
            }
            */

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void processDataSet(int fileNum) throws IOException {
        BufferedReader br = null;
        switch (fileNum) {
            case 1:
                br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_HI2)));
                break;
            case 2:
                br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_LIT)));
                break;
            case 3:
                br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_VENKA)));
                break;
            case 4:
                br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_YU)));
                break;
        }

        assert br != null;
        br.readLine();
        String line;
        while((line = br.readLine()) != null) {
            StringTokenizer tok = new StringTokenizer(line);

            int entrezGeneIdA = -1, entrezGeneIdB = -1;
            // Parse line differently depending on file
            switch (fileNum) {
                case 1:
                    entrezGeneIdA = Integer.parseInt(tok.nextToken());
                    tok.nextToken();    // Symbol A
                    entrezGeneIdB = Integer.parseInt(tok.nextToken());
                    break;
                case 2:
                    entrezGeneIdA = Integer.parseInt(tok.nextToken());
                    tok.nextToken();    // Symbol A
                    entrezGeneIdB = Integer.parseInt(tok.nextToken());
                    break;
                case 3:
                    entrezGeneIdA = Integer.parseInt(tok.nextToken());
                    tok.nextToken();    // DB_CCSB_ORF_ID
                    tok.nextToken();    // DB_Accession
                    entrezGeneIdB = Integer.parseInt(tok.nextToken());
                    break;
                case 4:
                    entrezGeneIdA = Integer.parseInt(tok.nextToken());
                    tok.nextToken();    // Symbol A
                    entrezGeneIdB = Integer.parseInt(tok.nextToken());
                    break;
            }

            String locusIdA = entrezMapInv.get(entrezGeneIdA);
            if(locusIdA == null) continue;  // Skip if no STRING id found in map
            String locusIdB = entrezMapInv.get(entrezGeneIdB);
            if(locusIdB == null) continue;  // Skip if no STRING id found in map

            StringPair keyPair1 = new StringPair(locusIdA, locusIdB);
            StringPair keyPair2 = new StringPair(locusIdB, locusIdA);
            ProteinInteraction pi1 = proteinInteractionsMap.get(keyPair1);
            ProteinInteraction pi2 = proteinInteractionsMap.get(keyPair2);

            // Only do this if the protein interaction exists
            if(pi1 != null && pi2 != null) {
                switch (fileNum) {
                    case 1:
                        pi1.hi_score = 950;
                        pi1.calculateCombinedScore();   // Update score
                        if(pi2 != pi1) {    // Same reference
                            pi2.hi_score = 950;
                            pi2.calculateCombinedScore();
                        }
                        break;
                    case 2:
                        pi1.lit_score = 900;
                        pi1.calculateCombinedScore();
                        if(pi2 != pi1) {    // Same reference
                            pi2.lit_score = 900;
                            pi2.calculateCombinedScore();
                        }
                        break;
                    case 3:
                        pi1.venka_score = 850;
                        pi1.calculateCombinedScore();
                        if(pi2 != pi1) {    // Same reference
                            pi2.venka_score = 850;
                            pi2.calculateCombinedScore();
                        }
                        break;
                    case 4:
                        pi1.yu_score = 850;
                        pi1.calculateCombinedScore();
                        if(pi2 != pi1) {    // Same reference
                            pi2.yu_score = 850;
                            pi2.calculateCombinedScore();
                        }
                        break;
                }
            }
        }
        br.close();
    }

    private void buildEntrezMap() throws IOException {
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
    private void processInteractions() throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(DATA_DIR + FILE_PROTEIN_LINKS)));
        br.readLine(); // The column header line

        String line;
        // protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
        HashSet<String> allProteins = new HashSet<String>(20780);
        proteinInteractionsMap = new HashMap<StringPair, ProteinInteraction>(4850630);
        proteinAssociations = new HashMap<String, HashSet<String>>();

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
            StringPair stringKey1 = new StringPair(protein1, protein2);
            StringPair stringKey2 = new StringPair(protein2, protein1);


            allProteins.add(protein1);
            allProteins.add(protein2);

            HashSet<String> p1 = proteinAssociations.get(protein1);
            if(p1 == null) {
                p1 = new HashSet<String>();
                proteinAssociations.put(protein1, p1);
            }
            p1.add(protein2);

            HashSet<String> p2 = proteinAssociations.get(protein2);
            if(p2 == null) {
                p2 = new HashSet<String>();
                proteinAssociations.put(protein2, p2);
            }
            p2.add(protein1);

            // Very slow
            proteinInteractionsMap.put(stringKey1, proteinInteraction);
            proteinInteractionsMap.put(stringKey2, proteinInteraction);
        }
        br.close();
        proteinIndexMap = new HashMap<String, Integer>(allProteins.size());
        sortedProteins = allProteins.toArray(new String[allProteins.size()]);
        Arrays.sort(sortedProteins);
        for (int i = 0; i < sortedProteins.length; i++) {
            proteinIndexMap.put(sortedProteins[i], i);
        }
    }
}

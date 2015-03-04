package mk.ukim.finki.oi.proteins;

/**
 * Created by robert on 4/3/15.
 */
public class ProteinInteraction {
    String protein1, protein2;
    int neighborhood, fusion, cooccurence,coexpression,experimental,database, textmining,combinedScore;
    private int hashcode;

    public ProteinInteraction(String protein1, String protein2, int neighborhood, int fusion, int cooccurence, int coexpression, int experimental, int database, int textmining, int combinedScore) {
        this.protein1 = protein1;
        this.protein2 = protein2;
        this.neighborhood = neighborhood;
        this.fusion = fusion;
        this.cooccurence = cooccurence;
        this.coexpression = coexpression;
        this.experimental = experimental;
        this.database = database;
        this.textmining = textmining;
        this.combinedScore = combinedScore;

        int result = protein1 != null ? protein1.hashCode() : 0;
        result = 31 * result + (protein2 != null ? protein2.hashCode() : 0);
        result = 31 * result + neighborhood;
        result = 31 * result + fusion;
        result = 31 * result + cooccurence;
        result = 31 * result + coexpression;
        result = 31 * result + experimental;
        result = 31 * result + database;
        result = 31 * result + textmining;
        result = 31 * result + combinedScore;
        hashcode = result;

    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ProteinInteraction that = (ProteinInteraction) o;

        if (coexpression != that.coexpression) return false;
        if (combinedScore != that.combinedScore) return false;
        if (cooccurence != that.cooccurence) return false;
        if (database != that.database) return false;
        if (experimental != that.experimental) return false;
        if (fusion != that.fusion) return false;
        if (neighborhood != that.neighborhood) return false;
        if (textmining != that.textmining) return false;
        if (protein1 != null ? !protein1.equals(that.protein1) : that.protein1 != null) return false;
        if (protein2 != null ? !protein2.equals(that.protein2) : that.protein2 != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return hashcode;
    }
}

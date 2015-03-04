package mk.ukim.finki.oi.proteins.p1;

import java.util.Comparator;
import java.util.Map;

/**
 * Created by Robert on 4/3/15.
 *
 */
public class StringPair {
    public String val1, val2;
    private int hashcode;

    public StringPair(String val1, String val2) {
        this.val1 = val1;
        this.val2 = val2;

        hashcode = val1 != null ? val1.hashCode() : 0;
        hashcode = 31 * hashcode  * (val2 != null ? val2.hashCode() : 0);


    }

    @Override
    public String toString() {
        return "(" + val1 + ", " + val2 + ')';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        StringPair that = (StringPair) o;

        if(!val1.equals(that.val1)) return false;
        if(!val2.equals(that.val2)) return false;
        return true;
    }

    @Override
    public int hashCode() {
        return hashcode;
    }

    public static Comparator<java.util.Map.Entry<StringPair, ProteinInteraction>> nameComparator = new Comparator<Map.Entry<StringPair, ProteinInteraction>>() {
        @Override
        public int compare(Map.Entry<StringPair, ProteinInteraction> o1, Map.Entry<StringPair, ProteinInteraction> o2) {
            int i = -o1.getKey().val1.compareTo(o2.getKey().val1);
            if (i != 0) return i;

            i = -o1.getKey().val2.compareTo(o2.getKey().val2);
            return i;
        }
    };

    public static Comparator<java.util.Map.Entry<StringPair, ProteinInteraction>> scoreComparator = new Comparator<Map.Entry<StringPair, ProteinInteraction>>() {
        @Override
        public int compare(Map.Entry<StringPair, ProteinInteraction> o1, Map.Entry<StringPair, ProteinInteraction> o2) {
            return o2.getValue().combinedScore - o1.getValue().combinedScore;
        }
    };
}

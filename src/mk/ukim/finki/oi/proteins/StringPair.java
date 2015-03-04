package mk.ukim.finki.oi.proteins;

/**
 * Created by robert on 4/3/15.
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
}

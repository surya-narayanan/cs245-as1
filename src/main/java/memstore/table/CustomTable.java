package memstore.table;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import memstore.data.ByteFormat;
import memstore.data.DataLoader;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Custom table implementation to adapt to provided query mix.
 *
 * Current implementation is column store
 */
public class CustomTable implements Table {
    int numCols;
    int numRows;
    HashMap<Integer, IntArrayList> index_col_zero;
    HashMap<Integer, IntArrayList> index_col_one;
    HashMap<Integer, IntArrayList> index_col_two;

    TreeMap<Integer, IntArrayList> index_col_zero_tree;
    TreeMap<Integer, IntArrayList> index_col_one_tree;
    TreeMap<Integer, IntArrayList> index_col_two_tree;
    HashMap<Integer, Long> index_precomputed_sum;

    long zero_column_sum = 0;
    long one_column_sum = 0;
    long two_column_sum = 0;
    long three_column_sum = 0;
    long four_column_sum = 0;

    private TreeMap<Integer, IntArrayList> temp;
    private ByteBuffer rows;
    ByteBuffer columns;
    ByteBuffer row_store_rows;

    public CustomTable() {
    }

    public void load(DataLoader loader) throws IOException {

        this.numCols = loader.getNumCols();
        List<ByteBuffer> rows = loader.getRows();
        numRows = rows.size();

        /* COLUMN TABLE ------------------------------------------------------------------*/
        this.columns = ByteBuffer.allocate(ByteFormat.FIELD_LEN * numRows * (numCols + 2));

        for (int rowId = 0; rowId < numRows; rowId++) {
            ByteBuffer curRow = rows.get(rowId);
            for (int colId = 0; colId < numCols; colId++) {
                int offset = ByteFormat.FIELD_LEN * ((colId * numRows) + rowId);
                this.columns.putInt(offset,
                        curRow.getInt(ByteFormat.FIELD_LEN * colId));
            }
        }

        /* COLUMN TABLE ------------------------------------------------------------------*/

        /* ROW TABLE ------------------------------------------------------------------*/
        this.rows = ByteBuffer.allocate(ByteFormat.FIELD_LEN * numRows * (numCols + 2));

        for (int rowId = 0; rowId < numRows; rowId++) {
            ByteBuffer curRow = rows.get(rowId);
            for (int colId = 0; colId < numCols; colId++) {
                int offset = ByteFormat.FIELD_LEN * ((rowId * (numCols + 2)) + colId);
                int data = curRow.getInt(ByteFormat.FIELD_LEN * colId);
                this.rows.putInt(offset, data);
                //System.out.println("offset" + offset + "data " + data);
            }

            //precomputed sum for all columns except 0th column
            int offset = ByteFormat.FIELD_LEN * ((rowId * (numCols + 2) ) + numCols);
            int data = 0;
            this.rows.putInt(offset, data);

            //precomputed sum for 0th column
            offset = ByteFormat.FIELD_LEN * ((rowId * (numCols + 2) ) + numCols + 1);
            data = 0;
            this.rows.putInt(offset, data);

            //System.out.println("offset" + offset + "data " + data);
        }
        /* ROW TABLE ------------------------------------------------------------------*/


        /* Indexes ------------------------------------------------------------------*/
        // Also creating an index on column zero since it gets used so
        // frequently in updates bench
        index_col_zero = new HashMap<Integer, IntArrayList>();

        for (int rowId = 0; rowId < numRows; rowId++) {

            int elem = getIntField(rowId, 0);

            if (index_col_zero.containsKey(elem)) {
                IntArrayList list = index_col_zero.get(elem);
                list.add(rowId);
            } else {
                IntArrayList newlist = new IntArrayList();
                newlist.add(rowId);
                index_col_zero.put(elem, newlist);
            }
        }

        /* Indexes ------------------------------------------------*/
        index_col_zero_tree = new TreeMap<Integer, IntArrayList>();
        index_col_zero_tree.putAll(index_col_zero);

        index_col_two = new HashMap<Integer, IntArrayList>();

        for (int rowId = 0; rowId < numRows; rowId++) {
            int elem = getIntField(rowId, 2);
            if (index_col_two.containsKey(elem)) {
                IntArrayList list = index_col_two.get(elem);
                list.add(rowId);
            } else {
                IntArrayList newlist = new IntArrayList();
                newlist.add(rowId);
                index_col_two.put(elem, newlist);
            }
        }
        index_col_two_tree = new TreeMap<Integer, IntArrayList>();
        index_col_two_tree.putAll(index_col_two);

        /* Indexes ------------------------------------------------*/

        index_col_one = new HashMap<Integer, IntArrayList>();

        for (int rowId = 0; rowId < numRows; rowId++) {
            int elem = getIntField(rowId, 1);
            if (index_col_one.containsKey(elem)) {
                IntArrayList list = index_col_one.get(elem);
                list.add(rowId);
            } else {
                IntArrayList newlist = new IntArrayList();
                newlist.add(rowId);
                index_col_one.put(elem, newlist);
            }
        }
        index_col_one_tree = new TreeMap<Integer, IntArrayList>();
        index_col_one_tree.putAll(index_col_one);

        /* Indexes ------------------------------------------------*/

        // Also creating an index on column zero since it gets used so
        // frequently in updates bench
        temp = new TreeMap<Integer, IntArrayList>();

        for (int rowId = 0; rowId < numRows; rowId++) {
            int elem = getIntField(rowId, 0);
            if (temp.containsKey(elem)) {
                IntArrayList list = temp.get(elem);
                list.add(rowId);
            } else {
                IntArrayList newlist = new IntArrayList();
                newlist.add(rowId);
                temp.put(elem, newlist);
            }
        }

        index_precomputed_sum = new HashMap<Integer, Long>();


        for (int rowId = 0; rowId < numRows; rowId++) {
            long row_sum = 0;
            for (int colId = 5; colId < numCols; colId++) {
                row_sum = row_sum + getIntField(rowId, colId);

            }
            index_precomputed_sum.put(rowId, row_sum);
        }


        /* Indexes ------------------------------------------------------------------*/

        /* Precomputed sums
         * =================================================Precomputed sums==========================================
         */

        //Every column after column 5
        for (int rowId = 0; rowId < numRows; rowId++) {
            int row_sum = 0;
            for (int colId = 5; colId < numCols; colId++) {
                row_sum = row_sum + getIntField(rowId, colId);
            }

            putIntField(rowId, numCols, row_sum);
            //System.out.println("Should be " + row_sum + "is" + getIntField(rowId, numCols));
        }

        //First 5 columns
        for (int rowId = 0; rowId < numRows; rowId++) {
            int row_sum = 0;
            for (int colId = 0; colId < 5; colId++) {
                row_sum = row_sum + getIntField(rowId, colId);
            }

            putIntField(rowId, numCols + 1, row_sum);
            //System.out.println("Should be " + row_sum + "is" + getIntField(rowId, numCols));
        }

        for (int rowId = 0; rowId < numRows; rowId++) {
            this.zero_column_sum = this.zero_column_sum + getIntField(rowId, 0);
            this.one_column_sum = this.zero_column_sum + getIntField(rowId, 1);
            this.two_column_sum = this.zero_column_sum + getIntField(rowId, 2);
            this.three_column_sum = this.zero_column_sum + getIntField(rowId, 3);
            this.four_column_sum = this.zero_column_sum + getIntField(rowId, 4);
        }
        /* Precomputed sums
         * =================================================Precomputed sums==========================================
         */

    }
    /**
     *
     * Returns the int field at row `rowId` and column `colId`.
     *
     * Column major format :
     *
     * int offset = ByteFormat.FIELD_LEN * ((colId * numRows) + rowId); return
     * column_store_columns.getInt(offset);
     *
     * Really slow method int offset = ByteFormat.FIELD_LEN * ((rowId * numCols)
     * + colId); return row_store_rows.getInt(offset);
     */
    @Override
    public int getIntField(int rowId, int colId) {
        /*
         * int offset = ByteFormat.FIELD_LEN * ((colId * numRows) + rowId);
         * return columns.getInt(offset);
         */
        int offset = ByteFormat.FIELD_LEN * ((rowId * (numCols + 2) ) + colId);
        //System.out.println(offset);
        //System.out.println("colId" + colId + "returned data" + rows.getInt(offset));
        return rows.getInt(offset);

    }

    /**
     * Inserts the passed-in int field at row `rowId` and column `colId`.
     *
     * Working copy int offset = ByteFormat.FIELD_LEN * ((colId * numRows) +
     * rowId); column_store_columns.putInt(offset, field);
     *
     *
     * Really slow method
     *
     * {
     *
     * //Gets the initial value before changes int was = getIntField(rowId,
     * colId); //Makes changes int offset = ByteFormat.FIELD_LEN * ((rowId *
     * numCols) + colId); row_store_rows.putInt(offset, field); //Updates index
     * to reflect chagnes if(colId == 0){
     *
     * IntArrayList old_index_arr = index_col_zero.get(was);
     * old_index_arr.rem(rowId);
     *
     * if(index_col_zero.containsKey(field)){ IntArrayList list =
     * index_col_zero.get(field); list.add(rowId); } else{ IntArrayList newlist
     * = new IntArrayList(); newlist.add(rowId); index_col_zero.put(field,
     * newlist); } }
     *
     * if(colId == 1){
     *
     * IntArrayList old_index_arr = index_col_one.get(was);
     * old_index_arr.rem(rowId);
     *
     * if(index_col_one.containsKey(field)){ IntArrayList list =
     * index_col_one.get(field); list.add(rowId); } else{ IntArrayList newlist =
     * new IntArrayList(); newlist.add(rowId); index_col_one.put(field,
     * newlist); } }
     *
     * if(colId == 2){
     *
     * IntArrayList old_index_arr = index_col_two.get(was);
     * old_index_arr.rem(rowId);
     *
     * if(index_col_two.containsKey(field)){ IntArrayList list =
     * index_col_two.get(field); list.add(rowId); } else{ IntArrayList newlist =
     * new IntArrayList(); newlist.add(rowId); index_col_two.put(field,
     * newlist); } }
     *
     *
     * }
     */

    @Override
    public void putIntField(int rowId, int colId, int field) {

        /*
         * int was = getIntField(rowId, colId);
         *
         * int offset = ByteFormat.FIELD_LEN * ((colId * numRows) + rowId);
         * columns.putInt(offset, field);
         *
         * if(colId == 0){ IntArrayList old_index_arr = index_col_zero.get(was);
         * old_index_arr.rem(rowId);
         *
         * if(index_col_zero.containsKey(field)){ IntArrayList list =
         * index_col_zero.get(field); list.add(rowId); } else{ IntArrayList
         * newlist = new IntArrayList(); newlist.add(rowId);
         * index_col_zero.put(field, newlist); }
         *
         * IntArrayList old_index_arr_tree = index_col_zero_tree.get(was);
         * old_index_arr_tree.rem(rowId);
         *
         * if(index_col_zero_tree.containsKey(field)){ IntArrayList list =
         * index_col_zero_tree.get(field); list.add(rowId); } else{ IntArrayList
         * newlist = new IntArrayList(); newlist.add(rowId);
         * index_col_zero_tree.put(field, newlist); }
         *
         * System.out.println(index_col_zero_tree.equals(index_col_zero) +
         * "comparing whether the tree map equals the hash map");
         *
         * }
         */

        // Gets the initial value before changes
        int was = getIntField(rowId, colId);
        // Makes changes
        int offset = ByteFormat.FIELD_LEN * ((rowId * (numCols + 2) ) + colId);
        rows.putInt(offset, field);

        /*
        if (1 == colId) {
            this.one_column_sum  = this.zero_column_sum - was  + field;
        }

        if (2 == colId) {
            this.two_column_sum  = this.zero_column_sum - was  + field;
        }

        if (3 == colId) {
            this.three_column_sum  = this.zero_column_sum - was  + field;
        }

        if (4 == colId) {
            this.four_column_sum  = this.zero_column_sum - was  + field;
        }
         */
        // Updates index to reflect chagnes


        if (colId < 5){
            int sum_was = getIntField(rowId, numCols + 1);
            offset = ByteFormat.FIELD_LEN * ((rowId * (numCols + 2)) + numCols + 1);
            rows.putInt(offset,sum_was - was + field );

            if (0 == colId) {

                //update precomputed sum

                this.zero_column_sum  = this.zero_column_sum - was  + field;

                //update index
                IntArrayList old_index_arr = temp.get(was);

                old_index_arr.rem(rowId);

                if (temp.containsKey(field)) {
                    IntArrayList list = temp.get(field);
                    list.add(rowId);
                } else {
                    IntArrayList newlist = new IntArrayList();
                    newlist.add(rowId);
                    temp.put(field, newlist);
                }

            }

            if (1 == colId) {

                //update index
                IntArrayList old_index_arr = index_col_one_tree.get(was);
                old_index_arr.rem(rowId);

                if (index_col_one_tree.containsKey(field)) {
                    IntArrayList list = index_col_one_tree.get(field);
                    list.add(rowId);
                } else {
                    IntArrayList newlist = new IntArrayList();
                    newlist.add(rowId);
                    index_col_one_tree.put(field, newlist);
                }
            }

            if (2 == colId) {

                //update index
                IntArrayList old_index_arr = index_col_two_tree.get(was);
                old_index_arr.rem(rowId);

                if (index_col_two_tree.containsKey(field)) {
                    IntArrayList list = index_col_two_tree.get(field);
                    list.add(rowId);
                } else {
                    IntArrayList newlist = new IntArrayList();
                    newlist.add(rowId);
                    index_col_two_tree.put(field, newlist);
                }
            }

        }

    }

    /**
     * Implements the query SELECT SUM(col0) FROM table;
     *
     * Returns the sum of all elements in the first column of the table.
     *
     *
     * Working, slow method
     *
     * { // TODO: Implement this! long sum = 0;
     *
     * for (int rowId = 0; rowId < numRows; rowId++) { sum = sum +
     * getIntField(rowId,0); } // TODO : Done
     *
     * return sum; }
     *
     * Method using indexes { long sum = 0; Set<Integer> keys =
     * index_col_zero.keySet(); for(int elem : keys){ sum = sum +
     * index_col_zero.get(elem).size() * elem ; } return sum; }
     *
     * { long sum = 0; Set<Integer> keys = index_col_zero.keySet(); for(int elem
     * : keys){ IntArrayList arr = index_col_zero.get(elem); long len =
     * arr.size(); sum = sum + (len * elem); } return sum; }
     *
     * 26th Jan 7:47. Using indexes, not sure if right.
     *
     * { long sum = 0;
     *
     * Set<Integer> keys = index_col_zero.keySet();
     *
     * for(int elem : keys){ IntArrayList arr = index_col_zero.get(elem); long
     * len = arr.size(); sum = sum + (len * elem); } return sum; }
     *
     *Chakresh's method
     *
     * long sum = 0;
     *         int indexColumn = 0;
     *         if( indexColumn == 0){
     *             Set<Integer> keys = temp.keySet();
     *             for(int elem : keys){
     *                 IntArrayList arr = temp.get(elem);
     *                 long len = arr.size();
     *                 sum = sum + (len * elem);
     *             }
     *         }
     *
     *         else{
     *             for (int rowId = 0; rowId < numRows; rowId++) {
     *                 sum = sum + getIntField(rowId,0);
     *             }
     *         }
     *         return sum;
     */
    @Override
    public long columnSum() {
        // TODO: Implement this!
		/*long sum = 0;

		for (int rowId = 0; rowId < numRows; rowId++) {
			sum = sum + getIntField(rowId, 0);
		}
		// TODO : Done

		return sum;
*/      return this.zero_column_sum;

    }

    /**
     * Implements the query SELECT SUM(col0) FROM table WHERE col1 > threshold1
     * AND col2 < threshold2;
     *
     * Returns the sum of all elements in the first column of the table, subject
     * to the passed-in predicates.
     */
    @Override
    public long predicatedColumnSum(int threshold1, int threshold2) {

        /**
         * Set<Integer> keys = index_col_one.keySet(); //comparing against
         * threshold 1 for(int elem : keys){ if(elem > threshold1){ // if so,
         * iterate over the row indices for which that holds and compare col2
         * with threshold 2 IntArrayList arr = index_col_one.get(elem); int len
         * = arr.size(); for(int i = 0; i < len; i++){
         * if(getIntField(arr.getInt(i), 2) < threshold2){ sum = sum +
         * getIntField(arr.getInt(i), 0) ; } } } } return sum; }
         *
         *
         * Method using intersection of two key value lists; too slow.
         *
         * long sum = 0;
         *
         * System.out.println(index_col_one.keySet() + "Threshold 1 " +
         * threshold1); System.out.println(index_col_two.keySet() +
         * "Threshold 2 " + threshold2);
         *
         * Set<Integer> one_keys = index_col_one.tailMap(threshold1 +
         * 1).keySet(); Set<Integer> two_keys = index_col_two.headMap(threshold2
         * - 1).keySet();
         *
         * List<Integer> one_values = new ArrayList<>(); List<Integer>
         * two_values = new ArrayList<>();
         *
         * for(int elem : one_keys){ one_values.addAll(index_col_one.get(elem));
         * }
         *
         * for(int elem : two_keys){ two_values.addAll(index_col_two.get(elem));
         * }
         *
         * System.out.println("two_values" + two_values);
         * System.out.println("one_values" + one_values); //get the intersection
         * of the two sets
         *
         *
         * one_values.retainAll(two_values);
         *
         * System.out.println("one_values" + one_values);
         *
         * //comparing against threshold 1 for(int elem : one_values){ sum = sum
         * + getIntField(elem, 0) ; }
         *
         * return sum;
         *
         *
         * Working, too slow
         *
         * long startTime = System.currentTimeMillis();
         *         Set<Integer> one_keys = index_col_one_tree.tailMap(threshold1, false).keySet();
         *         long endTime = System.currentTimeMillis();
         *         System.out.println("One keys: " + (endTime - startTime));
         *
         *         startTime = System.currentTimeMillis();
         *         Set<Integer> two_keys = index_col_two_tree.headMap(threshold2, false).keySet();
         *         endTime = System.currentTimeMillis();
         *         System.out.println("Two keys: " + (endTime - startTime));
         *
         *         Set<Integer> one_values = new TreeSet<>();
         *         Set<Integer> two_values = new TreeSet<>();
         *
         *         startTime = System.currentTimeMillis();
         *         for(int elem : one_keys){
         *             one_values.addAll(index_col_one_tree.get(elem));
         *         }
         *         endTime = System.currentTimeMillis();
         *         System.out.println("One values: " + (endTime - startTime));
         *
         *         startTime = System.currentTimeMillis();
         *         for(int elem : two_keys){
         *             two_values.addAll(index_col_two_tree.get(elem));
         *         }
         *         endTime = System.currentTimeMillis();
         *         System.out.println("Two values: " + (endTime - startTime));
         *
         *         //System.out.println("two_values" + two_values);
         *         //System.out.println("one_values" + one_values);
         *
         *         // get the intersection of the two sets
         *
         *         startTime = System.currentTimeMillis();
         *         one_values.retainAll(two_values);
         *         endTime = System.currentTimeMillis();
         *         System.out.println("Intersecting the two arrays: " + (endTime - startTime));
         *
         *         //System.out.println("one_values" + one_values);
         *
         *         long sum = 0;
         *         //comparing against threshold 1
         *
         *         startTime = System.currentTimeMillis();
         *         for(int elem : one_values){
         *
         *             sum = sum + getIntField(elem, 0) ;
         *
         *         }
         *         endTime = System.currentTimeMillis();
         *         System.out.println("Doing the sum: " + (endTime - startTime));
         *
         *         return sum;
         *
         * Working method
         *
         *
         * { // TODO: Implement this! long sum = 0;
         *
         * for (int rowId = 0; rowId < numRows; rowId++) { if
         * (getIntField(rowId,1) > threshold1) { if (getIntField(rowId,2) <
         * threshold2){ sum = sum + getIntField(rowId,0); } }
         *
         * } // TODO : Done
         *
         * return sum; } Returns the sum of all elements in the rows which pass
         * the predicate.
         *
         * Working iterative method
         * long sum = 0;
         *         for (int rowId = 0; rowId < numRows; rowId++) {
         *             if (getIntField(rowId, 1) > threshold1) {
         *                 if (getIntField(rowId, 2) < threshold2) {
         *                     sum = sum + getIntField(rowId, 0);
         *                 }
         *             }
         *         }
         *
         *       ###Method I am currently working on
         *
         *
         */



        /*
        System.out.println(index_col_one.keySet() + "Threshold 1 " +
                threshold1); System.out.println(index_col_two.keySet() +
                "Threshold 2 " + threshold2);
        */
        long sum = 0;
        for (int rowId = 0; rowId < numRows; rowId++) {
            if (getIntField(rowId, 1) > threshold1) {
                if (getIntField(rowId, 2) < threshold2) {
                    sum = sum + getIntField(rowId, 0);
                }
            }
        }

        return sum;
    }

    /**
     * Implements the query SELECT SUM(col0) + SUM(col1) + ... + SUM(coln) FROM
     * table WHERE col0 > threshold;
     *
     *
     * Buggy method;
     *
     * NavigableMap<Integer, IntArrayList> tailMap =
     * index_col_zero_tree.tailMap(threshold,false); Set<Integer> keys =
     * tailMap.keySet(); for(int key : keys){ // if so, iterate over the row
     * indices for which that holds IntArrayList arr = tailMap.get(key); for(int
     * i = 0; i < arr.size(); i++){ int rowId = arr.getInt(i); for(int colId =
     * 0; colId < numCols; colId++){ sum = sum + getIntField(rowId, colId); } }
     * }
     *
     * Working method: for (int rowId = 0; rowId < numRows; rowId++) { for (int
     * colId = 0; colId < numCols; colId++) { if (getIntField(rowId, 0) >
     * threshold){ sum = sum + getIntField(rowId,colId); } } } Returns the sum
     * of all elements in the rows which pass the predicate.
     *
     * Chakresh's method
     *
     * Set<Integer> keySet = temp.keySet();
     *         for (Integer elem : keySet) {
     *             if(elem > threshold){
     *                 IntArrayList arr = temp.get(elem);
     *
     *                 for (int i = 0; i < arr.size(); i++) {
     *                     for (int colId = 0; colId < numCols; colId++) {
     *                         sum = sum + getIntField(arr.getInt(i), colId);
     *                     }
     *                 }
     *             }
     *
     *         }
     *
     *  Surya Precompute methods
     *
     *   Set<Integer> keySet = temp.keySet();
     *         for (Integer elem : keySet) {
     *             if(elem > threshold){
     *                 IntArrayList arr = temp.get(elem);
     *
     *                 for (int i = 0; i < arr.size(); i++) {
     *                     for (int colId = 0; colId <= 4; colId++) {
     *                         sum = sum + getIntField(arr.getInt(i), colId);
     *                     }
     *                     sum = sum + index_precomputed_sum.get(i);
     *                 }
     *             }
     *
     *         }
     *
     *         return sum + perma_col_sum ;
     *
     *         ** Using tail map
     *
     *         Set<Integer> keySet = temp.tailMap(threshold, false).keySet();
     *
     *         for (Integer elem : keySet) {
     *                 IntArrayList arr = temp.get(elem);
     *
     *                 for (int i = 0; i < arr.size(); i++) {
     *                     for (int colId = 0; colId <= 4; colId++) {
     *                         sum = sum + getIntField(arr.getInt(i), colId);
     *                     }
     *                     sum = sum + getIntField(arr.getInt(i), numCols);
     *                 }
     *         }
     */
    @Override
    public long predicatedAllColumnsSum(int threshold) {

        long sum = 0;

        Set<Integer> keySet = temp.tailMap(threshold, false).keySet();

        for (Integer elem : keySet) {
            IntArrayList arr = temp.get(elem);

            for (int i = 0; i < arr.size(); i++) {
                sum = sum + getIntField(arr.getInt(i), numCols) + getIntField(arr.getInt(i), numCols + 1);
            }
        }

        return sum;


        /*
         * for (int rowId = 0; rowId < numRows; rowId++) { for (int colId = 0;
         * colId < numCols; colId++) { if (getIntField(rowId, 0) > threshold){
         * sum = sum + getIntField(rowId,colId); } } }
         */

        //System.out.println(index_precomputed_sum);

    }

    // Trying the get row sum function

    // Trying the get row sum function

    public void get_row_sum(int rowId) {
        System.out.println("rowId" + rowId);

        int offset = ByteFormat.FIELD_LEN * ((rowId * numCols));
        System.out.println("offset" + offset);
        int length = numCols * ByteFormat.FIELD_LEN;
        System.out.println("length" + length);
        byte[] dest = new byte[length];
        System.out.println("dest init" + dest);
        System.out.println("row_store_rows" + row_store_rows);

        byte[] temp = new byte[length];
        System.out.println("temp init" + temp);
        row_store_rows.get(temp, offset, 1);
        System.out.println("temp fill" + temp);

        row_store_rows.get(dest, offset, length);

        System.out.println("dest fill" + dest);

        System.out.println("ByteBuffer.wrap(dest)");
        System.out.println(ByteBuffer.wrap(dest));
        System.out.println("ByteBuffer.wrap(dest)" + ByteBuffer.wrap(dest));
        System.out.println("ByteBuffer.wrap(dest).order(ByteOrder.BIG_ENDIAN)"
                + ByteBuffer.wrap(dest).order(ByteOrder.BIG_ENDIAN));
        IntBuffer intBuf = ByteBuffer.wrap(dest).order(ByteOrder.BIG_ENDIAN)
                .asIntBuffer();

        System.out.println("intBuf" + intBuf);
        int[] array = new int[intBuf.remaining()];
        System.out.println("array" + IntStream.of(array));
        intBuf.get(array);
        System.out.println("array" + IntStream.of(array));
        System.out.println("sum" + IntStream.of(array).sum());
        // return rows.getInt(offset);
    }

    /**
     * Implements the query UPDATE(col3 = col3 + col2) WHERE col0 < threshold;
     *
     *
     * Returns the number of rows updated.
     *
     * Slow method using index { int count = 0; Set<Integer> keys =
     * index_col_zero.keySet(); //comparing against threshold for(int elem :
     * keys){ if(elem < threshold){ // if so, iterate over the row indices for
     * which that holds IntArrayList arr = index_col_zero.get(elem); for(int i =
     * 0; i < arr.size(); i++){ int rowId = arr.getInt(i); int field =
     * getIntField(rowId, 3) + getIntField(rowId, 2); putIntField(rowId, 3,
     * field); count = count + 1; } } } return count; }
     *
     * Working method: {
     *
     * int count = 0; for (int rowId = 0; rowId < numRows; rowId++) {
     *
     * if(getIntField(rowId, 0) < threshold){ int field = getIntField(rowId, 3)
     * + getIntField(rowId, 2); putIntField(rowId, 3, field); count = count + 1;
     * } } return count; }
     */
    @Override
    public int predicatedUpdate(int threshold)

    {
        int count = 0;
        for (int rowId = 0; rowId < numRows; rowId++) {

            if (getIntField(rowId, 0) < threshold) {
                int field = getIntField(rowId, 3) + getIntField(rowId, 2);
                putIntField(rowId, 3, field);
                count = count + 1;
            }
        }
        return count;
    }
}

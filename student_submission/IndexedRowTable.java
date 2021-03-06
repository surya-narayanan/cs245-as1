package memstore.table;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import memstore.data.ByteFormat;
import memstore.data.DataLoader;
import sun.jvm.hotspot.utilities.IntArray;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.*;

/**
 * IndexedRowTable, which stores data in row-major format.
 * That is, data is laid out like
 *   row 1 | row 2 | ... | row n.
 *
 * Also has a tree index on column `indexColumn`, which points
 * to all row indices with the given value.
 */
public class IndexedRowTable implements Table {

    int numCols;
    int numRows;
    private TreeMap<Integer, IntArrayList> index;
    private TreeMap<Integer, IntArrayList> index_col_zero;
    private ByteBuffer rows;
    private int indexColumn;

    public IndexedRowTable(int indexColumn) {
        this.indexColumn = indexColumn;
    }

    /**
     * Loads data into the table through passed-in data loader. Is not timed.
     *
     * @param loader Loader to load data from.
     * @throws IOException
     */
    @Override
    public void load(DataLoader loader) throws IOException{
        /*
        Loading the data
        ================================================================================================================
        ================================================================================================================
        */
        this.numCols = loader.getNumCols();
        List<ByteBuffer> rows = loader.getRows();
        numRows = rows.size();
        this.rows = ByteBuffer.allocate(ByteFormat.FIELD_LEN * numRows * numCols);

        for (int rowId = 0; rowId < numRows; rowId++) {
            ByteBuffer curRow = rows.get(rowId);
            for (int colId = 0; colId < numCols; colId++) {
                int offset = ByteFormat.FIELD_LEN * ((rowId * numCols) + colId);
                this.rows.putInt(offset, curRow.getInt(ByteFormat.FIELD_LEN * colId));
            }
        }
        /*
        ================================================================================================================
        ================================================================================================================
        End Loading the Data
        */

        /*
        Creating an index
        ========================================================
        ========================================================
        */
        index = new TreeMap<Integer, IntArrayList>();

        for (int rowId = 0; rowId < numRows; rowId++){
            int elem = getIntField(rowId, indexColumn);
            if(index.containsKey(elem)){
                IntArrayList list = index.get(elem);
                list.add(rowId);
            } else{
                IntArrayList newlist = new IntArrayList();
                newlist.add(rowId);
                index.put(elem, newlist);
            }
        }

        // Also creating an index on column zero since it gets used so frequently in updates bench
        index_col_zero = new TreeMap<Integer, IntArrayList>();

        for (int rowId = 0; rowId < numRows; rowId++){
            int elem = getIntField(rowId, 0);
            if(index_col_zero.containsKey(elem)){
                IntArrayList list = index_col_zero.get(elem);
                list.add(rowId);
            } else{
                IntArrayList newlist = new IntArrayList();
                newlist.add(rowId);
                index_col_zero.put(elem, newlist);
            }
        }
        /*
        ========================================================
        ========================================================
        Creating the index
        */

    }

    /**
     * Returns the int field at row `rowId` and column `colId`.
     */
    @Override
    public int getIntField(int rowId, int colId) {
        int offset = ByteFormat.FIELD_LEN * ((rowId * numCols) + colId);
        return rows.getInt(offset);
    }

    /**
     * Inserts the passed-in int field at row `rowId` and column `colId`.
     */
    @Override

    public void putIntField(int rowId, int colId, int field) {

        //Gets the initial value before changes
        int was = getIntField(rowId, colId);
        //Makes changes
        int offset = ByteFormat.FIELD_LEN * ((rowId * numCols) + colId);
        rows.putInt(offset, field);
        //Updates index to reflect chagnes
        if(indexColumn == colId){

            IntArrayList old_index_arr = index.get(was);
            old_index_arr.rem(rowId);

            if(index.containsKey(field)){
                IntArrayList list = index.get(field);
                list.add(rowId);
            } else{
                IntArrayList newlist = new IntArrayList();
                newlist.add(rowId);
                index.put(field, newlist);
            }
        }
    }

    /**
     * Implements the query
     *  SELECT SUM(col0) FROM table;
     *
     *  Returns the sum of all elements in the first column of the table.
     */
    @Override
    public long columnSum() {
        long sum = 0;

        if( indexColumn == 0){
            Set<Integer> keys = index.keySet();
            for(int elem : keys){
                IntArrayList arr = index.get(elem);
                long len = arr.size();
                sum = sum + (len * elem);
            }
        }

        else{
            for (int rowId = 0; rowId < numRows; rowId++) {
                sum = sum + getIntField(rowId,0);
            }
        }
        return sum;

    }


    /**
     * Implements the query
     *  SELECT SUM(col0) FROM table WHERE col1 > threshold1 AND col2 < threshold2;
     *
     *  Returns the sum of all elements in the first column of the table,
     *  subject to the passed-in predicates.
     */
    @Override
    public long predicatedColumnSum(int threshold1, int threshold2) {

        //updateindex(); //works but is too slow

        long sum = 0;

        // Making use of the int array list
        if (indexColumn == 1){
             /*
            / Checking if any of the values in the column are less than threshold
          */
            //getting the keys of the index to compare against threshold 1
            Set<Integer> keys = index.keySet();
            //comparing against threshold 1
            for(int elem : keys){
                if(elem > threshold1){
                    // if so, iterate over the row indices for which that holds and compare col2 with threshold 2
                    IntArrayList arr = index.get(elem);
                    int len = arr.size();
                    for(int i = 0; i < len; i++){
                        if(getIntField(arr.getInt(i), 2) < threshold2){
                            sum  = sum + getIntField(arr.getInt(i), 0) ;
                        }
                    }
                }
            }
        } else if (indexColumn == 2){
             /*
            / Checking if any of the values in the column are less than threshold
          */
            //getting the keys of the index to compare against threshold 1
            Set<Integer> keys = index.keySet();
            //comparing against threshold 1
            for(int elem : keys){
                if(elem < threshold2){
                    // if so, iterate over the row indices for which that holds and compare col2 with threshold 2
                    IntArrayList arr = index.get(elem);
                    int len  = arr.size();
                    for(int i = 0; i < len; i++){
                        if(getIntField(arr.getInt(i), 1) > threshold1){
                            sum  = sum + getIntField(arr.getInt(i), 0) ;
                        }

                    }
                }
            }
        } else  {
                for (int rowId = 0; rowId < numRows; rowId++) {
                    if (getIntField(rowId,1) > threshold1)
                    {
                        if (getIntField(rowId,2) < threshold2){
                            sum = sum + getIntField(rowId,0);
                        }
                    }
                }
        }

        return sum;
    }

    /**
     * Implements the query
     *  SELECT SUM(col0) + SUM(col1) + ... + SUM(coln) FROM table WHERE col0 > threshold;
     *
     *  Returns the sum of all elements in the rows which pass the predicate.
     */
    @Override
    public long predicatedAllColumnsSum(int threshold) {
        // TODO: Implement this!

        long sum = 0;
        for (int rowId = 0; rowId < numRows; rowId++) {
            for (int colId = 0; colId < numCols; colId++) {
                if (getIntField(rowId, 0) > threshold){
                    sum = sum + getIntField(rowId,colId);
                }
            }
        }
        // TODO : Done
        return sum;
    }

    /**
     * Implements the query
     *   UPDATE(col3 = col3 + col2) WHERE col0 < threshold;
     *
     *   Returns the number of rows updated.
     */
    @Override
     public int predicatedUpdate(int threshold) {
        int count = 0;
        Set<Integer> keys = index_col_zero.keySet();
        //comparing against threshold
        for(int elem : keys){
            if(elem < threshold){
                // if so, iterate over the row indices for which that holds
                IntArrayList arr = index_col_zero.get(elem);
                for(int i = 0; i < arr.size(); i++){
                    int rowId = arr.getInt(i);
                    int field = getIntField(rowId, 3) + getIntField(rowId, 2);
                    putIntField(rowId, 3, field);
                    count = count + 1;
                }
            }
        }
        return count;
     }

}

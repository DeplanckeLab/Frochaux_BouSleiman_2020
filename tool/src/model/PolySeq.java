package model;

public class PolySeq 
{
	public static final String[] POLY_DEL = {"", "-", "--", "---", "----", "-----", "------", "-------", "--------", "---------", "----------", "-----------", "------------", "-------------", "--------------", "---------------"};
	public static final String[] POLY_N = {"", "N", "NN", "NNN", "NNNN", "NNNNN", "NNNNNN", "NNNNNNN", "NNNNNNNN", "NNNNNNNNN", "NNNNNNNNNN", "NNNNNNNNNNN", "NNNNNNNNNNNN"};
	
	public static String getPolyDel(int nb)
	{
		if(nb < POLY_DEL.length) return POLY_DEL[nb];
		String del = POLY_DEL[POLY_DEL.length - 1];
		for(int i = POLY_DEL.length; i <= nb; i++) del += "-";
		return del;
	}
}

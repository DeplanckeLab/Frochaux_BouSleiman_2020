package model;

import java.util.HashMap;

public class EfficientHashMap // String/Integer HashMap
{
	private HashMap<String, HashMap<String, Integer>> model;
	
	public EfficientHashMap() 
	{
		model = new HashMap<String, HashMap<String, Integer>>();
	}
	
	public Integer get(String chr, String rsId)
	{
		HashMap<String, Integer> map = model.get(chr);
		if(map == null) return null;
		return map.get(rsId);
	}
	
	public void put(String chr, String rsId, int count)
	{
		HashMap<String, Integer> map = model.get(chr);
		if(map == null) map = new HashMap<String, Integer>();
		map.put(rsId, count);
		model.put(chr, map);
	}

	public void remove(String chr, String rsId)
	{
		HashMap<String, Integer> map = model.get(chr);
		if(map != null) map.remove(rsId);
	}
	
	public int size()
	{
		int sum = 0;
		for(String key:model.keySet()) sum += model.get(key).size();
		return sum;
	}
}

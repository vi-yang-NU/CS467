"use client";

import { useState, useEffect, useRef } from "react";
import { FoodEntry } from "@/lib/types";

interface Props {
  onSelect: (food: Omit<FoodEntry, "id">) => void;
  onCancel: () => void;
}

export default function FoodSearchInput({ onSelect, onCancel }: Props) {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<Omit<FoodEntry, "id">[]>([]);
  const [loading, setLoading] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    inputRef.current?.focus();
  }, []);

  useEffect(() => {
    if (!query.trim()) {
      setResults([]);
      return;
    }

    const timeout = setTimeout(async () => {
      setLoading(true);
      try {
        const res = await fetch(`/api/foods?q=${encodeURIComponent(query.trim())}`);
        const data = await res.json();
        setResults(data);
      } catch {
        setResults([]);
      } finally {
        setLoading(false);
      }
    }, 300);

    return () => clearTimeout(timeout);
  }, [query]);

  return (
    <div className="space-y-2">
      {/* Search input row */}
      <div className="flex gap-2">
        <input
          ref={inputRef}
          type="text"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Search for a food..."
          className="flex-1 border border-amber-300 rounded-xl px-4 py-2.5 text-sm focus:outline-none focus:ring-2 focus:ring-amber-400 focus:border-transparent"
        />
        <button
          onClick={onCancel}
          className="px-4 py-2.5 text-sm text-gray-500 hover:text-gray-700 border border-gray-200 rounded-xl active:bg-gray-50 transition-colors"
        >
          Cancel
        </button>
      </div>

      {/* Loading */}
      {loading && (
        <p className="text-center text-gray-400 text-sm py-3">Searching...</p>
      )}

      {/* Results */}
      {!loading && results.length > 0 && (
        <div className="border border-amber-100 rounded-xl overflow-hidden shadow-sm">
          {results.map((food, i) => (
            <button
              key={i}
              onClick={() => onSelect(food)}
              className="w-full flex items-center justify-between px-4 py-3 hover:bg-amber-50 active:bg-amber-100 transition-colors border-b border-amber-50 last:border-0 text-left"
            >
              <div className="flex-1 min-w-0">
                <p className="text-gray-800 font-medium text-sm">{food.name}</p>
                <p className="text-gray-400 text-xs">{food.servingSize}</p>
              </div>
              <span className="text-amber-600 font-semibold text-sm ml-4 shrink-0">
                {food.vitaminD} IU
              </span>
            </button>
          ))}
        </div>
      )}

      {/* No results */}
      {!loading && query.trim() && results.length === 0 && (
        <p className="text-center text-gray-400 text-sm py-3">No foods found</p>
      )}
    </div>
  );
}

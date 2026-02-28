"use client";

import { useState } from "react";
import FoodSearchInput from "./FoodSearchInput";
import { FoodEntry, MealType } from "@/lib/types";

const MEAL_LABELS: Record<MealType, string> = {
  breakfast: "Breakfast",
  lunch: "Lunch",
  dinner: "Dinner",
};

interface Props {
  meal: MealType;
  foods: FoodEntry[];
  onAddFood: (food: Omit<FoodEntry, "id">) => void;
  onRemoveFood: (id: string) => void;
}

export default function MealCard({ meal, foods, onAddFood, onRemoveFood }: Props) {
  const [searching, setSearching] = useState(false);
  const mealTotal = foods.reduce((sum, f) => sum + f.vitaminD, 0);

  return (
    <div className="bg-white rounded-2xl shadow-sm border border-amber-100 overflow-hidden">
      {/* Header */}
      <div className="flex justify-between items-center px-5 py-4 border-b border-gray-50">
        <h2 className="font-semibold text-gray-800 text-lg">{MEAL_LABELS[meal]}</h2>
        {mealTotal > 0 && (
          <span className="text-amber-600 font-semibold text-sm bg-amber-50 px-2.5 py-1 rounded-full">
            {mealTotal} IU
          </span>
        )}
      </div>

      {/* Food list + search */}
      <div className="px-5 py-3 space-y-1">
        {foods.length === 0 && !searching && (
          <p className="text-gray-400 text-sm py-2">No foods added yet</p>
        )}

        {foods.map((food) => (
          <div
            key={food.id}
            className="flex items-center justify-between py-2.5 border-b border-gray-50 last:border-0"
          >
            <div className="flex-1 min-w-0">
              <p className="text-gray-800 font-medium text-sm truncate">{food.name}</p>
              <p className="text-gray-400 text-xs">{food.servingSize}</p>
            </div>
            <div className="flex items-center gap-3 ml-3 shrink-0">
              <span className="text-amber-600 font-semibold text-sm">{food.vitaminD} IU</span>
              <button
                onClick={() => onRemoveFood(food.id)}
                className="w-6 h-6 flex items-center justify-center text-gray-300 hover:text-red-400 active:text-red-500 transition-colors rounded-full"
                aria-label={`Remove ${food.name}`}
              >
                &times;
              </button>
            </div>
          </div>
        ))}

        <div className="pt-1 pb-1">
          {searching ? (
            <FoodSearchInput
              onSelect={(food) => {
                onAddFood(food);
                setSearching(false);
              }}
              onCancel={() => setSearching(false)}
            />
          ) : (
            <button
              onClick={() => setSearching(true)}
              className="w-full py-3 text-amber-600 font-medium text-sm border-2 border-dashed border-amber-200 rounded-xl hover:border-amber-400 hover:bg-amber-50 active:bg-amber-100 transition-all"
            >
              + Add food
            </button>
          )}
        </div>
      </div>
    </div>
  );
}

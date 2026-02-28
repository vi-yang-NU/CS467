"use client";

import { useState } from "react";
import MealCard from "@/components/MealCard";
import VitaminDTotal from "@/components/VitaminDTotal";
import { FoodEntry, MealType } from "@/lib/types";

const MEALS: MealType[] = ["breakfast", "lunch", "dinner"];

export default function Home() {
  const [meals, setMeals] = useState<Record<MealType, FoodEntry[]>>({
    breakfast: [],
    lunch: [],
    dinner: [],
  });

  const addFood = (meal: MealType, food: Omit<FoodEntry, "id">) => {
    setMeals((prev) => ({
      ...prev,
      [meal]: [...prev[meal], { ...food, id: crypto.randomUUID() }],
    }));
  };

  const removeFood = (meal: MealType, id: string) => {
    setMeals((prev) => ({
      ...prev,
      [meal]: prev[meal].filter((f) => f.id !== id),
    }));
  };

  const totalVitD = Object.values(meals)
    .flat()
    .reduce((sum, f) => sum + f.vitaminD, 0);

  const today = new Date().toLocaleDateString("en-US", {
    weekday: "long",
    month: "long",
    day: "numeric",
  });

  return (
    <main className="min-h-screen bg-amber-50 pb-10">
      <header className="bg-amber-500 text-white px-4 py-5 shadow-md">
        <h1 className="text-2xl font-bold">UV Monitor</h1>
        <p className="text-amber-100 text-sm mt-0.5">{today}</p>
      </header>

      <div className="max-w-lg mx-auto px-4 pt-5 space-y-4">
        <VitaminDTotal total={totalVitD} />

        {MEALS.map((meal) => (
          <MealCard
            key={meal}
            meal={meal}
            foods={meals[meal]}
            onAddFood={(food) => addFood(meal, food)}
            onRemoveFood={(id) => removeFood(meal, id)}
          />
        ))}
      </div>
    </main>
  );
}

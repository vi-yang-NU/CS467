import { NextRequest, NextResponse } from "next/server";

// Placeholder food data — replace with CSV-backed backend lookup
const FOODS = [
  { name: "Salmon", vitaminD: 447, servingSize: "3 oz (85g)" },
  { name: "Swordfish", vitaminD: 566, servingSize: "3 oz (85g)" },
  { name: "Rainbow trout", vitaminD: 645, servingSize: "3 oz (85g)" },
  { name: "Mackerel", vitaminD: 388, servingSize: "3 oz (85g)" },
  { name: "Tuna, canned", vitaminD: 154, servingSize: "3 oz (85g)" },
  { name: "Herring", vitaminD: 216, servingSize: "3 oz (85g)" },
  { name: "Halibut", vitaminD: 196, servingSize: "3 oz (85g)" },
  { name: "Sardines, canned", vitaminD: 46, servingSize: "2 sardines (24g)" },
  { name: "Cod liver oil", vitaminD: 448, servingSize: "1 teaspoon (4.5ml)" },
  { name: "Egg", vitaminD: 44, servingSize: "1 large (50g)" },
  { name: "Egg yolk", vitaminD: 41, servingSize: "1 yolk (17g)" },
  { name: "Beef liver", vitaminD: 42, servingSize: "3 oz (85g)" },
  { name: "Milk, fortified", vitaminD: 115, servingSize: "1 cup (240ml)" },
  { name: "Soy milk, fortified", vitaminD: 107, servingSize: "1 cup (240ml)" },
  { name: "Orange juice, fortified", vitaminD: 100, servingSize: "1 cup (240ml)" },
  { name: "Yogurt, fortified", vitaminD: 80, servingSize: "6 oz (170g)" },
  { name: "Cheddar cheese", vitaminD: 6, servingSize: "1 oz (28g)" },
  { name: "Cereal, fortified", vitaminD: 80, servingSize: "1 cup (30g)" },
  { name: "Mushrooms, UV-exposed", vitaminD: 366, servingSize: "1/2 cup (48g)" },
  { name: "Mushrooms, white", vitaminD: 1, servingSize: "1/2 cup (48g)" },
];

export async function GET(request: NextRequest) {
  const q = request.nextUrl.searchParams.get("q")?.toLowerCase().trim() ?? "";

  if (!q) {
    return NextResponse.json([]);
  }

  const results = FOODS.filter((f) => f.name.toLowerCase().includes(q)).slice(0, 8);

  return NextResponse.json(results);
}

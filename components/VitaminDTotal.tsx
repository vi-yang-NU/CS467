const RECOMMENDED_IU = 600;

export default function VitaminDTotal({ total }: { total: number }) {
  const percentage = Math.min((total / RECOMMENDED_IU) * 100, 100);
  const reached = total >= RECOMMENDED_IU;

  return (
    <div className="bg-white rounded-2xl shadow-sm p-5 border border-amber-100">
      <div className="flex justify-between items-baseline mb-2">
        <span className="text-gray-600 font-semibold">Daily Vitamin D</span>
        <span className="text-xs text-gray-400">{RECOMMENDED_IU} IU recommended</span>
      </div>

      <div className="flex items-baseline gap-1.5 mb-3">
        <span className="text-4xl font-bold text-amber-600">{total}</span>
        <span className="text-gray-400 text-sm">IU</span>
      </div>

      <div className="w-full bg-amber-100 rounded-full h-2.5">
        <div
          className={`h-2.5 rounded-full transition-all duration-500 ${
            reached ? "bg-green-500" : "bg-amber-500"
          }`}
          style={{ width: `${percentage}%` }}
        />
      </div>

      <p className="text-sm mt-2 font-medium" style={{ color: reached ? "#16a34a" : "#9ca3af" }}>
        {reached
          ? "Daily goal reached!"
          : `${RECOMMENDED_IU - total} IU to reach daily goal`}
      </p>
    </div>
  );
}

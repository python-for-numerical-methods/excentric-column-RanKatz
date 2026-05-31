def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """  
import math
from scipy.optimize import bisect

def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """
    
    # חישוב עומס קריסה של אוילר (Euler buckling load)
    # עומס זה מהווה את הגבול העליון התיאורטי שלנו, שכן העמוד יקרוס בו בכל מקרה
    P_euler = (math.pi*2 * E * A * (r2)) / (L*2)
    
    # פונקציית המטרה: ההפרש בין המאמץ המחושב למאמץ המותר לפי נוסחת הסקנט
    # אנחנו מחפשים את הערך P שעבורו הפונקציה הזו שווה ל-0
    def stress_diff(P):
        if P <= 0:
            return -sigma_allow
            
        # חישוב הערך שבתוך פונקציית הסקנט (ברדיאנים)
        theta = (L / (2 * r)) * math.sqrt(P / (E * A))
        
        # Secant(x) = 1 / cos(x)
        secant_term = 1.0 / math.cos(theta)
        
        # חישוב המאמץ המקסימלי לעומס P הנתון
        sigma_max = (P / A) * (1.0 + (e * c / (r**2)) * secant_term)
        
        return sigma_max - sigma_allow

    # הגדרת טווח החיפוש עבור שיטת החצייה
    # הגבול התחתון הוא קרוב מאוד ל-0, והגבול העליון הוא כמעט עומס אוילר (כדי למנוע חלוקה ב-0 בקוסינוס)
    P_lower = 1e-5
    P_upper = P_euler * 0.9999
    
    try:
        # שימוש בשיטת החצייה (Bisection method) למציאת השורש
        P_crit = bisect(stress_diff, P_lower, P_upper)
        return float(P_crit)
    except ValueError:
        # שגיאה זו תקפוץ אם הנתונים שהוזנו לא מאפשרים מציאת שורש בטווח
        raise ValueError("Could not find a valid load. Check if inputs are physically reasonable.")

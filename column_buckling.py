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
    P_euler = (np.pi**2 * E * A * (r**2)) / (L**2)
    
    # פונקציית המטרה: ההפרש בין המאמץ המחושב למאמץ המותר לפי נוסחת הסקנט
    def stress_diff(P):
        if P <= 0:
            return -sigma_allow
            
        # חישוב הערך שבתוך פונקציית הסקנט (ברדיאנים)
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        
        # הגנה מפני קוסינוס שלילי או קרוב מדי לאפס (חציית האסימפטוטה של אוילר)
        if theta >= np.pi / 2:
            return float('inf') # מאמץ אינסופי מעבר לנקודת הקריסה
            
        # Secant(x) = 1 / cos(x)
        secant_term = 1.0 / np.cos(theta)
        
        # חישוב המאמץ המקסימלי לעומס P הנתון
        sigma_max = (P / A) * (1.0 + (e * c / (r**2)) * secant_term)
        
        return sigma_max - sigma_allow

    # הגדרת טווח החיפוש עבור שיטת החצייה
    P_lower = 0.0
    
    # מציאת חסם עליון אפקטיבי שמבטיח ערך חיובי לפונקציה
    # נתחיל קרוב לעומס אוילר, ואם הפונקציה לא חיובית, נתקרב אליו עוד יותר
    factor = 0.99
    P_upper = P_euler * factor
    
    # לולאה קטנה שמוודאת שהחסם העליון אכן נותן ערך חיובי לפני ששולחים ל-bisect
    iterations = 0
    while stress_diff(P_upper) <= 0 and iterations < 10:
        factor = factor + (1 - factor) * 0.9
        P_upper = P_euler * factor
        iterations += 1
    
    try:
        # שימוש בשיטת החצייה (Bisection method) למציאת השורש
        P_crit = bisect(stress_diff, P_lower, P_upper, xtol=1e-5)
        return float(P_crit)
    except ValueError:
        # אם בכל זאת אין היפוך סימן, נבדוק אם העמוד נכשל כבר בגלל אקסצנטריות ב-P כמעט אפסי
        if stress_diff(1e-5) > 0:
            return 0.0
        raise ValueError("Could not find a valid load. Check if inputs are physically reasonable.")

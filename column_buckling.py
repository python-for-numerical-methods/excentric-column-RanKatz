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
    
    # חישוב עומס קריסה של אוילר (Euler buckling load) עם תיקון החזקות
    # עומס זה מהווה את הגבול העליון התיאורטי שלנו
    P_euler = (np.pi**2 * E * A * (r**2)) / (L**2)
    
    # פונקציית המטרה: ההפרש בין המאמץ המחושב למאמץ המותר לפי נוסחת הסקנט
    def stress_diff(P):
        if P <= 0:
            return -sigma_allow
            
        # חישוב הערך שבתוך פונקציית הסקנט (ברדיאנים) - שימוש ב-np.sqrt
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        
        # Secant(x) = 1 / cos(x) - שימוש ב-np.cos
        secant_term = 1.0 / np.cos(theta)
        
        # חישוב המאמץ המקסימלי לעומס P הנתון
        sigma_max = (P / A) * (1.0 + (e * c / (r**2)) * secant_term)
        
        return sigma_max - sigma_allow

    # הגדרת טווח החיפוש עבור שיטת החצייה
    # נתחיל מ-0 ממש (או ערך קרוב מאוד אליו) ועד כמעט עומס אוילר
    P_lower = 0.0
    P_upper = P_euler * (1 - 1e-11) # חסם בטוח ומדויק יותר למניעת חלוקה באפס
    
    try:
        # שימוש בשיטת החצייה (Bisection method) למציאת השורש
        P_crit = bisect(stress_diff, P_lower, P_upper, xtol=1e-5)
        return float(P_crit)
    except ValueError:
        # שגיאה זו תקפוץ אם הנתונים שהוזנו לא מאפשרים מציאת שורש בטווח
        raise ValueError("Could not find a valid load. Check if inputs are physically reasonable.")

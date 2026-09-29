---
title: Παραδείγματα
# Στο HTML ο τίτλος και οι επικεφαλίδες εμφανίζονται χωρίς αριθμούς, ώστε η σελίδα να μη μοιάζει
# με ξεχωριστή ενότητα. Τα σχήματα και οι εξισώσεις αριθμούνται σε συνέχεια του κεφαλαίου, όπως
# στο PDF. Και τα δύο στηρίζονται στο tools/patch_myst_greek.py.
numbering:
  title: false
  figure:
    continue: true
  equation:
    continue: true
  table:
    continue: true
kernelspec:
  name: python3
  display_name: Python 3
---

## Γραμμική εκτίμηση κατάστασης (μοντέλο ΣΡ)

Τα παραδείγματα εφαρμόζουν τη μεθοδολογία της Ενότητας [«Εκτίμηση κατάστασης»](#sec:se) σε σύστημα τριών ζυγών. Το πρώτο παράδειγμα χρησιμοποιεί το γραμμικό μοντέλο ΣΡ και εξετάζει την εκτίμηση με σταθμισμένα ελάχιστα
τετράγωνα, τον έλεγχο $\chi^2$, την αναγνώριση εσφαλμένης μέτρησης με το κριτήριο του μεγίστου
κανονικοποιημένου υπολοίπου (LNR), την εύρωστη εκτίμηση WLAV και τη χρήση ψευδομετρήσεων. Το δεύτερο εφαρμόζει τη γενική μη γραμμική διατύπωση με πλήρες μοντέλο εναλλασσόμενου
ρεύματος: συναρτήσεις μετρήσεων και αναλυτικό ιακωβιανό, επαναληπτική επίλυση Gauss–Newton,
ακρίβεια της εκτίμησης, ανίχνευση χονδροειδούς σφάλματος και σύγκριση με το μοντέλο ΣΡ.

:::{tip} Διαδραστική εκτέλεση
Στην ηλεκτρονική έκδοση των σημειώσεων τα κελιά εκτελούνται στον φυλλομετρητή, χωρίς
εγκατάσταση λογισμικού, ενώ στο PDF εμφανίζονται τα αποτελέσματα με τις προεπιλεγμένες τιμές. Ο πυρήνας Python
ενεργοποιείται με το κουμπί **έναρξης περιβάλλοντος υπολογισμού** (η πρώτη φόρτωση διαρκεί
μερικά δευτερόλεπτα) και κάθε κελί εκτελείται με το αντίστοιχο κουμπί εκτέλεσης. Τα κελιά
κάθε ενότητας εκτελούνται με τη σειρά, επειδή χρησιμοποιούν συναρτήσεις και μεταβλητές που
ορίζονται στα προηγούμενα κελιά.

Οι παράμετροι μεταβάλλονται με δύο τρόπους: με τους **ολισθητές** των ενοτήτων «Παραμετρική διερεύνηση με ολισθητές» και «Σύγκριση με το γραμμικό μοντέλο ΣΡ», χωρίς
τροποποίηση του κώδικα, ή με το κουμπί **ανοίγματος σε Jupyter**, το οποίο φορτώνει τη
σελίδα ως σημειωματάριο με δυνατότητα πλήρους επεξεργασίας του κώδικα.
:::

### Το δίκτυο και το γραμμικό μοντέλο

Το [](#fig-diktyo-dc) απεικονίζει το δίκτυο και τις θέσεις των μετρητών. Ο ζυγός 1 τροφοδοτεί
τα φορτία των ζυγών 2 και 3 μέσω τριών κλάδων με αντιδράσεις $x_{12}$ = 0,20, $x_{13}$ = 0,40
και $x_{23}$ = 0,25 p.u. Μετρώνται οι ροές ενεργής ισχύος των τριών κλάδων, στο άκρο του ζυγού 1
για τους κλάδους 1–2 και 1–3 και στο άκρο του ζυγού 2 για τον κλάδο 2–3, καθώς και οι εγχύσεις
ενεργής ισχύος στους ζυγούς 2 και 3.

:::{figure} figures/diktyo-dc.svg
:label: fig-diktyo-dc
:width: 100%

Δίκτυο τριών ζυγών του γραμμικού παραδείγματος, με τις αντιδράσεις των κλάδων σε p.u. και τις
θέσεις των πέντε μετρητών ενεργής ισχύος. Ο ζυγός 1 είναι ο ζυγός αναφοράς.
:::

Χρησιμοποιείται το γραμμικοποιημένο μοντέλο συνεχούς ρεύματος (μοντέλο ΣΡ, DC model), με
τις παραδοχές: αμελητέες ωμικές αντιστάσεις, μοναδιαία μέτρα τάσεων και μικρές διαφορές
γωνιών. Οι ροές ενεργής ισχύος στους κλάδους και οι εγχύσεις στους ζυγούς είναι τότε
γραμμικές συναρτήσεις των γωνιών:

$$
P_{ij} = \frac{\theta_i - \theta_j}{x_{ij}}, \qquad
P_i = \sum_{j \in \mathcal{N}_i,\ j \neq i} \frac{\theta_i - \theta_j}{x_{ij}}
$$

Ο ζυγός 1 είναι ο ζυγός αναφοράς ($\theta_1 = 0$), οπότε το διάνυσμα κατάστασης είναι
$\mathbf{x} = [\theta_2, \theta_3]^{T}$ και $n = 2$. Ο ιακωβιανός πίνακας $\mathbf{H}$ είναι
σταθερός και η εκτίμηση υπολογίζεται χωρίς επαναλήψεις.

```{code-cell} python
import numpy as np

np.set_printoptions(precision=4, suppress=True)

# --- παράμετροι δικτύου (αντιδράσεις κλάδων σε p.u.) ---
x12, x13, x23 = 0.20, 0.40, 0.25

# Μετρήσεις: τρεις ροές κλάδων και δύο εγχύσεις ζυγών.
names = ["P12", "P13", "P23", "P2", "P3"]

H = np.array([
    [-1/x12,          0.0            ],   # P12 = (θ1 - θ2)/x12
    [ 0.0,           -1/x13          ],   # P13 = (θ1 - θ3)/x13
    [ 1/x23,         -1/x23          ],   # P23 = (θ2 - θ3)/x23
    [ 1/x12 + 1/x23, -1/x23          ],   # P2  = έγχυση στον ζυγό 2
    [-1/x23,          1/x13 + 1/x23  ],   # P3  = έγχυση στον ζυγό 3
])

n = H.shape[1]
m = H.shape[0]
print(f"m = {m} μετρήσεις, n = {n} άγνωστοι, "
      f"λόγος πλεονασμού η = {m/n:.1f}")
print("H =", H, sep="\n")
```

### Παραγωγή των μετρήσεων

Οι μετρήσεις παράγονται συνθετικά από γνωστή κατάσταση αναφοράς, με προσθήκη κανονικά
κατανεμημένου θορύβου μηδενικής μέσης τιμής και τυπικής απόκλισης $\sigma$. Η κατάσταση
αναφοράς χρησιμοποιείται μόνο για τον υπολογισμό του σφάλματος εκτίμησης.

```{code-cell} python
# --- παράμετροι: αληθής κατάσταση και ακρίβεια μετρητών ---
theta_true = np.array([-0.05, -0.08])   # γωνίες ζυγών 2 και 3 σε rad
sigma = 0.01                            # τυπική απόκλιση μετρητών, p.u.
seed = 2026                             # για αναπαραγώγιμο θόρυβο

R = np.diag(np.full(m, sigma**2))
Rinv = np.linalg.inv(R)

rng = np.random.default_rng(seed)
z_clean = H @ theta_true
z = z_clean + rng.normal(0.0, sigma, m)

for k, nm in enumerate(names):
    print(f"{nm:>4}:  αληθής {z_clean[k]:+.4f}   μετρούμενη {z[k]:+.4f}")
```

### Επίλυση με σταθμισμένα ελάχιστα τετράγωνα

```{code-cell} python
def wls(z, H, Rinv):
    """Εκτίμηση WLS με τα σχετικά μεγέθη ελέγχου."""
    G = H.T @ Rinv @ H                     # πίνακας κέρδους
    x_hat = np.linalg.solve(G, H.T @ Rinv @ z)
    r = z - H @ x_hat                      # υπόλοιπα
    J = float(r @ Rinv @ r)                # αντικειμενική συνάρτηση
    Omega = np.linalg.inv(Rinv) - H @ np.linalg.inv(G) @ H.T
    with np.errstate(divide="ignore", invalid="ignore"):
        rN = np.abs(r) / np.sqrt(np.diag(Omega))
    return x_hat, r, J, rN, Omega


x_hat, r, J, rN, Omega = wls(z, H, Rinv)

print("θ̂        =", x_hat)
print("θ αληθής  =", theta_true)
print("σφάλμα    =", (x_hat - theta_true) * 1e3, "mrad")
```

Το σφάλμα εκτίμησης είναι της τάξης λίγων mrad, καθώς η εκτίμηση WLS συνδυάζει και τις
πέντε μετρήσεις (λόγος πλεονασμού 2,5) και περιορίζει την επίδραση του θορύβου καθεμίας.

### Έλεγχος $\chi^2$

Οι βαθμοί ελευθερίας είναι $m - n = 3$ και το κατώφλι για πιθανότητα 95\%
είναι 7,815.

```{code-cell} python
dof = m - n
chi2_crit = 7.815      # χ² για 3 βαθμούς ελευθερίας, 95%

print(f"J = {J:.3f}   βαθμοί ελευθερίας = {dof}   κατώφλι = {chi2_crit}")
found = J > chi2_crit
print("Αποτέλεσμα:", "ανιχνεύονται" if found else "δεν ανιχνεύονται",
      "εσφαλμένα δεδομένα")
print()
for k, nm in enumerate(names):
    print(f"{nm:>4}:  r = {r[k]:+.4f}   r^N = {rN[k]:.2f}")
```

### Ανίχνευση και αναγνώριση χονδροειδούς σφάλματος

Στη μέτρηση `P13` προστίθεται χονδροειδές σφάλμα $8\sigma$. Η αλλοιωμένη τιμή παραμένει
εντός του εύλογου εύρους τιμών ροής, οπότε δεν εντοπίζεται με έλεγχο ορίων της μεμονωμένης
μέτρησης.

```{code-cell} python
# --- παράμετροι: ποια μέτρηση αλλοιώνεται και κατά πόσο ---
bad_index = names.index("P13")
bad_size = 8 * sigma

z_bad = z.copy()
z_bad[bad_index] += bad_size

x_b, r_b, J_b, rN_b, _ = wls(z_bad, H, Rinv)
suspect = int(np.argmax(rN_b))

print(f"J = {J_b:.1f}  έναντι κατωφλίου {chi2_crit}  ->  ανιχνεύεται")
print(f"μέγιστο r^N: {names[suspect]}  ({rN_b[suspect]:.2f})")
print()

keep = [i for i in range(m) if i != suspect]
R_keep = R[np.ix_(keep, keep)]
x_c, r_c, J_c, _, _ = wls(z_bad[keep], H[keep], np.linalg.inv(R_keep))
print(f"μετά την απόρριψη της {names[suspect]}:")
print(f"  θ̂ = {x_c}   J = {J_c:.3f}   βαθμοί ελευθερίας = {len(keep) - n}")
print("  σφάλμα ως προς την αληθή κατάσταση =",
      (x_c - theta_true) * 1e3, "mrad")
```

```{code-cell} python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6.2, 3.0))
xpos = np.arange(m)
ax.bar(xpos - 0.2, rN, width=0.4, color="#5b6b78",
       label="χωρίς χονδροειδές σφάλμα")
ax.bar(xpos + 0.2, rN_b, width=0.4, color="#8f2622",
       label="με χονδροειδές σφάλμα")
ax.axhline(3.0, ls="--", lw=1.1, color="#1f2a33")
ax.text(-0.45, 3.2, "κατώφλι 3σ", ha="left", fontsize=9)
ax.set_ylim(0, max(rN_b.max(), 3.0) * 1.18)
ax.set_xticks(xpos, names)
ax.set_ylabel("κανονικοποιημένο υπόλοιπο")
ax.legend(fontsize=9, frameon=False)
fig.tight_layout()
```

Η αλλοιωμένη μέτρηση αναγνωρίζεται από την ασυμφωνία της με τις υπόλοιπες μετρήσεις, η
οποία αποτυπώνεται στο κανονικοποιημένο υπόλοιπό της. Η αναγνώριση είναι δυνατή μόνο λόγω
του πλεονασμού των μετρήσεων. Το σφάλμα επηρεάζει όμως και τα υπόλοιπα άλλων μετρήσεων
(διάχυση του σφάλματος, smearing): στο γράφημα υπερβαίνει το κατώφλι και το κανονικοποιημένο
υπόλοιπο της `P3`. Για τον λόγο αυτόν απορρίπτεται κάθε φορά μόνο η μέτρηση με το μέγιστο
κανονικοποιημένο υπόλοιπο και η εκτίμηση επαναλαμβάνεται.

### Παραμετρική διερεύνηση με ολισθητές

Το κελί που ακολουθεί εκτελεί την ίδια διαδικασία για επιλεγόμενες τιμές τεσσάρων
παραμέτρων: τυπική απόκλιση θορύβου, αλλοιωμένη μέτρηση, μέγεθος χονδροειδούς σφάλματος και
σπόρος της γεννήτριας τυχαίων αριθμών. Σε κάθε μεταβολή παράγονται νέες μετρήσεις,
εκτελείται η εκτίμηση WLS και ελέγχεται αν ο έλεγχος $\chi^2$ ανιχνεύει το σφάλμα και αν
το κριτήριο LNR αναγνωρίζει την αλλοιωμένη μέτρηση.

```{code-cell} python
:tags: [hide-input]
import os

# Κατά την παραγωγή του βιβλίου (BOOK_BUILD=1) το κελί εκτελείται χωρίς
# ολισθητές, ώστε το PDF να περιέχει το γράφημα. Στον φυλλομετρητή το
# ipywidgets εγκαθίσταται κατά την εκτέλεση.
OLISTHITES = os.environ.get("BOOK_BUILD") != "1"
if OLISTHITES:
    try:
        import ipywidgets as widgets
    except ModuleNotFoundError:
        try:
            import micropip
            await micropip.install("ipywidgets")
            import ipywidgets as widgets
        except Exception as exc:
            OLISTHITES = False
            print("Οι ολισθητές δεν είναι διαθέσιμοι:", exc)


def peirama(sigma_pu=0.01, metrisi="P13", megethos=8.0, seed=7):
    R_ = np.diag(np.full(m, sigma_pu**2))
    Rinv_ = np.linalg.inv(R_)
    rng_ = np.random.default_rng(int(seed))

    z_ = H @ theta_true + rng_.normal(0.0, sigma_pu, m)
    k = names.index(metrisi)
    z_[k] += megethos * sigma_pu            # το χονδροειδές σφάλμα

    x_, r_, J_, rN_, _ = wls(z_, H, Rinv_)
    entopismeni = names[int(np.argmax(rN_))]
    anichneuetai = J_ > chi2_crit

    fig, ax = plt.subplots(figsize=(6.2, 2.8))
    colors = ["#8f2622" if nm == metrisi else "#5b6b78" for nm in names]
    ax.bar(np.arange(m), rN_, color=colors, width=0.55)
    ax.axhline(3.0, ls="--", lw=1.1, color="#1f2a33")
    ax.set_xticks(np.arange(m), names)
    ax.set_ylabel("κανονικοποιημένο υπόλοιπο")
    ax.set_ylim(0, max(rN_.max(), 3.5) * 1.15)
    ax.set_title(
        f"J = {J_:.1f}   (κατώφλι {chi2_crit})   ->   "
        + ("ανιχνεύεται" if anichneuetai else "ΔΕΝ ανιχνεύεται"),
        fontsize=10,
    )
    plt.show()

    print(f"αλλοιωμένη: {metrisi}   "
          f"ύποπτη κατά το κριτήριο: {entopismeni}", end="   ")
    sosti = entopismeni == metrisi and anichneuetai
    print("[σωστή]" if sosti else "[αστοχία]")
    print(f"σφάλμα εκτίμησης = {(x_ - theta_true) * 1e3} mrad")


if OLISTHITES:
    # Οι ελληνικές ετικέτες είναι μακρύτερες από τις αγγλικές, οπότε
    # το πλάτος περιγραφής ορίζεται ρητά για να μην αποκόπτονται.
    stil = {"description_width": "185px"}
    diataxi = widgets.Layout(width="520px")

    _ = widgets.interact(      # η ανάθεση αποτρέπει την εκτύπωση του repr
        peirama,
        sigma_pu=widgets.FloatSlider(
            value=0.01, min=0.002, max=0.05, step=0.002,
            description="θόρυβος σ (p.u.)", readout_format=".3f",
            style=stil, layout=diataxi),
        metrisi=widgets.Dropdown(
            options=names, value="P13", description="αλλοιωμένη μέτρηση",
            style=stil, layout=diataxi),
        megethos=widgets.FloatSlider(
            value=8.0, min=0.0, max=20.0, step=0.5,
            description="μέγεθος σφάλματος (×σ)", readout_format=".1f",
            style=stil, layout=diataxi),
        seed=widgets.IntSlider(
            value=7, min=1, max=50, description="τυχαίο δείγμα",
            style=stil, layout=diataxi),
    )
else:
    peirama()
```

Επαναλαμβάνοντας τη διαδικασία για 2000 τυχαία δείγματα θορύβου με $\sigma$ = 0,01 p.u.,
προκύπτουν δύο παρατηρήσεις:

- Για σφάλμα $2\sigma$ ο έλεγχος $\chi^2$ ανιχνεύει το σφάλμα μόνο στο 14–30\% των δειγμάτων,
  ανάλογα με τη μέτρηση, επειδή το σφάλμα είναι συγκρίσιμο με τον θόρυβο μέτρησης. Για
  σφάλμα $8\sigma$ η ανίχνευση υπερβαίνει το 99\%.
- Η αναγνώριση δυσχεραίνεται για μετρήσεις με μικρή διασπορά υπολοίπου και ισχυρά
  συσχετισμένα υπόλοιπα. Οι εγχύσεις έχουν $\Omega_{ii}$ ≈ 0,38–0,39 $\sigma^2$
  έναντι 0,60–0,83 $\sigma^2$ για τις ροές, ενώ τα υπόλοιπα των `P2` και `P12` έχουν
  συντελεστή συσχέτισης 0,86. Ως σωστή αναγνώριση μετρά η περίπτωση όπου ο έλεγχος $\chi^2$
  ανιχνεύει το σφάλμα και η αλλοιωμένη μέτρηση έχει το μέγιστο κανονικοποιημένο υπόλοιπο, όπως
  στην ένδειξη «[σωστή]» του κελιού με τους ολισθητές. Για σφάλμα $4\sigma$ αναγνωρίζονται
  έτσι σωστά οι `P13` και `P23` σε περισσότερο από το 80\% των δειγμάτων, οι εγχύσεις `P2` και
  `P3` όμως σε λιγότερο από το 45\%· για σφάλμα $8\sigma$ η `P2` αναγνωρίζεται στο 88\% των
  δειγμάτων.

### Εύρωστη εκτίμηση WLAV

Ο εκτιμητής WLAV της Ενότητας [«Παραλλαγές της εκτίμησης WLS»](#sec:wls-variants)
ελαχιστοποιεί το άθροισμα $\sum_i |r_i|/\sigma_i$ αντί για το άθροισμα των τετραγώνων των
υπολοίπων. Στη γενική περίπτωση το πρόβλημα επιλύεται με γραμμικό προγραμματισμό και η λύση
του έχει μηδενικό υπόλοιπο σε τουλάχιστον $n$ μετρήσεις. Στο παράδειγμα, με $n = 2$, η λύση
βρίσκεται επομένως απλούστερα: για κάθε ζεύγος μετρήσεων που προσδιορίζει και τις δύο γωνίες
υπολογίζεται η κατάσταση που ικανοποιεί ακριβώς τις δύο μετρήσεις, και επιλέγεται εκείνη με
το μικρότερο σταθμισμένο άθροισμα απόλυτων υπολοίπων. Η συνάρτηση `wlav` υλοποιεί αυτή την
απαρίθμηση, η οποία είναι πρακτική μόνο για μικρά δίκτυα.

Οι δύο εκτιμητές συγκρίνονται σε τέσσερις περιπτώσεις χονδροειδών σφαλμάτων. Η εκτίμηση WLS
υπολογίζεται χωρίς απόρριψη και με διαδοχική απόρριψη της μέτρησης με το μέγιστο
κανονικοποιημένο υπόλοιπο, έως ότου ο έλεγχος $\chi^2$ δεν ανιχνεύει πλέον εσφαλμένα δεδομένα
(εξωτερικός βρόχος, βλ. [](#fig:se-flow)).

```{code-cell} python
from itertools import combinations


def wlav(z, H, sig):
    """Εκτίμηση WLAV με απαρίθμηση: η λύση ικανοποιεί ακριβώς n μετρήσεις."""
    m_, n_ = H.shape
    best = None
    for S in combinations(range(m_), n_):
        HS = H[list(S)]
        if np.linalg.matrix_rank(HS) < n_:
            continue            # το ζεύγος δεν προσδιορίζει την κατάσταση
        x_S = np.linalg.solve(HS, z[list(S)])
        cost = float(np.sum(np.abs(z - H @ x_S) / sig))
        if best is None or cost < best[0]:
            best = (cost, x_S, S)
    _, x_hat, S = best
    return x_hat, z - H @ x_hat, S


CHI2_95 = {1: 3.841, 2: 5.991, 3: 7.815}      # κατώφλια χ² για πιθανότητα 95%


def wls_lnr(z, H, sig):
    """WLS με διαδοχική απόρριψη της μέτρησης με το μέγιστο r^N."""
    keep = list(range(len(z)))
    while True:
        Rinv_ = np.diag(1.0 / sig[keep]**2)
        x_, r_, J_, rN_, _ = wls(z[keep], H[keep], Rinv_)
        dof_ = len(keep) - H.shape[1]
        if dof_ == 0 or J_ <= CHI2_95[dof_]:
            return x_, [k for k in range(len(z)) if k not in keep]
        keep.pop(int(np.argmax(rN_)))


# --- παράμετροι: χονδροειδή σφάλματα (μέτρηση: μέγεθος σε σ) ---
cases = [{"P13": 8}, {"P13": 20}, {"P12": 8}, {"P13": 8, "P23": 8}]
sig = np.full(m, sigma)


def sfalma(x_):
    e = (x_ - theta_true) * 1e3
    return f"({e[0]:+5.1f}, {e[1]:+5.1f}) mrad"


for bad in cases:
    z_ = z.copy()
    for nm, k in bad.items():
        z_[names.index(nm)] += k * sigma
    x_w = wls(z_, H, Rinv)[0]
    x_n, rej = wls_lnr(z_, H, sig)
    x_l, _, S = wlav(z_, H, sig)
    perigrafi = ", ".join(f"{nm} +{k}σ" for nm, k in bad.items())
    print("χονδροειδές σφάλμα:", perigrafi)
    print(f"  WLS           σφάλμα {sfalma(x_w)}")
    print(f"  WLS και LNR   σφάλμα {sfalma(x_n)}   απορρίπτονται: "
          + ", ".join(names[k] for k in rej))
    print(f"  WLAV          σφάλμα {sfalma(x_l)}   ικανοποιούνται ακριβώς: "
          + ", ".join(names[k] for k in S))
    print()
```

Με σφάλμα $8\sigma$ στην `P13` και οι τρεις εκτιμήσεις έχουν σφάλμα λίγων mrad. Η WLAV
ικανοποιεί ακριβώς τις `P12` και `P3` και στην ουσία αγνοεί την αλλοιωμένη μέτρηση. Όταν το
χονδροειδές σφάλμα αυξάνεται σε $20\sigma$, το σφάλμα της WLS χωρίς απόρριψη αυξάνεται
ανάλογα, ενώ η εκτίμηση WLAV δεν μεταβάλλεται καθόλου, επειδή η λύση παραμένει στο ίδιο ζεύγος
μετρήσεων. Όταν αλλοιώνεται η `P12`, το κριτήριο LNR απορρίπτει την ορθή `P2`, της οποίας το
υπόλοιπο είναι ισχυρά συσχετισμένο με εκείνο της `P12`, και η εκτίμηση μετά την απόρριψη
είναι χειρότερη από την αρχική. Με δύο ταυτόχρονα χονδροειδή σφάλματα, στις `P13` και `P23`, η
διαδοχική απόρριψη αφαιρεί τις ορθές `P12` και `P3` και διατηρεί τις αλλοιωμένες. Αυτή είναι
η αστοχία του κριτηρίου LNR σε πολλαπλά εσφαλμένα δεδομένα που αναφέρεται στην Ενότητα
[«Ανίχνευση και αναγνώριση εσφαλμένων μετρήσεων»](#sec:baddata). Η WLAV δεν χρειάζεται
αναγνώριση και απόρριψη και έχει σφάλμα της τάξης του θορύβου σε όλες τις περιπτώσεις.

Το γράφημα συγκρίνει τα υπόλοιπα των δύο εκτιμητών για σφάλμα $8\sigma$ στην `P13`, σε
μονάδες τυπικής απόκλισης.

```{code-cell} python
_, r_l, S_l = wlav(z_bad, H, sig)
xpos = np.arange(m)

fig, ax = plt.subplots(figsize=(6.2, 3.0))
ax.bar(xpos - 0.2, np.abs(r_b) / sigma, width=0.4, color="#5b6b78",
       label="WLS")
ax.bar(xpos + 0.2, np.abs(r_l) / sigma, width=0.4, color="#8f2622",
       label="WLAV")
ax.set_xticks(xpos, names)
ax.set_ylabel("|υπόλοιπο| / σ")
ax.legend(fontsize=9, frameon=False)
fig.tight_layout()
```

Στη WLS μέρος του σφάλματος διαχέεται στα υπόλοιπα των άλλων μετρήσεων, ενώ η WLAV το
συγκεντρώνει στην αλλοιωμένη μέτρηση, με μηδενικά υπόλοιπα στις δύο μετρήσεις που ικανοποιεί
ακριβώς.

Η ευρωστία έχει όμως κόστος σε ακρίβεια. Χωρίς χονδροειδή σφάλματα, η εκτίμηση WLAV
στηρίζεται ουσιαστικά σε δύο μόνο μετρήσεις και δεν μεσοτιμεί τον θόρυβο όλων, όπως η WLS. Το
κελί που ακολουθεί συγκρίνει την ενεργό τιμή του σφάλματος των δύο εκτιμητών σε 1000 τυχαία
δείγματα θορύβου.

```{code-cell} python
rng_mc = np.random.default_rng(1)
err_w, err_l = [], []
for _ in range(1000):
    z_ = z_clean + rng_mc.normal(0.0, sigma, m)
    err_w.append(wls(z_, H, Rinv)[0] - theta_true)
    err_l.append(wlav(z_, H, sig)[0] - theta_true)

rms_w = np.sqrt(np.mean(np.square(err_w), axis=0)) * 1e3
rms_l = np.sqrt(np.mean(np.square(err_l), axis=0)) * 1e3
print("ενεργός τιμή του σφάλματος (θ2, θ3) χωρίς χονδροειδή σφάλματα:")
print(f"  WLS   ({rms_w[0]:.2f}, {rms_w[1]:.2f}) mrad")
print(f"  WLAV  ({rms_l[0]:.2f}, {rms_l[1]:.2f}) mrad")
```

Το σφάλμα της WLAV είναι περίπου 20–25\% μεγαλύτερο. Η WLS είναι επομένως ακριβέστερη όταν
οι μετρήσεις περιέχουν μόνο θόρυβο, ενώ η WLAV είναι ασφαλέστερη όταν υπάρχουν χονδροειδή
σφάλματα.

### Ψευδομετρήσεις

Έστω ότι, λόγω βλαβών στην τηλεμετρία, διατίθεται μόνο η μέτρηση `P12`. Η γωνία $\theta_3$
δεν προσδιορίζεται και το σύστημα δεν είναι παρατηρήσιμο (Ενότητα
[«Παρατηρησιμότητα»](#sec:observability)). Η παρατηρησιμότητα αποκαθίσταται με ψευδομετρήσεις
των εγχύσεων `P2` και `P3`, δηλαδή με προβλέψεις των φορτίων των ζυγών 2 και 3, με σφάλμα
πρόβλεψης +15\% και −10\% αντίστοιχα. Στις ψευδομετρήσεις αποδίδεται τυπική απόκλιση
0,04 p.u., της τάξης του σφάλματος πρόβλεψης και τετραπλάσια από εκείνη των μετρητών. Για
σύγκριση, η εκτίμηση επαναλαμβάνεται με τις ψευδομετρήσεις να θεωρούνται εξίσου ακριβείς με
τους μετρητές. Η ακρίβεια κάθε εκτίμησης υπολογίζεται από τα διαγώνια στοιχεία του
$\mathbf{G}^{-1}$.

```{code-cell} python
# --- παράμετροι: διαθέσιμη μέτρηση και ψευδομετρήσεις ---
diathesimes = ["P12"]                       # διαθέσιμες μετρήσεις
psevdo = ["P2", "P3"]                       # εγχύσεις από πρόβλεψη φορτίου
forecast_error = np.array([+0.15, -0.10])   # σχετικό σφάλμα πρόβλεψης
sigma_ps = 0.04                             # τυπ. απόκλιση ψευδομετρήσεων

iM = [names.index(nm) for nm in diathesimes]
iP = [names.index(nm) for nm in psevdo]
taxi = np.linalg.matrix_rank(H[iM])
print(f"τάξη του H με τις διαθέσιμες μετρήσεις: {taxi}, n = {n}")
print()

H_ps = H[iM + iP]
z_ps = np.concatenate([z[iM], z_clean[iP] * (1 + forecast_error)])


def akriveia(H_, sig_):
    """Τυπικές αποκλίσεις (mrad) από τα διαγώνια στοιχεία του G^-1."""
    G_ = H_.T @ np.diag(1.0 / sig_**2) @ H_
    return np.sqrt(np.diag(np.linalg.inv(G_))) * 1e3


periptoseis = [("όλες οι μετρήσεις (αναφορά)", z, H, np.full(m, sigma))]
for s in (sigma_ps, sigma):
    sig_ = np.array([sigma] * len(iM) + [s] * len(iP))
    etiketa = f"P12 και ψευδομετρήσεις, σ = {s:.2f}"
    periptoseis.append((etiketa, z_ps, H_ps, sig_))

print(f"{'':34}{'σφάλμα (mrad)':>16}{'τυπ. απόκλιση (mrad)':>22}")
print(f"{'':34}{'θ2':>8}{'θ3':>8}{'θ2':>14}{'θ3':>8}")
for label, z_, H_, sig_ in periptoseis:
    x_ = wls(z_, H_, np.diag(1.0 / sig_**2))[0]
    e = (x_ - theta_true) * 1e3
    s = akriveia(H_, sig_)
    print(f"{label:<34}{e[0]:+8.2f}{e[1]:+8.2f}{s[0]:14.2f}{s[1]:8.2f}")
```

Με μόνη την `P12` ο $\mathbf{H}$ έχει τάξη 1 και ο $\mathbf{G}$ δεν αντιστρέφεται. Με τις
ψευδομετρήσεις η εκτίμηση υπολογίζεται, αλλά η ακρίβειά της διαφέρει ανά μεταβλητή: η
$\theta_2$ προσδιορίζεται κυρίως από την ακριβή μέτρηση `P12`, ενώ η $\theta_3$ μόνο από τις
προβλέψεις, οπότε το σφάλμα της είναι περίπου 2,5 φορές μεγαλύτερο από εκείνο της εκτίμησης
με όλες τις μετρήσεις. Με τυπική απόκλιση 0,04 p.u. το σφάλμα της $\theta_3$ είναι συμβατό με
την τυπική απόκλιση που προκύπτει από τον $\mathbf{G}^{-1}$. Αν οι ψευδομετρήσεις θεωρηθούν
εξίσου ακριβείς με τους μετρητές, το σφάλμα της $\theta_3$ παραμένει της ίδιας τάξης, αλλά η
τυπική απόκλιση που αναφέρει ο εκτιμητής είναι σχεδόν τρεις φορές μικρότερη από το πραγματικό
σφάλμα. Η μεγάλη τυπική απόκλιση των ψευδομετρήσεων είναι επομένως αναγκαία, ώστε η
αβεβαιότητα των προβλέψεων να αποτυπώνεται στην ακρίβεια της εκτίμησης.

### Προτεινόμενες ασκήσεις (μοντέλο ΣΡ)

Προτείνονται οι ακόλουθες μεταβολές, με επανεκτέλεση των αντίστοιχων κελιών:

- **`sigma`**: αύξηση σε `0.03`. Πόσο μεγαλώνει το σφάλμα της εκτίμησης και ποιο μέγεθος
  χονδροειδούς σφάλματος παραμένει πλέον ανιχνεύσιμο;
- **`bad_size`**: μείωση σε `2 * sigma`. Ο έλεγχος $\chi^2$ εξακολουθεί να ανιχνεύει το σφάλμα;
- **`bad_index`**: αλλοίωση της `P2` αντί της `P13`. Το κριτήριο LNR αναγνωρίζει
  πάντοτε την αλλοιωμένη μέτρηση;
- **`cases`**: προσθήκη περιπτώσεων με τρία ταυτόχρονα χονδροειδή σφάλματα. Δίνει και τότε η
  WLAV σωστή εκτίμηση;
- **`forecast_error`** και **`sigma_ps`**: μεταβολή του σφάλματος πρόβλεψης και της τυπικής
  απόκλισης των ψευδομετρήσεων. Πώς μεταβάλλονται το σφάλμα και η τυπική απόκλιση της
  $\hat{\theta}_3$; Τι αλλάζει αν στις διαθέσιμες μετρήσεις προστεθεί η `P13`;

## Μη γραμμική εκτίμηση κατάστασης (μοντέλο ΕΡ)

Το παράδειγμα εφαρμόζει τη γενική μη γραμμική διατύπωση της Ενότητας [«Εκτίμηση κατάστασης»](#sec:se) στο ίδιο δίκτυο τριών
ζυγών, με πλήρες μοντέλο κλάδων (ωμική αντίσταση, αντίδραση σειράς και εγκάρσια χωρητική
επιδεκτικότητα) και μετρήσεις μέτρων τάσης, ροών και εγχύσεων ενεργής και άεργης ισχύος. Οι
συναρτήσεις μετρήσεων και ο ιακωβιανός πίνακας υλοποιούνται αναλυτικά σύμφωνα με τις σχέσεις
της ίδιας ενότητας και η εκτίμηση υπολογίζεται με τη μέθοδο Gauss–Newton από επίπεδο προφίλ.

### Δίκτυο και πίνακας αγωγιμοτήτων ζυγών

Τα δεδομένα των κλάδων σε τιμές p.u. (βάση 100 MVA) δίνει ο [](#tab-kladoi-ac).

:::{table} Δεδομένα των κλάδων του μη γραμμικού παραδείγματος σε p.u. (βάση 100 MVA).
:label: tab-kladoi-ac

| Κλάδος | $r$ (p.u.) | $x$ (p.u.) | $b$ (p.u.) |
|:------:|:----------:|:----------:|:----------:|
| 1–2    | 0,020      | 0,20       | 0,04       |
| 1–3    | 0,040      | 0,40       | 0,06       |
| 2–3    | 0,025      | 0,25       | 0,05       |
:::

Οι αντιδράσεις είναι ίδιες με εκείνες του γραμμικού παραδείγματος και ο λόγος $r/x$ είναι 0,1
σε όλους τους κλάδους. Η $b$ είναι η συνολική χωρητική επιδεκτικότητα της γραμμής και
κατανέμεται εξίσου στα δύο άκρα του ισοδύναμου κυκλώματος π.

Η τοπολογία του δικτύου είναι ίδια με εκείνη του γραμμικού παραδείγματος. Το [](#fig-diktyo-ac)
απεικονίζει επιπλέον τις θέσεις των μετρητών του μοντέλου ΕΡ, οι οποίοι περιγράφονται στην
ενότητα «Μετρήσεις».

:::{figure} figures/diktyo-ac.svg
:label: fig-diktyo-ac
:width: 100%

Δίκτυο τριών ζυγών του μη γραμμικού παραδείγματος και θέσεις των μετρητών: ροές ενεργής και
άεργης ισχύος στους τρεις κλάδους, εγχύσεις ενεργής και άεργης ισχύος στους ζυγούς 2 και 3 και
μέτρα τάσεων στους ζυγούς 1 και 3 (συνολικά 12 μετρήσεις).
:::

```{code-cell} python
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=4, suppress=True)

# --- δεδομένα κλάδων: (ζυγός i, ζυγός j, r, x, b ολικό) σε p.u. ---
KLADOI = [
    (1, 2, 0.020, 0.20, 0.04),
    (1, 3, 0.040, 0.40, 0.06),
    (2, 3, 0.025, 0.25, 0.05),
]
NB = 3                                    # πλήθος ζυγών


def build_ybus(kladoi):
    """Πίνακας αγωγιμοτήτων ζυγών και παράμετροι π κάθε κλάδου."""
    Y = np.zeros((NB, NB), dtype=complex)
    pi = {}
    for i, j, r, x, b in kladoi:
        i, j = i - 1, j - 1
        y = 1 / complex(r, x)          # σειριακή αγωγιμότητα g + jb
        Y[i, i] += y + 1j * b / 2
        Y[j, j] += y + 1j * b / 2
        Y[i, j] -= y
        Y[j, i] -= y
        # παράμετροι π: (g_ij, b_ij, b_si)
        pi[(i, j)] = pi[(j, i)] = (y.real, y.imag, b / 2)
    return Y, pi


Y_ac, PI_ac = build_ybus(KLADOI)
print("Ybus =", Y_ac, sep="\n")
```

### Συναρτήσεις μετρήσεων και ιακωβιανός πίνακας

Κάθε μέτρηση περιγράφεται από το είδος της (`V`, `Pi`, `Qi`, `Pf`, `Qf`) και τους ζυγούς στους
οποίους αναφέρεται. Η συνάρτηση `h_ac` υλοποιεί τις σχέσεις εγχύσεων και ροών της Ενότητας [«Συναρτήσεις μετρήσεων»](#sec:meas-func) και η
`H_ac` τις αναλυτικές παραγώγους τους. Το διάνυσμα κατάστασης είναι
$\mathbf{x} = [\theta_2, \theta_3, V_1, V_2, V_3]^{T}$ με $\theta_1 = 0$.

```{code-cell} python
def unpack(x):
    """Γωνίες (με θ1 = 0) και μέτρα τάσεων από το διάνυσμα κατάστασης."""
    return np.r_[0.0, x[:NB - 1]], x[NB - 1:]


def h_ac(x, meas, Y, pi):
    th, V = unpack(x)
    G, B = Y.real, Y.imag
    h = np.zeros(len(meas))
    for k, (kind, i, j) in enumerate(meas):
        if kind == "V":
            h[k] = V[i]
        elif kind in ("Pi", "Qi"):
            c, s = np.cos(th[i] - th), np.sin(th[i] - th)
            if kind == "Pi":
                h[k] = V[i] * np.sum(V * (G[i] * c + B[i] * s))
            else:
                h[k] = V[i] * np.sum(V * (G[i] * s - B[i] * c))
        else:
            g, b, bs = pi[(i, j)]
            c, s = np.cos(th[i] - th[j]), np.sin(th[i] - th[j])
            if kind == "Pf":
                h[k] = V[i]**2 * g - V[i] * V[j] * (g * c + b * s)
            else:
                h[k] = -V[i]**2 * (b + bs) - V[i] * V[j] * (g * s - b * c)
    return h


def H_ac(x, meas, Y, pi):
    th, V = unpack(x)
    G, B = Y.real, Y.imag
    Hf = np.zeros((len(meas), 2 * NB))        # στήλες: θ1..θN, V1..VN
    for k, (kind, i, j) in enumerate(meas):
        if kind == "V":
            Hf[k, NB + i] = 1.0
        elif kind in ("Pi", "Qi"):
            cv, sv = np.cos(th[i] - th), np.sin(th[i] - th)
            Pi = V[i] * np.sum(V * (G[i] * cv + B[i] * sv))
            Qi = V[i] * np.sum(V * (G[i] * sv - B[i] * cv))
            Gii, Bii = G[i, i], B[i, i]
            for q in range(NB):
                c, s = cv[q], sv[q]
                if q == i and kind == "Pi":
                    Hf[k, i] = -Qi - Bii * V[i]**2
                    Hf[k, NB + i] = Pi / V[i] + Gii * V[i]
                elif q == i:
                    Hf[k, i] = Pi - Gii * V[i]**2
                    Hf[k, NB + i] = Qi / V[i] - Bii * V[i]
                elif kind == "Pi":
                    Hf[k, q] = V[i] * V[q] * (G[i, q] * s - B[i, q] * c)
                    Hf[k, NB + q] = V[i] * (G[i, q] * c + B[i, q] * s)
                else:
                    Hf[k, q] = -V[i] * V[q] * (G[i, q] * c + B[i, q] * s)
                    Hf[k, NB + q] = V[i] * (G[i, q] * s - B[i, q] * c)
        else:
            g, b, bs = pi[(i, j)]
            c, s = np.cos(th[i] - th[j]), np.sin(th[i] - th[j])
            if kind == "Pf":
                dth = V[i] * V[j] * (g * s - b * c)
                Hf[k, i], Hf[k, j] = dth, -dth
                Hf[k, NB + i] = -V[j] * (g * c + b * s) + 2 * g * V[i]
                Hf[k, NB + j] = -V[i] * (g * c + b * s)
            else:
                dth = -V[i] * V[j] * (g * c + b * s)
                Hf[k, i], Hf[k, j] = dth, -dth
                Hf[k, NB + i] = (-V[j] * (g * s - b * c)
                                 - 2 * V[i] * (b + bs))
                Hf[k, NB + j] = -V[i] * (g * s - b * c)
    return np.delete(Hf, 0, axis=1)       # χωρίς τη στήλη της θ1
```

### Μετρήσεις

Το σύνολο μετρήσεων περιλαμβάνει τα μέτρα των τάσεων των ζυγών 1 και 3, τις ροές ενεργής και
άεργης ισχύος των τριών κλάδων, στο άκρο του ζυγού 1 για τους κλάδους 1–2 και 1–3 και στο άκρο
του ζυγού 2 για τον κλάδο 2–3, και τις εγχύσεις ενεργής και άεργης ισχύος στους ζυγούς 2 και 3
([](#fig-diktyo-ac)): $m = 12$ μετρήσεις για $n = 5$ μεταβλητές κατάστασης, με λόγο
πλεονασμού 2,4. Οι τυπικές αποκλίσεις είναι 0,004 p.u. για τις τάσεις, 0,008 p.u. για τις ροές
και 0,010 p.u. για τις εγχύσεις. Οι μετρήσεις παράγονται από γνωστή κατάσταση αναφοράς με
προσθήκη θορύβου. Πριν από την εκτίμηση, ο αναλυτικός ιακωβιανός ελέγχεται με κεντρικές
πεπερασμένες διαφορές.

```{code-cell} python
# --- μετρήσεις: (είδος, ζυγός i, ζυγός j), αρίθμηση ζυγών από το 0 ---
MEAS_AC = [("V", 0, None), ("V", 2, None),
           ("Pf", 0, 1), ("Qf", 0, 1), ("Pf", 0, 2), ("Qf", 0, 2),
           ("Pf", 1, 2), ("Qf", 1, 2),
           ("Pi", 1, None), ("Qi", 1, None),
           ("Pi", 2, None), ("Qi", 2, None)]
SIGMA_TYPE = {"V": 0.004, "Pf": 0.008, "Qf": 0.008,
              "Pi": 0.010, "Qi": 0.010}


def label(kind, i, j):
    base = {"V": "V", "Pi": "P", "Qi": "Q", "Pf": "P", "Qf": "Q"}[kind]
    return f"{base}{i + 1}" + ("" if j is None else f"{j + 1}")


names_ac = [label(*mm) for mm in MEAS_AC]
sig_ac = np.array([SIGMA_TYPE[kind] for kind, _, _ in MEAS_AC])
m_ac, n_ac = len(MEAS_AC), 2 * NB - 1

# --- κατάσταση αναφοράς ---
th_ref = np.array([0.0, -0.045, -0.075])          # rad
V_ref = np.array([1.020, 0.990, 0.975])           # p.u.
x_ref = np.r_[th_ref[1:], V_ref]

rng_ac = np.random.default_rng(2026)
z_ref = h_ac(x_ref, MEAS_AC, Y_ac, PI_ac)
z_ac = z_ref + rng_ac.normal(0.0, sig_ac)

# --- έλεγχος του αναλυτικού ιακωβιανού με κεντρικές διαφορές ---
def h_ref(x):
    return h_ac(x, MEAS_AC, Y_ac, PI_ac)


x_test, d = x_ref + 0.01, 1e-6
H_fd = np.column_stack(
    [(h_ref(x_test + d * e) - h_ref(x_test - d * e)) / (2 * d)
     for e in np.eye(n_ac)])
print(f"μέγιστη απόκλιση αναλυτικού και αριθμητικού ιακωβιανού: "
      f"{np.abs(H_ac(x_test, MEAS_AC, Y_ac, PI_ac) - H_fd).max():.1e}")
print()
print("μέτρηση    σ     αναφοράς   μετρούμενη")
for nm, s, zr, zm in zip(names_ac, sig_ac, z_ref, z_ac):
    print(f"{nm:>5}   {s:.3f}   {zr:+.4f}    {zm:+.4f}")
```

### Επίλυση με τη μέθοδο Gauss–Newton

Η συνάρτηση `wls_ac` υλοποιεί τον αλγόριθμο της Ενότητας [«Εκτίμηση σταθμισμένων ελαχίστων τετραγώνων»](#sec:wls): σε κάθε επανάληψη υπολογίζει τα
$\mathbf{h}(\mathbf{x}^{k})$, $\mathbf{H}(\mathbf{x}^{k})$ και $\mathbf{G}(\mathbf{x}^{k})$,
επιλύει τις κανονικές εξισώσεις και ελέγχει τη σύγκλιση με ανοχή $10^{-6}$.

```{code-cell} python
def wls_ac(z, sig, meas, Y, pi, x0, tol=1e-6, itmax=20, verbose=False):
    """Εκτίμηση WLS με τη μέθοδο Gauss–Newton."""
    Rinv = np.diag(1 / sig**2)
    x, history = x0.astype(float).copy(), []
    for k in range(itmax):
        r = z - h_ac(x, meas, Y, pi)
        H = H_ac(x, meas, Y, pi)
        G = H.T @ Rinv @ H
        dx = np.linalg.solve(G, H.T @ Rinv @ r)
        history.append((k, float(r @ Rinv @ r), np.abs(dx).max()))
        if verbose:
            _, Jk, step = history[-1]
            print(f"k = {k}   J(x^k) = {Jk:10.3f}   max|Δx| = {step:.2e}")
        x = x + dx
        if np.abs(dx).max() < tol:
            break
    r = z - h_ac(x, meas, Y, pi)
    H = H_ac(x, meas, Y, pi)
    G = H.T @ Rinv @ H
    Omega = np.diag(sig**2) - H @ np.linalg.solve(G, H.T)
    with np.errstate(divide="ignore", invalid="ignore"):
        rN = np.abs(r) / np.sqrt(np.diag(Omega))
    return dict(x=x, r=r, J=float(r @ Rinv @ r), rN=rN, G=G,
                history=history)


x_flat = np.r_[np.zeros(NB - 1), np.ones(NB)]          # επίπεδο προφίλ
est = wls_ac(z_ac, sig_ac, MEAS_AC, Y_ac, PI_ac, x_flat, verbose=True)

th_hat, V_hat = unpack(est["x"])
std = np.sqrt(np.diag(np.linalg.inv(est["G"])))
std_th, std_V = np.r_[0.0, std[:NB - 1]], std[NB - 1:]
print()
deg = np.degrees
print("ζυγός   θ εκτ.   θ αναφ.   σ_θ     V εκτ.   V αναφ.   σ_V")
print("         (°)      (°)      (°)     (p.u.)   (p.u.)    (p.u.)")
for b in range(NB):
    print(f"  {b + 1}    {deg(th_hat[b]):+7.3f}  {deg(th_ref[b]):+7.3f}   "
          f"{deg(std_th[b]):.3f}   {V_hat[b]:.4f}   {V_ref[b]:.4f}    "
          f"{std_V[b]:.4f}")
```

```{code-cell} python
it = [hh[0] for hh in est["history"]]
step = [hh[2] for hh in est["history"]]

fig, ax = plt.subplots(figsize=(6.2, 2.8))
ax.semilogy(it, step, "o-", color="#8f2622", lw=1.6)
ax.axhline(1e-6, ls="--", lw=1.0, color="#1f2a33")
ax.text(it[0], 1.6e-6, "ανοχή", ha="left", fontsize=9)
ax.set_xticks(it)
ax.set_xlabel("επανάληψη k")
ax.set_ylabel("max |Δx| (p.u., rad)")
ax.grid(True, which="both", lw=0.3, alpha=0.5)
fig.tight_layout()
```

Η μέθοδος συγκλίνει σε τρεις επαναλήψεις από επίπεδο προφίλ. Το μέγιστο βήμα κάθε επανάληψης
είναι περίπου ανάλογο του τετραγώνου του προηγούμενου, όπως αναμένεται για τη μέθοδο
Gauss–Newton όταν τα υπόλοιπα στη λύση είναι μικρά. Οι αποκλίσεις της εκτίμησης από την
κατάσταση αναφοράς είναι της τάξης των τυπικών αποκλίσεων που προκύπτουν από τα διαγώνια
στοιχεία του $\mathbf{G}^{-1}$.

### Έλεγχος χ² και κανονικοποιημένα υπόλοιπα

Για $m - n = 7$ βαθμούς ελευθερίας το κατώφλι του ελέγχου για πιθανότητα 95% είναι 14,067.

```{code-cell} python
CHI2_95 = {7: 14.067, 6: 12.592}       # κατώφλια χ² για πιθανότητα 95%

dof_ac = m_ac - n_ac
print(f"J(x̂) = {est['J']:.3f}   βαθμοί ελευθερίας = {dof_ac}   "
      f"κατώφλι = {CHI2_95[dof_ac]}")
detected = est["J"] > CHI2_95[dof_ac]
print("Αποτέλεσμα:", "ανιχνεύονται" if detected else "δεν ανιχνεύονται",
      "εσφαλμένα δεδομένα")
print()
for nm, rr, rn in zip(names_ac, est["r"], est["rN"]):
    print(f"{nm:>5}:  r = {rr:+.4f}   r^N = {rn:.2f}")
```

### Χονδροειδές σφάλμα σε μέτρηση άεργης ισχύος

Στη μέτρηση ροής άεργης ισχύος $Q_{12}$ προστίθεται σφάλμα $12\sigma$. Η εκτίμηση επαναλαμβάνεται,
η μέτρηση με το μέγιστο κανονικοποιημένο υπόλοιπο απορρίπτεται και η εκτίμηση υπολογίζεται εκ
νέου με τις υπόλοιπες έντεκα μετρήσεις.

```{code-cell} python
k_bad = names_ac.index("Q12")
z_bad_ac = z_ac.copy()
z_bad_ac[k_bad] += 12 * sig_ac[k_bad]

est_b = wls_ac(z_bad_ac, sig_ac, MEAS_AC, Y_ac, PI_ac, x_flat)
order = np.argsort(est_b["rN"])[::-1]
print(f"J(x̂) = {est_b['J']:.1f}  έναντι κατωφλίου {CHI2_95[dof_ac]}")
print("μεγαλύτερα κανονικοποιημένα υπόλοιπα:",
      ", ".join(f"{names_ac[q]} ({est_b['rN'][q]:.2f})" for q in order[:3]))

keep = [q for q in range(m_ac) if q != order[0]]
meas_keep = [MEAS_AC[q] for q in keep]
est_c = wls_ac(z_bad_ac[keep], sig_ac[keep], meas_keep, Y_ac, PI_ac, x_flat)
dof_c = len(keep) - n_ac
print()
print(f"μετά την απόρριψη της {names_ac[order[0]]}:")
print(f"  J(x̂) = {est_c['J']:.3f}   βαθμοί ελευθερίας = {dof_c}   "
      f"κατώφλι = {CHI2_95[dof_c]}")
e_c = np.abs(est_c["x"] - x_ref)
print(f"  μέγιστο σφάλμα εκτίμησης: {np.degrees(e_c[:NB - 1].max()):.3f}° "
      f"στις γωνίες, {e_c[NB - 1:].max():.4f} p.u. στις τάσεις")

fig, ax = plt.subplots(figsize=(6.2, 3.0))
xpos = np.arange(m_ac)
ax.bar(xpos - 0.2, est["rN"], width=0.4, color="#5b6b78",
       label="χωρίς χονδροειδές σφάλμα")
ax.bar(xpos + 0.2, est_b["rN"], width=0.4, color="#8f2622",
       label="σφάλμα 12σ στη Q12")
ax.axhline(3.0, ls="--", lw=1.1, color="#1f2a33")
ax.set_xticks(xpos, names_ac)
ax.set_ylabel("κανονικοποιημένο υπόλοιπο")
ax.legend(fontsize=9, frameon=False)
fig.tight_layout()
```

Ο έλεγχος $\chi^2$ ανιχνεύει το σφάλμα και το μέγιστο κανονικοποιημένο υπόλοιπο αντιστοιχεί
στη $Q_{12}$. Το δεύτερο μεγαλύτερο κανονικοποιημένο υπόλοιπο εμφανίζεται στην έγχυση $Q_2$,
της οποίας το υπόλοιπο έχει συντελεστή συσχέτισης 0,81 με εκείνο της $Q_{12}$, επειδή η έγχυση στον
ζυγό 2 περιλαμβάνει τη ροή άεργης ισχύος του κλάδου 1–2. Μετά την απόρριψη της $Q_{12}$ ο
έλεγχος $\chi^2$ δεν ανιχνεύει πλέον εσφαλμένα δεδομένα.

### Σύγκριση με το γραμμικό μοντέλο ΣΡ

Το κελί που ακολουθεί συγκρίνει τις γωνίες που εκτιμώνται με το γραμμικό μοντέλο ΣΡ, από τις
μετρήσεις ενεργής ισχύος, με εκείνες του μη γραμμικού εκτιμητή. Οι μετρήσεις παράγονται χωρίς
θόρυβο, ώστε η απόκλιση να οφείλεται αποκλειστικά στο σφάλμα του μοντέλου. Οι παράμετροι είναι
ο λόγος $r/x$ όλων των κλάδων, ο συντελεστής φόρτισης, με τον οποίο πολλαπλασιάζονται οι γωνίες
της κατάστασης αναφοράς, και η επιλογή μοναδιαίων μέτρων τάσεων, η οποία απομονώνει την
επίδραση του λόγου $r/x$ και της φόρτισης από εκείνη του προφίλ τάσεων.

```{code-cell} python
:tags: [hide-input]
import os

# Κατά την παραγωγή του βιβλίου (BOOK_BUILD=1) το κελί εκτελείται χωρίς
# ολισθητές, ώστε το PDF να περιέχει το γράφημα. Στον φυλλομετρητή το
# ipywidgets εγκαθίσταται κατά την εκτέλεση.
OLISTHITES = os.environ.get("BOOK_BUILD") != "1"
if OLISTHITES:
    try:
        import ipywidgets as widgets
    except ModuleNotFoundError:
        try:
            import micropip
            await micropip.install("ipywidgets")
            import ipywidgets as widgets
        except Exception as exc:
            OLISTHITES = False
            print("Οι ολισθητές δεν είναι διαθέσιμοι:", exc)


def dc_angles(z, sig, meas, kladoi):
    """Εκτίμηση γωνιών με το μοντέλο ΣΡ από τις μετρήσεις ενεργής ισχύος."""
    rows = [k for k, (kind, _, _) in enumerate(meas)
            if kind in ("Pf", "Pi")]
    Hd = np.zeros((len(rows), NB))
    for row, k in enumerate(rows):
        kind, i, j = meas[k]
        for a, b, _, x, _ in kladoi:
            a, b = a - 1, b - 1
            if kind == "Pf" and {a, b} == {i, j}:
                Hd[row, i] += 1 / x
                Hd[row, j] -= 1 / x
            elif kind == "Pi" and i in (a, b):
                other = b if a == i else a
                Hd[row, i] += 1 / x
                Hd[row, other] -= 1 / x
    Hd, W = Hd[:, 1:], np.diag(1 / sig[rows]**2)
    return np.linalg.solve(Hd.T @ W @ Hd, Hd.T @ W @ z[rows])


def sygkrisi_dc(r_pros_x=0.1, fortisi=1.0, monadiaies_taseis=False):
    kladoi = [(i, j, r_pros_x * x, x, b) for i, j, _, x, b in KLADOI]
    Y, pi = build_ybus(kladoi)
    V = np.ones(NB) if monadiaies_taseis else V_ref
    x_true = np.r_[th_ref[1:] * fortisi, V]
    z = h_ac(x_true, MEAS_AC, Y, pi)                          # χωρίς θόρυβο

    th_dc = dc_angles(z, sig_ac, MEAS_AC, kladoi)
    res = wls_ac(z, sig_ac, MEAS_AC, Y, pi, x_flat)
    err_dc = (th_dc - x_true[:NB - 1]) * 1e3                  # mrad
    err_ac = (res["x"][:NB - 1] - x_true[:NB - 1]) * 1e3

    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    pos = np.arange(NB - 1)
    ax.bar(pos - 0.18, err_dc, width=0.36, color="#5b6b78",
           label="μοντέλο ΣΡ")
    ax.bar(pos + 0.18, err_ac, width=0.36, color="#8f2622",
           label="μη γραμμικός εκτιμητής")
    ax.axhline(0.0, lw=0.8, color="#1f2a33")
    ax.set_xticks(pos, ["θ2", "θ3"])
    ax.set_ylabel("σφάλμα γωνίας (mrad)")
    ax.set_title(f"r/x = {r_pros_x:.2f}   φόρτιση × {fortisi:.2f}",
                 fontsize=10)
    ax.legend(fontsize=9, frameon=False, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    plt.show()

    e_dc, e_ac = np.abs(err_dc).max(), np.abs(err_ac).max()
    print(f"μέγιστο σφάλμα γωνίας, μοντέλο ΣΡ:            {e_dc:.2f} mrad "
          f"({np.degrees(e_dc / 1e3):.3f}°)")
    print(f"μέγιστο σφάλμα γωνίας, μη γραμμικός εκτιμητής: {e_ac:.1e} mrad "
          f"({len(res['history'])} επαναλήψεις)")


if OLISTHITES:
    stil = {"description_width": "185px"}
    diataxi = widgets.Layout(width="520px")
    _ = widgets.interact(
        sygkrisi_dc,
        r_pros_x=widgets.FloatSlider(
            value=0.1, min=0.0, max=0.5, step=0.05,
            description="λόγος r/x",
            readout_format=".2f", style=stil, layout=diataxi),
        fortisi=widgets.FloatSlider(
            value=1.0, min=0.5, max=3.0, step=0.25,
            description="συντελεστής φόρτισης",
            readout_format=".2f", style=stil, layout=diataxi),
        monadiaies_taseis=widgets.Checkbox(
            value=False, description="μοναδιαία μέτρα τάσεων",
            style=stil, layout=diataxi),
    )
else:
    sygkrisi_dc()
```

Ο μη γραμμικός εκτιμητής αναπαράγει την κατάσταση αναφοράς με σφάλμα στρογγυλοποίησης,
ανεξάρτητα από τον λόγο $r/x$ και τη φόρτιση. Με μοναδιαία μέτρα τάσεων, το σφάλμα του
μοντέλου ΣΡ αυξάνεται με τον λόγο $r/x$ και με τη φόρτιση: για $r/x$ = 0,3 και διπλάσια
φόρτιση υπερβαίνει τα 12 mrad (0,7°). Με μη μοναδιαίο προφίλ τάσεων προστίθεται το σφάλμα της
παραδοχής $V_i = 1$, το οποίο, ανάλογα με τις συνθήκες, ενισχύει ή αντισταθμίζει εν μέρει το
σφάλμα λόγω του λόγου $r/x$.

### Προτεινόμενες ασκήσεις (μοντέλο ΕΡ)

- **Περιοχή σύγκλισης**: εκκίνηση του `wls_ac` από αρχικές γωνίες ±0,5 rad ή από μέτρα
  τάσεων 0,8 p.u. Πόσες επαναλήψεις απαιτούνται και σε ποιες περιπτώσεις η μέθοδος αποκλίνει;
- **Παρατηρησιμότητα**: αφαίρεση των δύο μετρήσεων τάσης. Παραμένει ο $\mathbf{G}$
  αντιστρέψιμος και πώς μεταβάλλεται η ακρίβεια των εκτιμώμενων μέτρων τάσεων;

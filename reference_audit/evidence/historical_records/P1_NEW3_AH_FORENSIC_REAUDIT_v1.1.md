# EvoFungi-AMR Paper 1 — NEW-3A–H bağımsız yeniden denetim (v1.1)

**Denetim tarihi:** 13 Ekim 2026  
**Durum:** Bağımsız kod + sayısal reproduksiyon + özgün makale/metot taraması; fenotip veya yeni Evo2 çıkarımı yok.  
**Eski kayıtlar:** `P1_NEW3_AH_INTEGRATED_SCIENTIFIC_REVIEW_v1.0.md` ve dondurulmuş kaynaklar değiştirilmemiştir. Bu metin **v1.1 düzeltme ve denetim eki**dir, sonuç seçerek yeniden analiz değildir.

## A. Denetimin neyi kapsadığı / kapsamadığı

**Doğrudan incelenen kaynaklar:** Google Drive'daki 824 dosyalık toplama envanteri; güncel `EVO2FungiCongress` notebook'u (286 hücre) ve toplama anındaki notebook snapshot'ı; orijinal NEW-3A/B/C/D/E/F/G/H logları; H0B dondurması, H2/H3/H4 dosyaları; G1A ve F2A teknik denetimleri; Durand ve ark. 2026 özgün DMS makalesi (`iyag055`, DOI 10.1093/genetics/iyag055); Brixi ve ark. 2026 Evo2 makalesi (DOI 10.1038/s41586-026-10176-5); dos Reis ve ark. (2004) tAI çalışması (DOI 10.1093/nar/gkh834); güncel CRAN tAI 0.2.2; DMS metodoloji kaynakları. Farklı, hakemden geçmemiş Evo2 benchmark ön baskısı ayrıca literatür sınırlaması olarak not edilmiştir.

**Denetimin kesin sınırı:** H2'de 48 LOPO sonucu bağımsız ağırlıklı en küçük kareler uygulamasıyla tekrar hesaplanmıştır. 518 hosted Evo2 API isteği **yeniden yapılmamıştır**; checkpoint doğrulaması mevcut arşiv/sentinel kayıtlarına dayanır. H3'ün bütün bootstrap regresyonları bu denetimde yeniden çalıştırılmamıştır; önceden NPZ yeniden özetleme ve seçilmiş fitler için yapılan kontrol mevcuttur. 824 dosyanın her birinin baytları bu oturumda yeniden alınmamıştır; toplama envanteri tam kapsama kanıtı değildir. Bu sınırlar bir bilimsel QC başarısızlığı olarak gizlenmemelidir.

## B. Hatalı hücrelerden doğru uygulamaya — doğrulanan kronoloji

| Sıra | Kaynak/hücre | Ne yanlış veya eksikti? | Geçerli sonuç / karar |
|---|---|---|---|
| NEW-3A, 29 Eylül | A modül giriş dondurması | 589/435/28 giriş evreni; 278 grup-içi analiz henüz yeniden türetilmemişti | Fenotip ve Evo2-atıf ilişkisi görülmeden biyolojik G0/G1/G2 hiyerarşisi planlandı; G3 ancak kaynak erişimi sağlanırsa ikincildi |
| NEW-3B | B1B/B2B/B3/B3R1 | İlk HS2 metadata/koordinat anomalileri; T0/FINAL589 ayrı seçilim aşamaları gibi okunma riski | Nükleotit dizisi temelli kanonik pozisyon; FINAL589 = T0-HC tam küme; biyolojik ölüm iddiası yok |
| NEW-3C | Notebook 237–245 | CoCoPUTs endpoint/HTTP 206 sorunu; G1 kaynak satırında 76/77 şema uyuşmazlığı | Son R5E kaynak şema düzeltmesi ve C2B byte freeze; başarısız girişimler tarihsel olarak korunur |
| NEW-3D | 246–250 | G1/G2 kaynak ve özellik tablosu ön çalışmaları | 278/589/814 kesin-genotip manifestleri; özellikler Evo2 sonucu ve deneysel fenotip kullanılmadan donduruldu |
| NEW-3E | 251–254 | 28 pozisyon varsayımı başarısız; sonra JSON bool serileştirme hatası | Hücre 253 ile geçerli 27 LOPO pozisyonu, 1360 yapısal olarak gözlenen 278'de değerlendirilemez; kategorik gösterge kodlamaları donduruldu |
| NEW-3F1 | 255–256 | İlk gözlenen Evo2 skoru ilişkilendirme | 278 genotip/124 AA grubu, 27 pozisyon; LOPO M0–M3 sonuçları; bu aşamada yalnızca nokta tahmini |
| NEW-3F2 | 257–259 | İlk F2 kodu, manifestte doğrudan bulunmayan dummy sütunlarını bekleyerek STOP verdi | Ham kategorilerden dondurulmuş dummy sütunları türetildi; ilk F1 sonuçları aynı çıktı |
| NEW-3F2A/B | 260–261 | Primer bootstrap'ta 2.479/10.000 tekrar yapısal olarak hesaplanamaz; yok sayma veya regularizasyon uygun değil | Orijinal CI'lar sadece 7.521 hesaplanabilir tekrar için geçerli; sonradan destek-koruyan Dirichlet duyarlılığı ayrıca donduruldu |
| NEW-3F3 | 262 | LOHO primer CV olarak yorumlanma riski | Yalnızca ikincil, üç-hotspot betimleyici stres testi |
| NEW-3G0/G1/G1A/G2 | 263–266 | Bileşen Q² değerlerini toplama ve biyolojik bilgi yüzdesine dönüştürme riski; G1/G2 örtüşmesi | Skor düzeyinde bileşen toplamsallığı doğrulandı; Q² toplamsal değil; varyans, kovaryans ve G0'a koşullu ilişkiler denetlendi |
| NEW-3H0–H1B2 | 267–279 | Yanlış tarihsel F639 rekonstrüksiyonu ile sahte PASS; backend tercihi ve base64 katı-dekod hatası | Sekans/hash sentinel; hosted Evo2-7B aynı API ortamı; 28 REF + 490 ALT; immutable checkpoint ağırlık SHA **UNVERIFIED** |
| NEW-3H2–H4B | 280–283 | H4A kodunda tüm aralıklar sıfır içermeli şeklinde sonucu koşullayan hatalı `assert` | H2/H3 sayıları tutuldu; bilimsel kalite koşulu H4B'de kaldırılıp tarihsel hata açıkça belgelendi |
| NEW-3 toplama | 284 | Dosyalar çeşitli klasörlere yayılmış; ilk paketleme adayları kapsam sorunu | 824 dosyalık kopya ve envanter üretildi; tüm olası kaynakların bulunduğunu tek başına kanıtlamaz |

**Notebook snapshot bulgusu:** `EVO2FungiCongress.ipynb` toplanan snapshot 13.660.809 bayt; daha güncel Drive `EVO2FungiCongress` 13.662.591 bayt. İki notebook'un **286/286 hücre kaynak kodu birebir aynı**; farklılaşan tek hücre 284'ün kayıtlı çıktı/çalıştırma durumudur. Bu fark hesaplama tanımının değişmesi değil, notebook çıktısı/metadata kaydıdır. Eski denemeler (251, 252, 258, 274 vb.) geçmişle birlikte mevcut.

## C. Kaynak evreni ve ascertainment — biyolojik değil, deneysel/kaynak destek koşulu

- Teorik hedef NNK missense evreni: **814 kesin DNA genotipi, 532 AA ikamesi, 28 fiziksel pozisyon**. Bunun **208 AA grubu** çok-kodonludur ve toplam **490 kesin genotipi** oluşturur; 490, 814'ün alternatifi değil alt evrendir.
- PRE yüksek güven kümesi: **811/814**. T0 yüksek güven kümesi: **589/814**, ve `FINAL589` ile kesin genotip anahtarları birebir aynı. HS2'de T0-HC desteği **125/232**, HS1'de **231/258**, HS3'te **233/324**.
- Gözlenen **124 çok-kodon AA grubu, 278 genotip, 27 değerlendirilebilir pozisyon**. Teorik 208 çok-kodon AA grubundan **61**'i T0-HC'de tek kodona düşmüş, **23**'ü tümüyle destek dışı kalmış, **124**'ü çok-kodon karşılaştırması olarak kalmış.
- Durand ve ark. orijinal DMS çalışmasında yalnızca bütün tekrarlarında >10 okuma ile gözlenen varyantlar `high confidence` kabul edilmiştir. Kütüphane/PCR/okuma derinliği/klon çoğalması gibi faktörler destek-kayıp örüntülerine katkıda bulunabilir; T0-HC kaybı canlılık/ölüm ile özdeşleştirilmemelidir.
- 1360'ın gözlenen 278'de olmaması Evo2 sonucuna göre seçilim değil, **aynı AA grubunda en az iki alternatif kodon bulunmaması** kaynaklı; tam teorik490'da 1360 yer alır.
- **Kritik nüans:** Her grup aynı **mutant amino asidi** kodlayan farklı **mutant DNA kodonları** ile ilgilidir. Bu bireysel WT→mutant değişimlerinin her biri genel olarak **nonsynonymous**'tır; çalışma, WT amino asidi değiştirmeyen klasik synonymous-mutasyon DMS'i değildir.

## D. Bağımsız 48/48 H2 hesaplama tekrar testi

Aşağıdaki yeniden uygulama, dondurulmuş H0B üyelik ve M0–M3 özellikleri ile ALT490 skorlarından bağımsız WLS + `np.linalg.lstsq` kullanmıştır. Her test fiziki pozisyonundaki bütün varyantlar eğitim dışında bırakılmıştır; AA-grup merkezleme, ağırlık `1/n_g`, OOF Q² = 1 − ağırlıklı SSE/ağırlıklı SST. M1–M3 için özellik tanımları değiştirilmemiştir.

**Bağımsız karşılaştırma:** 4 panel × 3 bileşen × 4 model = **48/48**; H2 kaydedilmiş değerleriyle **maksimum mutlak fark `6.245004513516506e-16`**. Bütün pozisyon-dışarıda-bırakmalı tasarım matrisleri tam ranklı. Panel büyüklükleri: OBSERVED278 (278,124,27), COMMON124_COMPLETE (309,124,27), THEORETICAL490 (490,208,28), THEORETICAL_HISTORICAL27 (474,201,27).

OBSERVED278 primer total: **M0 0.528669, M1 0.550580, M2 0.527755, M3 0.547950**. M0 bileşenleri: **codon-local 0.624575**, **downstream50 0.047716**. Bunlar modelin **kendi olasılık skorlarına ait AA-grup-içi kontrastları** için açıklama performansıdır, gözlenen antifungal direnç veya biyolojik kısıt R²'si değildir.

**Bileşen ölçeğini de bağımsız doğruladım.** OBSERVED278 için AA-içi ve aynı ağırlıklarla: SST_codon ≈ 125.07473550877, SST_downstream ≈ 19.04063555078, `2×cross-product` ≈ −2.64840339709 ve SST_total ≈ 141.46696766246. Kodon ve downstream hedefleri farklı varyans ölçeklerinde; toplam hedef SST aritmetik toplamlarından kovaryans nedeniyle ayrılıyor. Toplam **skor** codon-local+downstream50 olsa bile Q²'ler toplanamaz veya birbirinden çıkarılıp biyolojik yüzdeler elde edilemez.

## E. Uncertainty, negatif sonuçlar ve epistemik sınırlar

- Gözlenen F2 primer ΔQ²(M1−M0): **+0.0219109**, bootstrap **10.000 planlı, 7.521 geçerli, 2.479 yapısal rank kaybı**. Geçerli tekrarlar koşuluyla percentile aralığı **[−0.026040, +0.084471]**. Nominal 10.000 tekrarın kayıp olduğu bir durumda bu, eksiksiz bootstrap kapsamasını kanıtlamaz.
- F2B'nin sonradan yapılmış destek-koruyan Dirichlet duyarlılık aralığı: **[−0.017645, +0.076084]**, **klasik örnekleme güven aralığı değil**. F2 yerine geçirilmez.
- G1A: G1/G2'nin ağırlıklı grup-içi korelasyonu ≈ **0.81297**; G0'a koşullu korelasyon ≈ **0.52859**. G0'ın G1 proxy varyasyonuna yönelik betimleyici R²'si ≈ **0.87876**, G2 için ≈ **0.70875**. Dolayısıyla kullanım/tAI proxy'leri geometriden bağımsız değil.
- H3, 4 panelin her birinde G1 artı ΔQ² nokta tahmini: **+0.021911 / +0.022992 / +0.018094 / +0.021001**; karşılık gelen destek-koruyan duyarlılık aralıklarının **tümü sıfır içeriyor**. 60/60 karşıtlık ve 27/27 eşleştirilmiş karşıtlık sıfır içeriyor; **bu bir hesap sonucu, PASS koşulu veya eşdeğerlik kanıtı değil**.
- Aynı 124 AA grubunda kodon tamamlama: ΔΔQ² = **+0.001081**, eşleştirilmiş duyarlılık aralığı **[−0.023085,+0.024482]**. Etkisizlik/eşdeğerlik kanıtlanmadı. 278→490 karşılaştırması grup/hotspot/pozisyon desteğini de değiştirir; nedensel ascertainment katsayısı sayılamaz.
- F3 LOHO sonuçları **yalnızca üç hotspot folduna** dayanıyor; başka gen veya türlere dış doğrulama değil.
- H4A'daki `interval contains zero` kuralı bilimsel olarak yanlıştı: bir sonuç sıfırı dışladığı için QC başarısız sayılmaz. H4B düzeltmesi doğru; tarihsel hata korunmalı.

## F. Temel bilimsel metodoloji denetimi

1. **Neyin öğrenildiği:** G0, hedeflenen WT→mutant DNA düzenleme geometrisinin aynı genotip üzerinden çıkarılan özellikleri. Autoregressive Evo2'nin kodon-local likelihood'ıyla ilişkisinin bir kısmı bu skorun **mimarisinden ve tanımından** beklenir. Biyolojik/evolüsyonel kısıt miktarı olarak yorumlanamaz.
2. **Aynı AA içinde kıyaslama:** Bir grupta mutant protein çıktısı sabit tutuluyor; bu protein-AA etkisinin basit karşılaştırmaya karışmasını azaltır. Ancak edit sayısı, GC, başlangıç konumu ve alternatif kodonun nükleotitleri değişir. Ayrı mekanistik translasyon testi yapılmadı.
3. **Veri sızıntısı:** LOPO, aynı fiziksel pozisyondaki hedef skorların modelin eğitimine girmesini engeller. Grup-içi merkezleme test kümesindeki **gözlenen skorların kendi grup ortalamasından** uzaklığını hesaplayan *kontrast hedefi* tanımlar; bu mutlak genotip likelihood'unun doğrudan sahada tahmini değildir. NEW-3'ün hedefi bu karşılaştırma olduğundan kabul edilebilir fakat kapsam açık yazılmalıdır.
4. **Mutant kodonlar ve fenotip:** Durand'ın orijinal hesap hattı amino asit dizisi düzeyinde seçim katsayısı üretir. `490 Evo2 scores` ≠ `490 experimental fitness measurements`. Aynı AA fenotipini birden fazla kodona kopyalayıp bağımsız veri saymak **pseudoreplication** olur. NEW-4'te özellikle engellenmeli.
5. **G1/G2 biyolojisi:** CoCoPUTs S288C genomik kodon kullanım sıklığı ve GtRNAdb S288C tRNA gen-kopyasına dayanan tAI, gerçek BY4741 DMS kültüründe ölçülmüş tRNA miktarı, kodon çözülme süresi veya protein katlanma hızı değildir. Bu değişkenlerin küçük ek ilişkisi bile mekanistik translasyon öğrenme kanıtı değildir.
6. **Autoregressive yön:** 2053-nt coding orientation ve ALT−REF local53 dondurulmuştur; sign sonucu iyileştirmek için değiştirilmez. Hosted 7B immutable ağırlık hash'i doğrulanmadığı için checkpoint kimliği **UNVERIFIED** olarak kalır.
7. **Orijinal Evo2 makalesiyle ilişki:** Brixi ve ark. (Nature, 2026), Evo2'nin kodon verimliliği ve triplet periodisitesi ile ilişkili likelihood örüntülerini zaten bildirmiştir; bu genel fenomen yeniymiş gibi sunulamaz. Bizim sınamamız, *aynı mutant AA ve G0 geometrisi koşulunda*, Fks1 hedef-NNK genotiplerindeki ek G1/G2 açıklamasıdır. Nature'da marginal örüntü raporlanmış olması ile bizde kontrollü ek Q²'nin belirsiz olması birbirini mantıksal olarak dışlamaz.
8. **Yayındaki biyolojik yer:** Bu NEW-3 modülünü model-davranışı/teknik skor-denetimi olarak sunmak savunulabilir. Makalenin asıl biyolojik iddiası (kontrol fitness vs ilaç seçilim farkı) NEW-4 ve önceki dondurulmuş fenotip benchmark modülleriyle desteklenmelidir; NEW-3 tek başına direnç sonucu değildir.

## G. Kaynak bütünlüğü ve açık iş

- Envanter: **824 kayıt, 137.732.157 bayt**, tekrar eden hedef anahtar yok, `MISSING_CORE_PATHS`/`UNRESOLVED_REFERENCED_PATHS` boş. Bu, **koddaki mevcut tarama kuralları** içinde sıfır eksik demek; bütün olası kaynakların mutlak eksiksizliği demek değil.
- Bu çalışma oturumunda eldeki yerel dosyalar ile envanter kaydı karşılaştırıldığında **54 eşsiz isim** eşleşti; tek görünen `RUN_MANIFEST.json` isim çakışması, **iki farklı kaynak klasöründe aynı dosya ismi** bulunduğu içindi, veri bozulması olarak yorumlanamaz.
- Notebook: Toplanan 13.660.809-bayt snapshot ile güncel 13.662.591-bayt sürümün **286/286 kod hücresi birebir aynı**. Notebook çıktı farkı hücre 284; mevcut kaynak için kod revizyonu yok.
- Önceki entegrasyon v1.0'da sözel **"geometry explains a substantial fraction"** ifadesi yanlış anlaşılmaya yatkın. Öneri: *"predefined edit geometry achieves appreciable held-out-position Q² for within-AA Evo2 codon-score contrasts, without establishing independent biological constraint"*.
- Önceki NEW-3 giriş dondurmasında **NEW-3I: tam-genetik-kod ve WT-synonymous sentetik duyarlılık** fikirleri yer aldı; G0 bunların kapsamının ayrıca denetlenmesi gerektiğini belirtti. **Tamamlandığına dair kayıt görülmedi**. Bunun Paper 1 için zorunlu olup olmadığı NEW-4'ten önce belgeli bir *scope decision* olarak çözümlenmeli; sessizce bitmiş veya zorunluymuş gibi ilan edilmemeli.
- Orijinal tüm 518 hosted API isteği bu denetimde yeniden yapılmadı; H1B2 sentinel/checkpoint zinciri ve tarihsel özdeş skorlarla uyumlu, ama immutable ağırlık kimliği çözülemedi.

## H. Makalede kullanılabilecek daha dar Results özeti (taslak)

> Within 124 Fks1 amino-acid-substitution groups represented by multiple exact mutant codons, predefined nucleotide-edit geometry achieved a pooled leave-one-position-out Q² of 0.529 for within-group contrasts of Evo 2 local53 likelihood differences. The corresponding Q² was 0.625 for the codon-local component and 0.048 for the downstream50 component, which differ in target variance and should not be interpreted as additive biological contributions. Adding the genomic codon-usage proxy yielded a small incremental Q² (+0.0219), but both estimability-qualified conventional bootstrap results and subsequent support-preserving sensitivity intervals included zero. A phenotype-blind complete targeted-NNK interrogation produced a similar qualitative pattern without establishing robust independent codon-usage or tRNA-related attribution. These analyses describe the genomic model's own score preferences; they do not quantify experimentally measured fitness or echinocandin resistance.

## I. Referanslar — özgün/güncel kaynaklar

1. Durand R. ve ark. (2026). *Mutational landscape and molecular bases of echinocandin resistance in Saccharomyces cerevisiae*. Genetics 233(1), iyag055. DOI: https://doi.org/10.1093/genetics/iyag055 ; özgün analiz deposu https://github.com/Landrylab/Durand_et_al_2026
2. Brixi G. ve ark. (2026). *Genome modelling and design across all domains of life with Evo 2*. Nature. DOI: https://doi.org/10.1038/s41586-026-10176-5 ; resmi repo https://github.com/ArcInstitute/evo2
3. dos Reis M, Savva R, Wernisch L. (2004). *Solving the riddle of codon usage preferences: a test for translational selection*. Nucleic Acids Research. DOI: https://doi.org/10.1093/nar/gkh834
4. CRAN tAI 0.2.2: https://stat.ethz.ch/CRAN/web/packages/tAI/index.html
5. Liu Y, Yang Q, Zhao F. (2021). *Synonymous but Not Silent: The Codon Usage Code for Gene Expression and Protein Folding*. Annual Review of Biochemistry. DOI: https://doi.org/10.1146/annurev-biochem-071320-112701
6. *DiMSum: an error model and pipeline for analyzing deep mutational scanning data and diagnosing common experimental pathologies*. https://pubmed.ncbi.nlm.nih.gov/32799905/
7. *A systematic analysis of regression models for protein engineering*. PLOS Computational Biology (2024). https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012061
8. Mathur V, Sachidanandam R. (2026). *Benchmarking DNA Foundation Models: Biological Blind Spots in Evo2 Variant-Effect Prediction*. **bioRxiv ön baskı, hakemli sonuç kabul edilmemeli**. https://www.biorxiv.org/content/10.64898/2026.03.10.710786v1 — abstract tarandı; ayrıntılı yöntem/kanıt doğrulaması tamamlanmadı.

## J. Sonraki adımlar — karara bağlanacak, şimdilik fenotip açılmayacak

1. NEW-3I'nin tarihsel rolünü ve Paper 1 açısından gerekli/ikincil/ertelenebilir statüsünü özgün A/D/G sözleşmelerine dayanarak belgelemek; sonuçlara bakarak zorunluluk değiştirmemek.
2. NEW-3A–H'nin bu v1.1 iddia düzeltmelerini *frozen output* üzerine yazmadan sürümlü arşivlemek. İki notebook snapshot'ının aynı kodu taşıdığını ve orijinal başarısız hücrelerin tarihsel saklandığını kayda geçirmek.
3. NEW-4'ten önce deneysel fenotip kaynağının tam olarak AA/genotip/tekrar düzeyini, koşul tanımlarını ve yüksek güven filtrelerini tekrar kilitlemek; aynı AA fenotipini kodon satırlarına kopyalayıp bağımsız ölçüm saymayı yasaklamak.
4. NEW-4 protokolü önceden dondurulduktan sonra kontrol-genel fitness ile üç ilaç koşulundaki kontrol-ayarlı seçilim etkilerini ayrı ayrı değerlendirmek; tekrar/pozisyon/assay bağımlılığını korumak.

**Son karar:** NEW-3H teknik kapanışı korunabilir; bütün NEW-3'ün biyolojik iddiaları v1.1'deki daha dar çerçeveyle sunulmalı. NEW-4'e fenotip kaynak düzeyi ve NEW-3I kapsamı kararı verilmeden geçilmemeli.

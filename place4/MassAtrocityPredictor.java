/*
 * Copyright 2012-2013 TopCoder, Inc.

 *
 * This code was developed under U.S. government contract NNH10CD71C.

 *
 * Licensed under the Apache License, Version 2.0 (the "License");

 * You may not use this file except in compliance with the License.

 * You may obtain a copy of the License at:
 *     http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software

 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

 * See the License for the specific language governing permissions and

 * limitations under the License.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;


public class MassAtrocityPredictor {
    public static final int R = 3671;
    public static final int C = 254;
    private static final int WARMUP_FOR_STAT = 25*30;
    private static final int NEWS_BUFFER_SIZE = 300;
    public static int TUNE_W = 600;
    public static int TUNE_FREQ = 100;

    double wGlobalFix = 1;
    double wGlobalTuned = 1;
    double wCountryTuned = 1;
    double wRegionTuned = 1;
    double wMixSum = wGlobalFix + wGlobalTuned + wCountryTuned + wRegionTuned;

    boolean needsTuning = wGlobalTuned + wCountryTuned + wRegionTuned > 0;

    public static final int POISSON = 0;
    public static final int OWN = 1;
    public static final int STAT = 2;
    public static final int NEWS = 3;
    public static final int NEWS2 = 4;
    public static final int COUNTRIES = 5;
    public static final int POISSONWCOUNTRY = 6;
    public static final int FRIENDS = 7; // should be last!
    public static final int ALGO_CNT = FRIENDS+1;

    public static final int R_00 = 0;
    public static final int R_01 = 1;
    public static final int R_10 = 2;
    public static final int R_11 = 3;
    public static final int R_12 = 4;
    public static final int R_21 = 5;
    public static final int R_22 = 6;
    public static final int R_02 = 7;
    public static final int R_20 = 8;
    public static final int R_xx = 9;
    public static final int R_CNT = R_xx+1;

    public static final int N_AA = 0;
    public static final int N_AP = 1;
    public static final int N_AM = 2;
    public static final int N_PA = 3;
    public static final int N_PP = 4;
    public static final int N_PM = 5;
    public static final int N_MA = 6;
    public static final int N_MP = 7;
    public static final int N_MM = 8;
    public static final int N_H = 9;
    public static final int N_L = 10;
    public static final int N_xx = 11;
    public static final int N_CNT = N_xx+1;

    public double[] WS = new double[ALGO_CNT];
    public double WSUM;
    private double[] WS_TUNED;
    private double WSUM_TUNED;

    private double[][] WS_R = new double[R_CNT][ALGO_CNT];
    private double[] WSUM_R = new double[R_CNT];

    private int day0;
    public static int day0fakePredict = 13999 - TUNE_W;

    public Region[] regions;
    private Country[] countries;
    private List<Stat> stats;
    private boolean predictingStarted = false;
    private boolean fakePredictingStarted = false;

    public MassAtrocityPredictor() {
        setWeights();

        regions = new Region[R];
        for (int i = 0; i < R; i++) {
            Region r = new Region();
            regions[i] = r;
            r.id = i;
        }

        countries = new Country[C];
        for (int i = 0; i < C; i++) {
            Country c = new Country();
            countries[i] = c;
            c.id = i;
            c.ws = WS.clone();
            c.sumWeights();
        }
        stats = new Vector<>();

    }

    private void setWeights() {
        WS[POISSON] = 51.316;
        WS[OWN] = 38.554;
        WS[STAT] = 46.651;
        WS[NEWS] = 30;
        WS[NEWS2] = 35.049;
        WS[COUNTRIES] = 235.795;
        WS[POISSONWCOUNTRY] = 50;
        WS[FRIENDS] = 0;

//		WS[POISSON] = 0;
//		WS[OWN] = 0;
//		WS[STAT] = 0;
//		WS[NEWS] = 1;
//		WS[NEWS2] = 0;
//		WS[COUNTRIES] = 0;
//		WS[POISSONWCOUNTRY] = 0;
//		WS[FRIENDS] = 0;

        WSUM = 0;
        for (int i = 0; i < ALGO_CNT-1; i++) WSUM += WS[i];

        WS_TUNED = WS.clone();
        WSUM_TUNED = WSUM;
    }

    public int receiveData(int dataSourceId, int day, String[] data) {
        if (dataSourceId == 1) { // news
            if (data == null) return 0;
            int dayIndex = (day-day0) % NEWS_BUFFER_SIZE;
            for (Region r: regions) {
                for (int i = 0; i < 21; i++) r.newsPerDay[dayIndex][i] = 0;
            }
            for (Country c: countries) {
                for (int i = 0; i < 21; i++) c.newsPerDay[dayIndex][i] = 0;
            }
            for (String s: data) {
                String[] parts = s.split(" ");
                char type = parts[10].charAt(0);
                Set<Integer> rs = new HashSet<>();
                if (!parts[4].equals("_")) rs.add(Integer.parseInt(parts[4]));
                if (!parts[9].equals("_")) rs.add(Integer.parseInt(parts[9]));
                if (!parts[14].equals("_")) rs.add(Integer.parseInt(parts[14]));
                for (int ri: rs) {
                    Region r = regions[ri];
                    r.newsPerDay[dayIndex][20]++; //total
                    r.newsPerDay[dayIndex][type - 'a']++;
                    Country c = r.country;
                    c.newsPerDay[dayIndex][20]++;
                    c.newsPerDay[dayIndex][type - 'a']++;
                }
            }
            // check short overflows
            for (Region r: regions) {
                for (int i = 0; i < 21; i++) {
                    if (r.newsPerDay[dayIndex][i] < 0) r.newsPerDay[dayIndex][i] = Short.MAX_VALUE;
                }
            }
            return 0;
        }
        if (dataSourceId == 2) { // geo
            day0 = day;
            for (int i = 0; i < data.length; i++) {
                String line = data[i];
                if (line.length() > 100) line = line.substring(0, 100);
                String[] parts = line.split(" ");
                if (parts.length >= 2) {
                    try {
                        int ci = Integer.parseInt(parts[0]);
                        i++;
                        String line2 = data[i];
                        String[] parts2 = line2.split(" ");
                        int ri = Integer.parseInt(parts2[0]);
                        Region r = regions[ri];
                        Country c = countries[ci];
                        c.rs.add(r);
                        r.ci = ci;
                        r.country = c;
                    }
                    catch (Exception e) {
                        // ignore
                    }
                }
            }
            return 0;
        }

        // atrocity
        int dScaled = day - day0;
        if (data != null) {
            for (String s: data) {
                String[] parts = s.split(" ");
                int ci = Integer.parseInt(parts[2]);
                int ri = Integer.parseInt(parts[3]);
                Region r = regions[ri];
                int cnt = r.as.size();

                if (cnt != 0 && r.as.get(cnt-1).day + 3 > day) {
                    continue;
                }

                A a = new A();
                a.day = day;
                a.region = ri;
                r.as.add(a);
                cnt++;

                int prevD = 0;
                if (cnt > 1) prevD = r.as.get(cnt-2).day;
                short prevCnt = (short)(cnt-1);
                for (int d = Math.max(0, prevD-day0);  d < day-day0; d++) {
                    r.cCnts[d] = prevCnt;
                    r.wghTs[d] = (short) (d + day0 - prevD);
                }
                r.cCnts[day-day0] = (short)cnt;
                r.wghTs[day-day0] = (short)0;

                Country c = countries[ci];
                c.as.add(a);

                if (day > day0 + WARMUP_FOR_STAT) {
                    if (cnt > 1) {
                        Stat stat = new Stat(r, true);
                        stats.add(stat);
                    }

                    Stat stat = new Stat(r, false);
                    r.lastFact = stat;
                }

                for (int i = Math.max(0, dScaled - 30); i < dScaled; i++) r.goodGuessPerDay[i] = 1;
            }
        }

        if (dScaled > 0) {
            for (Region r: regions) {
                int cnt = r.as.size();
                if (cnt > 0 && r.as.get(cnt-1).day == day) {
                    //already bookkept
                }
                else {
                    r.cCnts[dScaled] = r.cCnts[dScaled-1];
                    r.wghTs[dScaled] = (short) (r.wghTs[dScaled-1] + 1);
                }
            }
        }

        if (!predictingStarted && day >= day0fakePredict && needsTuning) {
            predict(day);
        }

        return 0;
    }

    public void receiveNews(int day, short[][] data) {
        for (short[] row: data) {
            Region r = regions[row[0]];
            for (int i = 0; i < 21; i++) {
                r.newsPerDay[(day-day0) % NEWS_BUFFER_SIZE][i] = row[i+1];
            }
        }
    }


    public double[] predictAtrocities(int day) {
        predictingStarted = true;
        return predict(day);
    }

    private double[] predict(int day) {
        if (!fakePredictingStarted) {
            if (WS[FRIENDS] != 0) recalcFriends(day);
            fakePredictingStarted = true;
        }

        double[] ret = new double[R];
        double[] retForCountries = new double[R];
        for (int i = 0; i < R; i++) retForCountries[i] = -1;

        for (Country c: countries) c.fillCountryNewsStats(day);
        doCountryNewsStats(day);

        if (WS[FRIENDS] != 0 && day % 10 == 0) recalcFriends(day);

        if (day > day0fakePredict && (day - day0fakePredict) % TUNE_FREQ == 0 &&
                needsTuning) tuneWeights(day);

        for (int ri = 0; ri < R; ri++) {
            try {
                Region r = regions[ri];

                double[] ps = new double[ALGO_CNT];
                ps[POISSON] = WS[POISSON] == 0 ? 0 : guessFromPoisson(day, ri);
                ps[OWN] = WS[OWN] == 0 ? 0 : guessFromOwn(day, ri);
                ps[STAT] = WS[STAT] == 0 ? 0 : guessFromStats(day, ri);
                ps[NEWS] = WS[NEWS] == 0 ? 0 : guessFromNews(day, ri);
                ps[NEWS2] = WS[NEWS2] == 0 ? 0 : guessFromNews2(day, ri);
                ps[COUNTRIES] = WS[COUNTRIES] == 0 ? 0 : guessFromCountries(day, ri, retForCountries);
                ps[POISSONWCOUNTRY] = WS[POISSONWCOUNTRY] == 0 ? 0 : guessFromPoissonWithCountry(day, ri);
                ps[FRIENDS] = WS[FRIENDS] == 0 ? 0 : guessFromFriends(day, ri);
                for (int i = 0; i < ALGO_CNT-1; i++) {
                    if (Double.isNaN(ps[i])) ps[i] = 0;
                }

                if (needsTuning) {
                    for (int i = 0; i < ALGO_CNT-1; i++) r.guessHistory[day - day0][i] = (short)(10000 * ps[i]);
                }

                double pGlobalFix = 0;
                for (int i = 0; i < ALGO_CNT-1; i++) pGlobalFix += WS[i]*ps[i];
                pGlobalFix /= WSUM;

                double pGlobalTuned = 0;
                double pRegionTuned = 0;
                double pCountry = 0;

                if (needsTuning) {
                    for (int i = 0; i < ALGO_CNT-1; i++) pGlobalTuned += WS_TUNED[i]*ps[i];
                    pGlobalTuned /= WSUM_TUNED;

                    int t = r.rType;
                    for (int i = 0; i < ALGO_CNT-1; i++) pRegionTuned += WS_R[t][i]*ps[i];
                    pRegionTuned /= WSUM_R[t];

                    double[] cws = r.country.ws;
                    for (int i = 0; i < ALGO_CNT-1; i++) pCountry += cws[i]*ps[i];
                    pCountry /= r.country.wSum;
                }

                double p = (wGlobalFix * pGlobalFix + wGlobalTuned * pGlobalTuned + wCountryTuned * pCountry + wRegionTuned * pRegionTuned) / wMixSum;

                if (ps[FRIENDS] > 0.3) p = max(ps[POISSON], ps[OWN], ps[STAT], ps[COUNTRIES], ps[POISSONWCOUNTRY]);
                p = getOverride(p, ri, day);

                if (p < 0) p = 0; if (p > 1) p = 1;
                ret[ri] = p;
            }
            catch (Exception e) {
                e.printStackTrace();
            }

        }
        return ret;
    }

    private double getOverride(double p, int ri, int day) {
        try {
            Region r = regions[ri];
            int cnt = r.as.size();
            if (cnt < 8) return p;

            int D = 100;
            int startIndex = day - day0 - D;
            int good = 0;
            for (int i = 0; i < D; i++) {
                if (r.goodGuessPerDay[startIndex + i] == 1) good++;
            }
            double pg = (double)good / D;
            if (pg > 0.9) return pg;
            return p;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }

    private double guessFromCountries(int day, int ri, double[] retForCountries) {
        try {
            if (retForCountries[ri] != -1) return retForCountries[ri];

            Region r = regions[ri];
            Country c = countries[r.ci];
            int cnt = c.as.size();
            if (cnt < 2) {
                for (Region r2: c.rs) retForCountries[r2.id] = 0;
                return 0;
            }
            A a1 = c.as.get(cnt - 1);
            int dLast = a1.day;
            int T = day - dLast;
            if ((cnt == 2 && T > 200) ||
                    (cnt == 3 && T > 400) ||
                    (T > 900)) {
                for (Region r2: c.rs) retForCountries[r2.id] = 0;
                return 0;
            }

            double sum = 0;
            int previous = 0;
            for (A a: c.as) {
                if (day - a.day < POISSON_D) {
                    if (a.day > previous) {
                        double w = POISSON_D - (day - a.day);
                        sum += w;
                        previous = a.day;
                    }
                }
            }
            if (sum == 0) {
                for (Region r2: c.rs) retForCountries[r2.id] = 0;
                return 0;
            }

            double avg = sum / POISSON_SUMW * POISSON_M;
            double p0 = poisson(avg, 0);
            double p = 1 - p0;
            for (Region r2: c.rs) {
                double pR = p * r2.as.size() / cnt;
                retForCountries[r2.id] = pR;
            }
            return retForCountries[r.id];
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }

    private int[][] countryNewsCounts;
    private int[][] countryNewsPos;

    private void doCountryNewsStats(int day) {
        int W = 600;
        countryNewsCounts = new int[N_CNT][21];
        countryNewsPos = new int[N_CNT][21];
        int d0 = day-day0;
        for (Region r: regions) if (r.as.size() > 1) {
            for (int d = Math.max(0, d0-W); d < d0; d++) {
                int gg = r.goodGuessPerDay[d];
                for (int t = 0; t < 21; t++) {
                    int newsType = r.country.newsTypeHistory[d][t];
                    countryNewsCounts[newsType][t]++;
                    countryNewsPos[newsType][t] += gg;
                }
            }
        }
    }

    private double guessFromNews(int day, int ri) {
        try {
            Region r = regions[ri];
            if (r.as.size() < 2) return 0;

            double sum = 0;
            double sumW = 0;
            for (int t = 0; t < 21; t++) {
                int newsType = r.country.newsTypeHistory[day - day0][t];
                int cnt = countryNewsCounts[newsType][t];
                int pos = countryNewsPos[newsType][t];
                double p = 0;
                if (cnt > 0) p = (double)pos / cnt;
                double w = 1;
                if (t == 5 || t == 13 || t == 19) w = 2;
                if (t == 10) w = 5;
                sum += p*w;
                sumW += w;
            }
            if (sumW == 0) return 0;
            double p = sum / sumW;
            return p;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }


    private double guessFromNews2(int day, int ri) {
        try {
            Region r = regions[ri];
            int cnt = r.as.size();
            if (cnt < 2) return 0;

            double sumP = 0;
            for (int type = 0; type < 21; type++) {
                int cg00 = 0;
                int cg01 = 0;
                int cg10 = 0;
                int cg11 = 0;

                for (int i = 0; i < NEWS_BUFFER_SIZE; i++) {
                    int d = day - day0 - i;
                    int guess = r.goodGuessPerDay[d];
                    int cntD = r.newsPerDay[d % NEWS_BUFFER_SIZE][type];
                    if (guess == 0) {
                        if (cntD == 0) cg00++; else cg10++;
                    }
                    else {
                        if (cntD == 0) cg01++; else cg11++;
                    }
                }

                double pc0 = 0;
                if (cg00 + cg01 > 0) pc0 = (double)cg01 / (cg00 + cg01);
                double pc1 = 0;
                if (cg10 + cg11 > 0) pc1 = (double)cg11 / (cg10 + cg11);
                int cntToday = r.newsPerDay[(day-day0) % NEWS_BUFFER_SIZE][type];
                if (cntToday == 0) sumP += pc0;
                else sumP += pc1;
            }
            return sumP / 21;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }

    private final int POISSON_M = 24;
    private final int POISSON_D = POISSON_M*30;
    private final double POISSON_SUMW = POISSON_D * (POISSON_D-1) / 2;
    private final static int MAX_K = 50;
    private final static double[] FACTS = new double[MAX_K + 1];
    static {
        FACTS[0] = 1;
        for (int i = 1; i <= MAX_K; i++) {
            FACTS[i] = i * FACTS[i-1];
        }
    }

    private double guessFromPoisson(int day, int ri) {
        try {
            Region r = regions[ri];
            int cnt = r.as.size();
            if (cnt < 2) return 0;
            A a1 = r.as.get(cnt - 1);
            int dLast = a1.day;
            int T = day - dLast;
            if (cnt == 2 && T > 200) return 0;
            if (cnt == 3 && T > 400) return 0;
            if (T > 900) return 0;

            double sum = 0;
            int previous = 0;
            for (A a: r.as) {
                if (day - a.day < POISSON_D) {
                    if (a.day > previous + 7) {
                        sum += POISSON_D - (day - a.day);
                        previous = a.day;
                    }
                }
            }
            double avg = sum / POISSON_SUMW * POISSON_M;
            double p0 = poisson(avg, 0);
            double p = 1 - p0;
            return p;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }

    private double guessFromPoissonWithCountry(int day, int ri) {
        try {
            Region r = regions[ri];
            int cnt = r.as.size();
            if (cnt < 2) return 0;
            A a1 = r.as.get(cnt - 1);
            int dLast = a1.day;
            int T = day - dLast;
            if (cnt == 2 && T > 200) return 0;
            if (cnt == 3 && T > 400) return 0;
            if (T > 900) return 0;

            double sum = 0;
            int previous = 0;
            for (A a: r.as) {
                if (day - a.day < POISSON_D) {
                    if (a.day > previous + 7) {
                        sum += POISSON_D - (day - a.day);
                        previous = a.day;
                    }
                }
            }
            double avg = sum / POISSON_SUMW * POISSON_M;

            Country c = r.country;
            int d1 = day-day0;
            int d2 = Math.max(0, d1-30);
            boolean wasInC = false;
            for (Region r2: c.rs) {
                if (r2 != r && r2.cCnts[d1] != r2.cCnts[d2]) {
                    wasInC = true;
                }
            }

            double p0 = poisson(avg, 0);
            double p = 1 - p0;
            if (wasInC) p *= 1.1;
            return p;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }


    private static final int lookBack = 1000;
    private static final double[] weightsForLookBack;
    static {
        weightsForLookBack = new double[lookBack];
        for (int i = 0; i < lookBack; i++) {
            double h = (double)(i + 1) / lookBack;
            weightsForLookBack[i] = Math.pow(h, 1.5);
        }
    }

    private double guessFromOwn(int day, int ri) {
        try {
            Region r = regions[ri];
            int cnt = r.as.size();
            if (cnt < 2) return 0;

            int startIndex = Math.max(0, day - day0 - lookBack);
            double sum = 0;
            double sumW = 0;
            for (int i = 0; i < lookBack; i++) {
                double w = weightsForLookBack[i];
                if (r.goodGuessPerDay[startIndex + i] == 1) sum += w;
                sumW += w;
            }
            double p = sum / sumW;
            return p;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }

    private double guessFromStats(int day, int ri) {
        try {
            Region r = regions[ri];
            int cnt = r.as.size();
            if (cnt < 2) return 0;

            A a1 = r.as.get(cnt - 1);
            int dLast = a1.day;
            int T = day - dLast;
            if (T > 400) return 0;

            Stat lastFact = r.lastFact;
            if (lastFact == null) return 0;

            List<StatDistancePair> sdps = new Vector<>();
            for (Stat s: stats) {
                int N = 6;
                if (s.T < T) continue;
                double dist = 0;
                int n = Math.min(s.deltas.size(), lastFact.deltas.size());
                for (int i = 0; i < N; i++) {
                    double d = 300;
                    if (i < n) d = Math.abs(s.deltas.get(i) - lastFact.deltas.get(i));
                    dist += d / (i+1);
                }
                StatDistancePair sdp = new StatDistancePair();
                sdp.stat = s;
                sdp.d = dist;
                sdps.add(sdp);
            }
            if (sdps.size() == 0) return 0;
            StatDistancePair[] sdpArr = sdps.toArray(new StatDistancePair[0]);
            Arrays.sort(sdpArr);

            double sum = 0;
            double sumW = 0;
            for (int i = 0; i < Math.min(500, sdpArr.length); i++) {
                StatDistancePair sdp = sdpArr[i];
                double w = 1.0 / (20+sdp.stat.N);
                w *= 1 / (1 + sdp.d);
                double locFactor = 1;
                if (sdp.stat.country == r.ci) locFactor = 10;
                if (sdp.stat.region == r.id) locFactor = 20;
                w *= locFactor;
                if (sdp.stat.T < T + 30) {
                    sum += w;
                }
                sumW += w;
            }
            double p = sum / sumW;
            return p;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }

    private double guessFromFriends(int day, int ri) {
        try {
            Region r = regions[ri];
            if (r.friends == null || r.friends.length == 0) return 0;
            double sum = 0;
            double sumW = 0;
            double hit = 0;
            for (Friend f: r.friends) {
                Region r2 = regions[f.ri];
                A a = r2.as.get(r2.as.size() - 1);
                if (day - a.day < 30) {
                    hit++;
                    sum += f.w * f.p;
                    sumW += f.w;
                }
            }
            if (hit < 6) return 0;
            double p = sum / sumW;
            return p;
        }
        catch (Exception e) {
            e.printStackTrace(); return 0;
        }
    }

    private void recalcFriends(int day) {
        try {
            for (int ri = 0; ri < R; ri++) {
                Region r = regions[ri];
                if (r.as.size() < 2) continue;
                List<Friend> fs = new Vector<>();
                for (int ri2 = 0; ri2 < R; ri2++) {
                    if (ri == ri2) continue;
                    Region r2 = regions[ri2];
                    int cnt = 0;
                    int pos = 0;
                    for (A a: r2.as) {
                        if (day - a.day < 900) {
                            cnt++;
                            if (r.goodGuessPerDay[a.day - day0] == 1) pos++;
                        }
                    }
                    if (cnt < 2) continue;
                    if (pos == 0) continue;
                    Friend f = new Friend();
                    f.ri = ri2;
                    f.p = (double)pos / cnt;
                    if (cnt < 6) f.w = (double)1 / (10 - cnt);
                    else f.w = 1.0;
                    fs.add(f);
                }
                Friend[] fArr = fs.toArray(new Friend[0]);
                Arrays.sort(fArr);
                r.friends = fArr;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }

    //////////////////////////////////////////////////
    private double[] bestWs;
    private double[] baseWs;
    private double bestScore;
    private double allTimebestScore;
    private Map<String, Double> keyToScoreMap;
    private List<List<Region>> rsPerType;
    private Random rand = new Random(0);

    private void sumWeights() {
        WSUM_TUNED = 0;
        for (int a = 0; a < ALGO_CNT-1; a++) WSUM_TUNED += WS_TUNED[a];
    }

    private void tuneWeights(int day) {
        try {
            System.out.println("  Retuning at " + day);

            for (Country c: countries) {
                if (c.as.size() < 2) continue;
                tune(c.rs, day);
                c.ws = bestWs.clone();
                c.sumWeights();
            }

            List<Region> rs = new Vector<>();
            for (Region r: regions) if (r.as.size() > 1) rs.add(r);
            tune(rs, day);
            WS_TUNED = bestWs.clone();
            for (double w: WS_TUNED) System.out.println("    " + f(w));
            sumWeights();

            classifyRs(day);
            for (int t = 0; t < R_CNT; t++) {
                rs = rsPerType.get(t);
                tune(rs, day);
                WS_R[t] = bestWs.clone();
                WSUM_R[t] = 0;
                for (int a = 0; a < ALGO_CNT-1; a++) WSUM_R[t] += WS_R[t][a];
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    private void tune(List<Region> rs, int day) {
        baseWs = WS.clone();
        bestScore = -Double.MAX_VALUE;
        allTimebestScore = -Double.MAX_VALUE;
        int prevImprovedAt = 0;
        int prevParam = -1;
        keyToScoreMap = new HashMap<>();

        for (int iteration = 0; iteration < 500; iteration++) {
            if (iteration > prevImprovedAt + 20) {
                bestScore = -Double.MAX_VALUE;
                for (int i = 0; i < ALGO_CNT-1; i++) baseWs[i] = 100 * rand.nextDouble();
                prevImprovedAt = iteration;
            }
            int param = 0;
            while (true) {
                param = rand.nextInt(ALGO_CNT-1);
                if (param != prevParam) break;
                prevParam = param;
            }
            double[] wsTested = baseWs.clone();
            double old = wsTested[param];
            int i = 0;
            while (i < 5) {
                wsTested[param] = old * Math.pow(1.1, i+1);
                boolean improved = iterate(wsTested, rs, day);
                i++;
                if (!improved) break;
                prevImprovedAt = iteration;
            }
            i = 0;
            while (i < 5) {
                wsTested[param] = old * Math.pow(1.1, -(i+1));
                boolean improved = iterate(wsTested, rs, day);
                i++;
                if (!improved) break;
                prevImprovedAt = iteration;
            }
        }
    }

    private boolean iterate(double[] wsTested, List<Region> rs, int day) {

        String key = "";
        for (int i = 0; i < ALGO_CNT-1; i++) key += f(wsTested[i]) + "-";
        if (keyToScoreMap.containsKey(key)) {
            return false;
        }

        double wSum = 0;
        for (int i = 0; i < ALGO_CNT-1; i++) wSum += wsTested[i];
        double sumScore = 0;
        double sumDailyW = 0;
        int d0 = day - day0 - TUNE_W;
        for (Region r: rs) {
            if (r.as.size() == 0) continue;

            for (int di = 0; di < TUNE_W; di++) {
                if (d0 + di < 0) continue;
                double p = 0;
                for (int i = 0; i < ALGO_CNT-1; i++) p += wsTested[i] * r.guessHistory[d0+di][i] / 10000;
                p /= wSum;
                double score = 0;
                double wgh = wgh(r.wghTs[d0 + di]);
                if (r.goodGuessPerDay[d0 + di] == 1) score = wgh * (p - p*p/2);
                else score = -wgh * p*p/2;
                double w = 1 + (double)di / TUNE_W;
                sumScore += score * w;
                sumDailyW += w;
            }
        }
        sumScore /= sumDailyW;

        boolean ret = false;
        if (sumScore > bestScore && (sumScore - bestScore > sumScore / 1e4)) {
            if (sumScore > allTimebestScore) {
                allTimebestScore = sumScore;
                bestWs = wsTested.clone();
            }
            bestScore = sumScore;
            baseWs = wsTested.clone();
            ret = true;
        }
        keyToScoreMap.put(key, sumScore);
        return ret;
    }

    private void classifyRs(int day) {
        int d = day - day0;
        for (Region r: regions) {
            r.rType = R_xx; // default
            int cnt3y = r.cCnts[d] - r.cCnts[d-1000];
            if (cnt3y < 4) {
                r.rType = R_00; continue;
            }
            int y1 = r.cCnts[d] - r.cCnts[d-300];
            int y2 = r.cCnts[d-300] - r.cCnts[d-600];
            if (cnt3y < 10) {
                if (y1 > y2 + 2) {
                    r.rType = R_10; continue;
                }
                else if (y1 + 2 < y2) {
                    r.rType = R_01; continue;
                }
                else if (y1 + y2 > 5) {
                    r.rType = R_11; continue;
                }
                else {
                    r.rType = R_00; continue;
                }
            }
            int q1 = r.cCnts[d] - r.cCnts[d-90];
            int q2 = r.cCnts[d-90] - r.cCnts[d-180];
            if (q1 > 6) {
                if (q2 > 6) {
                    r.rType = R_22; continue;
                }
                else if (q2 > 3) {
                    r.rType = R_21; continue;
                }
                else {
                    r.rType = R_20; continue;
                }
            }
            if (q2 > 6) {
                if (q1 > 3) {
                    r.rType = R_12; continue;
                }
                r.rType = R_02; continue;
            }
        }

        rsPerType = new Vector<>();
        for (int i = 0; i < R_CNT; i++) rsPerType.add(new Vector<Region>());
        for (Region r: regions) {
            int type = r.rType;
            rsPerType.get(type).add(r);
        }
    }


    /////////////////////////////////////////////////////////////////////

    private class A {
        int day;
        int region;

        @Override
        public String toString() {
            return day + ": " + region;
        }
    }

    private class Region {
        public int ci;
        public Country country;
        public int id;
        public List<A> as = new Vector<>();
        public Stat lastFact;
        public byte[] goodGuessPerDay = new byte[6400];
        public short[][] newsPerDay = new short[NEWS_BUFFER_SIZE][21];
        public Friend[] friends;
        public short[][] guessHistory = new short[6400][ALGO_CNT-1];
        public short[] wghTs = new short[6400];
        public short[] cCnts = new short[6400];
        public int rType;

        @Override
        public String toString() {
            return id + " " + as.toString();
        }
    }

    private class Country {
        public int id;
        public List<A> as = new Vector<>();
        public List<Region> rs = new Vector<>();
        public double[] ws;
        public double wSum;

        public byte[][] newsTypeHistory = new byte[6400][21];
        public int[][] newsPerDay = new int[NEWS_BUFFER_SIZE][21];

        @Override
        public String toString() {
            return id + " " + rs.size() + " " + as.size();
        }

        public void fillCountryNewsStats(int day) {
            int w1 = 30;
            int w2 = 10;
            int d0 = day - day0;
            try {
                int[] sums = new int[21];
                for (int i = 0; i < NEWS_BUFFER_SIZE; i++) {
                    for (int t = 0; t < 21; t++) sums[t] += newsPerDay[i][t];
                }
                double[] avgs = new double[21];
                for (int t = 0; t < 21; t++) avgs[t] = (double)sums[t] / NEWS_BUFFER_SIZE;

                double[] r1s = new double[21];
                double[] r2s = new double[21];
                for (int i = w2; i < w1+w2; i++) {
                    int d = (d0 + NEWS_BUFFER_SIZE - i) % NEWS_BUFFER_SIZE;
                    for (int t = 0; t < 21; t++) r1s[t] += newsPerDay[d][t];
                }
                for (int i = 0; i < w2; i++) {
                    int d = (d0 + NEWS_BUFFER_SIZE - i) % NEWS_BUFFER_SIZE;
                    for (int t = 0; t < 21; t++) r2s[t] += newsPerDay[d][t];
                }
                for (int t = 0; t < 21; t++) {
                    r1s[t] /= w1;
                    r2s[t] /= w2;
                    double a = avgs[t];
                    double r1 = 0;
                    if (a > 0) r1 = r1s[t] / a;
                    double r2 = 0;
                    if (a > 0) r2 = r2s[t] / a;

                    if (r2 > a*3) {
                        newsTypeHistory[d0][t] = N_H; continue;
                    }
                    if (r2 < a/3) {
                        newsTypeHistory[d0][t] = N_L; continue;
                    }
                    byte c = N_xx;
                    if (r1 > a*1.5) {
                        if (r2 > a*1.5) c = N_PP;
                        else if (r2 < a/1.5) c = N_PM;
                        else c = N_PA;
                    }
                    else if (r1 < a*1.5) {
                        if (r2 > a*1.5) c = N_MP;
                        else if (r2 < a/1.5) c = N_MM;
                        else c = N_MA;
                    }
                    else {
                        if (r2 > a*1.5) c = N_AP;
                        else if (r2 < a/1.5) c = N_AM;
                        else c = N_AA;
                    }
                    newsTypeHistory[d0][t] = c;
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }

        public void sumWeights() {
            wSum = 0;
            for (int a = 0; a < ALGO_CNT-1; a++) wSum += ws[a];
        }
    }

    private class Friend implements Comparable<Friend> {
        public int ri;
        public double p;
        public double w;

        @Override
        public int compareTo(Friend o) {
            if (this.w < o.w) return 1;
            if (this.w > o.w) return -1;
            return 0;
        }

        @Override
        public String toString() {
            return ri + " " + f(p) + " (" + f(w) + ")";
        }
    }

    private class Stat {
        public int N = 0;
        public List<Integer> deltas;
        public int dLast;
        public int T;
        public int region;
        public int country;

        public Stat(Region r, boolean train) {
            region = r.id;
            country = r.ci;
            int cnt = r.as.size();
            N = train ? cnt - 1 : cnt;
            dLast = r.as.get(N-1).day;
            if (train) T = r.as.get(N).day - r.as.get(N-1).day;
            deltas = new Vector<>();
            int previous = dLast;
            for (int i = N-2; i >= 0; i--) {
                A a = r.as.get(i);
                if (a.day > previous - 7) continue;
                previous = a.day;
                deltas.add(dLast - a.day);
            }
        }

        @Override
        public String toString() {
            String ret = N + " [";
            for (int i: deltas) ret += i + " ";
            ret += "] -> " + T + " (@" + dLast + ", r=" + region + ")";
            return ret;
        }
    }


    private class StatDistancePair implements Comparable<StatDistancePair> {
        public Stat stat;
        public double d;

        @Override
        public int compareTo(StatDistancePair o) {
            if (this.d < o.d) return -1;
            if (this.d > o.d) return 1;
            return o.stat.dLast - this.stat.dLast;
        }

        @Override
        public String toString() {
            return f(d) + ": " + stat.toString();
        }
    }

    private static double poisson(double avg, int k) {
        if (k > MAX_K) return 0; // what else?
        return Math.pow(avg, k) * Math.exp(-avg) / FACTS[k];
    }

    private static final DecimalFormat df = new DecimalFormat("0.###");
    private static final DecimalFormatSymbols dfs = new DecimalFormatSymbols();
    static {
        dfs.setDecimalSeparator('.');
        df.setDecimalFormatSymbols(dfs);
    }
    public static String f(double d) {
        return df.format(d);
    }

    private double max(double... ds) {
        double max = -Double.MAX_VALUE;
        for (double d: ds) max = Math.max(max, d);
        return max;
    }

    private static final double[] WGHS = new double[6400];
    static {
        for (int i = 0; i < WGHS.length; i++)
            WGHS[i] = Math.tanh((double)(i + 10) / 180);
    }

    public static double wgh(int delta) {
        if (delta < 0 || delta > WGHS.length-1) return 1;
        return WGHS[delta];
    }

    // ----- main method ----- //
    // added by TopCoder staff //

    static void error(String s) {
        System.out.println("ERROR: " + s);
        System.exit(0);
    }

    static String[] readData(String fileName) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            int n = Integer.parseInt(br.readLine());
            String[] res = new String[n];
            for (int i=0; i < n; i++) {
                res[i] = br.readLine();
            }
            br.close();
            return res;
        } catch (IOException e) {
            error("Unable to read data from \"" + fileName + "\".");
        }
        return null;
    }

    static final int MIN_DAY = 5440;
    static final int MIN_REC_DAY = 11284;
    static final int MAX_REC_DAY = 17644;
    static final int MAX_DAY = 17858;

    public static void main(String[] args) throws Exception {
        if (args.length != 5) {
            error("This program is to be executed as follows:\nmain <input folder> <output folder> <first training day> <first testing day> <last testing day>");
        }

        String dataFolder = args[0];
        String outFolder = args[1];

        String[] eventsData = readData(dataFolder + "/events.txt");
        String[] regionsData = readData(dataFolder + "/regions.txt");

        List<List<String>> events = new ArrayList<List<String>>();
        for (int i=0; i <= MAX_DAY; i++) {
            events.add(new ArrayList<String>());
        }

        for (int i=0; i < eventsData.length; i++) {
            String[] items = eventsData[i].split(" ");
            int day = Integer.parseInt(items[0]);
            String lat = items[1], lon = items[2];
            int country = Integer.parseInt(items[3]);
            int region = Integer.parseInt(items[4]);
            events.get(day).add(lat + " " + lon + " " + country + " " + region);
        }

        int startTrain = Integer.parseInt(args[2]), startTest = Integer.parseInt(args[3]), endTest = Integer.parseInt(args[4]);

        if (startTrain < MIN_DAY) {
            error("The value of <first training day> parameter must be " + MIN_DAY + " or above.");
        }

        if (startTrain < MIN_REC_DAY) {
            System.out.println("WARNING: The value of <first training day> parameter is recommended to be " + MIN_REC_DAY  + " or above.");
        }

        if (startTest < startTrain) {
            error("The value of <first testing day> parameter must be greater than or equal to the value of <first training day> parameter.");
        }

        if (endTest < startTest) {
            error("The value of <last testing day> parameter must be greater than or equal to the value of <first testing day> parameter.");
        }

        if (endTest > MAX_REC_DAY) {
            System.out.println("WARNING: The value of <last testing day> parameter is recommended to be " + MAX_REC_DAY + " or below.");
        }

        if (endTest > MAX_DAY) {
            error("The value of <last testing day> parameter must be " + MAX_DAY + " or below.");
        }

        MassAtrocityPredictor obj = new MassAtrocityPredictor();
        obj.receiveData(2, startTrain, regionsData);

        for (int curDay = startTrain; curDay <= endTest; curDay++) {
            System.out.println("Day = " + curDay);

            String[] data = readData(dataFolder + "/data_" + curDay + ".txt");

            obj.receiveData(1, curDay, data);
            obj.receiveData(0, curDay, events.get(curDay).toArray(new String[0]));

            if (curDay >= startTest) {
                double[] res = obj.predictAtrocities(curDay);
                String resFileName = outFolder + "/res_" + curDay + ".txt";
                try {
                    PrintWriter pw = new PrintWriter(resFileName);
                    for (int i=0; i < res.length; i++) {
                        pw.println(res[i]);
                    }
                    pw.flush();
                    pw.close();
                } catch (IOException e) {
                    error("Unable to save data to \"" + resFileName + "\".");
                }
            }
        }
    }
}

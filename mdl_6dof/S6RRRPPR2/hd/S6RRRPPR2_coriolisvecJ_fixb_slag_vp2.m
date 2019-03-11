% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:58
% EndTime: 2019-03-09 15:25:16
% DurationCPUTime: 8.11s
% Computational Cost: add. (10925->555), mult. (28734->718), div. (0->0), fcn. (20958->8), ass. (0->264)
t252 = sin(qJ(3));
t253 = sin(qJ(2));
t255 = cos(qJ(3));
t256 = cos(qJ(2));
t223 = -t252 * t253 + t255 * t256;
t211 = t223 * qJD(1);
t224 = t252 * t256 + t255 * t253;
t212 = t224 * qJD(1);
t249 = sin(pkin(10));
t250 = cos(pkin(10));
t265 = t211 * t249 + t250 * t212;
t338 = t265 / 0.2e1;
t168 = t211 * t250 - t212 * t249;
t342 = t168 / 0.2e1;
t376 = -Ifges(6,2) * t338 - Ifges(6,6) * t342;
t248 = qJD(2) + qJD(3);
t334 = -t248 / 0.2e1;
t375 = mrSges(5,1) - mrSges(6,2);
t365 = mrSges(5,3) + mrSges(6,1);
t352 = -pkin(8) - pkin(7);
t236 = t352 * t256;
t229 = qJD(1) * t236;
t216 = t255 * t229;
t235 = t352 * t253;
t228 = qJD(1) * t235;
t182 = -t228 * t252 + t216;
t303 = qJ(4) * t211;
t150 = t182 - t303;
t213 = t252 * t229;
t183 = t255 * t228 + t213;
t206 = t212 * qJ(4);
t151 = -t206 + t183;
t238 = t249 * t252 * pkin(2);
t289 = qJD(3) * t255;
t358 = -qJD(3) * t238 - t150 * t249 + (pkin(2) * t289 - t151) * t250;
t251 = sin(qJ(6));
t254 = cos(qJ(6));
t219 = qJD(2) * pkin(2) + t228;
t177 = t255 * t219 + t213;
t144 = t177 - t206;
t136 = pkin(3) * t248 + t144;
t178 = t219 * t252 - t216;
t145 = t178 + t303;
t137 = t249 * t145;
t83 = t136 * t250 - t137;
t279 = qJD(5) - t83;
t328 = pkin(5) * t265;
t353 = pkin(4) + pkin(9);
t44 = -t248 * t353 + t279 + t328;
t244 = -pkin(2) * t256 - pkin(1);
t234 = qJD(1) * t244;
t186 = -t211 * pkin(3) + qJD(4) + t234;
t261 = -qJ(5) * t265 + t186;
t52 = -t168 * t353 + t261;
t15 = t251 * t44 + t254 * t52;
t14 = -t251 * t52 + t254 * t44;
t315 = t14 * t251;
t269 = -t15 * t254 + t315;
t341 = -t168 / 0.2e1;
t327 = pkin(5) * t168;
t374 = -Ifges(6,4) + Ifges(5,5);
t373 = Ifges(6,5) - Ifges(5,6);
t148 = -t168 * t254 - t248 * t251;
t149 = -t168 * t251 + t248 * t254;
t160 = qJD(6) + t265;
t371 = Ifges(5,1) * t265 + Ifges(5,4) * t168 + t248 * Ifges(5,5) + t149 * Ifges(7,5) + t148 * Ifges(7,6) + t160 * Ifges(7,3);
t361 = -qJD(5) - t358;
t278 = mrSges(7,1) * t254 - mrSges(7,2) * t251;
t332 = -t251 / 0.2e1;
t297 = t250 * t145;
t84 = t249 * t136 + t297;
t82 = -qJ(5) * t248 - t84;
t46 = -t82 + t327;
t313 = t149 * Ifges(7,4);
t56 = t148 * Ifges(7,2) + t160 * Ifges(7,6) + t313;
t146 = Ifges(7,4) * t148;
t57 = Ifges(7,1) * t149 + Ifges(7,5) * t160 + t146;
t370 = -t254 * t56 / 0.2e1 + t57 * t332 + t278 * t46;
t88 = -pkin(4) * t168 + t261;
t369 = t186 * mrSges(5,1) - t88 * mrSges(6,2);
t81 = -pkin(4) * t248 + t279;
t368 = t81 * mrSges(6,1) + t14 * mrSges(7,1) + t186 * mrSges(5,2) - t15 * mrSges(7,2) - t83 * mrSges(5,3) - t88 * mrSges(6,3) + Ifges(6,4) * t334 - t376;
t339 = -t265 / 0.2e1;
t364 = Ifges(5,4) * t265;
t363 = Ifges(6,6) * t265;
t360 = t328 - t361;
t296 = t250 * t252;
t202 = (t249 * t255 + t296) * qJD(3) * pkin(2);
t91 = -t250 * t150 + t151 * t249;
t359 = t202 - t91;
t294 = -t248 * t375 + t265 * t365;
t188 = t252 * t235 - t255 * t236;
t184 = t248 * t223;
t173 = t184 * qJD(1);
t185 = t248 * t224;
t174 = t185 * qJD(1);
t123 = t173 * t250 - t174 * t249;
t286 = qJD(2) * t352;
t282 = qJD(1) * t286;
t220 = t253 * t282;
t221 = t256 * t282;
t126 = -qJD(3) * t178 - t220 * t252 + t221 * t255;
t258 = -qJ(4) * t173 - qJD(4) * t212 + t126;
t290 = qJD(3) * t252;
t125 = t219 * t289 + t255 * t220 + t252 * t221 + t229 * t290;
t69 = -qJ(4) * t174 + qJD(4) * t211 + t125;
t27 = t249 * t69 - t250 * t258;
t11 = pkin(5) * t123 + t27;
t122 = t173 * t249 + t250 * t174;
t293 = qJD(1) * t253;
t246 = pkin(2) * t293;
t156 = pkin(3) * t174 + qJD(2) * t246;
t263 = -qJ(5) * t123 - qJD(5) * t265 + t156;
t13 = t122 * t353 + t263;
t300 = qJD(6) * t15;
t2 = t11 * t254 - t13 * t251 - t300;
t1 = qJD(6) * t14 + t11 * t251 + t13 * t254;
t326 = t1 * t251;
t259 = m(7) * (-qJD(6) * t269 + t2 * t254 + t326);
t96 = -mrSges(7,2) * t160 + mrSges(7,3) * t148;
t97 = mrSges(7,1) * t160 - mrSges(7,3) * t149;
t268 = -t251 * t97 + t254 * t96;
t64 = qJD(6) * t148 + t122 * t251;
t34 = mrSges(7,1) * t123 - mrSges(7,3) * t64;
t65 = -qJD(6) * t149 + t122 * t254;
t35 = -mrSges(7,2) * t123 + mrSges(7,3) * t65;
t357 = qJD(6) * t268 + t251 * t35 + t254 * t34 + t259;
t270 = t14 * t254 + t15 * t251;
t356 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t64 + Ifges(7,6) * t65;
t355 = t64 / 0.2e1;
t354 = t65 / 0.2e1;
t351 = pkin(1) * mrSges(3,1);
t350 = pkin(1) * mrSges(3,2);
t348 = -t148 / 0.2e1;
t347 = t148 / 0.2e1;
t346 = -t149 / 0.2e1;
t345 = t149 / 0.2e1;
t344 = -t160 / 0.2e1;
t343 = t160 / 0.2e1;
t336 = t211 / 0.2e1;
t335 = t212 / 0.2e1;
t331 = t254 / 0.2e1;
t330 = m(4) * t234;
t329 = pkin(3) * t212;
t324 = -Ifges(5,4) - Ifges(6,6);
t28 = t249 * t258 + t250 * t69;
t323 = mrSges(4,3) * t211;
t321 = Ifges(3,4) * t253;
t320 = Ifges(7,4) * t251;
t319 = Ifges(7,4) * t254;
t318 = Ifges(7,5) * t251;
t317 = Ifges(7,6) * t254;
t187 = t255 * t235 + t236 * t252;
t163 = -qJ(4) * t224 + t187;
t164 = qJ(4) * t223 + t188;
t100 = -t250 * t163 + t164 * t249;
t316 = t100 * t27;
t311 = t265 * t82;
t310 = t212 * mrSges(4,3);
t309 = t212 * Ifges(4,4);
t308 = t84 * t265;
t154 = -mrSges(6,1) * t168 - mrSges(6,3) * t248;
t94 = -mrSges(7,1) * t148 + mrSges(7,2) * t149;
t307 = t154 - t94;
t305 = Ifges(3,5) * qJD(2);
t304 = Ifges(3,6) * qJD(2);
t302 = qJD(2) * mrSges(3,1);
t301 = qJD(2) * mrSges(3,2);
t127 = t184 * t249 + t250 * t185;
t299 = t127 * t251;
t298 = t127 * t254;
t152 = -mrSges(5,2) * t248 + mrSges(5,3) * t168;
t295 = t154 - t152;
t292 = qJD(1) * t256;
t291 = qJD(2) * t253;
t242 = -pkin(3) * t250 - pkin(4);
t285 = t305 / 0.2e1;
t284 = -t304 / 0.2e1;
t170 = pkin(2) * t291 + pkin(3) * t185;
t283 = -t2 - t300;
t230 = t253 * t286;
t231 = t256 * t286;
t132 = t255 * t230 + t252 * t231 + t235 * t289 + t236 * t290;
t89 = -qJ(4) * t185 + qJD(4) * t223 + t132;
t133 = -qJD(3) * t188 - t230 * t252 + t255 * t231;
t90 = -qJ(4) * t184 - qJD(4) * t224 + t133;
t32 = t249 * t89 - t250 * t90;
t85 = t144 * t249 + t297;
t243 = pkin(2) * t255 + pkin(3);
t204 = t243 * t250 - t238;
t198 = -pkin(4) - t204;
t24 = -qJD(5) * t248 - t28;
t280 = t1 * t254 - t2 * t251;
t277 = mrSges(7,1) * t251 + mrSges(7,2) * t254;
t276 = Ifges(7,1) * t254 - t320;
t275 = Ifges(7,1) * t251 + t319;
t274 = -Ifges(7,2) * t251 + t319;
t273 = Ifges(7,2) * t254 + t320;
t272 = Ifges(7,5) * t254 - Ifges(7,6) * t251;
t271 = t317 + t318;
t33 = t249 * t90 + t250 * t89;
t180 = -t250 * t223 + t224 * t249;
t181 = t223 * t249 + t224 * t250;
t192 = -t223 * pkin(3) + t244;
t264 = -t181 * qJ(5) + t192;
t74 = t180 * t353 + t264;
t75 = pkin(5) * t181 + t100;
t30 = t251 * t75 + t254 * t74;
t29 = -t251 * t74 + t254 * t75;
t267 = -t251 * t96 - t254 * t97;
t86 = t144 * t250 - t137;
t101 = t163 * t249 + t164 * t250;
t205 = pkin(2) * t296 + t243 * t249;
t95 = pkin(4) * t265 - qJ(5) * t168 + t329;
t93 = t246 + t95;
t128 = t184 * t250 - t185 * t249;
t262 = -qJ(5) * t128 - qJD(5) * t181 + t170;
t10 = -pkin(5) * t122 - t24;
t104 = t248 * Ifges(6,5) - Ifges(6,3) * t168 - t363;
t106 = Ifges(5,2) * t168 + t248 * Ifges(5,6) + t364;
t161 = t211 * Ifges(4,2) + t248 * Ifges(4,6) + t309;
t207 = Ifges(4,4) * t211;
t162 = t212 * Ifges(4,1) + t248 * Ifges(4,5) + t207;
t8 = t64 * Ifges(7,4) + t65 * Ifges(7,2) + t123 * Ifges(7,6);
t9 = t64 * Ifges(7,1) + t65 * Ifges(7,4) + t123 * Ifges(7,5);
t257 = t371 * t341 + (mrSges(7,3) * t315 + t370) * qJD(6) + (mrSges(7,3) * t269 - Ifges(5,2) * t341 + Ifges(6,3) * t342 + t271 * t344 + t273 * t348 + t275 * t346 + t373 * t334 - t369 + t370) * t265 + (-t364 + t104) * t339 - t234 * (mrSges(4,1) * t212 + mrSges(4,2) * t211) + (t363 + t106) * t338 + (Ifges(5,1) * t339 + Ifges(5,4) * t341 + Ifges(7,5) * t346 + Ifges(7,6) * t348 + Ifges(7,3) * t344 + t334 * t374 - t368 + t376) * t168 + t10 * t277 - Ifges(4,6) * t174 + Ifges(4,5) * t173 - t125 * mrSges(4,2) + t126 * mrSges(4,1) - t28 * mrSges(5,2) - t24 * mrSges(6,3) + t274 * t354 + t276 * t355 + t373 * t122 + (t272 / 0.2e1 + t374) * t123 - t375 * t27 - (t148 * t273 + t149 * t275 + t160 * t271) * qJD(6) / 0.2e1 - (-Ifges(4,2) * t212 + t162 + t207) * t211 / 0.2e1 + t177 * t323 + t9 * t331 + t8 * t332 + (Ifges(4,5) * t211 - Ifges(4,6) * t212) * t334 + t161 * t335 - t212 * (Ifges(4,1) * t211 - t309) / 0.2e1;
t245 = Ifges(3,4) * t292;
t241 = pkin(3) * t249 + qJ(5);
t233 = mrSges(3,3) * t292 - t301;
t232 = -mrSges(3,3) * t293 + t302;
t210 = Ifges(3,1) * t293 + t245 + t305;
t209 = t304 + (Ifges(3,2) * t256 + t321) * qJD(1);
t197 = qJ(5) + t205;
t195 = -pkin(9) + t198;
t191 = mrSges(4,1) * t248 - t310;
t190 = -mrSges(4,2) * t248 + t323;
t189 = t246 + t329;
t176 = -mrSges(4,1) * t211 + mrSges(4,2) * t212;
t159 = t265 * pkin(9);
t117 = Ifges(7,3) * t123;
t116 = t123 * mrSges(6,3);
t115 = t123 * mrSges(5,2);
t114 = mrSges(6,2) * t168 - mrSges(6,3) * t265;
t113 = -mrSges(5,1) * t168 + mrSges(5,2) * t265;
t99 = t180 * pkin(4) + t264;
t76 = -pkin(5) * t180 + t101;
t61 = t159 + t95;
t60 = t159 + t93;
t58 = t91 + t327;
t51 = t86 - t328;
t50 = t85 + t327;
t36 = pkin(4) * t127 + t262;
t31 = pkin(4) * t122 + t263;
t23 = -mrSges(7,1) * t65 + mrSges(7,2) * t64;
t22 = t127 * t353 + t262;
t21 = -pkin(5) * t127 + t33;
t20 = pkin(5) * t128 + t32;
t19 = t251 * t58 + t254 * t60;
t18 = -t251 * t60 + t254 * t58;
t17 = t251 * t50 + t254 * t61;
t16 = -t251 * t61 + t254 * t50;
t4 = -qJD(6) * t30 + t20 * t254 - t22 * t251;
t3 = qJD(6) * t29 + t20 * t251 + t22 * t254;
t5 = [t244 * (mrSges(4,1) * t174 + mrSges(4,2) * t173) + (t125 * t223 - t126 * t224 - t173 * t187 - t174 * t188 - t177 * t184 - t178 * t185) * mrSges(4,3) + (-t223 * t174 - t185 * t336) * Ifges(4,2) + (t223 * t173 - t174 * t224 + t184 * t336 - t185 * t335) * Ifges(4,4) + (-t84 * mrSges(5,3) + t82 * mrSges(6,1) + t104 / 0.2e1 - t106 / 0.2e1 - Ifges(5,4) * t338 + Ifges(6,6) * t339 + Ifges(6,3) * t341 - Ifges(5,2) * t342 + t369) * t127 + (Ifges(5,1) * t338 + Ifges(5,4) * t342 + Ifges(7,5) * t345 - Ifges(6,2) * t339 - Ifges(6,6) * t341 + Ifges(7,6) * t347 + Ifges(7,3) * t343 + t368) * t128 + (Ifges(7,5) * t299 + Ifges(7,6) * t298) * t343 + t99 * (-t122 * mrSges(6,2) - t116) + (Ifges(7,4) * t299 + Ifges(7,2) * t298) * t347 + (Ifges(7,1) * t299 + Ifges(7,4) * t298) * t345 + m(7) * (t1 * t30 + t10 * t76 + t14 * t4 + t15 * t3 + t2 * t29 + t21 * t46) + (-t10 * t278 + t275 * t355 + t273 * t354 - t31 * mrSges(6,2) + t156 * mrSges(5,1) + t8 * t331 + t251 * t9 / 0.2e1 + t24 * mrSges(6,1) - t28 * mrSges(5,3) + (t318 / 0.2e1 + t317 / 0.2e1 + t324) * t123 + (Ifges(6,3) + Ifges(5,2)) * t122 + (t272 * t343 + t274 * t347 + t276 * t345 + t277 * t46 + t331 * t57 + t332 * t56) * qJD(6)) * t180 + (t15 * t298 - t14 * t299 + (-qJD(6) * t270 + t280) * t180) * mrSges(7,3) + m(4) * (t125 * t188 + t126 * t187 + t132 * t178 + t133 * t177) + t132 * t190 + t133 * t191 + t192 * (t122 * mrSges(5,1) + t115) + t184 * t162 / 0.2e1 - t185 * t161 / 0.2e1 + t170 * t113 + t36 * t114 + t21 * t94 + t3 * t96 + t4 * t97 + t76 * t23 + t29 * t34 + t30 * t35 + t234 * (mrSges(4,1) * t185 + mrSges(4,2) * t184) + (t210 / 0.2e1 - pkin(7) * t232 + t285 + (-0.2e1 * t350 + 0.3e1 / 0.2e1 * Ifges(3,4) * t256) * qJD(1)) * t256 * qJD(2) + t371 * t128 / 0.2e1 + (Ifges(4,5) * t184 - Ifges(4,6) * t185 + t127 * t373 + t128 * t374) * t248 / 0.2e1 + t365 * (t100 * t123 - t101 * t122) + (t117 / 0.2e1 - t31 * mrSges(6,3) + t156 * mrSges(5,2) + t365 * t27 + t324 * t122 + (Ifges(7,3) / 0.2e1 + Ifges(5,1) + Ifges(6,2)) * t123 + t356) * t181 + t294 * t32 - t295 * t33 + t56 * t298 / 0.2e1 + t57 * t299 / 0.2e1 + t46 * (-mrSges(7,1) * t298 + mrSges(7,2) * t299) + m(5) * (t101 * t28 + t156 * t192 + t170 * t186 - t32 * t83 + t33 * t84 + t316) + m(6) * (-t101 * t24 + t31 * t99 + t32 * t81 - t33 * t82 + t36 * t88 + t316) + (t173 * t224 + t184 * t335) * Ifges(4,1) + (-pkin(7) * t233 - t209 / 0.2e1 + t284 + (-0.2e1 * t351 - 0.3e1 / 0.2e1 * t321 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t256) * qJD(1) + (t176 + 0.2e1 * t330 + qJD(1) * (-mrSges(4,1) * t223 + mrSges(4,2) * t224)) * pkin(2)) * t291; t361 * t154 + t358 * t152 + ((t190 * t255 - t191 * t252) * qJD(3) + (-t173 * t255 - t174 * t252) * mrSges(4,3)) * pkin(2) + (-t1 * mrSges(7,3) + t202 * t96 + (-qJD(6) * t97 + t35) * t195) * t251 + t360 * t94 + t257 - t189 * t113 - t183 * t190 - t182 * t191 + t197 * t23 - t93 * t114 - t19 * t96 - t18 * t97 + t195 * t259 + (t202 * t97 + (qJD(6) * t96 + t34) * t195 + t283 * mrSges(7,3)) * t254 + t178 * t310 + (-t122 * t205 - t123 * t204 + t308) * mrSges(5,3) + (-t122 * t197 + t123 * t198 - t311) * mrSges(6,1) + ((t285 - t245 / 0.2e1 - t210 / 0.2e1 + qJD(1) * t350 + (t232 - t302) * pkin(7)) * t256 + (t284 + t209 / 0.2e1 + (t351 + t321 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t256) * qJD(1) + (t233 + t301) * pkin(7) + (-t176 - t330) * pkin(2)) * t253) * qJD(1) + t294 * t359 + (t10 * t197 - t14 * t18 - t15 * t19 + t202 * t270 + t360 * t46) * m(7) + (-t197 * t24 + t198 * t27 + t359 * t81 + t361 * t82 - t88 * t93) * m(6) + (-t186 * t189 - t204 * t27 + t205 * t28 + t358 * t84 - t359 * t83) * m(5) + ((t125 * t252 + t126 * t255 + (-t177 * t252 + t178 * t255) * qJD(3)) * pkin(2) - t177 * t182 - t178 * t183) * m(4); (t308 + (-t122 * t249 - t123 * t250) * pkin(3)) * mrSges(5,3) + (t254 * t283 - t326) * mrSges(7,3) + t295 * t86 - t294 * t85 + (-t122 * t241 + t123 * t242 - t311) * mrSges(6,1) + t257 - t307 * qJD(5) + t241 * t23 - t177 * t190 - t113 * t329 - t95 * t114 - t51 * t94 - t17 * t96 - t16 * t97 + (t191 + t310) * t178 + t357 * (-pkin(9) + t242) + (t10 * t241 - t14 * t16 - t15 * t17 + (qJD(5) - t51) * t46) * m(7) + (-t24 * t241 + t242 * t27 - t81 * t85 - t88 * t95 + (-qJD(5) + t86) * t82) * m(6) + ((t249 * t28 - t250 * t27) * pkin(3) - t186 * t329 + t83 * t85 - t84 * t86) * m(5); -t251 * t34 + t254 * t35 + t115 - t116 + t375 * t122 + t267 * qJD(6) + (-t94 + t295) * t168 + (t267 - t294) * t265 + (-t160 * t270 - t168 * t46 + t280) * m(7) + (t168 * t82 - t265 * t81 + t31) * m(6) + (-t84 * t168 + t265 * t83 + t156) * m(5); t123 * mrSges(6,1) + t307 * t248 + (t114 + t268) * t265 - m(7) * (t248 * t46 + t265 * t269) + (t248 * t82 + t265 * t88 + t27) * m(6) + t357; t117 - t46 * (mrSges(7,1) * t149 + mrSges(7,2) * t148) + (Ifges(7,1) * t148 - t313) * t346 + t56 * t345 + (Ifges(7,5) * t148 - Ifges(7,6) * t149) * t344 - t14 * t96 + t15 * t97 + (t14 * t148 + t149 * t15) * mrSges(7,3) + (-Ifges(7,2) * t149 + t146 + t57) * t348 + t356;];
tauc  = t5(:);

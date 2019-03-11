% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:24
% EndTime: 2019-03-08 23:05:48
% DurationCPUTime: 11.07s
% Computational Cost: add. (9975->596), mult. (25158->849), div. (0->0), fcn. (18931->12), ass. (0->290)
t284 = sin(qJ(4));
t285 = sin(qJ(3));
t288 = cos(qJ(4));
t289 = cos(qJ(3));
t248 = t284 * t285 - t288 * t289;
t278 = qJD(3) + qJD(4);
t206 = t278 * t248;
t249 = t284 * t289 + t285 * t288;
t207 = t278 * t249;
t286 = sin(qJ(2));
t280 = sin(pkin(6));
t331 = qJD(1) * t280;
t318 = t286 * t331;
t325 = qJD(3) * t285;
t321 = pkin(3) * t325;
t420 = pkin(4) * t207 + qJ(5) * t206 - qJD(5) * t249 - t318 + t321;
t387 = -pkin(9) - pkin(8);
t319 = qJD(3) * t387;
t252 = t285 * t319;
t253 = t289 * t319;
t290 = cos(qJ(2));
t317 = t290 * t331;
t262 = t387 * t285;
t263 = t387 * t289;
t397 = t288 * t262 + t263 * t284;
t411 = qJD(4) * t397 + t248 * t317 + t252 * t288 + t253 * t284;
t279 = sin(pkin(12));
t281 = cos(pkin(12));
t413 = -t411 * t279 + t281 * t420;
t412 = t279 * t420 + t411 * t281;
t277 = t281 * pkin(10);
t419 = pkin(5) * t207 + t206 * t277 + t413;
t344 = t206 * t279;
t418 = -pkin(10) * t344 - t412;
t275 = -pkin(3) * t289 - pkin(2);
t201 = pkin(4) * t248 - qJ(5) * t249 + t275;
t222 = t262 * t284 - t263 * t288;
t127 = t279 * t201 + t281 * t222;
t341 = t249 * t279;
t104 = -pkin(10) * t341 + t127;
t283 = sin(qJ(6));
t287 = cos(qJ(6));
t126 = t281 * t201 - t222 * t279;
t91 = pkin(5) * t248 - t249 * t277 + t126;
t28 = -t104 * t283 + t287 * t91;
t417 = qJD(6) * t28 + t283 * t419 - t418 * t287;
t29 = t104 * t287 + t283 * t91;
t416 = -qJD(6) * t29 + t418 * t283 + t287 * t419;
t241 = t249 * qJD(2);
t217 = t241 * t281 + t278 * t279;
t326 = qJD(2) * t289;
t328 = qJD(2) * t285;
t240 = t284 * t328 - t288 * t326;
t255 = qJD(2) * pkin(8) + t318;
t308 = pkin(9) * qJD(2) + t255;
t282 = cos(pkin(6));
t330 = qJD(1) * t285;
t316 = t282 * t330;
t212 = t289 * t308 + t316;
t204 = t288 * t212;
t334 = t282 * t289;
t267 = qJD(1) * t334;
t299 = t308 * t285;
t211 = t267 - t299;
t205 = qJD(3) * pkin(3) + t211;
t129 = t205 * t284 + t204;
t120 = qJ(5) * t278 + t129;
t230 = qJD(2) * t275 - t317;
t152 = pkin(4) * t240 - qJ(5) * t241 + t230;
t64 = -t120 * t279 + t281 * t152;
t37 = pkin(5) * t240 - pkin(10) * t217 + t64;
t307 = -t241 * t279 + t281 * t278;
t65 = t281 * t120 + t279 * t152;
t49 = pkin(10) * t307 + t65;
t12 = -t283 * t49 + t287 * t37;
t415 = t12 * mrSges(7,1);
t13 = t283 * t37 + t287 * t49;
t414 = t13 * mrSges(7,2);
t141 = qJD(4) * t222 + t252 * t284 - t288 * t253;
t219 = t249 * t317;
t410 = t141 - t219;
t247 = t279 * t287 + t281 * t283;
t166 = t247 * t240;
t235 = t247 * qJD(6);
t409 = t166 + t235;
t296 = t279 * t283 - t281 * t287;
t167 = t296 * t240;
t234 = t296 * qJD(6);
t408 = t167 + t234;
t200 = mrSges(5,1) * t240 + mrSges(5,2) * t241;
t407 = -m(5) * t230 - t200;
t406 = -t217 * t283 + t287 * t307;
t138 = t217 * t287 + t283 * t307;
t405 = -Ifges(4,1) / 0.2e1;
t196 = t206 * qJD(2);
t56 = qJD(6) * t406 + t196 * t296;
t391 = t56 / 0.2e1;
t57 = -qJD(6) * t138 + t196 * t247;
t390 = t57 / 0.2e1;
t197 = t207 * qJD(2);
t382 = t197 / 0.2e1;
t404 = -Ifges(4,4) * t326 / 0.2e1;
t345 = t196 * t281;
t346 = t196 * t279;
t110 = -mrSges(6,1) * t346 - mrSges(6,2) * t345;
t17 = -t57 * mrSges(7,1) + t56 * mrSges(7,2);
t403 = t110 + t17;
t258 = (-pkin(10) - qJ(5)) * t279;
t347 = qJ(5) * t281;
t259 = t277 + t347;
t214 = t258 * t283 + t259 * t287;
t342 = t240 * t281;
t304 = t241 * pkin(5) + pkin(10) * t342;
t203 = t284 * t212;
t128 = t205 * t288 - t203;
t198 = pkin(4) * t241 + qJ(5) * t240;
t83 = -t128 * t279 + t281 * t198;
t58 = t304 + t83;
t343 = t240 * t279;
t322 = pkin(10) * t343;
t84 = t281 * t128 + t279 * t198;
t66 = t322 + t84;
t402 = -qJD(5) * t247 - qJD(6) * t214 + t283 * t66 - t287 * t58;
t213 = t258 * t287 - t259 * t283;
t401 = -qJD(5) * t296 + qJD(6) * t213 - t283 * t58 - t287 * t66;
t119 = -pkin(4) * t278 + qJD(5) - t128;
t303 = mrSges(6,1) * t279 + mrSges(6,2) * t281;
t400 = t119 * t303;
t329 = qJD(2) * t280;
t309 = qJD(1) * t329;
t306 = t290 * t309;
t332 = qJD(3) * t267 + t289 * t306;
t163 = -qJD(3) * t299 + t332;
t293 = -qJD(3) * t212 - t285 * t306;
t324 = qJD(4) * t288;
t44 = -qJD(4) * t203 + t288 * t163 + t205 * t324 + t284 * t293;
t42 = qJD(5) * t278 + t44;
t233 = qJD(2) * t321 + t286 * t309;
t85 = pkin(4) * t197 + qJ(5) * t196 - qJD(5) * t241 + t233;
t22 = t279 * t85 + t281 * t42;
t16 = pkin(10) * t346 + t22;
t21 = -t279 * t42 + t281 * t85;
t6 = pkin(5) * t197 + pkin(10) * t345 + t21;
t2 = qJD(6) * t12 + t16 * t287 + t283 * t6;
t3 = -qJD(6) * t13 - t16 * t283 + t287 * t6;
t396 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t56 + Ifges(7,6) * t57;
t394 = m(5) / 0.2e1;
t393 = Ifges(7,4) * t391 + Ifges(7,2) * t390 + Ifges(7,6) * t382;
t392 = Ifges(7,1) * t391 + Ifges(7,4) * t390 + Ifges(7,5) * t382;
t232 = qJD(6) + t240;
t364 = Ifges(7,4) * t138;
t60 = Ifges(7,2) * t406 + Ifges(7,6) * t232 + t364;
t389 = t60 / 0.2e1;
t134 = Ifges(7,4) * t406;
t61 = Ifges(7,1) * t138 + Ifges(7,5) * t232 + t134;
t388 = t61 / 0.2e1;
t386 = -t406 / 0.2e1;
t385 = t406 / 0.2e1;
t384 = -t138 / 0.2e1;
t383 = t138 / 0.2e1;
t381 = -t232 / 0.2e1;
t380 = t232 / 0.2e1;
t379 = -t240 / 0.2e1;
t378 = t240 / 0.2e1;
t377 = -t241 / 0.2e1;
t376 = -t279 / 0.2e1;
t375 = t281 / 0.2e1;
t372 = pkin(3) * t288;
t371 = t21 * mrSges(6,3);
t370 = mrSges(5,3) * t240;
t369 = mrSges(5,3) * t241;
t368 = Ifges(4,4) * t285;
t367 = Ifges(5,4) * t241;
t366 = Ifges(6,4) * t279;
t365 = Ifges(6,4) * t281;
t363 = Ifges(6,2) * t279;
t362 = t406 * Ifges(7,6);
t361 = t138 * Ifges(7,5);
t336 = t280 * t286;
t236 = -t285 * t336 + t334;
t237 = t282 * t285 + t289 * t336;
t297 = t288 * t236 - t237 * t284;
t45 = qJD(4) * t129 + t163 * t284 - t288 * t293;
t360 = t297 * t45;
t359 = t21 * t279;
t358 = t307 * Ifges(6,6);
t357 = t217 * Ifges(6,5);
t356 = t22 * t281;
t355 = t397 * t45;
t354 = t232 * Ifges(7,3);
t353 = t241 * Ifges(5,1);
t352 = t278 * Ifges(5,5);
t351 = t278 * Ifges(5,6);
t350 = t279 * t64;
t349 = Ifges(4,5) * qJD(3);
t348 = Ifges(4,6) * qJD(3);
t340 = t255 * t289;
t266 = pkin(3) * t324 + qJD(5);
t338 = t266 * t281;
t271 = pkin(3) * t284 + qJ(5);
t337 = t271 * t281;
t335 = t280 * t290;
t131 = t211 * t288 - t203;
t172 = pkin(3) * t328 + t198;
t78 = t281 * t131 + t279 * t172;
t333 = -mrSges(5,1) * t278 - mrSges(6,1) * t307 + mrSges(6,2) * t217 + t369;
t327 = qJD(2) * t286;
t323 = qJD(2) * qJD(3);
t68 = -mrSges(7,1) * t406 + mrSges(7,2) * t138;
t320 = t68 + t333;
t272 = -pkin(5) * t281 - pkin(4);
t315 = t280 * t327;
t314 = t290 * t329;
t313 = t349 / 0.2e1;
t312 = -t348 / 0.2e1;
t311 = (t217 * Ifges(6,4) + Ifges(6,2) * t307 + t240 * Ifges(6,6)) * t376;
t310 = (t217 * Ifges(6,1) + Ifges(6,4) * t307 + t240 * Ifges(6,5)) * t375;
t77 = -t131 * t279 + t281 * t172;
t130 = t211 * t284 + t204;
t302 = Ifges(6,1) * t281 - t366;
t301 = -t363 + t365;
t300 = Ifges(6,5) * t281 - Ifges(6,6) * t279;
t179 = t236 * t284 + t237 * t288;
t148 = -t179 * t279 - t281 * t335;
t149 = t179 * t281 - t279 * t335;
t79 = t148 * t287 - t149 * t283;
t80 = t148 * t283 + t149 * t287;
t188 = -t255 * t325 + t332;
t189 = -qJD(3) * t340 + (-qJD(3) * t282 - t314) * t330;
t298 = t188 * t289 - t189 * t285;
t244 = (-pkin(10) - t271) * t279;
t245 = t277 + t337;
t190 = t244 * t287 - t245 * t283;
t191 = t244 * t283 + t245 * t287;
t223 = -t255 * t285 + t267;
t256 = -qJD(2) * pkin(2) - t317;
t295 = t223 * mrSges(4,3) + t328 * t405 + t404 - t349 / 0.2e1 - t256 * mrSges(4,2);
t224 = t316 + t340;
t294 = t224 * mrSges(4,3) + t348 / 0.2e1 + (t289 * Ifges(4,2) + t368) * qJD(2) / 0.2e1 - t256 * mrSges(4,1);
t116 = t240 * Ifges(6,3) + t357 + t358;
t176 = -t240 * Ifges(5,2) + t351 + t367;
t229 = Ifges(5,4) * t240;
t177 = -t229 + t352 + t353;
t27 = -pkin(5) * t346 + t45;
t59 = t354 + t361 + t362;
t73 = t197 * Ifges(6,6) - t196 * t301;
t74 = t197 * Ifges(6,5) - t196 * t302;
t88 = -pkin(5) * t307 + t119;
t292 = (-Ifges(5,1) * t240 + t116 - t367 + t59) * t377 - t241 * t415 + (mrSges(7,1) * t409 - mrSges(7,2) * t408) * t88 + (t12 * t408 - t13 * t409 - t2 * t296 - t247 * t3) * mrSges(7,3) + (-Ifges(5,2) * t241 + t177 - t229) * t378 + t240 * t310 + t240 * t311 - t217 * (Ifges(6,5) * t241 - t240 * t302) / 0.2e1 + (-mrSges(6,1) * t281 + mrSges(6,2) * t279 - mrSges(5,1)) * t45 + (Ifges(6,5) * t279 + Ifges(7,5) * t247 + Ifges(6,6) * t281 - Ifges(7,6) * t296) * t382 + t27 * (mrSges(7,1) * t296 + mrSges(7,2) * t247) + (Ifges(7,4) * t247 - Ifges(7,2) * t296) * t390 + (Ifges(7,1) * t247 - Ifges(7,4) * t296) * t391 - t296 * t393 - t65 * (-mrSges(6,2) * t241 + mrSges(6,3) * t343) - t64 * (mrSges(6,1) * t241 + mrSges(6,3) * t342) - t230 * (mrSges(5,1) * t241 - mrSges(5,2) * t240) - t307 * (Ifges(6,6) * t241 - t240 * t301) / 0.2e1 + t279 * t74 / 0.2e1 - t278 * (-Ifges(5,5) * t240 - Ifges(5,6) * t241) / 0.2e1 + t241 * t176 / 0.2e1 + (-Ifges(7,5) * t234 - Ifges(7,6) * t235) * t380 + (-Ifges(7,1) * t234 - Ifges(7,4) * t235) * t383 + (-Ifges(7,4) * t234 - Ifges(7,2) * t235) * t385 + t240 * t400 + mrSges(6,3) * t356 - t128 * t370 - t44 * mrSges(5,2) + t73 * t375 + (Ifges(6,3) * t241 - t240 * t300) * t379 + (Ifges(7,5) * t167 + Ifges(7,6) * t166 + Ifges(7,3) * t241) * t381 + (Ifges(7,1) * t167 + Ifges(7,4) * t166 + Ifges(7,5) * t241) * t384 + (Ifges(7,4) * t167 + Ifges(7,2) * t166 + Ifges(7,6) * t241) * t386 - t234 * t388 - t235 * t389 + t247 * t392 + t241 * t414 - t166 * t60 / 0.2e1 - t167 * t61 / 0.2e1 - Ifges(5,5) * t196 - Ifges(5,6) * t197 - (Ifges(6,1) * t279 + t365) * t345 / 0.2e1 + (Ifges(6,2) * t281 + t366) * t346 / 0.2e1;
t291 = qJD(2) ^ 2;
t274 = -pkin(4) - t372;
t261 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t326;
t260 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t328;
t257 = t272 - t372;
t243 = (mrSges(4,1) * t285 + mrSges(4,2) * t289) * t323;
t227 = pkin(5) * t343;
t225 = -mrSges(5,2) * t278 - t370;
t210 = qJD(3) * t236 + t289 * t314;
t209 = -qJD(3) * t237 - t285 * t314;
t192 = Ifges(7,3) * t197;
t183 = t296 * t249;
t182 = t247 * t249;
t168 = pkin(5) * t341 - t397;
t160 = mrSges(6,1) * t240 - mrSges(6,3) * t217;
t159 = -mrSges(6,2) * t240 + mrSges(6,3) * t307;
t125 = -qJD(6) * t191 - t247 * t266;
t124 = qJD(6) * t190 - t266 * t296;
t114 = mrSges(5,1) * t197 - mrSges(5,2) * t196;
t112 = mrSges(6,1) * t197 + mrSges(6,3) * t345;
t111 = -mrSges(6,2) * t197 + mrSges(6,3) * t346;
t109 = mrSges(7,1) * t232 - mrSges(7,3) * t138;
t108 = -mrSges(7,2) * t232 + mrSges(7,3) * t406;
t107 = t130 - t227;
t101 = t129 - t227;
t94 = -pkin(5) * t344 + t141;
t87 = qJD(4) * t179 - t288 * t209 + t210 * t284;
t86 = qJD(4) * t297 + t209 * t284 + t210 * t288;
t76 = t206 * t247 + t234 * t249;
t75 = t206 * t296 - t235 * t249;
t72 = t279 * t315 + t281 * t86;
t71 = -t279 * t86 + t281 * t315;
t63 = t322 + t78;
t52 = t304 + t77;
t33 = -mrSges(7,2) * t197 + mrSges(7,3) * t57;
t32 = mrSges(7,1) * t197 - mrSges(7,3) * t56;
t20 = t283 * t52 + t287 * t63;
t19 = -t283 * t63 + t287 * t52;
t8 = -qJD(6) * t80 - t283 * t72 + t287 * t71;
t7 = qJD(6) * t79 + t283 * t71 + t287 * t72;
t1 = [-t179 * t197 * mrSges(5,3) + t7 * t108 + t8 * t109 + t149 * t111 + t148 * t112 + t72 * t159 + t71 * t160 + t209 * t260 + t210 * t261 + t86 * t225 + t79 * t32 + t80 * t33 + (-t236 * t289 - t237 * t285) * mrSges(4,3) * t323 + t320 * t87 - (-t196 * mrSges(5,3) + t403) * t297 + ((-mrSges(3,2) * t291 - t114 - t243) * t290 + (-mrSges(3,1) * t291 + (t200 + qJD(2) * (-mrSges(4,1) * t289 + mrSges(4,2) * t285)) * qJD(2)) * t286) * t280 + m(7) * (t12 * t8 + t13 * t7 + t2 * t80 - t27 * t297 + t3 * t79 + t87 * t88) + m(6) * (t119 * t87 + t148 * t21 + t149 * t22 + t64 * t71 + t65 * t72 - t360) + m(5) * (-t128 * t87 + t129 * t86 - t360 + t179 * t44 + (t230 * t327 - t233 * t290) * t280) + m(4) * (t188 * t237 + t189 * t236 + t209 * t223 + t210 * t224 + (t256 - t317) * t315); -t397 * t110 + (t128 * t206 - t129 * t207 + t196 * t397 - t197 * t222 - t248 * t44 + t249 * t45) * mrSges(5,3) + (t64 * mrSges(6,1) + t358 / 0.2e1 + t357 / 0.2e1 - t65 * mrSges(6,2) - t351 / 0.2e1 + t59 / 0.2e1 + t415 - t414 + t116 / 0.2e1 + t362 / 0.2e1 + t361 / 0.2e1 - t176 / 0.2e1 + t230 * mrSges(5,1) + t354 / 0.2e1 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t240) * t207 + (-(pkin(2) * t327 + (-t223 * t285 + t224 * t289) * t290) * m(4) + (t260 * t285 - t261 * t289) * t290 + (-t256 * m(4) + t407) * t286) * t331 + m(4) * ((-t223 * t289 - t224 * t285) * qJD(3) + t298) * pkin(8) + (-0.3e1 / 0.2e1 * t285 ^ 2 + 0.3e1 / 0.2e1 * t289 ^ 2) * Ifges(4,4) * t323 + (t21 * mrSges(6,1) - t22 * mrSges(6,2) + t233 * mrSges(5,1) + t192 / 0.2e1 - t300 * t196 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3) / 0.2e1) * t197 + t396) * t248 + (t196 * t248 - t197 * t249 - t206 * t379 + t207 * t377) * Ifges(5,4) + (t300 * t382 + t45 * t303 + t233 * mrSges(5,2) + t73 * t376 + t74 * t375 + (-t21 * t281 - t22 * t279) * mrSges(6,3) - (Ifges(5,1) + Ifges(6,1) * t281 ^ 2 / 0.2e1 + (-t365 + t363 / 0.2e1) * t279) * t196) * t249 + (-t128 * t410 + t129 * t411 + t222 * t44 + t230 * t321 + t233 * t275 - t355) * m(5) + t411 * t225 - (t353 / 0.2e1 + t300 * t378 + t400 + t307 * t301 / 0.2e1 + t217 * t302 / 0.2e1 + t352 / 0.2e1 + t311 + t310 + t177 / 0.2e1 + t230 * mrSges(5,2) + (-t279 * t65 - t281 * t64) * mrSges(6,3)) * t206 + t416 * t109 + t417 * t108 + (t168 * t27 + t2 * t29 + t28 * t3 + (-t219 + t94) * t88 + t417 * t13 + t416 * t12) * m(7) + t298 * mrSges(4,3) + (-t12 * t75 + t13 * t76 - t182 * t2 + t183 * t3) * mrSges(7,3) + (-Ifges(7,5) * t183 - Ifges(7,6) * t182) * t382 + (-Ifges(7,4) * t183 - Ifges(7,2) * t182) * t390 + (-Ifges(7,1) * t183 - Ifges(7,4) * t182) * t391 + t27 * (mrSges(7,1) * t182 - mrSges(7,2) * t183) + t275 * t114 - pkin(2) * t243 + t28 * t32 + t29 * t33 + t412 * t159 + t413 * t160 + (t119 * t410 + t126 * t21 + t127 * t22 + t412 * t65 + t413 * t64 - t355) * m(6) - t320 * t219 + t88 * (-mrSges(7,1) * t76 + mrSges(7,2) * t75) + t94 * t68 + ((-pkin(8) * t260 - t295 + t313) * t289 + (pkin(3) * t200 - pkin(8) * t261 + t312 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t326 - t294) * t285) * qJD(3) + (Ifges(7,5) * t75 + Ifges(7,6) * t76) * t380 + (Ifges(7,1) * t75 + Ifges(7,4) * t76) * t383 + (Ifges(7,4) * t75 + Ifges(7,2) * t76) * t385 + t75 * t388 + t76 * t389 - t183 * t392 - t182 * t393 + t333 * t141 + t126 * t112 + t127 * t111 + t168 * t17; ((t196 * t288 - t197 * t284) * mrSges(5,3) + (t225 * t288 + t284 * t320) * qJD(4)) * pkin(3) + ((t313 + t404 + t295) * t289 + (t312 + (t368 / 0.2e1 + (t405 + Ifges(4,2) / 0.2e1) * t289) * qJD(2) + t407 * pkin(3) + t294) * t285) * qJD(2) + t292 - m(7) * (t107 * t88 + t12 * t19 + t13 * t20) - m(6) * (t119 * t130 + t64 * t77 + t65 * t78) - m(5) * (-t128 * t130 + t129 * t131) + 0.2e1 * ((t284 * t44 - t288 * t45) * t394 + ((-t128 * t284 + t129 * t288) * t394 + (m(6) * t119 + m(7) * t88) * t284 / 0.2e1) * qJD(4)) * pkin(3) + m(7) * (t12 * t125 + t124 * t13 + t190 * t3 + t191 * t2 + t257 * t27) + t274 * t110 + t224 * t260 - t223 * t261 + t257 * t17 + (-t20 + t124) * t108 + t111 * t337 + t129 * t369 + (-t19 + t125) * t109 - t107 * t68 - t333 * t130 + (-t78 + t338) * t159 - t77 * t160 - t188 * mrSges(4,2) + t189 * mrSges(4,1) + t190 * t32 + m(6) * (t22 * t337 - t266 * t350 - t271 * t359 + t274 * t45 + t338 * t65) + t191 * t33 + (-t112 * t271 - t160 * t266 - t371) * t279 - t131 * t225; t402 * t109 + t401 * t108 + t292 + (-qJ(5) * t112 - qJD(5) * t160 - t371) * t279 + (-t333 + t369) * t129 + t272 * t17 + (qJD(5) * t281 - t84) * t159 + t111 * t347 - t101 * t68 - pkin(4) * t110 - t83 * t160 + t213 * t32 + t214 * t33 - t128 * t225 + (-t101 * t88 + t12 * t402 + t13 * t401 + t2 * t214 + t213 * t3 + t27 * t272) * m(7) + (-t119 * t129 - t64 * t83 - t65 * t84 - pkin(4) * t45 + (t281 * t65 - t350) * qJD(5) + (t356 - t359) * qJ(5)) * m(6); -t406 * t108 + t138 * t109 - t307 * t159 + t217 * t160 + (t12 * t138 - t13 * t406 + t27) * m(7) + (t217 * t64 - t307 * t65 + t45) * m(6) + t403; t192 - t88 * (mrSges(7,1) * t138 + mrSges(7,2) * t406) + (Ifges(7,1) * t406 - t364) * t384 + t60 * t383 + (Ifges(7,5) * t406 - Ifges(7,6) * t138) * t381 - t12 * t108 + t13 * t109 + (t12 * t406 + t13 * t138) * mrSges(7,3) + (-Ifges(7,2) * t138 + t134 + t61) * t386 + t396;];
tauc  = t1(:);

% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:52
% EndTime: 2019-03-09 13:04:29
% DurationCPUTime: 19.86s
% Computational Cost: add. (10350->772), mult. (26626->1008), div. (0->0), fcn. (19034->8), ass. (0->337)
t444 = Ifges(6,1) + Ifges(7,1);
t443 = Ifges(7,4) + Ifges(6,5);
t442 = Ifges(7,5) - Ifges(6,4);
t456 = Ifges(6,6) - Ifges(7,6);
t252 = cos(pkin(6));
t348 = qJD(1) * t252;
t244 = qJD(2) + t348;
t257 = cos(qJ(4));
t254 = sin(qJ(4));
t258 = cos(qJ(2));
t251 = sin(pkin(6));
t349 = qJD(1) * t251;
t326 = t258 * t349;
t316 = t254 * t326;
t189 = t244 * t257 - t316;
t253 = sin(qJ(5));
t256 = cos(qJ(5));
t255 = sin(qJ(2));
t327 = t255 * t349;
t288 = qJD(4) + t327;
t126 = t253 * t189 - t256 * t288;
t127 = t256 * t189 + t253 * t288;
t335 = pkin(1) * t348;
t241 = t258 * t335;
t259 = -pkin(2) - pkin(9);
t412 = pkin(3) + pkin(8);
t119 = t244 * t259 + t327 * t412 + qJD(3) - t241;
t319 = -qJ(3) * t255 - pkin(1);
t166 = (t258 * t259 + t319) * t251;
t146 = qJD(1) * t166;
t76 = t257 * t119 - t254 * t146;
t61 = -pkin(4) * t288 - t76;
t24 = t126 * pkin(5) - t127 * qJ(6) + t61;
t77 = t254 * t119 + t257 * t146;
t62 = pkin(10) * t288 + t77;
t203 = pkin(8) * t326 + t255 * t335;
t175 = pkin(3) * t326 + t203;
t233 = t244 * qJ(3);
t140 = t233 + t175;
t188 = -t244 * t254 - t257 * t326;
t78 = -pkin(4) * t188 - pkin(10) * t189 + t140;
t21 = -t253 * t62 + t256 * t78;
t22 = t253 * t78 + t256 * t62;
t286 = t21 * t256 + t22 * t253;
t185 = qJD(5) - t188;
t435 = qJD(6) - t21;
t17 = -pkin(5) * t185 + t435;
t18 = qJ(6) * t185 + t22;
t287 = t17 * t256 - t18 * t253;
t305 = mrSges(7,1) * t253 - mrSges(7,3) * t256;
t307 = mrSges(6,1) * t253 + mrSges(6,2) * t256;
t389 = t256 / 0.2e1;
t391 = t253 / 0.2e1;
t392 = -t253 / 0.2e1;
t123 = Ifges(6,4) * t126;
t370 = Ifges(7,5) * t126;
t437 = t127 * t444 + t443 * t185 - t123 + t370;
t122 = Ifges(7,5) * t127;
t55 = Ifges(7,6) * t185 + Ifges(7,3) * t126 + t122;
t374 = Ifges(6,4) * t127;
t58 = -Ifges(6,2) * t126 + Ifges(6,6) * t185 + t374;
t457 = -t287 * mrSges(7,2) + t286 * mrSges(6,3) - t24 * t305 - t61 * t307 - t389 * t437 - t391 * t55 - t392 * t58;
t315 = t257 * t327;
t340 = qJD(5) * t253;
t321 = t259 * t340;
t311 = pkin(4) * t257 + pkin(10) * t254;
t214 = qJD(4) * t311 + qJD(3);
t219 = pkin(4) * t254 - pkin(10) * t257 + qJ(3);
t341 = qJD(4) * t257;
t322 = t259 * t341;
t339 = qJD(5) * t256;
t329 = t253 * t214 + t219 * t339 + t256 * t322;
t114 = t241 + (-t311 - t412) * t327;
t236 = pkin(2) * t327;
t289 = pkin(9) * t255 - qJ(3) * t258;
t173 = t289 * t349 + t236;
t102 = t257 * t173 + t254 * t175;
t90 = pkin(10) * t326 + t102;
t51 = t253 * t114 + t256 * t90;
t455 = -t51 + (qJD(6) - t321) * t254 + t329 + (t315 + t341) * qJ(6);
t318 = t253 * t259 - pkin(5);
t355 = t254 * t259;
t434 = t253 * t219 + t256 * t355;
t425 = -qJD(5) * t434 + t256 * t214;
t50 = t114 * t256 - t253 * t90;
t454 = -pkin(5) * t315 + t318 * t341 - t425 + t50;
t181 = t253 * t254 * t327 - t256 * t326;
t354 = t255 * t256;
t182 = (t253 * t258 + t254 * t354) * t349;
t290 = pkin(5) * t253 - qJ(6) * t256;
t280 = -t259 + t290;
t291 = pkin(5) * t256 + qJ(6) * t253;
t343 = qJD(4) * t254;
t101 = -t254 * t173 + t175 * t257;
t89 = -pkin(4) * t326 - t101;
t453 = -pkin(5) * t181 + qJ(6) * t182 + (qJD(5) * t291 - qJD(6) * t256) * t257 - t280 * t343 - t89;
t368 = Ifges(7,5) * t256;
t293 = Ifges(7,3) * t253 + t368;
t372 = Ifges(6,4) * t256;
t300 = -Ifges(6,2) * t253 + t372;
t396 = t185 / 0.2e1;
t403 = t127 / 0.2e1;
t405 = t126 / 0.2e1;
t406 = -t126 / 0.2e1;
t369 = Ifges(7,5) * t253;
t373 = Ifges(6,4) * t253;
t428 = t256 * t444 + t369 - t373;
t430 = -t253 * t456 + t256 * t443;
t452 = t293 * t405 + t300 * t406 + t396 * t430 + t403 * t428 - t457;
t202 = pkin(8) * t327 - t241;
t451 = -qJD(3) - t202;
t450 = t24 * mrSges(7,1) + t61 * mrSges(6,1) + t55 / 0.2e1 - t58 / 0.2e1 - t18 * mrSges(7,2) - t22 * mrSges(6,3);
t447 = -t244 / 0.2e1;
t449 = mrSges(4,1) + mrSges(3,3);
t448 = mrSges(4,2) - mrSges(3,1);
t320 = qJD(2) * t349;
t314 = t255 * t320;
t136 = -qJD(4) * t316 + t244 * t341 - t257 * t314;
t347 = qJD(2) * t255;
t323 = t254 * t347;
t135 = -t244 * t343 + (-t258 * t341 + t323) * t349;
t313 = t258 * t320;
t70 = -qJD(5) * t126 + t256 * t135 + t253 * t313;
t71 = qJD(5) * t127 + t253 * t135 - t256 * t313;
t438 = t443 * t136 + t442 * t71 + t444 * t70;
t446 = t327 / 0.2e1;
t445 = -t349 / 0.2e1;
t441 = Ifges(7,2) + Ifges(6,3);
t12 = Ifges(6,5) * t70 - Ifges(6,6) * t71 + Ifges(6,3) * t136;
t13 = Ifges(7,4) * t70 + Ifges(7,2) * t136 + Ifges(7,6) * t71;
t439 = t13 + t12;
t436 = -qJD(6) * t253 + t185 * t290 - t77;
t183 = Ifges(5,4) * t188;
t277 = Ifges(5,5) * t288;
t376 = Ifges(5,1) * t189;
t105 = t183 + t277 + t376;
t433 = t76 * mrSges(5,3) - t105 / 0.2e1 - t140 * mrSges(5,2) - t183 / 0.2e1;
t431 = t253 * t443 + t256 * t456;
t429 = t253 * t444 - t368 + t372;
t229 = pkin(2) * t314;
t345 = qJD(3) * t255;
t266 = (qJD(2) * t289 - t345) * t251;
t124 = qJD(1) * t266 + t229;
t248 = t252 * t255 * pkin(1);
t356 = t251 * t258;
t176 = (t356 * t412 + t248) * qJD(2);
t150 = qJD(1) * t176;
t27 = t119 * t341 + t257 * t124 - t146 * t343 + t254 * t150;
t25 = pkin(10) * t313 + t27;
t230 = qJD(2) * t241;
t232 = t244 * qJD(3);
t357 = t251 * t255;
t317 = t412 * t357;
t282 = qJD(2) * t317;
t125 = -qJD(1) * t282 + t230 + t232;
t54 = pkin(4) * t136 - pkin(10) * t135 + t125;
t3 = t256 * t25 + t253 * t54 + t78 * t339 - t340 * t62;
t4 = -qJD(5) * t22 - t25 * t253 + t256 * t54;
t309 = -t253 * t4 + t256 * t3;
t1 = qJ(6) * t136 + qJD(6) * t185 + t3;
t2 = -pkin(5) * t136 - t4;
t310 = t1 * t256 + t2 * t253;
t28 = -t119 * t343 - t254 * t124 - t146 * t341 + t150 * t257;
t427 = t28 * mrSges(5,1) - t27 * mrSges(5,2) + Ifges(5,5) * t135 - Ifges(5,6) * t136;
t426 = -Ifges(3,4) * t326 / 0.2e1 + Ifges(3,5) * t447;
t346 = qJD(2) * t258;
t324 = t251 * t346;
t325 = t251 * t347;
t238 = pkin(2) * t325;
t144 = t238 + t266;
t245 = pkin(8) * t357;
t388 = pkin(1) * t258;
t328 = -pkin(2) - t388;
t149 = pkin(3) * t357 + t245 + (-pkin(9) + t328) * t252;
t43 = t257 * t144 + t149 * t341 - t166 * t343 + t254 * t176;
t34 = pkin(10) * t324 + t43;
t211 = pkin(8) * t356 + t248;
t193 = -t252 * qJ(3) - t211;
t165 = pkin(3) * t356 - t193;
t208 = t252 * t254 + t257 * t356;
t330 = t254 * t356;
t209 = t252 * t257 - t330;
t100 = pkin(4) * t208 - pkin(10) * t209 + t165;
t93 = t254 * t149 + t257 * t166;
t87 = pkin(10) * t357 + t93;
t381 = t253 * t100 + t256 * t87;
t338 = t252 * t388;
t242 = qJD(2) * t338;
t249 = t252 * qJD(3);
t148 = t242 + t249 - t282;
t159 = -qJD(4) * t208 + t251 * t323;
t160 = -qJD(4) * t330 + t252 * t341 - t257 * t325;
t69 = pkin(4) * t160 - pkin(10) * t159 + t148;
t9 = -qJD(5) * t381 - t253 * t34 + t256 * t69;
t276 = Ifges(5,6) * t288;
t367 = Ifges(5,2) * t188;
t375 = Ifges(5,4) * t189;
t104 = t276 + t367 + t375;
t331 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t332 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t334 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t56 = t127 * Ifges(6,5) - t126 * Ifges(6,6) + t185 * Ifges(6,3);
t57 = t127 * Ifges(7,4) + t185 * Ifges(7,2) + t126 * Ifges(7,6);
t423 = -t332 * t126 - t334 * t127 - t331 * t185 + t17 * mrSges(7,1) + t22 * mrSges(6,2) + t77 * mrSges(5,3) + t104 / 0.2e1 - t56 / 0.2e1 - t57 / 0.2e1 + t375 / 0.2e1 - t140 * mrSges(5,1) - t18 * mrSges(7,3) - t21 * mrSges(6,1);
t422 = m(5) / 0.2e1;
t421 = m(6) / 0.2e1;
t420 = m(7) / 0.2e1;
t419 = Ifges(3,5) / 0.2e1;
t418 = Ifges(4,5) / 0.2e1;
t417 = Ifges(5,2) / 0.2e1;
t416 = t70 / 0.2e1;
t415 = -t71 / 0.2e1;
t414 = t71 / 0.2e1;
t411 = m(5) * t76;
t410 = m(5) * t77;
t409 = m(6) * t61;
t408 = pkin(1) * mrSges(3,1);
t407 = pkin(1) * mrSges(3,2);
t404 = -t127 / 0.2e1;
t402 = t136 / 0.2e1;
t397 = -t185 / 0.2e1;
t395 = -t208 / 0.2e1;
t393 = t209 / 0.2e1;
t390 = -t256 / 0.2e1;
t383 = -qJD(4) / 0.2e1;
t382 = qJD(4) / 0.2e1;
t94 = -mrSges(7,2) * t126 + mrSges(7,3) * t185;
t378 = mrSges(6,3) * t126;
t95 = -mrSges(6,2) * t185 - t378;
t380 = t94 + t95;
t377 = mrSges(6,3) * t127;
t96 = mrSges(6,1) * t185 - t377;
t97 = -mrSges(7,1) * t185 + mrSges(7,2) * t127;
t379 = -t96 + t97;
t366 = Ifges(4,6) * t258;
t364 = t135 * Ifges(5,1);
t363 = t135 * Ifges(5,4);
t362 = t136 * Ifges(5,4);
t360 = t244 * Ifges(4,5);
t115 = mrSges(5,1) * t313 - mrSges(5,3) * t135;
t20 = mrSges(6,1) * t71 + mrSges(6,2) * t70;
t359 = t115 - t20;
t118 = pkin(4) * t189 - pkin(10) * t188;
t39 = t253 * t118 + t256 * t76;
t142 = mrSges(5,1) * t288 - t189 * mrSges(5,3);
t81 = mrSges(6,1) * t126 + mrSges(6,2) * t127;
t358 = t81 - t142;
t352 = t256 * t219;
t117 = -mrSges(5,1) * t188 + mrSges(5,2) * t189;
t198 = -mrSges(4,1) * t326 - mrSges(4,3) * t244;
t351 = -t198 + t117;
t350 = t244 * t448 + t327 * t449;
t344 = qJD(4) * t253;
t342 = qJD(4) * t256;
t337 = t244 / 0.2e1 - qJD(2);
t333 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t47 = -t136 * mrSges(7,1) + t70 * mrSges(7,2);
t92 = t149 * t257 - t254 * t166;
t308 = mrSges(6,1) * t256 - mrSges(6,2) * t253;
t306 = mrSges(7,1) * t256 + mrSges(7,3) * t253;
t299 = Ifges(6,2) * t256 + t373;
t296 = Ifges(5,5) * t159 - Ifges(5,6) * t160;
t292 = -Ifges(7,3) * t256 + t369;
t36 = t100 * t256 - t253 * t87;
t38 = t118 * t256 - t253 * t76;
t281 = -qJD(4) - t327 / 0.2e1;
t204 = -pkin(8) * t325 + t242;
t44 = -t254 * t144 - t149 * t343 - t166 * t341 + t176 * t257;
t279 = -t209 * t253 + t251 * t354;
t162 = t209 * t256 + t253 * t357;
t8 = t100 * t339 + t253 * t69 + t256 * t34 - t340 * t87;
t190 = -pkin(8) * t314 + t230;
t86 = -pkin(4) * t357 - t92;
t194 = (-pkin(2) * t258 + t319) * t251;
t46 = mrSges(6,1) * t136 - mrSges(6,3) * t70;
t48 = -mrSges(6,2) * t136 - mrSges(6,3) * t71;
t49 = -mrSges(7,2) * t71 + mrSges(7,3) * t136;
t272 = (t48 + t49) * t256 + (-t46 + t47) * t253;
t271 = -t253 * t380 + t256 * t379;
t270 = (-qJ(3) * t346 - t345) * t251;
t178 = qJD(1) * t194;
t269 = Ifges(3,6) * t447 + (Ifges(3,4) * t255 + Ifges(3,2) * t258) * t445 + t360 / 0.2e1 + (-Ifges(4,6) * t255 - Ifges(4,3) * t258) * t349 / 0.2e1 - t178 * mrSges(4,2) - t203 * mrSges(3,3);
t205 = t211 * qJD(2);
t268 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t156 = -t190 - t232;
t191 = qJD(1) * t205;
t267 = -t190 * mrSges(3,2) - t156 * mrSges(4,3) + t191 * t448;
t35 = -pkin(4) * t324 - t44;
t26 = -pkin(4) * t313 - t28;
t264 = -t376 / 0.2e1 + t433;
t262 = t202 * mrSges(3,3) + t76 * mrSges(5,1) + Ifges(5,6) * t188 + Ifges(5,5) * t189 + Ifges(3,1) * t446 + Ifges(4,4) * t447 + (-Ifges(4,2) * t255 - t366) * t445 - t178 * mrSges(4,3) - t77 * mrSges(5,2) - t426 + (t288 / 0.2e1 + t382) * Ifges(5,3);
t261 = -t367 / 0.2e1 - t423;
t228 = Ifges(3,5) * t313;
t227 = Ifges(4,5) * t314;
t226 = Ifges(5,3) * t313;
t220 = -pkin(4) - t291;
t210 = -t245 + t338;
t201 = -qJ(3) * t326 + t236;
t200 = (mrSges(4,2) * t258 - mrSges(4,3) * t255) * t349;
t197 = -mrSges(3,2) * t244 + mrSges(3,3) * t326;
t195 = t252 * t328 + t245;
t192 = t280 * t257;
t186 = -t253 * t355 + t352;
t184 = -t204 - t249;
t180 = t244 * t253 - t256 * t315;
t179 = t244 * t256 + t253 * t315;
t177 = t238 + t270;
t174 = -qJD(1) * t317 + t241;
t172 = t254 * t318 - t352;
t171 = -t233 - t203;
t164 = qJ(6) * t254 + t434;
t163 = -pkin(2) * t244 - t451;
t153 = qJD(1) * t270 + t229;
t141 = -mrSges(5,2) * t288 + t188 * mrSges(5,3);
t116 = -mrSges(5,2) * t313 - mrSges(5,3) * t136;
t109 = -t253 * t322 + t425;
t108 = -t254 * t321 + t329;
t85 = qJD(5) * t279 + t159 * t256 + t253 * t324;
t84 = qJD(5) * t162 + t159 * t253 - t256 * t324;
t82 = mrSges(5,1) * t136 + mrSges(5,2) * t135;
t80 = mrSges(7,1) * t126 - mrSges(7,3) * t127;
t79 = pkin(5) * t127 + qJ(6) * t126;
t74 = Ifges(5,5) * t313 - t362 + t364;
t73 = -t136 * Ifges(5,2) + Ifges(5,6) * t313 + t363;
t45 = -pkin(5) * t279 - qJ(6) * t162 + t86;
t32 = -pkin(5) * t208 - t36;
t31 = qJ(6) * t208 + t381;
t30 = -pkin(5) * t189 - t38;
t29 = qJ(6) * t189 + t39;
t19 = mrSges(7,1) * t71 - mrSges(7,3) * t70;
t14 = Ifges(6,4) * t70 - Ifges(6,2) * t71 + Ifges(6,6) * t136;
t11 = Ifges(7,5) * t70 + Ifges(7,6) * t136 + Ifges(7,3) * t71;
t10 = pkin(5) * t84 - qJ(6) * t85 - qJD(6) * t162 + t35;
t7 = pkin(5) * t71 - qJ(6) * t70 - qJD(6) * t127 + t26;
t6 = -pkin(5) * t160 - t9;
t5 = qJ(6) * t160 + qJD(6) * t208 + t8;
t15 = [-(-t1 * mrSges(7,2) - t3 * mrSges(6,3) + Ifges(7,3) * t414 - Ifges(6,2) * t415 + t11 / 0.2e1 - t14 / 0.2e1 + t7 * mrSges(7,1) + t26 * mrSges(6,1) + t442 * t416 - t456 * t402) * t279 + (-Ifges(6,2) * t406 + Ifges(7,3) * t405 - t396 * t456 + t442 * t403 + t450) * t84 + t296 * t382 + m(7) * (t1 * t31 + t10 * t24 + t17 * t6 + t18 * t5 + t2 * t32 + t45 * t7) + m(5) * (t125 * t165 + t140 * t148 + t27 * t93 + t28 * t92 + t43 * t77 + t44 * t76) + m(4) * (t153 * t194 + t156 * t193 + t163 * t205 + t171 * t184 + t177 * t178 + t191 * t195) + m(3) * (t190 * t211 - t191 * t210 + t202 * t205 + t203 * t204) + t437 * t85 / 0.2e1 + t438 * t162 / 0.2e1 + t439 * t208 / 0.2e1 + t350 * t205 + (t1 * t208 + t160 * t18 - t162 * t7 - t24 * t85) * mrSges(7,3) + t74 * t393 + t73 * t395 + (t255 * t296 / 0.2e1 + ((-t194 * mrSges(4,3) + Ifges(5,5) * t393 + Ifges(5,6) * t395 + t195 * mrSges(4,1) - t210 * mrSges(3,3) + (-Ifges(4,4) + t419) * t252 + (t258 * t333 - 0.2e1 * t407) * t251) * t258 + (t193 * mrSges(4,1) - t194 * mrSges(4,2) - t211 * mrSges(3,3) + (-Ifges(3,6) + t418) * t252 + (-t255 * t333 - 0.2e1 * t408) * t251 + (Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t356) * t255) * qJD(2)) * t349 + (-t159 * t76 - t160 * t77 - t208 * t27 - t209 * t28) * mrSges(5,3) + (t57 + t56) * t160 / 0.2e1 + m(6) * (t21 * t9 + t22 * t8 + t26 * t86 + t3 * t381 + t35 * t61 + t36 * t4) + t381 * t48 + (Ifges(7,5) * t85 + Ifges(7,6) * t160) * t405 + (Ifges(7,5) * t162 + Ifges(7,6) * t208) * t414 + (-t160 * t22 + t162 * t26 - t208 * t3 + t61 * t85) * mrSges(6,2) + (t162 * t443 + t208 * t441) * t402 + (t160 * t441 + t443 * t85) * t396 + (t160 * t443 + t444 * t85) * t403 + (t162 * t444 + t208 * t443) * t416 + (Ifges(6,4) * t85 + Ifges(6,6) * t160) * t406 + (Ifges(6,4) * t162 + Ifges(6,6) * t208) * t415 + t188 * (Ifges(5,4) * t159 - Ifges(5,2) * t160) / 0.2e1 - t136 * (Ifges(5,4) * t209 - Ifges(5,2) * t208) / 0.2e1 + t189 * (Ifges(5,1) * t159 - Ifges(5,4) * t160) / 0.2e1 + t135 * (Ifges(5,1) * t209 - Ifges(5,4) * t208) / 0.2e1 + t125 * (mrSges(5,1) * t208 + mrSges(5,2) * t209) + t204 * t197 + t2 * (-mrSges(7,1) * t208 + mrSges(7,2) * t162) + t4 * (mrSges(6,1) * t208 - mrSges(6,3) * t162) + t184 * t198 + t177 * t200 + ((-mrSges(4,1) * t156 + mrSges(4,2) * t153 + mrSges(3,3) * t190) * t258 + (t226 / 0.2e1 - t153 * mrSges(4,3) + t449 * t191 + t427) * t255 + ((t171 * mrSges(4,1) + (t418 - Ifges(3,6) / 0.2e1) * t244 + t269) * t255 + (t163 * mrSges(4,1) + (-Ifges(4,4) / 0.2e1 + t419) * t244 + t262) * t258) * qJD(2)) * t251 + (t227 / 0.2e1 + t228 / 0.2e1 + t267) * t252 + t45 * t19 + t36 * t46 + t32 * t47 + t31 * t49 + t10 * t80 + t35 * t81 + t86 * t20 + t5 * t94 + t8 * t95 + t9 * t96 + t6 * t97 + t92 * t115 + t93 * t116 + t43 * t141 + t44 * t142 + t148 * t117 + t159 * t105 / 0.2e1 + t17 * (-mrSges(7,1) * t160 + mrSges(7,2) * t85) + t21 * (mrSges(6,1) * t160 - mrSges(6,3) * t85) - t160 * t104 / 0.2e1 + t140 * (mrSges(5,1) * t160 + mrSges(5,2) * t159) + t165 * t82; (-pkin(2) * t191 - qJ(3) * t156 - t163 * t203 + t171 * t451 - t178 * t201) * m(4) + (-t198 + t197) * t202 + t227 + t228 + ((t337 * Ifges(4,4) + (-pkin(2) * qJD(2) - t163) * mrSges(4,1) + qJD(2) * (Ifges(5,5) * t257 - Ifges(5,6) * t254) / 0.2e1 + (t407 - t366 / 0.2e1) * t349 - t262 + t426) * t258 + (-t360 / 0.2e1 + (t408 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t255) * t349 + t337 * Ifges(3,6) + (-qJ(3) * qJD(2) - t171) * mrSges(4,1) + (-Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t326 + (Ifges(5,5) * t281 + t264) * t254 + (Ifges(5,6) * t281 + t261) * t257 - t269) * t255) * t349 + m(6) * (t22 * t108 + t21 * t109 + t4 * t186 + t3 * t434) + t434 * t48 - t350 * t203 + t351 * qJD(3) - m(6) * (t21 * t50 + t22 * t51 + t61 * t89) - m(5) * (t101 * t76 + t102 * t77 + t140 * t174) + t453 * t80 + t454 * t97 + (t125 * mrSges(5,1) - t363 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1 - t73 / 0.2e1 - t27 * mrSges(5,3) + t332 * t71 + t334 * t70 + (m(5) * t27 + t116) * t259 + (t417 + t331) * t136 + t268) * t254 + t267 + (-Ifges(6,2) * t405 + Ifges(7,3) * t406 - t397 * t456 + t404 * t442 - t450) * t181 + (t1 * t164 + t17 * t454 + t172 * t2 + t18 * t455 + t192 * t7 + t24 * t453) * m(7) + t455 * t94 + m(5) * (t125 * qJ(3) + t140 * qJD(3)) + (-t50 + t109) * t96 + (-t51 + t108) * t95 - t201 * t200 + t192 * t19 + t186 * t46 + (((t141 + t410) * t259 + Ifges(5,6) * t383 + t261) * t257 + (t264 + Ifges(5,5) * t383 + (t358 + t409 - t411) * t259 + t293 * t406 + t300 * t405 + t428 * t404 + t430 * t397 + t457) * t254) * qJD(4) + (-mrSges(6,2) * t61 - mrSges(7,2) * t17 + mrSges(6,3) * t21 + mrSges(7,3) * t24 + Ifges(6,4) * t405 + Ifges(7,5) * t406 + t443 * t397 + t404 * t444 - t437 / 0.2e1) * t182 + (t293 * t414 + t300 * t415 + t7 * t305 + t26 * t307 + t11 * t391 + t14 * t392 + t74 / 0.2e1 + t364 / 0.2e1 - t362 / 0.2e1 - t28 * mrSges(5,3) + t125 * mrSges(5,2) + (-t253 * t3 - t256 * t4) * mrSges(6,3) + (-t1 * t253 + t2 * t256) * mrSges(7,2) + (m(5) * t28 - m(6) * t26 + t359) * t259 + (t292 * t406 + t299 * t405 + t24 * t306 + t61 * t308 + t58 * t390 + (t21 * t253 - t22 * t256) * mrSges(6,3) + (-t17 * t253 - t18 * t256) * mrSges(7,2) + t429 * t404 + t431 * t397 + t437 * t392) * qJD(5) + t428 * t416 + t430 * t402 + (qJD(5) * t55 + t438) * t389) * t257 + qJ(3) * t82 - t89 * t81 - t102 * t141 - t101 * t142 + t164 * t49 + t172 * t47 - t174 * t117; -t380 * t180 + t379 * t179 + (mrSges(4,1) * t346 + t200 * t255) * t349 + (t141 * t327 - t19 + (t253 * t379 + t256 * t380 + t141) * qJD(4) + t359) * t257 + (qJD(5) * t271 + t116 + t272 + t288 * (t80 + t358)) * t254 - m(6) * (t179 * t21 + t180 * t22) - m(7) * (-t17 * t179 + t18 * t180) + 0.2e1 * ((qJD(4) * t77 + t28) * t422 + t410 * t446 + (t17 * t344 + t18 * t342 - t7) * t420 + (-t21 * t344 + t22 * t342 - t26) * t421) * t257 + 0.2e1 * ((-qJD(4) * t76 + t27) * t422 + (qJD(4) * t24 + t17 * t339 - t18 * t340 + t310) * t420 + (qJD(4) * t61 - t21 * t339 - t22 * t340 + t309) * t421 + (-t411 / 0.2e1 + t24 * t420 + t409 / 0.2e1) * t327) * t254 + (-m(5) * t140 - t351) * t244 + (t171 * t244 + t178 * t327 + t191) * m(4); ((-Ifges(5,1) / 0.2e1 + t417) * t189 - t277 / 0.2e1 + t433 - t452) * t188 + t452 * qJD(5) - t358 * t77 + t226 + t436 * t80 + (-t17 * t30 - t18 * t29 + t220 * t7 + t436 * t24) * m(7) + t438 * t391 + (t276 / 0.2e1 + t423) * t189 + t14 * t389 + t11 * t390 + (-pkin(4) * t26 - t21 * t38 - t22 * t39 - t61 * t77) * m(6) + t310 * mrSges(7,2) + t309 * mrSges(6,3) + (m(6) * t309 + m(7) * t310 + (-m(6) * t286 + m(7) * t287 + t271) * qJD(5) + t272) * pkin(10) + t429 * t416 + t431 * t402 - t7 * t306 - t26 * t308 + t220 * t19 + t292 * t414 + t299 * t415 + t427 - pkin(4) * t20 - t29 * t94 - t39 * t95 - t38 * t96 - t30 * t97 - t76 * t141; (t377 - t379) * t22 + (-t378 - t380) * t21 + t268 + (t126 * t17 + t127 * t18) * mrSges(7,2) - pkin(5) * t47 + qJ(6) * t49 - t79 * t80 + qJD(6) * t94 - t24 * (mrSges(7,1) * t127 + mrSges(7,3) * t126) + (Ifges(7,3) * t127 - t370) * t406 + t58 * t403 - t61 * (mrSges(6,1) * t127 - mrSges(6,2) * t126) + (-t126 * t443 - t127 * t456) * t397 + (-pkin(5) * t2 + qJ(6) * t1 - t17 * t22 + t18 * t435 - t24 * t79) * m(7) + (-Ifges(6,2) * t127 - t123 + t437) * t405 + (-t126 * t444 + t122 - t374 + t55) * t404 + t439; t127 * t80 - t185 * t94 + 0.2e1 * (t2 / 0.2e1 + t24 * t403 + t18 * t397) * m(7) + t47;];
tauc  = t15(:);

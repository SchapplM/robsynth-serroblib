% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:18
% EndTime: 2019-03-09 16:35:36
% DurationCPUTime: 10.26s
% Computational Cost: add. (13983->613), mult. (36732->788), div. (0->0), fcn. (26933->8), ass. (0->290)
t456 = Ifges(6,1) + Ifges(7,1);
t452 = Ifges(6,5) + Ifges(7,4);
t455 = Ifges(6,6) - Ifges(7,6);
t354 = qJD(2) + qJD(3);
t457 = t354 / 0.2e1;
t288 = sin(qJ(5));
t291 = cos(qJ(5));
t446 = -t455 * t288 + t452 * t291;
t388 = Ifges(7,5) * t288;
t391 = Ifges(6,4) * t288;
t444 = t456 * t291 + t388 - t391;
t293 = cos(qJ(2));
t427 = -pkin(8) - pkin(7);
t273 = t427 * t293;
t265 = qJD(1) * t273;
t289 = sin(qJ(3));
t250 = t289 * t265;
t290 = sin(qJ(2));
t272 = t427 * t290;
t264 = qJD(1) * t272;
t257 = qJD(2) * pkin(2) + t264;
t292 = cos(qJ(3));
t213 = t292 * t257 + t250;
t260 = t289 * t293 + t292 * t290;
t249 = t260 * qJD(1);
t242 = t249 * qJ(4);
t179 = t213 - t242;
t170 = pkin(3) * t354 + t179;
t253 = t292 * t265;
t214 = t257 * t289 - t253;
t365 = t292 * t293;
t259 = -t289 * t290 + t365;
t248 = t259 * qJD(1);
t378 = qJ(4) * t248;
t180 = t214 + t378;
t286 = sin(pkin(10));
t171 = t286 * t180;
t287 = cos(pkin(10));
t106 = t287 * t170 - t171;
t281 = -pkin(2) * t293 - pkin(1);
t271 = qJD(1) * t281;
t222 = -t248 * pkin(3) + qJD(4) + t271;
t454 = t222 * mrSges(5,2) - t106 * mrSges(5,3);
t341 = t287 * t248 - t249 * t286;
t310 = t248 * t286 + t287 * t249;
t393 = Ifges(5,4) * t310;
t453 = Ifges(5,6) * t457 + Ifges(5,2) * t341 / 0.2e1 + t393 / 0.2e1;
t220 = t354 * t259;
t208 = t220 * qJD(1);
t221 = t354 * t260;
t209 = t221 * qJD(1);
t154 = t208 * t286 + t287 * t209;
t155 = t208 * t287 - t209 * t286;
t184 = t288 * t310 - t291 * t354;
t95 = -qJD(5) * t184 + t291 * t155;
t185 = t288 * t354 + t291 * t310;
t96 = qJD(5) * t185 + t288 * t155;
t451 = (-Ifges(6,4) + Ifges(7,5)) * t96 + t456 * t95 + t452 * t154;
t182 = Ifges(6,4) * t184;
t195 = qJD(5) - t341;
t389 = Ifges(7,5) * t184;
t438 = t456 * t185 + t452 * t195 - t182 + t389;
t218 = -t264 * t289 + t253;
t186 = t218 - t378;
t219 = t292 * t264 + t250;
t187 = -t242 + t219;
t368 = t287 * t289;
t385 = pkin(2) * qJD(3);
t434 = -t287 * t186 + t187 * t286 - (t286 * t292 + t368) * t385;
t355 = qJD(5) * t291;
t356 = qJD(5) * t288;
t373 = t341 * t291;
t374 = t341 * t288;
t449 = -qJD(6) * t288 + (-t355 + t373) * qJ(6) + (t356 - t374) * pkin(5);
t104 = -pkin(4) * t354 - t106;
t331 = mrSges(7,1) * t288 - mrSges(7,3) * t291;
t333 = mrSges(6,1) * t288 + mrSges(6,2) * t291;
t40 = t184 * pkin(5) - t185 * qJ(6) + t104;
t448 = t104 * t333 + t40 * t331;
t447 = t452 * t288 + t455 * t291;
t387 = Ifges(7,5) * t291;
t390 = Ifges(6,4) * t291;
t445 = t456 * t288 - t387 + t390;
t369 = t287 * t180;
t107 = t286 * t170 + t369;
t105 = pkin(9) * t354 + t107;
t113 = -pkin(4) * t341 - pkin(9) * t310 + t222;
t261 = t289 * t272;
t297 = (t365 * t427 - t261) * qJD(2) * qJD(1);
t309 = -t208 * qJ(4) - t249 * qJD(4);
t357 = qJD(3) * t292;
t358 = qJD(3) * t289;
t347 = qJD(2) * t427;
t157 = t249 * t347 + t257 * t357 + t265 * t358;
t99 = -qJ(4) * t209 + qJD(4) * t248 + t157;
t32 = t287 * t99 + (-t257 * t358 + t265 * t357 + t297 + t309) * t286;
t361 = qJD(1) * t290;
t283 = pkin(2) * t361;
t190 = pkin(3) * t209 + qJD(2) * t283;
t55 = pkin(4) * t154 - pkin(9) * t155 + t190;
t6 = -t105 * t356 + t113 * t355 + t288 * t55 + t291 * t32;
t2 = qJ(6) * t154 + qJD(6) * t195 + t6;
t37 = -t105 * t288 + t113 * t291;
t435 = qJD(6) - t37;
t24 = -pkin(5) * t195 + t435;
t38 = t105 * t291 + t113 * t288;
t7 = -qJD(5) * t38 - t288 * t32 + t291 * t55;
t4 = -pkin(5) * t154 - t7;
t443 = t2 * t291 + t24 * t355 + t288 * t4;
t320 = Ifges(7,3) * t288 + t387;
t326 = -Ifges(6,2) * t288 + t390;
t412 = t288 / 0.2e1;
t413 = -t288 / 0.2e1;
t420 = t185 / 0.2e1;
t422 = t184 / 0.2e1;
t423 = -t184 / 0.2e1;
t181 = Ifges(7,5) * t185;
t84 = t195 * Ifges(7,6) + t184 * Ifges(7,3) + t181;
t392 = Ifges(6,4) * t185;
t87 = -t184 * Ifges(6,2) + t195 * Ifges(6,6) + t392;
t442 = t320 * t422 + t326 * t423 + t84 * t412 + t87 * t413 + t448 + t444 * t420 + t446 * t195 / 0.2e1;
t441 = -t222 * mrSges(5,1) - t37 * mrSges(6,1) + t24 * mrSges(7,1) + t38 * mrSges(6,2) + t453;
t417 = -t310 / 0.2e1;
t440 = -t341 / 0.2e1;
t439 = -t354 / 0.2e1;
t394 = Ifges(5,4) * t341;
t110 = t179 * t286 + t369;
t437 = -t110 + t449;
t436 = t449 - t434;
t119 = t186 * t286 + t187 * t287;
t370 = t286 * t289;
t239 = (t287 * t292 - t370) * t385;
t433 = -t119 + t239;
t216 = -t287 * t259 + t260 * t286;
t217 = t259 * t286 + t260 * t287;
t229 = -t259 * pkin(3) + t281;
t135 = t216 * pkin(4) - t217 * pkin(9) + t229;
t224 = t292 * t272 + t273 * t289;
t198 = -qJ(4) * t260 + t224;
t225 = -t292 * t273 + t261;
t199 = qJ(4) * t259 + t225;
t137 = t198 * t286 + t199 * t287;
t432 = t288 * t135 + t291 * t137;
t25 = qJ(6) * t195 + t38;
t402 = t7 * t288;
t403 = t291 * t6;
t431 = m(7) * (-t25 * t356 + t443) + m(6) * (-t355 * t37 - t356 * t38 - t402 + t403);
t266 = t290 * t347;
t267 = t293 * t347;
t166 = t292 * t266 + t289 * t267 + t272 * t357 + t273 * t358;
t116 = -qJ(4) * t221 + qJD(4) * t259 + t166;
t167 = -qJD(3) * t225 - t266 * t289 + t292 * t267;
t117 = -qJ(4) * t220 - qJD(4) * t260 + t167;
t45 = t116 * t287 + t117 * t286;
t159 = t220 * t286 + t287 * t221;
t160 = t220 * t287 - t221 * t286;
t359 = qJD(2) * t290;
t205 = pkin(2) * t359 + pkin(3) * t221;
t71 = pkin(4) * t159 - pkin(9) * t160 + t205;
t12 = -qJD(5) * t432 - t288 * t45 + t291 * t71;
t430 = t95 / 0.2e1;
t429 = -t96 / 0.2e1;
t428 = t96 / 0.2e1;
t426 = pkin(1) * mrSges(3,1);
t425 = pkin(1) * mrSges(3,2);
t424 = t154 / 0.2e1;
t421 = -t185 / 0.2e1;
t419 = -t195 / 0.2e1;
t415 = t248 / 0.2e1;
t414 = t249 / 0.2e1;
t411 = -t291 / 0.2e1;
t410 = t291 / 0.2e1;
t409 = m(4) * t271;
t408 = pkin(3) * t249;
t407 = pkin(3) * t287;
t406 = pkin(5) * t310;
t46 = mrSges(6,1) * t154 - mrSges(6,3) * t95;
t47 = -t154 * mrSges(7,1) + t95 * mrSges(7,2);
t401 = -t46 + t47;
t48 = -mrSges(6,2) * t154 - mrSges(6,3) * t96;
t49 = -mrSges(7,2) * t96 + mrSges(7,3) * t154;
t400 = t48 + t49;
t399 = mrSges(6,3) * t184;
t398 = mrSges(6,3) * t185;
t397 = Ifges(5,1) * t310;
t396 = Ifges(3,4) * t290;
t395 = Ifges(4,4) * t249;
t136 = -t287 * t198 + t199 * t286;
t158 = -qJD(3) * t214 + t297;
t31 = t286 * t99 - t287 * (t158 + t309);
t383 = t136 * t31;
t382 = t248 * mrSges(4,3);
t381 = t249 * mrSges(4,3);
t380 = Ifges(3,5) * qJD(2);
t379 = Ifges(3,6) * qJD(2);
t377 = qJD(2) * mrSges(3,1);
t376 = qJD(2) * mrSges(3,2);
t375 = t107 * t310;
t372 = t239 * t288;
t371 = t239 * t291;
t111 = t179 * t287 - t171;
t125 = pkin(4) * t310 - pkin(9) * t341 + t408;
t51 = t291 * t111 + t288 * t125;
t120 = t125 + t283;
t53 = t291 * t119 + t288 * t120;
t364 = -mrSges(5,1) * t354 + mrSges(6,1) * t184 + mrSges(6,2) * t185 + mrSges(5,3) * t310;
t127 = -mrSges(7,2) * t184 + mrSges(7,3) * t195;
t128 = -mrSges(6,2) * t195 - t399;
t363 = t127 + t128;
t129 = mrSges(6,1) * t195 - t398;
t130 = -mrSges(7,1) * t195 + mrSges(7,2) * t185;
t362 = t129 - t130;
t280 = pkin(2) * t292 + pkin(3);
t241 = pkin(2) * t368 + t286 * t280;
t360 = qJD(1) * t293;
t353 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t352 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t351 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t346 = t380 / 0.2e1;
t345 = -t379 / 0.2e1;
t342 = t154 * mrSges(5,1) + t155 * mrSges(5,2);
t44 = t116 * t286 - t287 * t117;
t240 = -pkin(2) * t370 + t280 * t287;
t340 = Ifges(5,5) * t354;
t336 = -t2 * t288 + t291 * t4;
t335 = -t288 * t6 - t291 * t7;
t334 = mrSges(6,1) * t291 - mrSges(6,2) * t288;
t332 = mrSges(7,1) * t291 + mrSges(7,3) * t288;
t325 = Ifges(6,2) * t291 + t391;
t319 = -Ifges(7,3) * t291 + t388;
t318 = pkin(5) * t291 + qJ(6) * t288;
t317 = pkin(5) * t288 - qJ(6) * t291;
t316 = t24 * t288 + t25 * t291;
t314 = -t288 * t38 - t291 * t37;
t313 = t288 * t37 - t291 * t38;
t50 = -t111 * t288 + t125 * t291;
t52 = -t119 * t288 + t120 * t291;
t64 = t135 * t291 - t137 * t288;
t308 = -pkin(4) - t318;
t307 = -qJD(5) * t362 + t400;
t11 = t135 * t355 - t137 * t356 + t288 * t71 + t291 * t45;
t298 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t13 = pkin(5) * t96 - qJ(6) * t95 - qJD(6) * t185 + t31;
t140 = t340 + t394 + t397;
t19 = t95 * Ifges(7,5) + t154 * Ifges(7,6) + t96 * Ifges(7,3);
t196 = Ifges(4,2) * t248 + t354 * Ifges(4,6) + t395;
t243 = Ifges(4,4) * t248;
t197 = Ifges(4,1) * t249 + t354 * Ifges(4,5) + t243;
t20 = t95 * Ifges(6,4) - t96 * Ifges(6,2) + t154 * Ifges(6,6);
t85 = t185 * Ifges(6,5) - t184 * Ifges(6,6) + t195 * Ifges(6,3);
t86 = t185 * Ifges(7,4) + t195 * Ifges(7,2) + t184 * Ifges(7,6);
t294 = (t140 + t394) * t440 + (Ifges(5,1) * t417 + Ifges(5,5) * t439 + t320 * t423 + t326 * t422 + t419 * t446 + t421 * t444 - t448 - t454) * t341 - t249 * (Ifges(4,1) * t248 - t395) / 0.2e1 + t442 * qJD(5) + t213 * t382 + t445 * t430 + t447 * t424 + (Ifges(4,5) * t248 - Ifges(4,6) * t249) * t439 - (-Ifges(4,2) * t249 + t197 + t243) * t248 / 0.2e1 + (t37 * t373 + t374 * t38 + t403) * mrSges(6,3) - t13 * t332 + (t355 / 0.2e1 - t373 / 0.2e1) * t438 + (-t393 + t86 + t85) * t417 + (-t334 - mrSges(5,1)) * t31 + (-t24 * t373 + t443) * mrSges(7,2) - t271 * (mrSges(4,1) * t249 + mrSges(4,2) * t248) + t451 * t412 + (-Ifges(5,6) * t439 - Ifges(5,2) * t440 + Ifges(6,6) * t422 + Ifges(7,6) * t423 + t452 * t421 + (Ifges(6,3) + Ifges(7,2)) * t419 + t441) * t310 + t20 * t410 + t19 * t411 + t196 * t414 - t32 * mrSges(5,2) + (t87 / 0.2e1 - t84 / 0.2e1) * t374 + t319 * t428 + t325 * t429 - Ifges(5,6) * t154 + Ifges(5,5) * t155 - t157 * mrSges(4,2) + t158 * mrSges(4,1) + Ifges(4,5) * t208 - Ifges(4,6) * t209;
t282 = Ifges(3,4) * t360;
t279 = -pkin(4) - t407;
t269 = mrSges(3,3) * t360 - t376;
t268 = -mrSges(3,3) * t361 + t377;
t256 = t308 - t407;
t247 = Ifges(3,1) * t361 + t282 + t380;
t246 = t379 + (t293 * Ifges(3,2) + t396) * qJD(1);
t233 = -pkin(4) - t240;
t228 = mrSges(4,1) * t354 - t381;
t227 = -mrSges(4,2) * t354 + t382;
t226 = t283 + t408;
t223 = -t240 + t308;
t211 = -mrSges(4,1) * t248 + mrSges(4,2) * t249;
t193 = t310 * qJ(6);
t188 = -mrSges(5,2) * t354 + mrSges(5,3) * t341;
t152 = Ifges(7,2) * t154;
t150 = Ifges(6,3) * t154;
t148 = -mrSges(5,1) * t341 + mrSges(5,2) * t310;
t146 = -mrSges(7,2) * t374 + mrSges(7,3) * t310;
t122 = mrSges(7,1) * t184 - mrSges(7,3) * t185;
t121 = pkin(5) * t185 + qJ(6) * t184;
t94 = Ifges(7,4) * t95;
t93 = Ifges(6,5) * t95;
t92 = Ifges(6,6) * t96;
t91 = Ifges(7,6) * t96;
t74 = t217 * t317 + t136;
t43 = -pkin(5) * t216 - t64;
t42 = qJ(6) * t216 + t432;
t36 = -t52 - t406;
t35 = t193 + t53;
t34 = -t50 - t406;
t33 = t193 + t51;
t27 = mrSges(6,1) * t96 + mrSges(6,2) * t95;
t26 = mrSges(7,1) * t96 - mrSges(7,3) * t95;
t14 = t317 * t160 + (qJD(5) * t318 - qJD(6) * t291) * t217 + t44;
t9 = -pkin(5) * t159 - t12;
t8 = qJ(6) * t159 + qJD(6) * t216 + t11;
t1 = [(t140 / 0.2e1 + t442 + t438 * t410 + (t24 * t291 - t25 * t288) * mrSges(7,2) + t314 * mrSges(6,3) + t397 / 0.2e1 + t394 / 0.2e1 + t340 / 0.2e1 + t454) * t160 + (-t246 / 0.2e1 - pkin(7) * t269 + t345 + (-0.2e1 * t426 - 0.3e1 / 0.2e1 * t396 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t293) * qJD(1) + (0.2e1 * t409 + qJD(1) * (-mrSges(4,1) * t259 + mrSges(4,2) * t260) + t211) * pkin(2)) * t359 + m(5) * (-t106 * t44 + t107 * t45 + t137 * t32 + t190 * t229 + t205 * t222 + t383) + t364 * t44 + (t260 * t208 + t220 * t414) * Ifges(4,1) + (t93 / 0.2e1 - t92 / 0.2e1 + t150 / 0.2e1 + t94 / 0.2e1 + t152 / 0.2e1 + t91 / 0.2e1 - Ifges(5,4) * t155 + t190 * mrSges(5,1) - t32 * mrSges(5,3) + t351 * t96 + t353 * t95 + (Ifges(5,2) + t352) * t154 + t298) * t216 + (t136 * t155 - t137 * t154) * mrSges(5,3) + m(6) * (t104 * t44 + t11 * t38 + t12 * t37 + t432 * t6 + t64 * t7 + t383) + t432 * t48 + t271 * (mrSges(4,1) * t221 + mrSges(4,2) * t220) + (t247 / 0.2e1 - pkin(7) * t268 + t346 + (-0.2e1 * t425 + 0.3e1 / 0.2e1 * Ifges(3,4) * t293) * qJD(1)) * qJD(2) * t293 + t229 * t342 + m(4) * (t157 * t225 + t158 * t224 + t166 * t214 + t167 * t213) + m(7) * (t13 * t74 + t14 * t40 + t2 * t42 + t24 * t9 + t25 * t8 + t4 * t43) + (Ifges(5,1) * t155 - Ifges(5,4) * t154 + t19 * t412 + t190 * mrSges(5,2) + t13 * t331 + t320 * t428 + t326 * t429 + (mrSges(5,3) + t333) * t31 + t335 * mrSges(6,3) + t336 * mrSges(7,2) + (-mrSges(7,2) * t316 + mrSges(6,3) * t313 + t104 * t334 + t319 * t423 + t325 * t422 + t332 * t40 + t411 * t87 + t419 * t447 + t421 * t445) * qJD(5) + t444 * t430 + t446 * t424 + (qJD(5) * t438 + t20) * t413 + (qJD(5) * t84 + t451) * t410) * t217 + (t85 / 0.2e1 + t86 / 0.2e1 - t107 * mrSges(5,3) + t25 * mrSges(7,3) + t352 * t195 + t353 * t185 + t351 * t184 - t441 - t453) * t159 + t43 * t47 + t42 * t49 + t64 * t46 + (-t259 * t209 - t221 * t415) * Ifges(4,2) + (t259 * t208 - t260 * t209 + t220 * t415 - t221 * t414) * Ifges(4,4) + (t157 * t259 - t158 * t260 - t208 * t224 - t209 * t225 - t213 * t220 - t214 * t221) * mrSges(4,3) + t281 * (mrSges(4,1) * t209 + mrSges(4,2) * t208) + t74 * t26 + (Ifges(4,5) * t220 - Ifges(4,6) * t221) * t457 + t14 * t122 + t8 * t127 + t11 * t128 + t12 * t129 + t9 * t130 + t136 * t27 + t45 * t188 + t205 * t148 + t220 * t197 / 0.2e1 - t221 * t196 / 0.2e1 + t166 * t227 + t167 * t228; ((-qJD(5) * t363 + t401) * t288 + t307 * t291 + t431) * (pkin(9) + t241) + t433 * t188 + (t106 * t434 + t107 * t433 - t222 * t226 - t240 * t31 + t241 * t32) * m(5) + (t233 * t31 + (t371 - t53) * t38 + (-t372 - t52) * t37 - t434 * t104) * m(6) - t434 * t364 + ((-t282 / 0.2e1 - t247 / 0.2e1 + t346 + qJD(1) * t425 + (t268 - t377) * pkin(7)) * t293 + (t246 / 0.2e1 + t345 + (t426 + t396 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t293) * qJD(1) + (t269 + t376) * pkin(7) + (-t211 - t409) * pkin(2)) * t290) * qJD(1) + (-t7 * mrSges(6,3) - t362 * t239 + (-t25 * mrSges(7,2) - t38 * mrSges(6,3)) * qJD(5)) * t288 + (-t37 * qJD(5) * mrSges(6,3) + t239 * t363) * t291 + ((t157 * t289 + t158 * t292 + (-t213 * t289 + t214 * t292) * qJD(3)) * pkin(2) - t213 * t218 - t214 * t219) * m(4) + t214 * t381 + t294 + (-t154 * t241 - t155 * t240 + t375) * mrSges(5,3) + (t13 * t223 + t436 * t40 + (-t35 + t371) * t25 + (-t36 + t372) * t24) * m(7) + t436 * t122 + ((t227 * t292 - t228 * t289) * qJD(3) + (-t208 * t292 - t209 * t289) * mrSges(4,3)) * pkin(2) - t35 * t127 - t53 * t128 - t52 * t129 - t36 * t130 - t25 * t146 + t223 * t26 - t226 * t148 - t219 * t227 - t218 * t228 + t233 * t27; -t364 * t110 + (t228 + t381) * t214 + t294 + t437 * t122 + (-mrSges(7,2) * t356 - t146) * t25 - t148 * t408 + (t375 + (-t154 * t286 - t155 * t287) * pkin(3)) * mrSges(5,3) + (qJD(5) * t314 - t402) * mrSges(6,3) + (t400 * t291 + t401 * t288 + (-t288 * t363 - t291 * t362) * qJD(5) + t431) * (pkin(3) * t286 + pkin(9)) - t33 * t127 - t51 * t128 - t50 * t129 - t34 * t130 - t111 * t188 - t213 * t227 + t256 * t26 + t279 * t27 + (t13 * t256 - t24 * t34 - t25 * t33 + t40 * t437) * m(7) + (-t104 * t110 + t279 * t31 - t37 * t50 - t38 * t51) * m(6) + (t106 * t110 - t107 * t111 - t222 * t408 + (t286 * t32 - t287 * t31) * pkin(3)) * m(5); -t341 * t188 + (-t122 - t364) * t310 + (t195 * t363 - t401) * t291 + (t341 * t362 + t307) * t288 + t342 + (t195 * t316 - t310 * t40 - t336) * m(7) + (-t104 * t310 - t195 * t313 - t335) * m(6) + (t106 * t310 - t107 * t341 + t190) * m(5); (t362 + t398) * t38 + (-t363 - t399) * t37 + t152 + t150 + t298 + t94 - pkin(5) * t47 + qJ(6) * t49 + t87 * t420 + (Ifges(7,3) * t185 - t389) * t423 - t121 * t122 + qJD(6) * t127 + t91 - t92 + t93 - t40 * (mrSges(7,1) * t185 + mrSges(7,3) * t184) - t104 * (mrSges(6,1) * t185 - mrSges(6,2) * t184) + (t184 * t24 + t185 * t25) * mrSges(7,2) + (-t452 * t184 - t455 * t185) * t419 + (-pkin(5) * t4 + qJ(6) * t2 - t121 * t40 - t24 * t38 + t25 * t435) * m(7) + (-Ifges(6,2) * t185 - t182 + t438) * t422 + (-t456 * t184 + t181 - t392 + t84) * t421; t185 * t122 - t195 * t127 + 0.2e1 * (t4 / 0.2e1 + t40 * t420 + t25 * t419) * m(7) + t47;];
tauc  = t1(:);

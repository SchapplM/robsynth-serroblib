% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:32:52
% EndTime: 2018-11-23 17:33:01
% DurationCPUTime: 8.85s
% Computational Cost: add. (17249->634), mult. (45545->874), div. (0->0), fcn. (34383->10), ass. (0->304)
t330 = sin(qJ(3));
t331 = sin(qJ(2));
t333 = cos(qJ(3));
t334 = cos(qJ(2));
t365 = t333 * t334;
t296 = -t330 * t331 + t365;
t282 = t296 * qJD(1);
t297 = t330 * t334 + t333 * t331;
t283 = t297 * qJD(1);
t327 = sin(pkin(10));
t385 = cos(pkin(10));
t231 = -t385 * t282 + t283 * t327;
t326 = sin(pkin(11));
t328 = cos(pkin(11));
t329 = sin(qJ(6));
t332 = cos(qJ(6));
t338 = t326 * t329 - t328 * t332;
t164 = t338 * t231;
t278 = t338 * qJD(6);
t441 = -t278 - t164;
t418 = t231 / 0.2e1;
t460 = -t326 / 0.2e1;
t325 = qJD(2) + qJD(3);
t337 = t327 * t282 + t283 * t385;
t213 = t325 * t326 + t328 * t337;
t430 = -pkin(8) - pkin(7);
t311 = t430 * t334;
t302 = qJD(1) * t311;
t286 = t330 * t302;
t310 = t430 * t331;
t301 = qJD(1) * t310;
t292 = qJD(2) * pkin(2) + t301;
t243 = t333 * t292 + t286;
t275 = t283 * qJ(4);
t208 = t243 - t275;
t199 = pkin(3) * t325 + t208;
t289 = t333 * t302;
t244 = t292 * t330 - t289;
t382 = qJ(4) * t282;
t209 = t244 + t382;
t352 = t385 * t209;
t129 = t327 * t199 + t352;
t125 = qJ(5) * t325 + t129;
t320 = -pkin(2) * t334 - pkin(1);
t309 = qJD(1) * t320;
t252 = -t282 * pkin(3) + qJD(4) + t309;
t138 = t231 * pkin(4) - qJ(5) * t337 + t252;
t72 = -t125 * t326 + t328 * t138;
t37 = pkin(5) * t231 - pkin(9) * t213 + t72;
t349 = t328 * t325 - t326 * t337;
t73 = t328 * t125 + t326 * t138;
t44 = pkin(9) * t349 + t73;
t12 = -t329 * t44 + t332 * t37;
t459 = t12 * mrSges(7,1);
t13 = t329 * t37 + t332 * t44;
t458 = t13 * mrSges(7,2);
t457 = mrSges(5,3) * t231;
t456 = t231 * Ifges(5,4);
t295 = t326 * t332 + t328 * t329;
t163 = t295 * t231;
t279 = t295 * qJD(6);
t455 = t279 + t163;
t454 = -t213 * t329 + t332 * t349;
t142 = t213 * t332 + t329 * t349;
t250 = t325 * t296;
t239 = t250 * qJD(1);
t251 = t325 * t297;
t240 = t251 * qJD(1);
t185 = t239 * t385 - t327 * t240;
t48 = qJD(6) * t454 - t185 * t338;
t434 = t48 / 0.2e1;
t49 = -qJD(6) * t142 - t185 * t295;
t433 = t49 / 0.2e1;
t224 = qJD(6) + t231;
t400 = Ifges(7,4) * t142;
t60 = Ifges(7,2) * t454 + Ifges(7,6) * t224 + t400;
t432 = t60 / 0.2e1;
t136 = Ifges(7,4) * t454;
t61 = Ifges(7,1) * t142 + Ifges(7,5) * t224 + t136;
t431 = t61 / 0.2e1;
t184 = t239 * t327 + t240 * t385;
t423 = t184 / 0.2e1;
t453 = -t337 / 0.2e1;
t452 = (Ifges(6,4) * t213 + Ifges(6,2) * t349 + Ifges(6,6) * t231) * t460;
t451 = Ifges(5,4) * t337;
t408 = pkin(3) * t327;
t317 = qJ(5) + t408;
t284 = (-pkin(9) - t317) * t326;
t324 = t328 * pkin(9);
t285 = t317 * t328 + t324;
t235 = t284 * t329 + t285 * t332;
t374 = t231 * t328;
t346 = pkin(5) * t337 + pkin(9) * t374;
t202 = t327 * t209;
t135 = t208 * t385 - t202;
t409 = pkin(3) * t283;
t154 = pkin(4) * t337 + qJ(5) * t231 + t409;
t81 = -t135 * t326 + t328 * t154;
t41 = t346 + t81;
t375 = t231 * t326;
t358 = pkin(9) * t375;
t82 = t328 * t135 + t326 * t154;
t62 = t358 + t82;
t450 = -qJD(5) * t295 - qJD(6) * t235 + t329 * t62 - t332 * t41;
t234 = t284 * t332 - t285 * t329;
t449 = -qJD(5) * t338 + qJD(6) * t234 - t329 * t41 - t332 * t62;
t319 = pkin(2) * t333 + pkin(3);
t351 = t385 * t330;
t274 = pkin(2) * t351 + t327 * t319;
t266 = qJ(5) + t274;
t253 = (-pkin(9) - t266) * t326;
t254 = t266 * t328 + t324;
t193 = t253 * t329 + t254 * t332;
t313 = t327 * t330 * pkin(2);
t359 = qJD(3) * t333;
t272 = t385 * pkin(2) * t359 - qJD(3) * t313;
t263 = qJD(5) + t272;
t248 = -t301 * t330 + t289;
t215 = t248 - t382;
t249 = t333 * t301 + t286;
t216 = -t275 + t249;
t149 = t327 * t215 + t216 * t385;
t363 = qJD(1) * t331;
t322 = pkin(2) * t363;
t151 = t154 + t322;
t83 = -t149 * t326 + t328 * t151;
t42 = t346 + t83;
t84 = t328 * t149 + t326 * t151;
t63 = t358 + t84;
t448 = -qJD(6) * t193 - t263 * t295 + t329 * t63 - t332 * t42;
t192 = t253 * t332 - t254 * t329;
t447 = qJD(6) * t192 - t263 * t338 - t329 * t42 - t332 * t63;
t128 = t199 * t385 - t202;
t124 = -t325 * pkin(4) + qJD(5) - t128;
t345 = mrSges(6,1) * t326 + mrSges(6,2) * t328;
t446 = t124 * t345;
t148 = -t385 * t215 + t216 * t327;
t271 = (t327 * t333 + t351) * qJD(3) * pkin(2);
t443 = t148 - t271;
t442 = -t149 + t272;
t364 = -mrSges(5,1) * t325 - mrSges(6,1) * t349 + mrSges(6,2) * t213 + mrSges(5,3) * t337;
t298 = t330 * t310;
t256 = -t333 * t311 + t298;
t159 = -mrSges(6,2) * t231 + mrSges(6,3) * t349;
t160 = mrSges(6,1) * t231 - mrSges(6,3) * t213;
t439 = t328 * t159 - t326 * t160;
t356 = qJD(2) * t430;
t360 = qJD(3) * t330;
t188 = t283 * t356 + t292 * t359 + t302 * t360;
t117 = -qJ(4) * t240 + qJD(4) * t282 + t188;
t189 = -t244 * qJD(3) + (t365 * t430 - t298) * qJD(2) * qJD(1);
t336 = -t239 * qJ(4) - t283 * qJD(4) + t189;
t58 = t385 * t117 + t327 * t336;
t54 = qJD(5) * t325 + t58;
t219 = pkin(3) * t240 + qJD(2) * t322;
t71 = pkin(4) * t184 - qJ(5) * t185 - qJD(5) * t337 + t219;
t23 = t326 * t71 + t328 * t54;
t378 = t185 * t326;
t14 = -pkin(9) * t378 + t23;
t22 = -t326 * t54 + t328 * t71;
t377 = t185 * t328;
t6 = pkin(5) * t184 - pkin(9) * t377 + t22;
t2 = qJD(6) * t12 + t14 * t332 + t329 * t6;
t3 = -qJD(6) * t13 - t14 * t329 + t332 * t6;
t438 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t48 + Ifges(7,6) * t49;
t393 = t213 * Ifges(6,5);
t394 = t349 * Ifges(6,6);
t109 = t231 * Ifges(6,3) + t393 + t394;
t386 = t325 * Ifges(5,6);
t170 = -t231 * Ifges(5,2) + t386 + t451;
t387 = t325 * Ifges(5,5);
t389 = t337 * Ifges(5,1);
t171 = t387 + t389 - t456;
t403 = Ifges(4,4) * t283;
t225 = t282 * Ifges(4,2) + t325 * Ifges(4,6) + t403;
t276 = Ifges(4,4) * t282;
t226 = t283 * Ifges(4,1) + t325 * Ifges(4,5) + t276;
t57 = t117 * t327 - t385 * t336;
t34 = pkin(5) * t378 + t57;
t342 = Ifges(6,5) * t328 - Ifges(6,6) * t326;
t399 = Ifges(6,2) * t326;
t401 = Ifges(6,4) * t328;
t343 = -t399 + t401;
t402 = Ifges(6,4) * t326;
t344 = Ifges(6,1) * t328 - t402;
t369 = t328 * (t213 * Ifges(6,1) + Ifges(6,4) * t349 + t231 * Ifges(6,5));
t388 = t283 * mrSges(4,3);
t390 = t23 * t328;
t405 = mrSges(4,3) * t282;
t411 = t328 / 0.2e1;
t413 = t283 / 0.2e1;
t419 = -t231 / 0.2e1;
t421 = t224 / 0.2e1;
t422 = -t224 / 0.2e1;
t424 = t142 / 0.2e1;
t425 = -t142 / 0.2e1;
t426 = t454 / 0.2e1;
t427 = -t454 / 0.2e1;
t435 = Ifges(7,1) * t434 + Ifges(7,4) * t433 + Ifges(7,5) * t423;
t436 = Ifges(7,4) * t434 + Ifges(7,2) * t433 + Ifges(7,6) * t423;
t391 = t224 * Ifges(7,3);
t396 = t142 * Ifges(7,5);
t397 = t454 * Ifges(7,6);
t59 = t391 + t396 + t397;
t69 = t184 * Ifges(6,6) + t185 * t343;
t70 = t184 * Ifges(6,5) + t185 * t344;
t97 = -pkin(5) * t349 + t124;
t437 = t244 * t388 + mrSges(6,3) * t390 + t337 * t458 + t243 * t405 + (-t12 * t441 - t13 * t455 - t2 * t338 - t295 * t3) * mrSges(7,3) + t69 * t411 + t225 * t413 - t309 * (mrSges(4,1) * t283 + mrSges(4,2) * t282) + t295 * t435 - t337 * t459 + (-Ifges(5,1) * t231 + t109 - t451 + t59) * t453 + t337 * t170 / 0.2e1 - t72 * (mrSges(6,1) * t337 + mrSges(6,3) * t374) - t73 * (-mrSges(6,2) * t337 + mrSges(6,3) * t375) - t455 * t432 + (mrSges(7,1) * t455 + mrSges(7,2) * t441) * t97 + t441 * t431 + (-Ifges(7,5) * t278 - Ifges(7,6) * t279) * t421 + (-Ifges(7,1) * t278 - Ifges(7,4) * t279) * t424 + (-Ifges(7,4) * t278 - Ifges(7,2) * t279) * t426 - t128 * t457 + (-Ifges(5,2) * t337 + t171 + t369 - t456) * t418 - t58 * mrSges(5,2) + (Ifges(7,5) * t164 + Ifges(7,6) * t163 + Ifges(7,3) * t337) * t422 + (Ifges(7,1) * t164 + Ifges(7,4) * t163 + Ifges(7,5) * t337) * t425 + (Ifges(7,4) * t164 + Ifges(7,2) * t163 + Ifges(7,6) * t337) * t427 - (-Ifges(4,2) * t283 + t226 + t276) * t282 / 0.2e1 + (Ifges(6,3) * t337 - t231 * t342) * t419 - t252 * (mrSges(5,1) * t337 - mrSges(5,2) * t231) - t349 * (Ifges(6,6) * t337 - t231 * t343) / 0.2e1 - t213 * (Ifges(6,5) * t337 - t231 * t344) / 0.2e1 + t231 * t446 - (Ifges(4,5) * t282 - Ifges(5,5) * t231 - Ifges(4,6) * t283 - Ifges(5,6) * t337) * t325 / 0.2e1 + t231 * t452 - t338 * t436 + (Ifges(7,4) * t295 - Ifges(7,2) * t338) * t433 + (Ifges(7,1) * t295 - Ifges(7,4) * t338) * t434 + t34 * (mrSges(7,1) * t338 + mrSges(7,2) * t295) + (Ifges(6,5) * t326 + Ifges(7,5) * t295 + Ifges(6,6) * t328 - Ifges(7,6) * t338) * t423 - Ifges(5,6) * t184 + Ifges(5,5) * t185 - t188 * mrSges(4,2) + t189 * mrSges(4,1) + Ifges(4,5) * t239 - Ifges(4,6) * t240 + t326 * t70 / 0.2e1 + (-mrSges(6,1) * t328 + mrSges(6,2) * t326 - mrSges(5,1)) * t57 + (Ifges(6,1) * t326 + t401) * t377 / 0.2e1 - (Ifges(6,2) * t328 + t402) * t378 / 0.2e1 - t283 * (Ifges(4,1) * t282 - t403) / 0.2e1;
t429 = pkin(1) * mrSges(3,1);
t428 = pkin(1) * mrSges(3,2);
t414 = t282 / 0.2e1;
t410 = m(4) * t309;
t406 = t328 * pkin(5);
t303 = t331 * t356;
t304 = t334 * t356;
t195 = t333 * t303 + t330 * t304 + t310 * t359 + t311 * t360;
t146 = -qJ(4) * t251 + qJD(4) * t296 + t195;
t196 = -qJD(3) * t256 - t303 * t330 + t333 * t304;
t147 = -qJ(4) * t250 - qJD(4) * t297 + t196;
t80 = t146 * t385 + t327 * t147;
t190 = t250 * t327 + t251 * t385;
t191 = t250 * t385 - t327 * t251;
t361 = qJD(2) * t331;
t236 = pkin(2) * t361 + pkin(3) * t251;
t247 = t327 * t296 + t297 * t385;
t87 = pkin(4) * t190 - qJ(5) * t191 - qJD(5) * t247 + t236;
t28 = t326 * t87 + t328 * t80;
t404 = Ifges(3,4) * t331;
t255 = t333 * t310 + t311 * t330;
t227 = -qJ(4) * t297 + t255;
t228 = qJ(4) * t296 + t256;
t166 = -t385 * t227 + t228 * t327;
t395 = t166 * t57;
t392 = t22 * t326;
t384 = Ifges(3,5) * qJD(2);
t383 = Ifges(3,6) * qJD(2);
t381 = qJD(2) * mrSges(3,1);
t380 = qJD(2) * mrSges(3,2);
t379 = t129 * t337;
t376 = t191 * t326;
t373 = t247 * t326;
t100 = -mrSges(6,2) * t184 - mrSges(6,3) * t378;
t370 = t328 * t100;
t246 = -t296 * t385 + t297 * t327;
t261 = -t296 * pkin(3) + t320;
t165 = t246 * pkin(4) - t247 * qJ(5) + t261;
t167 = t327 * t227 + t228 * t385;
t93 = t326 * t165 + t328 * t167;
t99 = mrSges(6,1) * t378 + mrSges(6,2) * t377;
t362 = qJD(1) * t334;
t78 = -mrSges(7,1) * t454 + mrSges(7,2) * t142;
t357 = t78 + t364;
t355 = t385 * pkin(3);
t354 = t384 / 0.2e1;
t353 = -t383 / 0.2e1;
t20 = -t49 * mrSges(7,1) + t48 * mrSges(7,2);
t27 = -t326 * t80 + t328 * t87;
t350 = t184 * mrSges(5,1) + t185 * mrSges(5,2);
t79 = t146 * t327 - t385 * t147;
t92 = t328 * t165 - t167 * t326;
t134 = t208 * t327 + t352;
t318 = -t355 - pkin(4);
t341 = t22 * t328 + t23 * t326;
t340 = t390 - t392;
t339 = -t326 * t72 + t328 * t73;
t64 = pkin(5) * t246 - t247 * t324 + t92;
t76 = -pkin(9) * t373 + t93;
t25 = -t329 * t76 + t332 * t64;
t26 = t329 * t64 + t332 * t76;
t273 = t319 * t385 - t313;
t267 = -pkin(4) - t273;
t321 = Ifges(3,4) * t362;
t307 = mrSges(3,3) * t362 - t380;
t306 = -mrSges(3,3) * t363 + t381;
t305 = t318 - t406;
t281 = Ifges(3,1) * t363 + t321 + t384;
t280 = t383 + (t334 * Ifges(3,2) + t404) * qJD(1);
t260 = t267 - t406;
t259 = mrSges(4,1) * t325 - t388;
t258 = -mrSges(4,2) * t325 + t405;
t257 = t322 + t409;
t242 = -mrSges(4,1) * t282 + mrSges(4,2) * t283;
t221 = pkin(5) * t375;
t217 = -mrSges(5,2) * t325 - t457;
t183 = t338 * t247;
t182 = t295 * t247;
t179 = Ifges(7,3) * t184;
t177 = mrSges(5,1) * t231 + mrSges(5,2) * t337;
t120 = pkin(5) * t373 + t166;
t106 = t148 - t221;
t105 = mrSges(7,1) * t224 - mrSges(7,3) * t142;
t104 = -mrSges(7,2) * t224 + mrSges(7,3) * t454;
t103 = t134 - t221;
t101 = mrSges(6,1) * t184 - mrSges(6,3) * t377;
t75 = -t191 * t295 + t247 * t278;
t74 = -t191 * t338 - t247 * t279;
t43 = pkin(5) * t376 + t79;
t31 = -mrSges(7,2) * t184 + mrSges(7,3) * t49;
t30 = mrSges(7,1) * t184 - mrSges(7,3) * t48;
t24 = -pkin(9) * t376 + t28;
t15 = pkin(5) * t190 - t191 * t324 + t27;
t5 = -qJD(6) * t26 + t15 * t332 - t24 * t329;
t4 = qJD(6) * t25 + t15 * t329 + t24 * t332;
t1 = [(-Ifges(7,5) * t183 - Ifges(7,6) * t182) * t423 + (-Ifges(7,4) * t183 - Ifges(7,2) * t182) * t433 + (-Ifges(7,1) * t183 - Ifges(7,4) * t182) * t434 + (-t12 * t74 + t13 * t75 - t182 * t2 + t183 * t3) * mrSges(7,3) + t34 * (mrSges(7,1) * t182 - mrSges(7,2) * t183) + (Ifges(7,5) * t74 + Ifges(7,6) * t75) * t421 + (Ifges(7,1) * t74 + Ifges(7,4) * t75) * t424 + (Ifges(7,4) * t74 + Ifges(7,2) * t75) * t426 + t74 * t431 + t75 * t432 - t183 * t435 - t182 * t436 + (t459 - t458 + t59 / 0.2e1 + t109 / 0.2e1 + t397 / 0.2e1 + t396 / 0.2e1 - t170 / 0.2e1 + t391 / 0.2e1 + t252 * mrSges(5,1) - t386 / 0.2e1 - t73 * mrSges(6,2) + t72 * mrSges(6,1) + t394 / 0.2e1 + t393 / 0.2e1 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t231) * t190 + (t389 / 0.2e1 + t369 / 0.2e1 + t452 + t171 / 0.2e1 + t252 * mrSges(5,2) + t387 / 0.2e1 + t446 + t349 * t343 / 0.2e1 + t213 * t344 / 0.2e1 + t342 * t418 + (-t326 * t73 - t328 * t72) * mrSges(6,3)) * t191 + (t188 * t296 - t189 * t297 - t239 * t255 - t240 * t256 - t243 * t250 - t244 * t251) * mrSges(4,3) + (-t296 * t240 - t251 * t414) * Ifges(4,2) + (t296 * t239 - t297 * t240 + t250 * t414 - t251 * t413) * Ifges(4,4) + t320 * (mrSges(4,1) * t240 + mrSges(4,2) * t239) + m(4) * (t188 * t256 + t189 * t255 + t195 * t244 + t196 * t243) + t43 * t78 + (-t128 * t191 - t129 * t190 + t166 * t185 - t167 * t184 - t246 * t58 + t247 * t57) * mrSges(5,3) + (t179 / 0.2e1 + t219 * mrSges(5,1) - t23 * mrSges(6,2) + t22 * mrSges(6,1) + t342 * t185 + (Ifges(5,2) + Ifges(7,3) / 0.2e1 + Ifges(6,3)) * t184 + t438) * t246 + (t70 * t411 + t69 * t460 + t219 * mrSges(5,2) + t57 * t345 + t342 * t423 - t341 * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1) * t328 ^ 2 / 0.2e1 + (-t401 + t399 / 0.2e1) * t326) * t185) * t247 + m(7) * (t12 * t5 + t120 * t34 + t13 * t4 + t2 * t26 + t25 * t3 + t43 * t97) + t309 * (mrSges(4,1) * t251 + mrSges(4,2) * t250) + t325 * (Ifges(4,5) * t250 - Ifges(4,6) * t251) / 0.2e1 + (-t247 * t184 - t246 * t185 + t190 * t453 + t191 * t419) * Ifges(5,4) + t25 * t30 + t26 * t31 + t97 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + t93 * t100 + t92 * t101 + t4 * t104 + t5 * t105 + t120 * t20 + t28 * t159 + t27 * t160 + t166 * t99 + (t281 / 0.2e1 - pkin(7) * t306 + t354 + (-0.2e1 * t428 + 0.3e1 / 0.2e1 * Ifges(3,4) * t334) * qJD(1)) * t334 * qJD(2) + t80 * t217 + t236 * t177 + t261 * t350 + t250 * t226 / 0.2e1 - t251 * t225 / 0.2e1 + t195 * t258 + t196 * t259 + t364 * t79 + m(6) * (t124 * t79 + t22 * t92 + t23 * t93 + t27 * t72 + t28 * t73 + t395) + m(5) * (-t128 * t79 + t129 * t80 + t167 * t58 + t219 * t261 + t236 * t252 + t395) + (t297 * t239 + t250 * t413) * Ifges(4,1) + (-t280 / 0.2e1 - pkin(7) * t307 + t353 + (-0.2e1 * t429 - 0.3e1 / 0.2e1 * t404 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t334) * qJD(1) + (0.2e1 * t410 + qJD(1) * (-mrSges(4,1) * t296 + mrSges(4,2) * t297) + t242) * pkin(2)) * t361; t266 * t370 + (-t22 * mrSges(6,3) - t266 * t101 - t263 * t160) * t326 + (-t124 * t443 + t263 * t339 + t266 * t340 + t267 * t57 - t72 * t83 - t73 * t84) * m(6) + ((t258 * t333 - t259 * t330) * qJD(3) + (-t239 * t333 - t240 * t330) * mrSges(4,3)) * pkin(2) + t447 * t104 + (t192 * t3 + t193 * t2 + t260 * t34 + (-t106 + t271) * t97 + t447 * t13 + t448 * t12) * m(7) + t448 * t105 + (t263 * t328 - t84) * t159 + t442 * t217 + (t128 * t443 + t129 * t442 - t252 * t257 - t273 * t57 + t274 * t58) * m(5) + ((t188 * t330 + t189 * t333 + (-t243 * t330 + t244 * t333) * qJD(3)) * pkin(2) - t243 * t248 - t244 * t249) * m(4) + t437 - t106 * t78 - t83 * t160 + t192 * t30 + t193 * t31 - t257 * t177 - t249 * t258 - t248 * t259 + t260 * t20 + t267 * t99 + t357 * t271 - t364 * t148 + (-t184 * t274 - t185 * t273 + t379) * mrSges(5,3) + ((-t321 / 0.2e1 - t281 / 0.2e1 + t354 + qJD(1) * t428 + (t306 - t381) * pkin(7)) * t334 + (t280 / 0.2e1 + t353 + (t429 + t404 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t334) * qJD(1) + (t307 + t380) * pkin(7) + (-t242 - t410) * pkin(2)) * t331) * qJD(1); (qJD(5) * t339 - t124 * t134 + t317 * t340 + t318 * t57 - t72 * t81 - t73 * t82) * m(6) + t450 * t105 + ((t327 * t58 - t385 * t57) * pkin(3) + t128 * t134 - t129 * t135 - t252 * t409) * m(5) + (-t103 * t97 + t12 * t450 + t13 * t449 + t2 * t235 + t234 * t3 + t305 * t34) * m(7) - t364 * t134 + t449 * t104 + t439 * qJD(5) + (-t184 * t408 - t185 * t355 + t379) * mrSges(5,3) + t437 - t103 * t78 - t82 * t159 - t81 * t160 - t135 * t217 + t234 * t30 + t235 * t31 + (-t326 * t101 + t370) * t317 - t243 * t258 + t244 * t259 + t305 * t20 + t318 * t99 - mrSges(6,3) * t392 - t177 * t409; t326 * t100 + t328 * t101 - t338 * t30 + t295 * t31 - t455 * t105 + t441 * t104 - (-t217 - t439) * t231 - t357 * t337 + t350 + (-t12 * t455 + t13 * t441 + t2 * t295 - t3 * t338 - t337 * t97) * m(7) + (-t124 * t337 + t231 * t339 + t341) * m(6) + (t128 * t337 + t129 * t231 + t219) * m(5); -t454 * t104 + t142 * t105 - t349 * t159 + t213 * t160 + t20 + t99 + (t12 * t142 - t13 * t454 + t34) * m(7) + (t213 * t72 - t349 * t73 + t57) * m(6); t179 - t97 * (mrSges(7,1) * t142 + mrSges(7,2) * t454) + (Ifges(7,1) * t454 - t400) * t425 + t60 * t424 + (Ifges(7,5) * t454 - Ifges(7,6) * t142) * t422 - t12 * t104 + t13 * t105 + (t12 * t454 + t13 * t142) * mrSges(7,3) + (-Ifges(7,2) * t142 + t136 + t61) * t427 + t438;];
tauc  = t1(:);

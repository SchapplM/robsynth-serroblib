% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:52:56
% EndTime: 2019-03-09 21:53:22
% DurationCPUTime: 12.72s
% Computational Cost: add. (25541->630), mult. (67808->853), div. (0->0), fcn. (51565->10), ass. (0->312)
t298 = qJD(2) + qJD(3);
t297 = qJD(4) + t298;
t478 = t297 / 0.2e1;
t302 = sin(qJ(3));
t303 = sin(qJ(2));
t306 = cos(qJ(3));
t307 = cos(qJ(2));
t271 = -t302 * t303 + t306 * t307;
t261 = t271 * qJD(1);
t272 = t302 * t307 + t306 * t303;
t262 = t272 * qJD(1);
t301 = sin(qJ(4));
t305 = cos(qJ(4));
t218 = t261 * t301 + t262 * t305;
t299 = sin(pkin(11));
t341 = t305 * t261 - t262 * t301;
t379 = cos(pkin(11));
t469 = t218 * t379 + t299 * t341;
t481 = t469 / 0.2e1;
t470 = -t218 * t299 + t341 * t379;
t480 = t470 / 0.2e1;
t153 = qJD(6) - t470;
t432 = -t153 / 0.2e1;
t300 = sin(qJ(6));
t304 = cos(qJ(6));
t140 = t297 * t300 + t304 * t469;
t434 = -t140 / 0.2e1;
t139 = t297 * t304 - t300 * t469;
t435 = -t139 / 0.2e1;
t479 = Ifges(7,5) * t434 + Ifges(7,6) * t435 + Ifges(7,3) * t432;
t477 = Ifges(6,1) * t481 + Ifges(6,4) * t480 + Ifges(6,5) * t478;
t420 = -t300 / 0.2e1;
t419 = t304 / 0.2e1;
t293 = -pkin(2) * t307 - pkin(1);
t283 = qJD(1) * t293;
t237 = -t261 * pkin(3) + t283;
t169 = -pkin(4) * t341 + qJD(5) + t237;
t336 = mrSges(7,1) * t300 + mrSges(7,2) * t304;
t256 = t262 * pkin(9);
t439 = -pkin(8) - pkin(7);
t453 = t439 * t307;
t277 = qJD(1) * t453;
t263 = t302 * t277;
t284 = t439 * t303;
t276 = qJD(1) * t284;
t269 = qJD(2) * pkin(2) + t276;
t450 = t306 * t269 + t263;
t189 = t450 - t256;
t178 = pkin(3) * t298 + t189;
t266 = t306 * t277;
t230 = t269 * t302 - t266;
t414 = pkin(9) * t261;
t190 = t230 + t414;
t179 = t301 * t190;
t124 = t305 * t178 - t179;
t210 = qJ(5) * t218;
t105 = t124 - t210;
t102 = pkin(4) * t297 + t105;
t181 = t305 * t190;
t125 = t178 * t301 + t181;
t376 = qJ(5) * t341;
t106 = t125 + t376;
t370 = t299 * t106;
t54 = t102 * t379 - t370;
t52 = -t297 * pkin(5) - t54;
t321 = t52 * t336;
t476 = t169 * mrSges(6,2) + t321 + t477;
t474 = Ifges(6,4) * t481 + Ifges(6,2) * t480 + Ifges(6,6) * t478 + t479;
t421 = -t297 / 0.2e1;
t430 = -t469 / 0.2e1;
t431 = -t470 / 0.2e1;
t103 = t379 * t106;
t55 = t299 * t102 + t103;
t53 = pkin(10) * t297 + t55;
t82 = -pkin(5) * t470 - pkin(10) * t469 + t169;
t20 = -t300 * t53 + t304 * t82;
t21 = t300 * t82 + t304 * t53;
t446 = t169 * mrSges(6,1) + t20 * mrSges(7,1) - t21 * mrSges(7,2);
t475 = -Ifges(6,4) * t430 - Ifges(6,2) * t431 - Ifges(6,6) * t421 - t446 + t474 + t479;
t398 = t140 * Ifges(7,4);
t65 = t139 * Ifges(7,2) + t153 * Ifges(7,6) + t398;
t138 = Ifges(7,4) * t139;
t66 = t140 * Ifges(7,1) + t153 * Ifges(7,5) + t138;
t473 = t66 * t419 + t65 * t420;
t472 = pkin(5) * t469 - pkin(10) * t470;
t235 = t298 * t271;
t225 = t235 * qJD(1);
t236 = t298 * t272;
t226 = t236 * qJD(1);
t326 = -t225 * t305 + t226 * t301;
t122 = qJD(4) * t341 - t326;
t123 = -qJD(4) * t218 - t225 * t301 - t226 * t305;
t81 = t122 * t379 + t299 * t123;
t46 = qJD(6) * t139 + t304 * t81;
t47 = -qJD(6) * t140 - t300 * t81;
t80 = t122 * t299 - t123 * t379;
t12 = t46 * Ifges(7,4) + t47 * Ifges(7,2) + t80 * Ifges(7,6);
t13 = t46 * Ifges(7,1) + t47 * Ifges(7,4) + t80 * Ifges(7,5);
t211 = Ifges(5,4) * t341;
t150 = Ifges(5,1) * t218 + t297 * Ifges(5,5) + t211;
t366 = qJD(1) * t303;
t295 = pkin(2) * t366;
t195 = pkin(3) * t226 + qJD(2) * t295;
t100 = -pkin(4) * t123 + t195;
t17 = pkin(5) * t80 - pkin(10) * t81 + t100;
t238 = t306 * t284 + t302 * t453;
t357 = qJD(1) * qJD(2);
t167 = t450 * qJD(3) + t238 * t357;
t411 = t226 * pkin(9);
t310 = t167 - t411;
t468 = -t302 * t284 + t306 * t453;
t168 = -t230 * qJD(3) + t357 * t468;
t412 = t225 * pkin(9);
t311 = t168 - t412;
t360 = qJD(4) * t305;
t361 = qJD(4) * t301;
t49 = t178 * t360 - t190 * t361 + t301 * t311 + t305 * t310;
t30 = qJ(5) * t123 + qJD(5) * t341 + t49;
t355 = t306 * t439;
t356 = t302 * t439;
t309 = ((-t301 * t356 + t305 * t355) * t307 + (-t301 * t355 - t305 * t356) * t303) * t357;
t315 = -t122 * qJ(5) - t218 * qJD(5) - t178 * t361 - t190 * t360;
t362 = qJD(3) * t306;
t363 = qJD(3) * t302;
t9 = t379 * t30 + (-t301 * (t269 * t362 + t277 * t363 - t411) + t305 * (-t269 * t363 + t277 * t362 - t412) + t315 + t309) * t299;
t2 = qJD(6) * t20 + t17 * t300 + t304 * t9;
t330 = Ifges(7,5) * t300 + Ifges(7,6) * t304;
t331 = Ifges(7,5) * t304 - Ifges(7,6) * t300;
t402 = Ifges(7,4) * t300;
t332 = Ifges(7,2) * t304 + t402;
t401 = Ifges(7,4) * t304;
t333 = -Ifges(7,2) * t300 + t401;
t334 = Ifges(7,1) * t300 + t401;
t335 = Ifges(7,1) * t304 - t402;
t337 = mrSges(7,1) * t304 - mrSges(7,2) * t300;
t403 = Ifges(5,4) * t218;
t405 = mrSges(7,3) * t304;
t406 = mrSges(7,3) * t300;
t442 = t80 / 0.2e1;
t444 = t47 / 0.2e1;
t445 = t46 / 0.2e1;
t433 = t140 / 0.2e1;
t467 = t153 * t331 / 0.2e1 + t335 * t433 + t139 * t333 / 0.2e1 + t473;
t50 = -t125 * qJD(4) + t326 * pkin(9) + (-t230 * t305 - t301 * t450) * qJD(3) + t309;
t8 = t299 * t30 - t379 * (-t301 * t310 + t305 * t311 + t315);
t471 = -t49 * mrSges(5,2) - t9 * mrSges(6,2) + t2 * t405 + t300 * t13 / 0.2e1 + t12 * t419 + Ifges(5,6) * t123 + Ifges(5,5) * t122 + t334 * t445 + t332 * t444 + t50 * mrSges(5,1) + t330 * t442 - Ifges(6,6) * t80 + Ifges(6,5) * t81 + (-mrSges(6,1) - t337) * t8 + (t321 + t467) * qJD(6) + (Ifges(5,5) * t341 - Ifges(5,6) * t218) * t421 - (-Ifges(5,2) * t218 + t150 + t211) * t341 / 0.2e1 - t237 * (mrSges(5,1) * t218 + mrSges(5,2) * t341) - (Ifges(5,1) * t341 - t403) * t218 / 0.2e1 + (t54 * mrSges(6,3) + Ifges(6,1) * t430 + Ifges(6,4) * t431 + Ifges(6,5) * t421 + t20 * t405 + t21 * t406 + t331 * t432 + t333 * t435 + t335 * t434 - t473 - t476) * t470;
t416 = pkin(4) * t218;
t149 = Ifges(5,2) * t341 + t297 * Ifges(5,6) + t403;
t463 = t149 / 0.2e1;
t373 = qJD(6) * t21;
t3 = t17 * t304 - t300 * t9 - t373;
t329 = -t20 * t304 - t21 * t300;
t314 = m(7) * (qJD(6) * t329 + t2 * t304 - t3 * t300);
t408 = mrSges(5,3) * t341;
t131 = t305 * t189 - t179;
t109 = -t210 + t131;
t130 = -t189 * t301 - t181;
t324 = t130 - t376;
t344 = t379 * t301;
t400 = pkin(3) * qJD(4);
t459 = -t109 * t299 + t324 * t379 + (t299 * t305 + t344) * t400;
t369 = t299 * t301;
t251 = (t305 * t379 - t369) * t400;
t59 = t109 * t379 + t299 * t324;
t458 = t251 - t59;
t233 = -t276 * t302 + t266;
t191 = t233 - t414;
t234 = t306 * t276 + t263;
t192 = -t256 + t234;
t134 = t301 * t191 + t305 * t192;
t111 = -t210 + t134;
t292 = pkin(2) * t306 + pkin(3);
t368 = t301 * t302;
t223 = t292 * t360 + (-t302 * t361 + (t305 * t306 - t368) * qJD(3)) * pkin(2);
t367 = t302 * t305;
t224 = -t292 * t361 + (-t302 * t360 + (-t301 * t306 - t367) * qJD(3)) * pkin(2);
t133 = t305 * t191 - t192 * t301;
t323 = t133 - t376;
t457 = (t224 - t323) * t379 + (t111 - t223) * t299;
t163 = t223 * t379 + t299 * t224;
t61 = t111 * t379 + t299 * t323;
t456 = -t61 + t163;
t380 = -mrSges(6,1) * t297 - mrSges(7,1) * t139 + mrSges(7,2) * t140 + mrSges(6,3) * t469;
t212 = -pkin(9) * t272 + t238;
t213 = pkin(9) * t271 - t468;
t152 = t301 * t212 + t305 * t213;
t452 = t223 - t134;
t451 = t224 - t133;
t449 = -t20 * t300 + t21 * t304;
t448 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t46 + Ifges(7,6) * t47;
t438 = pkin(1) * mrSges(3,1);
t437 = pkin(1) * mrSges(3,2);
t231 = t271 * t305 - t272 * t301;
t117 = qJ(5) * t231 + t152;
t151 = t305 * t212 - t213 * t301;
t232 = t271 * t301 + t272 * t305;
t322 = -qJ(5) * t232 + t151;
t74 = t117 * t299 - t322 * t379;
t436 = t74 * t8;
t428 = t341 / 0.2e1;
t426 = t218 / 0.2e1;
t424 = t261 / 0.2e1;
t423 = -t262 / 0.2e1;
t422 = t262 / 0.2e1;
t418 = m(4) * t283;
t417 = pkin(3) * t262;
t415 = pkin(4) * t299;
t413 = t20 * mrSges(7,3);
t409 = mrSges(4,3) * t261;
t407 = mrSges(5,3) * t218;
t404 = Ifges(3,4) * t303;
t395 = t469 * t55;
t388 = t262 * mrSges(4,3);
t387 = t262 * Ifges(4,4);
t378 = Ifges(3,5) * qJD(2);
t377 = Ifges(3,6) * qJD(2);
t375 = qJD(2) * mrSges(3,1);
t374 = qJD(2) * mrSges(3,2);
t372 = t125 * t218;
t257 = -pkin(2) * t368 + t305 * t292;
t252 = pkin(4) + t257;
t258 = pkin(2) * t367 + t292 * t301;
t203 = t299 * t252 + t379 * t258;
t291 = pkin(3) * t305 + pkin(4);
t254 = pkin(3) * t344 + t299 * t291;
t365 = qJD(1) * t307;
t364 = qJD(2) * t303;
t359 = qJD(6) * t300;
t358 = qJD(6) * t304;
t349 = qJD(2) * t439;
t348 = t379 * pkin(4);
t347 = t378 / 0.2e1;
t346 = -t377 / 0.2e1;
t345 = t80 * mrSges(6,1) + t81 * mrSges(6,2);
t219 = pkin(2) * t364 + pkin(3) * t236;
t18 = mrSges(7,1) * t80 - mrSges(7,3) * t46;
t90 = -mrSges(7,2) * t153 + mrSges(7,3) * t139;
t343 = -qJD(6) * t90 - t18;
t177 = t417 + t416;
t339 = -t2 * t300 - t3 * t304;
t338 = (-t3 - t373) * mrSges(7,3);
t75 = t117 * t379 + t299 * t322;
t165 = -t231 * t379 + t232 * t299;
t166 = t299 * t231 + t232 * t379;
t245 = -t271 * pkin(3) + t293;
t182 = -t231 * pkin(4) + t245;
t88 = t165 * pkin(5) - t166 * pkin(10) + t182;
t32 = t300 * t88 + t304 * t75;
t31 = -t300 * t75 + t304 * t88;
t91 = mrSges(7,1) * t153 - mrSges(7,3) * t140;
t327 = -t300 * t91 + t304 * t90;
t136 = -qJD(4) * t232 - t235 * t301 - t236 * t305;
t112 = -pkin(4) * t136 + t219;
t202 = t252 * t379 - t299 * t258;
t278 = t303 * t349;
t279 = t307 * t349;
t173 = t306 * t278 + t302 * t279 + t284 * t362 + t363 * t453;
t146 = -pkin(9) * t236 + t173;
t174 = qJD(3) * t468 - t278 * t302 + t306 * t279;
t147 = -pkin(9) * t235 + t174;
t67 = t305 * t146 + t301 * t147 + t212 * t360 - t213 * t361;
t253 = -pkin(3) * t369 + t291 * t379;
t84 = t177 + t472;
t68 = -qJD(4) * t152 - t146 * t301 + t305 * t147;
t135 = qJD(4) * t231 + t235 * t305 - t236 * t301;
t313 = -qJ(5) * t135 - qJD(5) * t232 + t68;
t208 = t261 * Ifges(4,2) + t298 * Ifges(4,6) + t387;
t255 = Ifges(4,4) * t261;
t209 = Ifges(4,1) * t262 + t298 * Ifges(4,5) + t255;
t308 = t208 * t422 + (Ifges(4,1) * t261 - t387) * t423 + t124 * t408 + t450 * t409 - t283 * (mrSges(4,1) * t262 + mrSges(4,2) * t261) + t218 * t463 - t298 * (Ifges(4,5) * t261 - Ifges(4,6) * t262) / 0.2e1 + Ifges(4,5) * t225 - Ifges(4,6) * t226 - t167 * mrSges(4,2) + t168 * mrSges(4,1) - (-Ifges(4,2) * t262 + t209 + t255) * t261 / 0.2e1 + t475 * t469 + t471;
t294 = Ifges(3,4) * t365;
t290 = -t348 - pkin(5);
t281 = mrSges(3,3) * t365 - t374;
t280 = -mrSges(3,3) * t366 + t375;
t260 = Ifges(3,1) * t366 + t294 + t378;
t259 = t377 + (t307 * Ifges(3,2) + t404) * qJD(1);
t248 = pkin(10) + t254;
t247 = -pkin(5) - t253;
t242 = mrSges(4,1) * t298 - t388;
t241 = -mrSges(4,2) * t298 + t409;
t240 = t295 + t417;
t228 = -mrSges(4,1) * t261 + mrSges(4,2) * t262;
t198 = pkin(10) + t203;
t197 = -pkin(5) - t202;
t194 = mrSges(5,1) * t297 - t407;
t193 = -mrSges(5,2) * t297 + t408;
t170 = t177 + t295;
t161 = -mrSges(5,1) * t341 + mrSges(5,2) * t218;
t144 = -mrSges(6,2) * t297 + mrSges(6,3) * t470;
t99 = -mrSges(6,1) * t470 + mrSges(6,2) * t469;
t87 = t416 + t472;
t86 = t135 * t379 + t299 * t136;
t85 = t135 * t299 - t136 * t379;
t83 = t295 + t84;
t77 = Ifges(7,3) * t80;
t57 = t105 * t379 - t370;
t56 = t105 * t299 + t103;
t34 = qJ(5) * t136 + qJD(5) * t231 + t67;
t28 = t300 * t87 + t304 * t57;
t27 = -t300 * t57 + t304 * t87;
t26 = t300 * t83 + t304 * t61;
t25 = -t300 * t61 + t304 * t83;
t24 = t300 * t84 + t304 * t59;
t23 = -t300 * t59 + t304 * t84;
t22 = pkin(5) * t85 - pkin(10) * t86 + t112;
t19 = -mrSges(7,2) * t80 + mrSges(7,3) * t47;
t16 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t15 = t299 * t313 + t34 * t379;
t14 = t299 * t34 - t313 * t379;
t5 = -qJD(6) * t32 - t15 * t300 + t22 * t304;
t4 = qJD(6) * t31 + t15 * t304 + t22 * t300;
t1 = [(t329 * mrSges(7,3) + t467 + t476 + t477) * t86 + (Ifges(5,5) * t135 + Ifges(5,6) * t136) * t478 + m(5) * (t124 * t68 + t125 * t67 + t151 * t50 + t152 * t49 + t195 * t245 + t219 * t237) + t298 * (Ifges(4,5) * t235 - Ifges(4,6) * t236) / 0.2e1 + t283 * (mrSges(4,1) * t236 + mrSges(4,2) * t235) + (t335 * t445 + t333 * t444 + t331 * t442 + t12 * t420 + t13 * t419 + Ifges(6,1) * t81 - Ifges(6,4) * t80 + t100 * mrSges(6,2) + (mrSges(6,3) + t336) * t8 + t339 * mrSges(7,3) + (t66 * t420 + t52 * t337 + t332 * t435 + t334 * t434 + t330 * t432 - t304 * t65 / 0.2e1 - t449 * mrSges(7,3)) * qJD(6)) * t166 + (t446 - 0.2e1 * t474) * t85 + (-t9 * mrSges(6,3) + t77 / 0.2e1 - Ifges(6,4) * t81 + t100 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t80 + t448) * t165 + t136 * t463 + (t167 * t271 - t168 * t272 - t225 * t238 + t226 * t468 - t230 * t236 - t235 * t450) * mrSges(4,3) + m(4) * (-t167 * t468 + t168 * t238 + t173 * t230 + t174 * t450) + (-pkin(7) * t280 + t260 / 0.2e1 + t347 + (-0.2e1 * t437 + 0.3e1 / 0.2e1 * Ifges(3,4) * t307) * qJD(1)) * t307 * qJD(2) + (-t54 * t86 - t55 * t85 + t74 * t81 - t75 * t80) * mrSges(6,3) + t173 * t241 + t174 * t242 + t245 * (-mrSges(5,1) * t123 + mrSges(5,2) * t122) + t195 * (-mrSges(5,1) * t231 + mrSges(5,2) * t232) + t235 * t209 / 0.2e1 - t236 * t208 / 0.2e1 + t237 * (-mrSges(5,1) * t136 + mrSges(5,2) * t135) + t219 * t161 + t67 * t193 + t68 * t194 + t182 * t345 + t31 * t18 + t32 * t19 + (-t271 * t226 - t236 * t424) * Ifges(4,2) + (t271 * t225 - t226 * t272 + t235 * t424 - t236 * t422) * Ifges(4,4) + t293 * (mrSges(4,1) * t226 + mrSges(4,2) * t225) + (-t122 * t151 + t123 * t152 - t124 * t135 + t125 * t136 + t231 * t49 - t232 * t50) * mrSges(5,3) + t380 * t14 + t74 * t16 + t4 * t90 + t5 * t91 + (t225 * t272 + t235 * t422) * Ifges(4,1) + (t122 * t232 + t135 * t426) * Ifges(5,1) + (t123 * t231 + t136 * t428) * Ifges(5,2) + (t122 * t231 + t123 * t232 + t135 * t428 + t136 * t426) * Ifges(5,4) + m(7) * (t14 * t52 + t2 * t32 + t20 * t5 + t21 * t4 + t3 * t31 + t436) + m(6) * (t100 * t182 + t112 * t169 - t14 * t54 + t15 * t55 + t75 * t9 + t436) + (-pkin(7) * t281 - t259 / 0.2e1 + t346 + (-0.2e1 * t438 - 0.3e1 / 0.2e1 * t404 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t307) * qJD(1) + (t228 + qJD(1) * (-mrSges(4,1) * t271 + mrSges(4,2) * t272) + 0.2e1 * t418) * pkin(2)) * t364 + t112 * t99 + t15 * t144 + t135 * t150 / 0.2e1; t451 * t194 + t452 * t193 + t198 * t314 + t308 + t230 * t388 + t456 * t144 + ((t241 * t306 - t242 * t302) * qJD(3) + (-t225 * t306 - t226 * t302) * mrSges(4,3)) * pkin(2) - t240 * t161 - t234 * t241 - t233 * t242 + t197 * t16 + (-t163 * t91 + t198 * t343 + t338) * t300 + (-t122 * t257 + t123 * t258 + t372) * mrSges(5,3) + (-t202 * t81 - t203 * t80 + t395) * mrSges(6,3) + (t163 * t90 + t198 * t19 + (-t198 * t91 - t413) * qJD(6)) * t304 - t26 * t90 - t25 * t91 + ((t347 - t260 / 0.2e1 - t294 / 0.2e1 + qJD(1) * t437 + (t280 - t375) * pkin(7)) * t307 + (t346 + t259 / 0.2e1 + (t438 + t404 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t307) * qJD(1) + (t281 + t374) * pkin(7) + (-t228 - t418) * pkin(2)) * t303) * qJD(1) - t170 * t99 - t380 * t457 + (t163 * t449 + t197 * t8 - t20 * t25 - t21 * t26 - t457 * t52) * m(7) + (-t169 * t170 - t202 * t8 + t203 * t9 + t456 * t55 + t457 * t54) * m(6) + (t124 * t451 + t125 * t452 - t237 * t240 + t257 * t50 + t258 * t49) * m(5) + (-t450 * t233 - t230 * t234 + (t167 * t302 + t168 * t306 + (t230 * t306 - t302 * t450) * qJD(3)) * pkin(2)) * m(4); t248 * t314 - m(5) * (t124 * t130 + t125 * t131) + t308 + t458 * t144 + (t248 * t343 - t251 * t91 + t338) * t300 + (t248 * t19 + t251 * t90 + (-t248 * t91 - t413) * qJD(6)) * t304 + (-t253 * t81 - t254 * t80 + t395) * mrSges(6,3) - t450 * t241 + t247 * t16 - t131 * t193 - t130 * t194 + (t242 + t388) * t230 + mrSges(5,3) * t372 - t24 * t90 - t23 * t91 - t177 * t99 + t380 * t459 + (-t20 * t23 - t21 * t24 + t247 * t8 + t251 * t449 + t459 * t52) * m(7) + (-t169 * t177 - t253 * t8 + t254 * t9 + t458 * t55 - t459 * t54) * m(6) + (-t262 * t161 + (t193 * t305 - t194 * t301) * qJD(4) + (-t122 * t305 + t123 * t301) * mrSges(5,3) + (-t124 * t361 + t125 * t360 + 0.2e1 * t237 * t423 + t301 * t49 + t305 * t50) * m(5)) * pkin(3); -t380 * t56 + ((t299 * t9 - t379 * t8) * pkin(4) - t169 * t416 + t54 * t56 - t55 * t57) * m(6) + (-t18 * t300 + t19 * t304 - t358 * t91 - t359 * t90 + t314) * (pkin(10) + t415) + (t408 - t193) * t124 - t21 * mrSges(7,3) * t359 + t149 * t426 + (-t20 * t27 - t21 * t28 + t290 * t8 - t52 * t56) * m(7) + (-t348 * t81 - t415 * t80) * mrSges(6,3) + (t407 + t194) * t125 + (t55 * mrSges(6,3) + t475) * t469 + t290 * t16 - t3 * t406 - t358 * t413 - t99 * t416 - t28 * t90 - t27 * t91 - t57 * t144 + t471; t304 * t18 + t300 * t19 - t380 * t469 + t327 * qJD(6) + (-t144 - t327) * t470 + t345 + (t153 * t449 - t469 * t52 - t339) * m(7) + (t469 * t54 - t470 * t55 + t100) * m(6); t77 - t52 * (mrSges(7,1) * t140 + mrSges(7,2) * t139) + (Ifges(7,1) * t139 - t398) * t434 + t65 * t433 + (Ifges(7,5) * t139 - Ifges(7,6) * t140) * t432 - t20 * t90 + t21 * t91 + (t139 * t20 + t140 * t21) * mrSges(7,3) + (-Ifges(7,2) * t140 + t138 + t66) * t435 + t448;];
tauc  = t1(:);

% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:34
% EndTime: 2019-03-09 12:37:24
% DurationCPUTime: 25.32s
% Computational Cost: add. (17616->794), mult. (46931->1090), div. (0->0), fcn. (37654->10), ass. (0->342)
t296 = sin(pkin(11));
t298 = cos(pkin(11));
t301 = sin(qJ(4));
t304 = cos(qJ(4));
t267 = t296 * t304 + t298 * t301;
t297 = sin(pkin(6));
t305 = cos(qJ(2));
t385 = t297 * t305;
t316 = t267 * t385;
t223 = qJD(1) * t316;
t262 = t267 * qJD(4);
t501 = t223 - t262;
t302 = sin(qJ(2));
t336 = pkin(2) * t302 - qJ(3) * t305;
t380 = qJD(1) * t297;
t253 = t336 * t380;
t365 = t302 * t380;
t299 = cos(pkin(6));
t379 = qJD(1) * t299;
t371 = pkin(1) * t379;
t254 = -pkin(8) * t365 + t305 * t371;
t190 = t298 * t253 - t296 * t254;
t384 = t298 * t305;
t317 = (pkin(3) * t302 - pkin(9) * t384) * t297;
t163 = qJD(1) * t317 + t190;
t191 = t296 * t253 + t298 * t254;
t363 = t305 * t380;
t360 = t296 * t363;
t174 = -pkin(9) * t360 + t191;
t266 = t296 * t301 - t304 * t298;
t422 = pkin(9) + qJ(3);
t279 = t422 * t296;
t280 = t422 * t298;
t474 = -t304 * t279 - t280 * t301;
t495 = -t266 * qJD(3) + qJD(4) * t474 - t301 * t163 - t304 * t174;
t504 = -pkin(10) * t365 + t495;
t255 = pkin(8) * t363 + t302 * t371;
t216 = pkin(3) * t360 + t255;
t315 = t266 * t385;
t224 = qJD(1) * t315;
t261 = t266 * qJD(4);
t503 = -t216 + (-t224 + t261) * pkin(10) - t501 * pkin(4);
t486 = Ifges(6,1) + Ifges(7,1);
t485 = Ifges(6,5) + Ifges(7,4);
t502 = Ifges(6,6) - Ifges(7,6);
t222 = -t279 * t301 + t280 * t304;
t494 = -qJD(3) * t267 - qJD(4) * t222 - t163 * t304 + t301 * t174;
t484 = Ifges(7,5) - Ifges(6,4);
t295 = -pkin(3) * t298 - pkin(2);
t203 = pkin(4) * t266 - pkin(10) * t267 + t295;
t300 = sin(qJ(5));
t303 = cos(qJ(5));
t373 = qJD(5) * t303;
t374 = qJD(5) * t300;
t498 = t203 * t373 - t222 * t374 + t503 * t300 + t303 * t504;
t475 = t300 * t203 + t303 * t222;
t497 = -qJD(5) * t475 - t300 * t504 + t503 * t303;
t493 = pkin(4) * t365 - t494;
t289 = qJD(2) + t379;
t240 = t289 * t298 - t296 * t365;
t241 = t289 * t296 + t298 * t365;
t361 = t304 * t240 - t241 * t301;
t178 = qJD(5) - t361;
t483 = Ifges(7,2) + Ifges(6,3);
t322 = t240 * t301 + t304 * t241;
t333 = -qJD(4) + t363;
t154 = t300 * t322 + t303 * t333;
t149 = Ifges(6,4) * t154;
t155 = -t300 * t333 + t303 * t322;
t405 = Ifges(7,5) * t154;
t479 = t155 * t486 + t485 * t178 - t149 + t405;
t500 = -qJ(6) * t501 + qJD(6) * t266 + t498;
t499 = pkin(5) * t501 - t497;
t192 = -t224 * t300 - t303 * t365;
t193 = -t224 * t303 + t300 * t365;
t334 = pkin(5) * t300 - qJ(6) * t303;
t335 = pkin(5) * t303 + qJ(6) * t300;
t496 = -pkin(5) * t192 + qJ(6) * t193 - t334 * t261 + (qJD(5) * t335 - qJD(6) * t303) * t267 + t493;
t177 = Ifges(5,4) * t361;
t218 = -t289 * pkin(2) + qJD(3) - t254;
t179 = -t240 * pkin(3) + t218;
t227 = qJ(3) * t289 + t255;
t248 = (-pkin(2) * t305 - qJ(3) * t302 - pkin(1)) * t297;
t233 = qJD(1) * t248;
t165 = -t296 * t227 + t298 * t233;
t118 = -pkin(3) * t363 - t241 * pkin(9) + t165;
t166 = t298 * t227 + t296 * t233;
t135 = pkin(9) * t240 + t166;
t67 = t304 * t118 - t301 * t135;
t62 = pkin(4) * t333 - t67;
t30 = t154 * pkin(5) - t155 * qJ(6) + t62;
t68 = t301 * t118 + t304 * t135;
t63 = -pkin(10) * t333 + t68;
t83 = -pkin(4) * t361 - pkin(10) * t322 + t179;
t24 = -t300 * t63 + t303 * t83;
t25 = t300 * t83 + t303 * t63;
t329 = t24 * t303 + t25 * t300;
t476 = qJD(6) - t24;
t17 = -pkin(5) * t178 + t476;
t18 = qJ(6) * t178 + t25;
t332 = t17 * t303 - t18 * t300;
t403 = Ifges(7,5) * t303;
t338 = Ifges(7,3) * t300 + t403;
t409 = Ifges(6,4) * t303;
t344 = -Ifges(6,2) * t300 + t409;
t351 = mrSges(7,1) * t300 - mrSges(7,3) * t303;
t353 = mrSges(6,1) * t300 + mrSges(6,2) * t303;
t434 = t303 / 0.2e1;
t436 = t300 / 0.2e1;
t437 = -t300 / 0.2e1;
t448 = t178 / 0.2e1;
t451 = t155 / 0.2e1;
t453 = t154 / 0.2e1;
t454 = -t154 / 0.2e1;
t404 = Ifges(7,5) * t300;
t410 = Ifges(6,4) * t300;
t472 = t303 * t486 + t404 - t410;
t473 = -t300 * t502 + t303 * t485;
t148 = Ifges(7,5) * t155;
t55 = t178 * Ifges(7,6) + t154 * Ifges(7,3) + t148;
t411 = Ifges(6,4) * t155;
t58 = -t154 * Ifges(6,2) + t178 * Ifges(6,6) + t411;
t306 = t332 * mrSges(7,2) - t329 * mrSges(6,3) + t30 * t351 + t338 * t453 + t344 * t454 + t353 * t62 + t434 * t479 + t436 * t55 + t437 * t58 + t448 * t473 + t451 * t472;
t492 = -t306 - t179 * mrSges(5,2) - t177 / 0.2e1 + t67 * mrSges(5,3);
t491 = t300 * t485 + t303 * t502;
t490 = t300 * t486 - t403 + t409;
t489 = qJD(5) * t267;
t488 = t30 * mrSges(7,1) + t62 * mrSges(6,1) + t55 / 0.2e1 - t58 / 0.2e1 - t18 * mrSges(7,2) - t25 * mrSges(6,3);
t487 = -Ifges(3,6) * t289 / 0.2e1;
t312 = qJD(2) * t316;
t134 = qJD(1) * t312 + qJD(4) * t322;
t311 = qJD(2) * t315;
t133 = -qJD(1) * t311 + qJD(4) * t361;
t378 = qJD(2) * t297;
t362 = qJD(1) * t378;
t358 = t302 * t362;
t78 = -qJD(5) * t154 + t303 * t133 + t300 * t358;
t79 = qJD(5) * t155 + t300 * t133 - t303 * t358;
t12 = Ifges(6,5) * t78 - Ifges(6,6) * t79 + Ifges(6,3) * t134;
t13 = Ifges(7,4) * t78 + Ifges(7,2) * t134 + Ifges(7,6) * t79;
t481 = t13 + t12;
t480 = t485 * t134 + t484 * t79 + t486 * t78;
t478 = -qJD(6) * t300 + t178 * t334 - t68;
t386 = t297 * t302;
t259 = -t296 * t386 + t298 * t299;
t260 = t296 * t299 + t298 * t386;
t195 = t259 * t301 + t260 * t304;
t290 = pkin(8) * t386;
t433 = pkin(1) * t305;
t250 = t290 + (-pkin(2) - t433) * t299;
t199 = -t259 * pkin(3) + t250;
t321 = t304 * t259 - t260 * t301;
t106 = -pkin(4) * t321 - t195 * pkin(10) + t199;
t264 = t299 * t302 * pkin(1) + pkin(8) * t385;
t247 = qJ(3) * t299 + t264;
t186 = -t296 * t247 + t298 * t248;
t142 = -pkin(3) * t385 - t260 * pkin(9) + t186;
t187 = t298 * t247 + t296 * t248;
t160 = pkin(9) * t259 + t187;
t88 = t301 * t142 + t304 * t160;
t86 = -pkin(10) * t385 + t88;
t477 = t300 * t106 + t303 * t86;
t230 = (qJD(2) * t336 - qJD(3) * t302) * t297;
t212 = qJD(1) * t230;
t372 = t299 * t433;
t288 = qJD(2) * t372;
t245 = -pkin(8) * t358 + qJD(1) * t288;
t213 = qJD(3) * t289 + t245;
t157 = t298 * t212 - t296 * t213;
t314 = qJD(2) * t317;
t119 = qJD(1) * t314 + t157;
t158 = t296 * t212 + t298 * t213;
t357 = t305 * t362;
t320 = t296 * t357;
t136 = -pkin(9) * t320 + t158;
t375 = qJD(4) * t304;
t376 = qJD(4) * t301;
t22 = t118 * t375 + t301 * t119 - t135 * t376 + t304 * t136;
t20 = pkin(10) * t358 + t22;
t257 = t264 * qJD(2);
t246 = qJD(1) * t257;
t206 = pkin(3) * t320 + t246;
t54 = t134 * pkin(4) - t133 * pkin(10) + t206;
t3 = t303 * t20 + t300 * t54 + t83 * t373 - t374 * t63;
t4 = -qJD(5) * t25 - t20 * t300 + t303 * t54;
t471 = t3 * t303 - t300 * t4;
t1 = qJ(6) * t134 + qJD(6) * t178 + t3;
t2 = -pkin(5) * t134 - t4;
t470 = t1 * t303 + t2 * t300;
t284 = Ifges(3,4) * t363;
t413 = Ifges(4,4) * t298;
t345 = -Ifges(4,2) * t296 + t413;
t414 = Ifges(4,4) * t296;
t350 = Ifges(4,1) * t298 - t414;
t419 = mrSges(4,2) * t298;
t439 = t298 / 0.2e1;
t440 = -t296 / 0.2e1;
t469 = -(t165 * t298 + t166 * t296) * mrSges(4,3) + (Ifges(4,1) * t241 + Ifges(4,4) * t240 - Ifges(4,5) * t363) * t439 + (Ifges(4,4) * t241 + Ifges(4,2) * t240 - Ifges(4,6) * t363) * t440 - t254 * mrSges(3,3) + Ifges(3,1) * t365 / 0.2e1 + t284 / 0.2e1 + t289 * Ifges(3,5) + t218 * (mrSges(4,1) * t296 + t419) + t240 * t345 / 0.2e1 + t241 * t350 / 0.2e1;
t377 = qJD(2) * t302;
t364 = t297 * t377;
t256 = -pkin(8) * t364 + t288;
t239 = qJD(3) * t299 + t256;
t172 = t298 * t230 - t296 * t239;
t140 = t172 + t314;
t173 = t296 * t230 + t298 * t239;
t387 = t296 * t305;
t359 = t378 * t387;
t156 = -pkin(9) * t359 + t173;
t36 = t301 * t140 + t142 * t375 + t304 * t156 - t160 * t376;
t34 = pkin(10) * t364 + t36;
t150 = qJD(4) * t321 - t311;
t151 = qJD(4) * t195 + t312;
t217 = pkin(3) * t359 + t257;
t71 = t151 * pkin(4) - t150 * pkin(10) + t217;
t9 = -qJD(5) * t477 - t300 * t34 + t303 * t71;
t368 = -Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1;
t369 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t370 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t412 = Ifges(5,4) * t322;
t318 = Ifges(5,6) * t333;
t402 = Ifges(5,2) * t361;
t108 = -t318 + t402 + t412;
t459 = t108 / 0.2e1;
t56 = t155 * Ifges(6,5) - t154 * Ifges(6,6) + t178 * Ifges(6,3);
t57 = t155 * Ifges(7,4) + t178 * Ifges(7,2) + t154 * Ifges(7,6);
t467 = t369 * t154 - t370 * t155 + t368 * t178 + t17 * mrSges(7,1) + t25 * mrSges(6,2) + t68 * mrSges(5,3) + t459 - t56 / 0.2e1 - t57 / 0.2e1 + t412 / 0.2e1 - t179 * mrSges(5,1) - t18 * mrSges(7,3) - t24 * mrSges(6,1);
t466 = Ifges(5,2) / 0.2e1;
t465 = t78 / 0.2e1;
t464 = -t79 / 0.2e1;
t463 = t79 / 0.2e1;
t461 = pkin(1) * mrSges(3,1);
t460 = pkin(1) * mrSges(3,2);
t319 = Ifges(5,5) * t333;
t416 = Ifges(5,1) * t322;
t109 = t177 - t319 + t416;
t458 = -t109 / 0.2e1;
t457 = t109 / 0.2e1;
t456 = t134 / 0.2e1;
t452 = -t155 / 0.2e1;
t449 = -t178 / 0.2e1;
t446 = t321 / 0.2e1;
t444 = t195 / 0.2e1;
t443 = -t223 / 0.2e1;
t442 = t259 / 0.2e1;
t441 = t260 / 0.2e1;
t438 = t299 / 0.2e1;
t435 = -t303 / 0.2e1;
t430 = t22 * mrSges(5,2);
t23 = -t118 * t376 + t119 * t304 - t135 * t375 - t301 * t136;
t429 = t23 * mrSges(5,1);
t424 = -qJD(4) / 0.2e1;
t423 = qJD(4) / 0.2e1;
t46 = mrSges(6,1) * t134 - mrSges(6,3) * t78;
t47 = -t134 * mrSges(7,1) + t78 * mrSges(7,2);
t421 = -t46 + t47;
t48 = -mrSges(6,2) * t134 - mrSges(6,3) * t79;
t49 = -mrSges(7,2) * t79 + mrSges(7,3) * t134;
t420 = t48 + t49;
t418 = mrSges(6,3) * t154;
t417 = mrSges(6,3) * t155;
t415 = Ifges(3,4) * t302;
t407 = Ifges(5,5) * t150;
t399 = Ifges(5,6) * t151;
t397 = t133 * Ifges(5,1);
t396 = t133 * Ifges(5,4);
t395 = t134 * Ifges(5,4);
t392 = t296 * Ifges(4,6);
t391 = t298 * Ifges(4,5);
t162 = -mrSges(5,1) * t333 - mrSges(5,3) * t322;
t92 = mrSges(6,1) * t154 + mrSges(6,2) * t155;
t390 = t162 - t92;
t112 = pkin(4) * t322 - pkin(10) * t361;
t39 = t300 * t112 + t303 * t67;
t101 = -mrSges(7,2) * t154 + mrSges(7,3) * t178;
t102 = -mrSges(6,2) * t178 - t418;
t383 = t101 + t102;
t103 = mrSges(6,1) * t178 - t417;
t104 = -mrSges(7,1) * t178 + mrSges(7,2) * t155;
t382 = t103 - t104;
t381 = -mrSges(3,1) * t289 - mrSges(4,1) * t240 + mrSges(4,2) * t241 + mrSges(3,3) * t365;
t228 = mrSges(4,1) * t320 + t357 * t419;
t367 = t300 * t385;
t366 = Ifges(5,5) * t133 - Ifges(5,6) * t134 + Ifges(5,3) * t358;
t77 = t134 * mrSges(5,1) + t133 * mrSges(5,2);
t87 = t142 * t304 - t301 * t160;
t356 = -t1 * t300 + t2 * t303;
t355 = -t3 * t300 - t4 * t303;
t85 = pkin(4) * t385 - t87;
t354 = mrSges(6,1) * t303 - mrSges(6,2) * t300;
t352 = mrSges(7,1) * t303 + mrSges(7,3) * t300;
t343 = Ifges(6,2) * t303 + t410;
t337 = -Ifges(7,3) * t303 + t404;
t331 = t17 * t300 + t18 * t303;
t328 = t24 * t300 - t25 * t303;
t41 = t106 * t303 - t300 * t86;
t38 = t112 * t303 - t300 * t67;
t145 = t203 * t303 - t222 * t300;
t37 = t140 * t304 - t142 * t376 - t301 * t156 - t160 * t375;
t175 = t300 * t195 + t303 * t385;
t8 = t106 * t373 + t300 * t71 + t303 * t34 - t374 * t86;
t313 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t35 = -pkin(4) * t364 - t37;
t21 = -pkin(4) * t358 - t23;
t309 = t165 * mrSges(4,1) + t67 * mrSges(5,1) + Ifges(5,6) * t361 + Ifges(5,5) * t322 - Ifges(4,3) * t363 / 0.2e1 + Ifges(4,6) * t240 + Ifges(4,5) * t241 + t487 - (t305 * Ifges(3,2) + t415) * t380 / 0.2e1 - t166 * mrSges(4,2) - t255 * mrSges(3,3) - t68 * mrSges(5,2) + (-t333 / 0.2e1 + t423) * Ifges(5,3);
t282 = Ifges(3,5) * t357;
t275 = -pkin(4) - t335;
t263 = -t290 + t372;
t252 = -t289 * mrSges(3,2) + mrSges(3,3) * t363;
t235 = (mrSges(4,1) * t302 - mrSges(4,3) * t384) * t362;
t234 = (-mrSges(4,2) * t302 - mrSges(4,3) * t387) * t362;
t209 = -mrSges(4,1) * t363 - t241 * mrSges(4,3);
t208 = mrSges(4,2) * t363 + t240 * mrSges(4,3);
t201 = (Ifges(4,5) * t302 + t305 * t350) * t362;
t200 = (Ifges(4,6) * t302 + t305 * t345) * t362;
t176 = t303 * t195 - t367;
t164 = t267 * t334 - t474;
t161 = mrSges(5,2) * t333 + mrSges(5,3) * t361;
t125 = -pkin(5) * t266 - t145;
t124 = qJ(6) * t266 + t475;
t114 = -mrSges(5,2) * t358 - mrSges(5,3) * t134;
t113 = mrSges(5,1) * t358 - mrSges(5,3) * t133;
t111 = -mrSges(5,1) * t361 + mrSges(5,2) * t322;
t94 = -qJD(5) * t367 + t150 * t300 + t195 * t373 - t303 * t364;
t93 = -qJD(5) * t175 + t303 * t150 + t300 * t364;
t91 = mrSges(7,1) * t154 - mrSges(7,3) * t155;
t90 = pkin(5) * t155 + qJ(6) * t154;
t66 = Ifges(5,5) * t358 - t395 + t397;
t65 = -t134 * Ifges(5,2) + Ifges(5,6) * t358 + t396;
t43 = pkin(5) * t175 - qJ(6) * t176 + t85;
t32 = pkin(5) * t321 - t41;
t31 = -qJ(6) * t321 + t477;
t29 = mrSges(6,1) * t79 + mrSges(6,2) * t78;
t28 = mrSges(7,1) * t79 - mrSges(7,3) * t78;
t27 = -pkin(5) * t322 - t38;
t26 = qJ(6) * t322 + t39;
t14 = t78 * Ifges(6,4) - t79 * Ifges(6,2) + t134 * Ifges(6,6);
t11 = t78 * Ifges(7,5) + t134 * Ifges(7,6) + t79 * Ifges(7,3);
t10 = pkin(5) * t94 - qJ(6) * t93 - qJD(6) * t176 + t35;
t7 = pkin(5) * t79 - qJ(6) * t78 - qJD(6) * t155 + t21;
t6 = -pkin(5) * t151 - t9;
t5 = qJ(6) * t151 - qJD(6) * t321 + t8;
t15 = [(Ifges(6,4) * t93 + Ifges(6,6) * t151) * t454 - t385 * t429 + ((-Ifges(3,6) * t299 + Ifges(5,5) * t444 + Ifges(5,6) * t446 + Ifges(4,5) * t441 + Ifges(4,6) * t442 - t264 * mrSges(3,3) + (-0.2e1 * t461 - 0.3e1 / 0.2e1 * t415) * t297) * t377 + (-t407 / 0.2e1 + t399 / 0.2e1 + (Ifges(3,5) * t438 - t263 * mrSges(3,3) + (Ifges(4,1) * t260 + Ifges(4,4) * t259) * t439 + (Ifges(4,4) * t260 + Ifges(4,2) * t259) * t440 + (-0.2e1 * t460 + (0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * t391 + 0.3e1 / 0.2e1 * t392) * t305) * t297 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,3)) * t386) * qJD(2)) * t305) * t380 + m(7) * (t1 * t31 + t10 * t30 + t17 * t6 + t18 * t5 + t2 * t32 + t43 * t7) + m(5) * (t179 * t217 + t199 * t206 + t22 * t88 + t23 * t87 + t36 * t68 + t37 * t67) + m(4) * (t157 * t186 + t158 * t187 + t165 * t172 + t166 * t173 + t218 * t257 + t246 * t250) + m(3) * (t245 * t264 - t246 * t263 - t254 * t257 + t255 * t256) + t256 * t252 + t250 * t228 + t322 * (Ifges(5,1) * t150 - Ifges(5,4) * t151) / 0.2e1 + (Ifges(7,5) * t93 + Ifges(7,6) * t151) * t453 + (t57 + t56) * t151 / 0.2e1 + t361 * (Ifges(5,4) * t150 - Ifges(5,2) * t151) / 0.2e1 + m(6) * (t21 * t85 + t24 * t9 + t25 * t8 + t3 * t477 + t35 * t62 + t4 * t41) + t477 * t48 + (Ifges(6,4) * t176 - Ifges(6,6) * t321) * t464 + (-t151 * t25 + t176 * t21 + t3 * t321 + t62 * t93) * mrSges(6,2) + (Ifges(7,5) * t176 - Ifges(7,6) * t321) * t463 + (-t1 * t321 + t151 * t18 - t176 * t7 - t30 * t93) * mrSges(7,3) + (-t150 * t67 - t151 * t68 - t195 * t23 + t22 * t321) * mrSges(5,3) - t134 * (Ifges(5,4) * t195 + Ifges(5,2) * t321 - Ifges(5,6) * t385) / 0.2e1 + t133 * (Ifges(5,1) * t195 + Ifges(5,4) * t321 - Ifges(5,5) * t385) / 0.2e1 + t206 * (-mrSges(5,1) * t321 + mrSges(5,2) * t195) + t2 * (mrSges(7,1) * t321 + mrSges(7,2) * t176) + t4 * (-mrSges(6,1) * t321 - mrSges(6,3) * t176) - t481 * t321 / 0.2e1 + (t176 * t485 - t321 * t483) * t456 + (t176 * t486 - t321 * t485) * t465 + t157 * (-mrSges(4,1) * t385 - t260 * mrSges(4,3)) + t245 * (-t299 * mrSges(3,2) + mrSges(3,3) * t385) + t158 * (mrSges(4,2) * t385 + t259 * mrSges(4,3)) + (-mrSges(3,1) * t299 - mrSges(4,1) * t259 + mrSges(4,2) * t260 + mrSges(3,3) * t386) * t246 - t366 * t385 / 0.2e1 + (-Ifges(6,2) * t454 + Ifges(7,3) * t453 - t448 * t502 + t484 * t451 + t488) * t94 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) + t7 * mrSges(7,1) + t21 * mrSges(6,1) + t11 / 0.2e1 - t14 / 0.2e1 + Ifges(7,3) * t463 - Ifges(6,2) * t464 + t484 * t465 - t502 * t456) * t175 + t187 * t234 + t186 * t235 + t173 * t208 + t172 * t209 + t217 * t111 + t199 * t77 + t179 * (mrSges(5,1) * t151 + mrSges(5,2) * t150) + t36 * t161 + t37 * t162 - t151 * t108 / 0.2e1 + t24 * (mrSges(6,1) * t151 - mrSges(6,3) * t93) + t17 * (-mrSges(7,1) * t151 + mrSges(7,2) * t93) + t87 * t113 + t88 * t114 + t8 * t102 + t9 * t103 + t6 * t104 + t5 * t101 + t35 * t92 + t10 * t91 + t85 * t29 + t31 * t49 + t43 * t28 + t41 * t46 + t32 * t47 + t381 * t257 + (-t399 + t407) * t423 + t385 * t430 + t282 * t438 + t201 * t441 + t200 * t442 + t66 * t444 + t65 * t446 + t150 * t457 + t479 * t93 / 0.2e1 + t480 * t176 / 0.2e1 + (t151 * t483 + t485 * t93) * t448 + (t151 * t485 + t486 * t93) * t451 + (t469 * t305 + (t309 + t487) * t302) * t378; (Ifges(5,6) * t424 - t402 / 0.2e1 - t467) * t262 + t282 + (-t157 * mrSges(4,3) - qJD(3) * t209 - qJ(3) * t235 + t201 / 0.2e1 + t246 * mrSges(4,2)) * t296 + (t12 / 0.2e1 + t13 / 0.2e1 - t65 / 0.2e1 + t206 * mrSges(5,1) - t396 / 0.2e1 - t22 * mrSges(5,3) - t369 * t79 + t370 * t78 + (t466 - t368) * t134 + t313) * t266 + (-t165 * t190 - t166 * t191 - t218 * t255 - pkin(2) * t246 + (-t165 * t296 + t166 * t298) * qJD(3) + (-t157 * t296 + t158 * t298) * qJ(3)) * m(4) + (-t193 * t62 + t223 * t25) * mrSges(6,2) + (t158 * mrSges(4,3) + qJD(3) * t208 + qJ(3) * t234 + t200 / 0.2e1 - t246 * mrSges(4,1)) * t298 - t254 * t252 - t245 * mrSges(3,2) - t246 * mrSges(3,1) - (Ifges(5,5) * t423 + t457 + t416 / 0.2e1 - t492) * t261 + t475 * t48 + (t223 * t68 - t224 * t67) * mrSges(5,3) - t322 * (-Ifges(5,1) * t224 - Ifges(5,4) * t223) / 0.2e1 + (t11 * t436 + t338 * t463 + t344 * t464 + t21 * t353 + t7 * t351 + t206 * mrSges(5,2) + t397 / 0.2e1 - t23 * mrSges(5,3) - t395 / 0.2e1 + t66 / 0.2e1 + t355 * mrSges(6,3) + t356 * mrSges(7,2) + (-mrSges(7,2) * t331 + mrSges(6,3) * t328 + t30 * t352 + t337 * t454 + t343 * t453 + t354 * t62 + t435 * t58) * qJD(5) + t472 * t465 + t473 * t456 + (qJD(5) * t479 + t14) * t437 + (qJD(5) * t55 + t480) * t434) * t267 + t493 * t92 + t494 * t162 + t495 * t161 + (-Ifges(6,2) * t453 + Ifges(7,3) * t454 - t449 * t502 + t452 * t484 - t488) * t192 - pkin(2) * t228 + t222 * t114 - t24 * (mrSges(6,1) * t223 - mrSges(6,3) * t193) - t17 * (-mrSges(7,1) * t223 + mrSges(7,2) * t193) + (-t179 * t216 + t206 * t295 + t22 * t222 + t23 * t474 + t494 * t67 + t495 * t68) * m(5) - (t29 - t113) * t474 + (t145 * t4 - t21 * t474 + t24 * t497 + t25 * t498 + t475 * t3 + t493 * t62) * m(6) + (-t18 * t223 + t193 * t30) * mrSges(7,3) - t191 * t208 - t190 * t209 - t216 * t111 + t164 * t28 + t145 * t46 + t124 * t49 + t125 * t47 - t361 * (-Ifges(5,4) * t224 - Ifges(5,2) * t223) / 0.2e1 - t179 * (mrSges(5,1) * t223 - mrSges(5,2) * t224) + ((-t284 / 0.2e1 + ((Ifges(4,1) * t296 + t413) * t439 + (Ifges(4,2) * t298 + t414) * t440) * qJD(2) + (t460 + (t391 / 0.2e1 - t392 / 0.2e1) * t305) * t380 + (t262 / 0.2e1 + t443) * Ifges(5,6) + (t261 / 0.2e1 - t224 / 0.2e1) * Ifges(5,5) - t469) * t305 + ((t289 / 0.2e1 - qJD(2)) * Ifges(3,6) + ((t461 + t415 / 0.2e1) * t297 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t385) * qJD(1) - t309 + (Ifges(4,5) * t296 + Ifges(5,5) * t267 + Ifges(4,6) * t298 - Ifges(5,6) * t266) * qJD(2) / 0.2e1) * t302) * t380 + (-Ifges(5,5) * t224 - Ifges(5,6) * t223) * t424 - t381 * t255 + t295 * t77 + t496 * t91 + t497 * t103 + (Ifges(6,4) * t193 + Ifges(6,6) * t223) * t453 + (Ifges(7,5) * t193 + Ifges(7,6) * t223) * t454 + t498 * t102 + t499 * t104 + (t1 * t124 + t125 * t2 + t164 * t7 + t17 * t499 + t18 * t500 + t30 * t496) * m(7) + t500 * t101 - t479 * t193 / 0.2e1 + (t193 * t485 + t223 * t483 + t489 * t491) * t449 + (t193 * t486 + t485 * t223 + t490 * t489) * t452 + t57 * t443 + t56 * t443 - t224 * t458 + t223 * t459; -t361 * t161 - t240 * t208 + t241 * t209 + (-t91 + t390) * t322 + (t178 * t383 - t421) * t303 + (-t178 * t382 + t420) * t300 + t77 + t228 + (t178 * t331 - t322 * t30 - t356) * m(7) + (-t178 * t328 - t322 * t62 - t355) * m(6) + (t322 * t67 - t361 * t68 + t206) * m(5) + (t165 * t241 - t166 * t240 + t246) * m(4); (-t318 / 0.2e1 + t467) * t322 + t366 + t490 * t465 + t491 * t456 + (t458 + (t466 - Ifges(5,1) / 0.2e1) * t322 + t319 / 0.2e1 + t492) * t361 - t7 * t352 - t21 * t354 + (-pkin(4) * t21 - t24 * t38 - t25 * t39 - t62 * t68) * m(6) - t67 * t161 - t39 * t102 - t38 * t103 - t27 * t104 - t26 * t101 + t390 * t68 - pkin(4) * t29 + t275 * t28 + t429 - t430 + t470 * mrSges(7,2) + t471 * mrSges(6,3) + (m(6) * t471 + t421 * t300 + t420 * t303 + m(7) * t470 + (-m(6) * t329 + m(7) * t332 - t300 * t383 - t303 * t382) * qJD(5)) * pkin(10) + t478 * t91 + (-t17 * t27 - t18 * t26 + t275 * t7 + t30 * t478) * m(7) + t14 * t434 + t11 * t435 + t337 * t463 + t343 * t464 + t480 * t436 + t306 * qJD(5); t313 + (t382 + t417) * t25 + (-t383 - t418) * t24 - t30 * (mrSges(7,1) * t155 + mrSges(7,3) * t154) - t62 * (mrSges(6,1) * t155 - mrSges(6,2) * t154) + t58 * t451 + (Ifges(7,3) * t155 - t405) * t454 + qJD(6) * t101 - t90 * t91 + qJ(6) * t49 - pkin(5) * t47 + (t154 * t17 + t155 * t18) * mrSges(7,2) + (-t154 * t485 - t155 * t502) * t449 + (-pkin(5) * t2 + qJ(6) * t1 - t17 * t25 + t18 * t476 - t30 * t90) * m(7) + (-Ifges(6,2) * t155 - t149 + t479) * t453 + (-t154 * t486 + t148 - t411 + t55) * t452 + t481; -t178 * t101 + t155 * t91 + 0.2e1 * (t2 / 0.2e1 + t30 * t451 + t18 * t449) * m(7) + t47;];
tauc  = t15(:);

% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRP5
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:46
% EndTime: 2019-03-09 16:48:33
% DurationCPUTime: 23.74s
% Computational Cost: add. (14203->741), mult. (35648->991), div. (0->0), fcn. (25390->8), ass. (0->321)
t304 = sin(qJ(2));
t306 = cos(qJ(2));
t333 = pkin(2) * t304 - pkin(8) * t306;
t264 = t333 * qJD(1);
t305 = cos(qJ(3));
t303 = sin(qJ(3));
t367 = qJD(1) * t304;
t346 = t303 * t367;
t210 = pkin(7) * t346 + t305 * t264;
t376 = t305 * t306;
t322 = pkin(3) * t304 - qJ(4) * t376;
t402 = -qJ(4) - pkin(8);
t336 = qJD(3) * t402;
t495 = -qJD(1) * t322 - qJD(4) * t303 + t305 * t336 - t210;
t246 = t303 * t264;
t359 = qJD(4) * t305;
t377 = t304 * t305;
t378 = t303 * t306;
t494 = t246 + (-pkin(7) * t377 - qJ(4) * t378) * qJD(1) - t303 * t336 - t359;
t300 = sin(pkin(10));
t301 = cos(pkin(10));
t256 = t300 * t305 + t301 * t303;
t312 = t256 * t306;
t218 = qJD(1) * t312;
t242 = t256 * qJD(3);
t493 = t218 - t242;
t323 = t300 * t303 - t301 * t305;
t311 = t323 * t306;
t219 = qJD(1) * t311;
t243 = t323 * qJD(3);
t370 = -t219 + t243;
t477 = t494 * t300 + t301 * t495;
t476 = t300 * t495 - t494 * t301;
t492 = pkin(9) * t493 + t476;
t491 = -pkin(4) * t367 + pkin(9) * t370 + t477;
t364 = qJD(2) * t305;
t262 = -t346 + t364;
t345 = t305 * t367;
t263 = qJD(2) * t303 + t345;
t192 = t262 * t301 - t263 * t300;
t193 = t262 * t300 + t263 * t301;
t302 = sin(qJ(5));
t411 = cos(qJ(5));
t124 = -t411 * t192 + t302 * t193;
t357 = qJD(2) * qJD(3);
t361 = qJD(3) * t304;
t363 = qJD(2) * t306;
t208 = t305 * t357 + (-t303 * t361 + t305 * t363) * qJD(1);
t360 = qJD(3) * t305;
t474 = t303 * t363 + t304 * t360;
t209 = -qJD(1) * t474 - t303 * t357;
t147 = -t208 * t300 + t209 * t301;
t148 = t208 * t301 + t209 * t300;
t55 = -qJD(5) * t124 + t302 * t147 + t148 * t411;
t444 = t55 / 0.2e1;
t471 = t192 * t302 + t193 * t411;
t56 = qJD(5) * t471 - t147 * t411 + t302 * t148;
t442 = t56 / 0.2e1;
t465 = Ifges(6,1) + Ifges(7,1);
t464 = Ifges(7,4) + Ifges(6,5);
t463 = Ifges(7,5) - Ifges(6,4);
t314 = -t302 * t256 - t323 * t411;
t120 = qJD(5) * t314 - t302 * t242 - t243 * t411;
t153 = -t302 * t218 - t219 * t411;
t373 = t120 - t153;
t187 = t256 * t411 - t302 * t323;
t121 = qJD(5) * t187 + t242 * t411 - t302 * t243;
t152 = t218 * t411 - t219 * t302;
t372 = t121 - t152;
t366 = qJD(1) * t306;
t296 = pkin(7) * t366;
t408 = pkin(3) * t303;
t253 = t366 * t408 + t296;
t362 = qJD(3) * t303;
t490 = pkin(3) * t362 - t253;
t70 = pkin(5) * t471 + qJ(6) * t124;
t342 = Ifges(3,5) * qJD(2) / 0.2e1;
t365 = qJD(2) * t304;
t337 = qJD(1) * t365;
t489 = t463 * t442 + t465 * t444 + t464 * t337 / 0.2e1;
t462 = -Ifges(6,6) + Ifges(7,6);
t487 = Ifges(5,3) + Ifges(4,3);
t486 = Ifges(6,3) + Ifges(7,2);
t272 = t402 * t303;
t273 = t402 * t305;
t202 = t301 * t272 + t273 * t300;
t174 = -pkin(9) * t256 + t202;
t203 = t300 * t272 - t301 * t273;
t175 = -pkin(9) * t323 + t203;
t316 = t174 * t411 - t302 * t175;
t481 = qJD(5) * t316 + t491 * t302 + t411 * t492;
t104 = t302 * t174 + t175 * t411;
t480 = -qJD(5) * t104 - t302 * t492 + t491 * t411;
t475 = -pkin(4) * t493 + t490;
t485 = pkin(9) * t193;
t483 = -qJ(6) * t367 + t481;
t482 = pkin(5) * t367 - t480;
t479 = pkin(5) * t372 - qJ(6) * t373 - qJD(6) * t187 + t475;
t119 = Ifges(6,4) * t124;
t289 = qJD(3) - t366;
t277 = qJD(5) + t289;
t392 = Ifges(7,5) * t124;
t461 = t277 * t464 + t465 * t471 - t119 + t392;
t478 = Ifges(5,4) * t193;
t292 = pkin(3) * t301 + pkin(4);
t406 = pkin(9) * t192;
t268 = -pkin(2) * t306 - t304 * pkin(8) - pkin(1);
t250 = t268 * qJD(1);
t276 = qJD(2) * pkin(8) + t296;
t198 = t305 * t250 - t276 * t303;
t160 = -qJ(4) * t263 + t198;
t199 = t250 * t303 + t276 * t305;
t161 = qJ(4) * t262 + t199;
t380 = t301 * t161;
t95 = -t160 * t300 - t380;
t321 = t95 - t406;
t338 = qJD(5) * t411;
t409 = pkin(3) * t300;
t356 = t302 * t409;
t154 = t300 * t161;
t96 = t301 * t160 - t154;
t79 = t96 - t485;
t457 = -qJD(5) * t356 + t292 * t338 - t302 * t321 - t411 * t79;
t294 = Ifges(3,4) * t366;
t387 = t263 * Ifges(4,4);
t178 = t262 * Ifges(4,2) + t289 * Ifges(4,6) + t387;
t254 = Ifges(4,4) * t262;
t179 = t263 * Ifges(4,1) + t289 * Ifges(4,5) + t254;
t275 = -qJD(2) * pkin(2) + pkin(7) * t367;
t324 = t198 * t305 + t199 * t303;
t397 = Ifges(4,4) * t305;
t328 = -Ifges(4,2) * t303 + t397;
t398 = Ifges(4,4) * t303;
t330 = Ifges(4,1) * t305 - t398;
t331 = mrSges(4,1) * t303 + mrSges(4,2) * t305;
t391 = Ifges(4,6) * t303;
t394 = Ifges(4,5) * t305;
t413 = t305 / 0.2e1;
t414 = -t303 / 0.2e1;
t415 = t289 / 0.2e1;
t419 = t263 / 0.2e1;
t307 = -t324 * mrSges(4,3) + t275 * t331 + t262 * t328 / 0.2e1 + t330 * t419 + (-t391 + t394) * t415 + t178 * t414 + t179 * t413;
t473 = t307 + Ifges(3,1) * t367 / 0.2e1 + t294 / 0.2e1 + t342;
t206 = -t262 * pkin(3) + qJD(4) + t275;
t137 = -pkin(4) * t192 + t206;
t149 = pkin(3) * t289 + t160;
t88 = t301 * t149 - t154;
t73 = pkin(4) * t289 - t485 + t88;
t89 = t300 * t149 + t380;
t74 = t89 + t406;
t21 = -t302 * t74 + t411 * t73;
t456 = qJD(6) - t21;
t19 = -t277 * pkin(5) + t456;
t38 = t124 * pkin(5) - qJ(6) * t471 + t137;
t472 = -mrSges(6,2) * t137 - mrSges(7,2) * t19 + mrSges(6,3) * t21 + t38 * mrSges(7,3) - t461 / 0.2e1;
t470 = -Ifges(6,6) / 0.2e1;
t469 = Ifges(7,6) / 0.2e1;
t113 = t192 * Ifges(5,2) + Ifges(5,6) * t289 + t478;
t468 = t113 / 0.2e1;
t436 = -t124 / 0.2e1;
t435 = t124 / 0.2e1;
t432 = t471 / 0.2e1;
t467 = -t192 / 0.2e1;
t466 = t192 / 0.2e1;
t417 = t277 / 0.2e1;
t341 = -Ifges(3,6) * qJD(2) / 0.2e1;
t258 = t305 * t268;
t407 = pkin(7) * t303;
t195 = -qJ(4) * t377 + t258 + (-pkin(3) - t407) * t306;
t291 = pkin(7) * t376;
t221 = t303 * t268 + t291;
t379 = t303 * t304;
t201 = -qJ(4) * t379 + t221;
t130 = t300 * t195 + t301 * t201;
t231 = t256 * t304;
t105 = -pkin(9) * t231 + t130;
t129 = t301 * t195 - t300 * t201;
t232 = t323 * t304;
t99 = -pkin(4) * t306 + t232 * pkin(9) + t129;
t460 = t411 * t105 + t302 * t99;
t238 = t302 * t292 + t411 * t409;
t459 = -t238 * qJD(5) + t302 * t79 - t321 * t411;
t458 = qJD(6) + t457;
t400 = mrSges(6,3) * t471;
t108 = mrSges(6,1) * t277 - t400;
t109 = -mrSges(7,1) * t277 + mrSges(7,2) * t471;
t374 = t108 - t109;
t455 = Ifges(4,5) * t208 + Ifges(4,6) * t209;
t454 = Ifges(5,5) * t148 + Ifges(5,6) * t147 + t337 * t487 + t455;
t453 = t337 * t486 + t462 * t56 + t464 * t55;
t452 = qJD(1) * pkin(1) * mrSges(3,2);
t265 = t333 * qJD(2);
t251 = qJD(1) * t265;
t335 = pkin(7) * t337;
t136 = -qJD(3) * t199 + t305 * t251 + t303 * t335;
t87 = pkin(3) * t337 - qJ(4) * t208 - qJD(4) * t263 + t136;
t135 = t250 * t360 + t303 * t251 - t276 * t362 - t305 * t335;
t93 = qJ(4) * t209 + qJD(4) * t262 + t135;
t31 = -t300 * t93 + t301 * t87;
t24 = pkin(4) * t337 - pkin(9) * t148 + t31;
t32 = t300 * t87 + t301 * t93;
t26 = pkin(9) * t147 + t32;
t358 = qJD(5) * t302;
t5 = t302 * t24 + t411 * t26 + t73 * t338 - t358 * t74;
t2 = qJ(6) * t337 + qJD(6) * t277 + t5;
t22 = t302 * t73 + t411 * t74;
t6 = -qJD(5) * t22 + t24 * t411 - t302 * t26;
t3 = -pkin(5) * t337 - t6;
t451 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t169 = -qJD(2) * t311 - t242 * t304;
t368 = t305 * t265 + t365 * t407;
t117 = -t304 * t359 + t322 * qJD(2) + (-t291 + (qJ(4) * t304 - t268) * t303) * qJD(3) + t368;
t369 = t303 * t265 + t268 * t360;
t127 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t377 + (-qJD(4) * t304 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t306) * t303 + t369;
t68 = t301 * t117 - t127 * t300;
t44 = pkin(4) * t365 - pkin(9) * t169 + t68;
t168 = -qJD(2) * t312 + t323 * t361;
t69 = t300 * t117 + t301 * t127;
t47 = pkin(9) * t168 + t69;
t10 = -qJD(5) * t460 - t302 * t47 + t411 * t44;
t450 = -t136 * mrSges(4,1) - t31 * mrSges(5,1) + t135 * mrSges(4,2) + t32 * mrSges(5,2);
t20 = t277 * qJ(6) + t22;
t350 = Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t351 = t469 + t470;
t352 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t353 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t399 = Ifges(3,4) * t304;
t448 = t351 * t124 + t353 * t471 + t352 * t277 + t350 * t289 + t198 * mrSges(4,1) + t20 * mrSges(7,3) + t21 * mrSges(6,1) + t88 * mrSges(5,1) + Ifges(5,6) * t192 + Ifges(5,5) * t193 + t263 * Ifges(4,5) + t262 * Ifges(4,6) + t341 - (t306 * Ifges(3,2) + t399) * qJD(1) / 0.2e1 + Ifges(6,6) * t436 + Ifges(7,6) * t435 - t19 * mrSges(7,1) - t199 * mrSges(4,2) - t22 * mrSges(6,2) - t89 * mrSges(5,2) + t464 * t432 + t486 * t417 + t487 * t415;
t118 = Ifges(7,5) * t471;
t61 = t277 * Ifges(7,6) + t124 * Ifges(7,3) + t118;
t395 = Ifges(6,4) * t471;
t64 = -t124 * Ifges(6,2) + t277 * Ifges(6,6) + t395;
t447 = t137 * mrSges(6,1) + t38 * mrSges(7,1) + t61 / 0.2e1 - t64 / 0.2e1 - t20 * mrSges(7,2) - t22 * mrSges(6,3);
t446 = Ifges(7,5) * t444 + Ifges(7,3) * t442 + t337 * t469;
t445 = -t55 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t442 + t337 * t470;
t443 = -t56 / 0.2e1;
t440 = pkin(1) * mrSges(3,1);
t433 = -t471 / 0.2e1;
t431 = t147 / 0.2e1;
t430 = t148 / 0.2e1;
t427 = -t193 / 0.2e1;
t426 = t193 / 0.2e1;
t425 = t208 / 0.2e1;
t424 = t209 / 0.2e1;
t423 = -t231 / 0.2e1;
t422 = -t232 / 0.2e1;
t421 = -t262 / 0.2e1;
t420 = -t263 / 0.2e1;
t418 = -t277 / 0.2e1;
t416 = -t289 / 0.2e1;
t410 = pkin(3) * t263;
t401 = mrSges(6,3) * t124;
t396 = Ifges(5,4) * t192;
t382 = qJD(2) * mrSges(3,2);
t106 = -mrSges(7,2) * t124 + mrSges(7,3) * t277;
t107 = -mrSges(6,2) * t277 - t401;
t375 = -t106 - t107;
t266 = pkin(3) * t379 + t304 * pkin(7);
t215 = pkin(3) * t474 + pkin(7) * t363;
t293 = -pkin(3) * t305 - pkin(2);
t17 = t56 * mrSges(6,1) + t55 * mrSges(6,2);
t16 = t56 * mrSges(7,1) - t55 * mrSges(7,3);
t190 = -pkin(3) * t209 + qJD(2) * t296;
t83 = -t147 * mrSges(5,1) + t148 * mrSges(5,2);
t334 = m(4) * t275 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t262 + mrSges(4,2) * t263 + mrSges(3,3) * t367;
t196 = pkin(4) * t231 + t266;
t151 = pkin(4) * t193 + t410;
t332 = mrSges(4,1) * t305 - mrSges(4,2) * t303;
t329 = Ifges(4,1) * t303 + t397;
t327 = Ifges(4,2) * t305 + t398;
t326 = Ifges(4,5) * t303 + Ifges(4,6) * t305;
t325 = t135 * t305 - t136 * t303;
t222 = pkin(4) * t323 + t293;
t131 = -pkin(4) * t168 + t215;
t42 = -mrSges(7,1) * t337 + t55 * mrSges(7,2);
t101 = -pkin(4) * t147 + t190;
t59 = -t302 * t105 + t411 * t99;
t315 = -t231 * t411 + t302 * t232;
t167 = -t302 * t231 - t232 * t411;
t9 = -t105 * t358 + t302 * t44 + t99 * t338 + t411 * t47;
t237 = t292 * t411 - t356;
t310 = t451 + t453;
t270 = mrSges(3,3) * t366 - t382;
t233 = -pkin(5) - t237;
t230 = qJ(6) + t238;
t220 = -pkin(7) * t378 + t258;
t214 = mrSges(4,1) * t289 - mrSges(4,3) * t263;
t213 = -mrSges(4,2) * t289 + mrSges(4,3) * t262;
t211 = -pkin(7) * t345 + t246;
t189 = -mrSges(4,2) * t337 + mrSges(4,3) * t209;
t188 = mrSges(4,1) * t337 - mrSges(4,3) * t208;
t165 = mrSges(5,1) * t289 - mrSges(5,3) * t193;
t164 = -t289 * mrSges(5,2) + mrSges(5,3) * t192;
t159 = -qJD(3) * t221 + t368;
t158 = (-t304 * t364 - t306 * t362) * pkin(7) + t369;
t150 = -mrSges(4,1) * t209 + mrSges(4,2) * t208;
t139 = t208 * Ifges(4,1) + t209 * Ifges(4,4) + Ifges(4,5) * t337;
t138 = t208 * Ifges(4,4) + t209 * Ifges(4,2) + Ifges(4,6) * t337;
t133 = mrSges(5,1) * t337 - mrSges(5,3) * t148;
t132 = -mrSges(5,2) * t337 + mrSges(5,3) * t147;
t128 = -mrSges(5,1) * t192 + t193 * mrSges(5,2);
t114 = Ifges(5,1) * t193 + Ifges(5,5) * t289 + t396;
t100 = -pkin(5) * t314 - qJ(6) * t187 + t222;
t82 = -pkin(5) * t315 - qJ(6) * t167 + t196;
t81 = t148 * Ifges(5,1) + t147 * Ifges(5,4) + Ifges(5,5) * t337;
t80 = t148 * Ifges(5,4) + t147 * Ifges(5,2) + Ifges(5,6) * t337;
t78 = qJD(5) * t167 - t168 * t411 + t302 * t169;
t77 = qJD(5) * t315 + t302 * t168 + t169 * t411;
t72 = mrSges(6,1) * t124 + mrSges(6,2) * t471;
t71 = mrSges(7,1) * t124 - mrSges(7,3) * t471;
t58 = t306 * pkin(5) - t59;
t57 = -qJ(6) * t306 + t460;
t46 = t151 + t70;
t43 = -mrSges(6,2) * t337 - mrSges(6,3) * t56;
t41 = mrSges(6,1) * t337 - mrSges(6,3) * t55;
t40 = -mrSges(7,2) * t56 + mrSges(7,3) * t337;
t18 = pkin(5) * t78 - qJ(6) * t77 - qJD(6) * t167 + t131;
t11 = pkin(5) * t56 - qJ(6) * t55 - qJD(6) * t471 + t101;
t8 = -pkin(5) * t365 - t10;
t7 = qJ(6) * t365 - qJD(6) * t306 + t9;
t1 = [(t334 * pkin(7) + t342 - 0.2e1 * t452 + t473) * t363 + (Ifges(5,4) * t169 + Ifges(5,2) * t168) * t466 + t168 * t468 + (-Ifges(6,2) * t436 + Ifges(7,3) * t435 + t417 * t462 + t432 * t463 + t447) * t78 + (-Ifges(5,1) * t232 - Ifges(5,4) * t231) * t430 + (-Ifges(5,4) * t232 - Ifges(5,2) * t231) * t431 + (t168 * t89 - t169 * t88 - t231 * t32 + t232 * t31) * mrSges(5,3) + t190 * (mrSges(5,1) * t231 - mrSges(5,2) * t232) + t167 * t489 + m(6) * (t10 * t21 + t101 * t196 + t131 * t137 + t22 * t9 + t460 * t5 + t59 * t6) + t460 * t43 + (t167 * t465 - t315 * t463) * t444 + (t139 * t413 + t138 * t414 + pkin(7) * t150 + t330 * t425 + t328 * t424 + (-t135 * t303 - t136 * t305) * mrSges(4,3) + (t327 * t421 + t329 * t420 + t275 * t332 + t326 * t416 + t179 * t414 - t305 * t178 / 0.2e1 + (t198 * t303 - t199 * t305) * mrSges(4,3)) * qJD(3) + (-pkin(7) * t270 + ((-0.3e1 / 0.2e1 * Ifges(3,4) + t394 / 0.2e1 - t391 / 0.2e1) * t304 + Ifges(5,5) * t422 + Ifges(5,6) * t423 - 0.2e1 * t440 + t353 * t167 - t351 * t315) * qJD(1) + t341 + t448) * qJD(2)) * t304 - t315 * t445 - t315 * t446 + (t167 * t3 + t2 * t315) * mrSges(7,2) + (-t167 * t6 + t315 * t5) * mrSges(6,3) + (Ifges(6,4) * t167 + Ifges(6,2) * t315) * t443 + (Ifges(7,5) * t167 - Ifges(7,3) * t315) * t442 + t11 * (-mrSges(7,1) * t315 - mrSges(7,3) * t167) + t101 * (-mrSges(6,1) * t315 + mrSges(6,2) * t167) + m(5) * (t129 * t31 + t130 * t32 + t190 * t266 + t206 * t215 + t68 * t88 + t69 * t89) + m(7) * (t11 * t82 + t18 * t38 + t19 * t8 + t2 * t57 + t20 * t7 + t3 * t58) + (-Ifges(5,5) * t430 - Ifges(5,6) * t431 - Ifges(6,6) * t443 - Ifges(7,6) * t442 - t464 * t444 + (0.3e1 / 0.2e1 * Ifges(3,4) * t363 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + (m(4) * pkin(7) + t331) * pkin(7) - t350 - t352) * t365) * qJD(1) + t450 - t451) * t306 - (t453 + t454 + t455) * t306 / 0.2e1 + (Ifges(5,1) * t169 + Ifges(5,4) * t168) * t426 + t81 * t422 + t80 * t423 + (Ifges(5,5) * t169 + Ifges(5,6) * t168) * t415 + (Ifges(6,4) * t436 + Ifges(7,5) * t435 + t464 * t417 + t465 * t432 - t472) * t77 + m(4) * (t135 * t221 + t136 * t220 + t199 * t158 + t198 * t159) + t57 * t40 + t58 * t42 + t59 * t41 + t18 * t71 + t82 * t16 + t7 * t106 + t9 * t107 + t10 * t108 + t8 * t109 + t131 * t72 + t130 * t132 + t129 * t133 + t69 * t164 + t68 * t165 + t169 * t114 / 0.2e1 + t196 * t17 + t206 * (-mrSges(5,1) * t168 + mrSges(5,2) * t169) + t158 * t213 + t159 * t214 + t215 * t128 + t220 * t188 + t221 * t189 + t266 * t83; (-mrSges(5,1) * t493 - mrSges(5,2) * t370) * t206 + (-t256 * t31 - t32 * t323 + t370 * t88 + t493 * t89) * mrSges(5,3) + (Ifges(5,4) * t256 - Ifges(5,2) * t323) * t431 + t190 * (mrSges(5,1) * t323 + mrSges(5,2) * t256) - t323 * t80 / 0.2e1 + (t218 / 0.2e1 - t242 / 0.2e1) * t113 + (-Ifges(5,4) * t219 - Ifges(5,2) * t218) * t467 + (t152 * t462 + t153 * t464) * t418 + (-Ifges(5,1) * t219 - Ifges(5,4) * t218) * t427 + (t120 * t464 + t121 * t462) * t417 + (-Ifges(5,5) * t219 - Ifges(5,6) * t218) * t416 + (t219 / 0.2e1 - t243 / 0.2e1) * t114 + (t152 * t463 + t153 * t465) * t433 + (t120 * t465 + t121 * t463) * t432 + (t187 * t465 - t314 * t463) * t444 + (Ifges(6,4) * t120 + Ifges(7,5) * t153 - Ifges(6,2) * t121 + Ifges(7,3) * t152) * t436 + ((t452 - t294 / 0.2e1 + t342 + ((-m(4) * pkin(2) - mrSges(3,1) - t332) * qJD(2) - t334) * pkin(7) - t473) * t306 + ((t270 + t382) * pkin(7) + (t440 + t399 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t306) * qJD(1) + t341 - t448) * t304 + (Ifges(5,5) * t256 - Ifges(5,6) * t323 + t187 * t464 - t314 * t462 + t326) * t365 / 0.2e1) * qJD(1) + (-Ifges(5,4) * t243 - Ifges(5,2) * t242) * t466 + t461 * (t120 / 0.2e1 - t153 / 0.2e1) + t479 * t71 + t480 * t108 + (-t188 * t303 + t189 * t305 + (-t213 * t303 - t214 * t305) * qJD(3) + m(4) * (-qJD(3) * t324 + t325)) * pkin(8) + (mrSges(6,1) * t372 + mrSges(6,2) * t373) * t137 + (mrSges(7,1) * t372 - mrSges(7,3) * t373) * t38 + (-Ifges(5,1) * t243 - Ifges(5,4) * t242) * t426 + (-Ifges(5,5) * t243 - Ifges(5,6) * t242) * t415 + t187 * t489 + (t101 * t222 + t104 * t5 + t137 * t475 + t21 * t480 + t22 * t481 + t316 * t6) * m(6) + t481 * t107 + t482 * t109 + (t61 - t64) * (t121 / 0.2e1 - t152 / 0.2e1) - (t42 - t41) * t316 + t483 * t106 + (t100 * t11 + t104 * t2 + t19 * t482 + t20 * t483 - t316 * t3 + t38 * t479) * m(7) + t11 * (-mrSges(7,1) * t314 - mrSges(7,3) * t187) + t101 * (-mrSges(6,1) * t314 + mrSges(6,2) * t187) - t314 * t445 - t314 * t446 + (t187 * t3 + t19 * t373 + t2 * t314 - t20 * t372) * mrSges(7,2) + (-t187 * t6 - t21 * t373 - t22 * t372 + t314 * t5) * mrSges(6,3) + (Ifges(7,5) * t187 - Ifges(7,3) * t314) * t442 + (Ifges(6,4) * t187 + Ifges(6,2) * t314) * t443 + (Ifges(5,1) * t256 - Ifges(5,4) * t323) * t430 + (t128 * t408 + t307) * qJD(3) - m(4) * (t198 * t210 + t199 * t211) + t475 * t72 + t476 * t164 + t477 * t165 + (Ifges(6,4) * t153 + Ifges(7,5) * t120 - Ifges(6,2) * t152 + Ifges(7,3) * t121) * t435 + t327 * t424 + t329 * t425 + t138 * t413 + (t40 + t43) * t104 + (t190 * t293 + t202 * t31 + t203 * t32 + t206 * t490 + t476 * t89 + t477 * t88) * m(5) + t325 * mrSges(4,3) + t100 * t16 - pkin(2) * t150 + t202 * t133 + t203 * t132 - t211 * t213 - t210 * t214 + t222 * t17 - t253 * t128 + t256 * t81 / 0.2e1 + t293 * t83 + t303 * t139 / 0.2e1; t310 + t454 + (t198 * t262 + t199 * t263) * mrSges(4,3) + (-Ifges(4,2) * t263 + t179 + t254) * t421 + (t192 * t88 + t193 * t89) * mrSges(5,3) + (Ifges(4,5) * t262 + Ifges(5,5) * t192 - Ifges(4,6) * t263 - Ifges(5,6) * t193) * t416 - t206 * (mrSges(5,1) * t193 + mrSges(5,2) * t192) + t193 * t468 + (-Ifges(5,2) * t193 + t114 + t396) * t467 + (Ifges(5,1) * t192 - t478) * t427 + ((t300 * t32 + t301 * t31) * pkin(3) - t206 * t410 - t88 * t95 - t89 * t96) * m(5) + t457 * t107 + t458 * t106 + (-t137 * t151 + t21 * t459 + t22 * t457 + t237 * t6 + t238 * t5) * m(6) + t374 * t459 + (-t19 * t459 + t2 * t230 + t20 * t458 + t233 * t3 - t38 * t46) * m(7) - (Ifges(6,4) * t435 + Ifges(7,5) * t436 + t464 * t418 + t465 * t433 + t472) * t124 + t178 * t419 + (Ifges(4,1) * t262 - t387) * t420 - t450 + (-t128 * t263 + t132 * t300 + t133 * t301) * pkin(3) - t46 * t71 + (-Ifges(6,2) * t435 + Ifges(7,3) * t436 + t418 * t462 + t433 * t463 - t447) * t471 - t151 * t72 - t96 * t164 - t95 * t165 - t198 * t213 + t199 * t214 + t230 * t40 + t233 * t42 + t237 * t41 + t238 * t43 - t275 * (mrSges(4,1) * t263 + mrSges(4,2) * t262); t374 * t471 - t375 * t124 - t192 * t164 + t193 * t165 + t16 + t17 + t83 + (t124 * t20 - t19 * t471 + t11) * m(7) + (t124 * t22 + t21 * t471 + t101) * m(6) + (-t192 * t89 + t193 * t88 + t190) * m(5); t310 + (t124 * t19 + t20 * t471) * mrSges(7,2) + (t374 + t400) * t22 + qJ(6) * t40 - pkin(5) * t42 - t70 * t71 + qJD(6) * t106 - t38 * (mrSges(7,1) * t471 + mrSges(7,3) * t124) + (Ifges(7,3) * t471 - t392) * t436 + t64 * t432 - t137 * (mrSges(6,1) * t471 - mrSges(6,2) * t124) + (t375 - t401) * t21 + (-t124 * t464 + t462 * t471) * t418 + (-pkin(5) * t3 + qJ(6) * t2 - t19 * t22 + t20 * t456 - t38 * t70) * m(7) + (-Ifges(6,2) * t471 - t119 + t461) * t435 + (-t124 * t465 + t118 - t395 + t61) * t433; -t277 * t106 + t471 * t71 + 0.2e1 * (t3 / 0.2e1 + t38 * t432 + t20 * t418) * m(7) + t42;];
tauc  = t1(:);

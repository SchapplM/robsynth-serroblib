% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:47
% EndTime: 2019-03-09 12:17:36
% DurationCPUTime: 31.70s
% Computational Cost: add. (9642->748), mult. (20051->952), div. (0->0), fcn. (13042->8), ass. (0->324)
t269 = sin(qJ(5));
t273 = cos(qJ(5));
t208 = t273 * mrSges(7,1) + t269 * mrSges(7,3);
t312 = t273 * pkin(5) + t269 * qJ(6);
t326 = mrSges(6,1) * t273 - mrSges(6,2) * t269;
t522 = m(7) * t312 + t208 + t326;
t581 = -mrSges(5,1) - t522;
t571 = Ifges(6,1) + Ifges(7,1);
t570 = Ifges(4,4) + Ifges(3,5);
t580 = -Ifges(6,4) + Ifges(7,5);
t569 = Ifges(7,4) + Ifges(6,5);
t568 = Ifges(7,2) + Ifges(6,3);
t567 = Ifges(4,6) - Ifges(3,6);
t566 = Ifges(6,6) - Ifges(7,6);
t271 = sin(qJ(2));
t394 = qJD(1) * t271;
t249 = pkin(7) * t394;
t192 = pkin(8) * t394 - t249;
t579 = -qJD(3) + t192;
t270 = sin(qJ(4));
t274 = cos(qJ(4));
t275 = cos(qJ(2));
t393 = qJD(1) * t275;
t161 = -t270 * t394 - t274 * t393;
t155 = qJD(5) - t161;
t272 = sin(qJ(1));
t305 = t275 * t270 - t271 * t274;
t145 = t305 * t272;
t177 = t270 * t271 + t274 * t275;
t146 = t177 * t272;
t575 = mrSges(7,2) + mrSges(6,3);
t545 = mrSges(5,2) - t575;
t578 = t581 * t145 - t146 * t545;
t276 = cos(qJ(1));
t398 = t275 * t276;
t402 = t271 * t276;
t147 = -t270 * t402 - t274 * t398;
t148 = t270 * t398 - t274 * t402;
t577 = t147 * t545 + t581 * t148;
t337 = t177 * mrSges(5,1) - mrSges(5,2) * t305;
t576 = t522 * t177 + t305 * t575 + t337;
t163 = -t270 * t393 + t274 * t394;
t377 = qJD(2) - qJD(4);
t121 = t273 * t163 - t269 * t377;
t263 = -qJDD(2) + qJDD(4);
t382 = qJD(1) * qJD(2);
t196 = -t275 * qJDD(1) + t271 * t382;
t197 = qJDD(1) * t271 + t275 * t382;
t289 = t177 * qJD(4);
t79 = -qJD(1) * t289 + t196 * t270 + t197 * t274;
t48 = qJD(5) * t121 - t273 * t263 + t269 * t79;
t491 = -t48 / 0.2e1;
t574 = Ifges(6,2) * t491;
t165 = -qJD(1) * pkin(1) - pkin(2) * t393 - qJ(3) * t394;
t129 = pkin(3) * t393 - t165;
t151 = Ifges(5,4) * t161;
t152 = Ifges(5,4) * t163;
t73 = -pkin(4) * t161 - pkin(9) * t163 + t129;
t277 = -pkin(2) - pkin(3);
t359 = t277 * qJD(2);
t139 = t359 - t579;
t250 = pkin(7) * t393;
t193 = -pkin(8) * t393 + t250;
t266 = qJD(2) * qJ(3);
t164 = t193 + t266;
t92 = t270 * t139 + t274 * t164;
t87 = -pkin(9) * t377 + t92;
t30 = -t269 * t87 + t273 * t73;
t23 = -pkin(5) * t155 + qJD(6) - t30;
t31 = t269 * t73 + t273 * t87;
t25 = qJ(6) * t155 + t31;
t174 = t197 * pkin(7);
t345 = qJDD(3) + t174;
t102 = -pkin(8) * t197 + qJDD(2) * t277 + t345;
t173 = t196 * pkin(7);
t131 = qJDD(2) * qJ(3) + qJD(2) * qJD(3) - t173;
t106 = pkin(8) * t196 + t131;
t386 = qJD(4) * t274;
t388 = qJD(4) * t270;
t29 = t102 * t274 - t270 * t106 - t139 * t388 - t164 * t386;
t27 = -pkin(4) * t263 - t29;
t28 = t270 * t102 + t274 * t106 + t139 * t386 - t164 * t388;
t439 = Ifges(7,5) * t273;
t313 = Ifges(7,3) * t269 + t439;
t445 = Ifges(6,4) * t273;
t318 = -Ifges(6,2) * t269 + t445;
t384 = qJD(5) * t273;
t346 = t384 / 0.2e1;
t385 = qJD(5) * t269;
t347 = -t385 / 0.2e1;
t474 = t269 / 0.2e1;
t112 = Ifges(7,5) * t121;
t120 = t269 * t163 + t273 * t377;
t55 = t155 * Ifges(7,6) + t120 * Ifges(7,3) + t112;
t361 = t55 * t474;
t401 = t273 * t161;
t411 = t161 * t269;
t424 = t163 * mrSges(5,3);
t427 = t161 * mrSges(5,3);
t440 = Ifges(7,5) * t269;
t446 = Ifges(6,4) * t269;
t477 = t163 / 0.2e1;
t480 = t155 / 0.2e1;
t481 = -t155 / 0.2e1;
t482 = t121 / 0.2e1;
t483 = -t121 / 0.2e1;
t484 = t120 / 0.2e1;
t485 = -t120 / 0.2e1;
t80 = qJD(1) * qJD(4) * t305 + t196 * t274 - t197 * t270;
t78 = qJDD(5) - t80;
t488 = t78 / 0.2e1;
t490 = t48 / 0.2e1;
t47 = -qJD(5) * t120 + t269 * t263 + t273 * t79;
t492 = t47 / 0.2e1;
t5 = pkin(5) * t48 - qJ(6) * t47 - qJD(6) * t121 + t27;
t519 = t273 * t571 + t440 - t446;
t520 = -t269 * t566 + t273 * t569;
t113 = Ifges(6,4) * t120;
t433 = t120 * Ifges(7,5);
t533 = t121 * t571 + t569 * t155 - t113 + t433;
t546 = Ifges(6,4) * t492 + t574 + Ifges(6,6) * t488 - t47 * Ifges(7,5) / 0.2e1 - t78 * Ifges(7,6) / 0.2e1 + Ifges(7,3) * t491;
t254 = t271 * qJD(3);
t267 = qJDD(1) * pkin(1);
t99 = t196 * pkin(2) - t197 * qJ(3) - qJD(1) * t254 - t267;
t77 = -pkin(3) * t196 - t99;
t18 = -pkin(4) * t80 - pkin(9) * t79 + t77;
t26 = pkin(9) * t263 + t28;
t3 = t269 * t18 + t273 * t26 + t73 * t384 - t385 * t87;
t1 = qJ(6) * t78 + qJD(6) * t155 + t3;
t4 = -qJD(5) * t31 + t18 * t273 - t26 * t269;
t2 = -pkin(5) * t78 + qJDD(6) - t4;
t334 = t1 * t273 + t2 * t269;
t553 = t23 * t384 - t25 * t385 + t334;
t333 = -t269 * t4 + t273 * t3;
t554 = -t30 * t384 - t31 * t385 + t333;
t325 = mrSges(6,1) * t269 + mrSges(6,2) * t273;
t91 = t274 * t139 - t270 * t164;
t86 = pkin(4) * t377 - t91;
t561 = t86 * t325;
t324 = t269 * mrSges(7,1) - t273 * mrSges(7,3);
t34 = t120 * pkin(5) - t121 * qJ(6) + t86;
t562 = t34 * t324;
t564 = -t120 * t566 + t121 * t569 + t155 * t568;
t565 = t47 * t571 + t580 * t48 + t569 * t78;
t432 = t121 * Ifges(6,4);
t58 = -t120 * Ifges(6,2) + t155 * Ifges(6,6) + t432;
t93 = Ifges(5,2) * t161 - t377 * Ifges(5,6) + t152;
t94 = Ifges(5,1) * t163 - t377 * Ifges(5,5) + t151;
t573 = (-t439 + t445) * t492 + t25 * (mrSges(7,2) * t411 - mrSges(7,3) * t163) + t31 * (mrSges(6,2) * t163 + mrSges(6,3) * t411) - t23 * (-mrSges(7,1) * t163 + mrSges(7,2) * t401) + t30 * (-mrSges(6,1) * t163 + mrSges(6,3) * t401) + t377 * (Ifges(5,5) * t161 - Ifges(5,6) * t163) / 0.2e1 - t129 * (mrSges(5,1) * t163 + mrSges(5,2) * t161) + t554 * mrSges(6,3) - t161 * t562 + (t313 * t484 + t318 * t485 + t480 * t520 + t482 * t519 + t361 + t561 + t562) * qJD(5) - t28 * mrSges(5,2) + t29 * mrSges(5,1) + (t58 / 0.2e1 - t55 / 0.2e1) * t411 - t161 * t561 - t27 * t326 + t440 * t490 + t446 * t491 + (-Ifges(7,3) * t490 + t488 * t566 + t546 + t574) * t273 + t553 * mrSges(7,2) - (Ifges(5,1) * t161 - t152 + t564) * t163 / 0.2e1 - (-Ifges(5,2) * t163 + t151 + t94) * t161 / 0.2e1 + (t161 * t520 + t163 * t568) * t481 + (t161 * t519 + t163 * t569) * t483 + (-t401 / 0.2e1 + t346) * t533 + (Ifges(6,6) * t163 + t161 * t318) * t484 + (Ifges(7,6) * t163 + t161 * t313) * t485 + Ifges(5,3) * t263 - t5 * t208 + t565 * t474 + (t488 * t569 + t492 * t571) * t269 + Ifges(5,5) * t79 + Ifges(5,6) * t80 + t58 * t347 + t92 * t424 + t91 * t427 + t93 * t477;
t534 = -m(6) - m(7);
t572 = m(6) * t86;
t16 = mrSges(6,1) * t48 + mrSges(6,2) * t47;
t70 = mrSges(5,1) * t263 - mrSges(5,3) * t79;
t563 = t70 - t16;
t529 = mrSges(5,1) * t377 + mrSges(6,1) * t120 + mrSges(6,2) * t121 + t424;
t199 = t274 * qJ(3) + t270 * t277;
t527 = -qJD(4) * t199 - t274 * t193 + t270 * t579;
t311 = pkin(5) * t269 - qJ(6) * t273;
t560 = -t269 * qJD(6) + t155 * t311;
t247 = Ifges(3,4) * t393;
t442 = Ifges(4,5) * t275;
t323 = t271 * Ifges(4,1) - t442;
t559 = Ifges(3,1) * t394 + qJD(1) * t323 + qJD(2) * t570 + t247;
t367 = mrSges(4,2) * t394;
t558 = -mrSges(3,3) * t394 - t367 + (mrSges(3,1) + mrSges(4,1)) * qJD(2);
t261 = t276 * pkin(7);
t256 = t271 * qJ(3);
t344 = -pkin(1) - t256;
t557 = t272 * (t275 * t277 + t344) - pkin(8) * t276 + t261;
t556 = t271 * t567 + t275 * t570;
t330 = t275 * mrSges(4,1) + t271 * mrSges(4,3);
t332 = mrSges(3,1) * t275 - mrSges(3,2) * t271;
t555 = t330 + t332;
t110 = qJD(2) * t177 - t289;
t296 = t110 * t269 - t305 * t384;
t551 = t47 * t569 - t48 * t566 + t568 * t78;
t550 = -t173 * t275 + t174 * t271;
t143 = -qJDD(2) * pkin(2) + t345;
t549 = t131 * t275 + t143 * t271;
t548 = -t177 * pkin(4) - pkin(9) * t305;
t19 = -mrSges(7,2) * t48 + mrSges(7,3) * t78;
t20 = mrSges(6,1) * t78 - mrSges(6,3) * t47;
t21 = -t78 * mrSges(7,1) + t47 * mrSges(7,2);
t22 = -mrSges(6,2) * t78 - mrSges(6,3) * t48;
t547 = (-t20 + t21) * t269 + (t19 + t22) * t273;
t391 = qJD(2) * t274;
t544 = -mrSges(2,1) - t555;
t543 = mrSges(5,3) - mrSges(3,3) - mrSges(4,2) + mrSges(2,2);
t499 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t540 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t539 = m(6) * ((-t269 * t31 - t273 * t30) * qJD(5) + t333) + m(7) * ((t23 * t273 - t25 * t269) * qJD(5) + t334) + t547;
t538 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t453 = mrSges(7,2) * t120;
t82 = mrSges(7,3) * t155 - t453;
t451 = mrSges(6,3) * t120;
t83 = -mrSges(6,2) * t155 - t451;
t455 = t82 + t83;
t450 = mrSges(6,3) * t121;
t84 = mrSges(6,1) * t155 - t450;
t452 = mrSges(7,2) * t121;
t85 = -mrSges(7,1) * t155 + t452;
t454 = -t85 + t84;
t487 = pkin(7) - pkin(8);
t217 = t487 * t271;
t218 = t487 * t275;
t123 = t217 * t270 + t218 * t274;
t260 = t275 * pkin(2);
t396 = t260 + t256;
t360 = t275 * pkin(3) + t396;
t510 = t360 - t548;
t90 = pkin(1) + t510;
t532 = t273 * t123 + t269 * t90;
t531 = -t92 + t560;
t530 = -t560 - t527;
t108 = t192 * t274 + t193 * t270;
t198 = -t270 * qJ(3) + t274 * t277;
t140 = t274 * qJD(3) + qJD(4) * t198;
t526 = -t108 + t140;
t524 = t274 * t217 - t218 * t270;
t191 = -pkin(9) + t199;
t518 = -t140 * t269 - t191 * t384;
t515 = t140 * t273 - t191 * t385;
t509 = -t148 * pkin(4) - pkin(9) * t147;
t508 = -t145 * pkin(4) + pkin(9) * t146;
t507 = g(1) * t276 + g(2) * t272;
t506 = qJD(3) + t249;
t390 = qJD(2) * t275;
t109 = t270 * t390 - t275 * t388 + (t386 - t391) * t271;
t397 = qJ(3) * t390 + t254;
t128 = t271 * t359 + t397;
t49 = pkin(4) * t109 - pkin(9) * t110 + t128;
t392 = qJD(2) * t271;
t194 = t487 * t392;
t195 = qJD(2) * t218;
t65 = qJD(4) * t524 - t194 * t274 + t195 * t270;
t13 = -qJD(5) * t532 - t269 * t65 + t273 * t49;
t495 = m(5) / 0.2e1;
t494 = m(6) / 0.2e1;
t493 = m(7) / 0.2e1;
t473 = t271 / 0.2e1;
t472 = pkin(5) * t163;
t471 = pkin(7) * t271;
t470 = pkin(7) * t275;
t104 = pkin(4) * t163 - pkin(9) * t161;
t51 = t269 * t104 + t273 * t91;
t448 = Ifges(3,4) * t271;
t447 = Ifges(3,4) * t275;
t443 = Ifges(4,5) * t271;
t420 = t275 * mrSges(4,3);
t239 = qJ(3) * t393;
t368 = t271 * t277;
t144 = qJD(1) * t368 + t239;
t74 = -t104 + t144;
t39 = t273 * t108 + t269 * t74;
t418 = qJ(6) * t163;
t416 = t110 * t273;
t408 = t305 * t269;
t395 = t276 * pkin(1) + t272 * pkin(7);
t389 = qJD(4) * t269;
t387 = qJD(4) * t273;
t369 = m(4) + m(5) - t534;
t366 = mrSges(4,2) * t393;
t339 = -t382 / 0.2e1;
t114 = t146 * t269 - t276 * t273;
t149 = -qJDD(2) * mrSges(4,1) + t197 * mrSges(4,2);
t336 = pkin(2) * t398 + qJ(3) * t402 + t395;
t331 = mrSges(3,1) * t271 + mrSges(3,2) * t275;
t320 = t275 * Ifges(3,2) + t448;
t50 = t104 * t273 - t269 * t91;
t38 = -t108 * t269 + t273 * t74;
t53 = -t123 * t269 + t273 * t90;
t115 = t146 * t273 + t269 * t276;
t200 = -qJD(2) * pkin(2) + t506;
t205 = t250 + t266;
t306 = t200 * t275 - t205 * t271;
t229 = t275 * t272 * qJ(3);
t303 = t272 * t368 + t229;
t230 = qJ(3) * t398;
t302 = t276 * t368 + t230;
t201 = -pkin(4) - t312;
t297 = pkin(1) * t331;
t295 = -t305 * t385 - t416;
t12 = -t123 * t385 + t269 * t49 + t273 * t65 + t90 * t384;
t293 = t165 * (t271 * mrSges(4,1) - t420);
t292 = t271 * (Ifges(3,1) * t275 - t448);
t291 = t275 * (Ifges(4,3) * t271 + t442);
t288 = pkin(3) * t398 - pkin(8) * t272 + t336;
t287 = t420 + (-m(4) * pkin(2) - mrSges(4,1)) * t271;
t66 = qJD(4) * t123 - t194 * t270 - t274 * t195;
t246 = Ifges(4,5) * t394;
t207 = qJD(2) * mrSges(4,3) + t366;
t206 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t393;
t202 = -pkin(1) - t396;
t190 = pkin(4) - t198;
t181 = pkin(2) * t394 - t239;
t180 = t330 * qJD(1);
t168 = pkin(1) + t360;
t162 = t269 * t394 + t273 * t391;
t160 = t269 * t391 - t273 * t394;
t157 = Ifges(3,6) * qJD(2) + qJD(1) * t320;
t156 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t393 + t246;
t154 = pkin(2) * t392 - t397;
t150 = -mrSges(4,2) * t196 + qJDD(2) * mrSges(4,3);
t130 = -t198 - t201;
t126 = mrSges(5,2) * t377 + t427;
t119 = -t147 * t273 - t272 * t269;
t118 = -t147 * t269 + t272 * t273;
t103 = t325 * t305;
t100 = -mrSges(5,1) * t161 + mrSges(5,2) * t163;
t71 = -mrSges(5,2) * t263 + mrSges(5,3) * t80;
t68 = mrSges(7,1) * t120 - mrSges(7,3) * t121;
t67 = pkin(5) * t121 + qJ(6) * t120;
t63 = -t305 * t311 - t524;
t40 = -pkin(5) * t177 - t53;
t37 = qJ(6) * t177 + t532;
t36 = -t50 - t472;
t35 = t418 + t51;
t33 = -t38 + t472;
t32 = -t418 + t39;
t15 = mrSges(7,1) * t48 - mrSges(7,3) * t47;
t14 = t311 * t110 - (qJD(5) * t312 - qJD(6) * t273) * t305 + t66;
t7 = -pkin(5) * t109 - t13;
t6 = qJ(6) * t109 + qJD(6) * t177 + t12;
t8 = [-t202 * mrSges(4,3) * t197 + ((-t2 * mrSges(7,2) + t4 * mrSges(6,3) - t565 / 0.2e1) * t273 + t29 * mrSges(5,3) - Ifges(5,1) * t79 - Ifges(5,4) * t80 - Ifges(5,5) * t263 - t313 * t490 - t318 * t491 - t324 * t5 - t346 * t55 - t347 * t533 - t488 * t520 - t492 * t519) * t305 + (t1 * mrSges(7,2) + t3 * mrSges(6,3) + t546) * t408 + (m(3) * pkin(1) ^ 2 + Ifges(2,3)) * qJDD(1) + (-t109 * t92 - t110 * t91) * mrSges(5,3) + (t551 / 0.2e1 - t28 * mrSges(5,3) - Ifges(5,4) * t79 - Ifges(5,2) * t80 - Ifges(5,6) * t263 + Ifges(6,6) * t491 + Ifges(7,6) * t490 + t488 * t568 + t569 * t492 + t538) * t177 + (-t157 / 0.2e1 + t156 / 0.2e1) * t392 + (t271 * (Ifges(4,1) * t275 + t443) + t275 * (-Ifges(3,2) * t271 + t447) + t292) * t382 / 0.2e1 + (t271 * Ifges(3,1) + t323 + t447) * t197 / 0.2e1 + m(6) * (t12 * t31 + t13 * t30 + t3 * t532 + t4 * t53) + t532 * t22 + t332 * t267 + t559 * t390 / 0.2e1 + t37 * t19 - t297 * t382 - t377 * (Ifges(5,5) * t110 - Ifges(5,6) * t109) / 0.2e1 + (t443 / 0.2e1 - pkin(1) * mrSges(3,1) + t202 * mrSges(4,1) - mrSges(3,3) * t470 - t320 / 0.2e1 + (Ifges(4,5) - Ifges(3,4)) * t473 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t275) * t196 + m(5) * (t123 * t28 + t128 * t129 + t168 * t77 + t65 * t92) - t275 * (Ifges(4,5) * t197 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (t293 + t556 * qJD(2) / 0.2e1) * qJD(2) + m(4) * (t154 * t165 + t202 * t99) + (-pkin(1) * t197 - qJDD(2) * t470) * mrSges(3,2) + t275 * (Ifges(3,4) * t197 + Ifges(3,6) * qJDD(2)) / 0.2e1 + m(7) * (t1 * t37 + t14 * t34 + t2 * t40 + t23 * t7 + t25 * t6 + t5 * t63) + (t200 * t390 - t205 * t392 + t549) * mrSges(4,2) + (t197 * t471 + t550) * mrSges(3,3) - t296 * t58 / 0.2e1 + t77 * t337 - (-m(5) * t29 + m(6) * t27 - t563) * t524 + (-qJDD(2) * mrSges(3,1) + t149) * t471 + ((-t207 - t206) * t392 - t558 * t390 + m(4) * (qJD(2) * t306 + t549) + m(3) * t550) * pkin(7) + (-Ifges(7,5) * t295 + Ifges(7,6) * t109 + Ifges(7,3) * t296) * t484 + (-Ifges(6,4) * t295 - Ifges(6,2) * t296 + Ifges(6,6) * t109) * t485 - t154 * t180 + t40 * t21 + (t569 * t109 - t571 * t295 + t296 * t580) * t482 + t533 * t416 / 0.2e1 + t564 * t109 / 0.2e1 + t53 * t20 + (t109 * t568 - t295 * t569 - t296 * t566) * t480 + ((Ifges(3,1) + Ifges(4,1)) * t197 + t570 * qJDD(2)) * t473 + (t271 * t570 - t275 * t567) * qJDD(2) / 0.2e1 + t63 * t15 + t30 * (mrSges(6,1) * t109 + mrSges(6,3) * t295) + t23 * (-mrSges(7,1) * t109 - mrSges(7,2) * t295) + t14 * t68 + t34 * (mrSges(7,1) * t296 + mrSges(7,3) * t295) + t86 * (mrSges(6,1) * t296 - mrSges(6,2) * t295) + t31 * (-mrSges(6,2) * t109 - mrSges(6,3) * t296) + t25 * (-mrSges(7,2) * t296 + mrSges(7,3) * t109) + t6 * t82 + t12 * t83 + t13 * t84 + t7 * t85 + (-m(5) * t91 + t529 + t572) * t66 + (t146 * mrSges(5,1) - t557 * m(5) + t534 * (-t146 * pkin(4) - pkin(9) * t145 + t557) + (-m(4) - m(3)) * t261 + t499 * t115 + t540 * t114 + (-m(4) * (t344 - t260) + m(3) * pkin(1) - t544) * t272 - t545 * t145 + t543 * t276) * g(1) + (-m(3) * t395 - m(4) * t336 - m(5) * t288 + t147 * mrSges(5,1) + t534 * (-t147 * pkin(4) + pkin(9) * t148 + t288) - t499 * t119 - t540 * t118 + t544 * t276 + t545 * t148 + t543 * t272) * g(2) - t27 * t103 - t109 * t93 / 0.2e1 + t110 * t94 / 0.2e1 + t291 * t339 + t110 * t361 + t150 * t470 + (Ifges(5,1) * t110 - Ifges(5,4) * t109) * t477 + t123 * t71 + t65 * t126 + t128 * t100 + t129 * (mrSges(5,1) * t109 + mrSges(5,2) * t110) - t99 * t330 + t161 * (Ifges(5,4) * t110 - Ifges(5,2) * t109) / 0.2e1 + t168 * (-mrSges(5,1) * t80 + mrSges(5,2) * t79); t556 * t339 + t558 * t250 - (-Ifges(3,2) * t394 + t247 + t559) * t393 / 0.2e1 - t527 * t529 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) + t157 * t394 / 0.2e1 - t573 - (Ifges(4,1) * t393 + t156 + t246) * t394 / 0.2e1 + t205 * t367 + t206 * t249 + (-m(4) * t396 - m(5) * t360 + t510 * t534 - t555 - t576) * g(3) + (-m(4) * t230 - m(5) * t302 - t276 * t287 + t534 * (t302 - t509) + t577) * g(1) + (-m(4) * t229 - m(5) * t303 - t272 * t287 + t534 * (t303 - t508) + t578) * g(2) - t200 * t366 + t539 * t191 + (t190 * t27 + (-t269 * t30 + t273 * t31) * t140 - t30 * t38 - t31 * t39 - t527 * t86) * m(6) + (t130 * t5 + (t23 * t269 + t25 * t273) * t140 - t23 * t33 - t25 * t32 + t530 * t34) * m(7) + (-pkin(2) * t143 + qJ(3) * t131 + qJD(3) * t205 - t165 * t181) * m(4) + (-pkin(7) * t306 * m(4) - t293 + (t297 + t291 / 0.2e1 - t292 / 0.2e1) * qJD(1)) * qJD(1) + t198 * t70 + t199 * t71 + t190 * t16 + t181 * t180 + t173 * mrSges(3,2) - t174 * mrSges(3,1) + t506 * t207 + t507 * t331 + (-t39 + t515) * t83 + (-t32 + t515) * t82 + (-t33 - t518) * t85 + (-t38 + t518) * t84 + t567 * t196 + t570 * t197 + t526 * t126 + (-t129 * t144 + t198 * t29 + t199 * t28 + t526 * t92 + t527 * t91) * m(5) + t530 * t68 + t130 * t15 + t131 * mrSges(4,3) - t143 * mrSges(4,1) - t144 * t100 - pkin(2) * t149 + qJ(3) * t150; -qJD(2) * t207 - t455 * t162 + t454 * t160 + t369 * t275 * g(3) + (-qJD(2) * t126 - t15 + (-t269 * t454 + t273 * t455 + t126) * qJD(4) + t563) * t274 + ((-t100 - t180) * qJD(1) - t507 * t369) * t271 + (t71 + (-t269 * t455 - t273 * t454) * qJD(5) + t377 * (-t68 - t529) + t547) * t270 - m(7) * (t160 * t23 + t162 * t25) - m(6) * (-t160 * t30 + t162 * t31) + 0.2e1 * ((qJD(4) * t92 + t29) * t495 + (-t30 * t389 + t31 * t387 - t27) * t494 + (t23 * t389 + t25 * t387 - t5) * t493) * t274 + 0.2e1 * ((-qJD(4) * t91 + t28) * t495 + (qJD(4) * t86 + t554) * t494 + (qJD(4) * t34 + t553) * t493 + (-m(7) * t34 / 0.2e1 - t572 / 0.2e1 + t91 * t495) * qJD(2)) * t270 + t149 + (-t129 * t394 - t391 * t92) * m(5) + (-qJD(2) * t205 + t165 * t394 + t143) * m(4); (-pkin(4) * t27 - t30 * t50 - t31 * t51 - t86 * t92) * m(6) + (t201 * t5 - t23 * t36 - t25 * t35 + t34 * t531) * m(7) + (-t384 * t454 - t385 * t455 + t539) * pkin(9) - pkin(4) * t16 + t201 * t15 - t35 * t82 - t51 * t83 - t50 * t84 - t36 * t85 - t529 * t92 + t531 * t68 - t91 * t126 + (t508 * t534 - t578) * g(2) + (t509 * t534 - t577) * g(1) + (t534 * t548 + t576) * g(3) + t573; (t118 * t499 - t119 * t540) * g(1) + (t114 * t499 - t115 * t540) * g(2) + ((-m(7) * t311 - t324) * t305 - t103) * g(3) + t538 + qJ(6) * t19 - pkin(5) * t21 + (-t120 * t569 - t121 * t566) * t481 + (-t120 * t571 + t112 - t432 + t55) * t483 + (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t25 - t34 * t67) * m(7) + (Ifges(7,3) * t121 - t433) * t485 - t67 * t68 + qJD(6) * t82 + t551 + t25 * t452 + t23 * t453 + t58 * t482 + (-m(7) * t23 + t450 + t454) * t31 + (-m(7) * t25 - t451 - t455) * t30 - t34 * (mrSges(7,1) * t121 + mrSges(7,3) * t120) + (-Ifges(6,2) * t121 - t113 + t533) * t484 - t86 * (mrSges(6,1) * t121 - mrSges(6,2) * t120); t121 * t68 - t155 * t82 + (-g(1) * t118 - g(2) * t114 + g(3) * t408 + t34 * t121 - t25 * t155 + t2) * m(7) + t21;];
tau  = t8;

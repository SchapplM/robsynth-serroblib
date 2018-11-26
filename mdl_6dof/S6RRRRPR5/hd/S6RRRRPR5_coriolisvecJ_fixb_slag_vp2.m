% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:14:36
% EndTime: 2018-11-23 18:14:49
% DurationCPUTime: 12.69s
% Computational Cost: add. (14198->713), mult. (34124->932), div. (0->0), fcn. (23960->8), ass. (0->346)
t563 = Ifges(5,1) + Ifges(6,1);
t549 = Ifges(5,5) + Ifges(6,4);
t562 = Ifges(5,6) - Ifges(6,6);
t343 = sin(qJ(4));
t417 = qJD(4) * t343;
t344 = sin(qJ(3));
t348 = cos(qJ(3));
t349 = cos(qJ(2));
t419 = qJD(1) * t349;
t345 = sin(qJ(2));
t420 = qJD(1) * t345;
t294 = -t344 * t420 + t348 * t419;
t431 = t294 * t343;
t564 = t417 - t431;
t471 = t343 / 0.2e1;
t561 = t564 * pkin(10);
t347 = cos(qJ(4));
t540 = -t343 * t562 + t347 * t549;
t452 = Ifges(6,5) * t343;
t456 = Ifges(5,4) * t343;
t539 = t347 * t563 + t452 - t456;
t310 = t344 * t349 + t345 * t348;
t295 = t310 * qJD(1);
t334 = -pkin(2) * t349 - pkin(1);
t323 = qJD(1) * t334;
t205 = -pkin(3) * t294 - pkin(9) * t295 + t323;
t496 = -pkin(8) - pkin(7);
t327 = t496 * t349;
t313 = qJD(1) * t327;
t297 = t348 * t313;
t325 = t496 * t345;
t312 = qJD(1) * t325;
t299 = qJD(2) * pkin(2) + t312;
t249 = t344 * t299 - t297;
t413 = qJD(2) + qJD(3);
t219 = pkin(9) * t413 + t249;
t113 = t347 * t205 - t343 * t219;
t509 = qJD(5) - t113;
t332 = pkin(2) * t344 + pkin(9);
t449 = pkin(2) * qJD(3);
t408 = t348 * t449;
t244 = pkin(3) * t295 - pkin(9) * t294;
t213 = pkin(2) * t420 + t244;
t296 = t344 * t313;
t252 = t312 * t348 + t296;
t137 = t343 * t213 + t347 * t252;
t288 = t295 * qJ(5);
t97 = t288 + t137;
t559 = -t332 * t417 + t347 * t408 + t561 - t97;
t245 = t343 * t252;
t463 = -pkin(10) + t332;
t306 = t463 * t347;
t350 = -pkin(4) - pkin(5);
t404 = t350 * t295;
t465 = pkin(10) * t294;
t558 = qJD(4) * t306 + t343 * t408 - t245 - (-t213 - t465) * t347 - t404;
t248 = t348 * t299 + t296;
t144 = t343 * t244 + t347 * t248;
t106 = t288 + t144;
t557 = -pkin(9) * t417 - t106 + t561;
t231 = t343 * t248;
t495 = pkin(9) - pkin(10);
t326 = t495 * t347;
t556 = qJD(4) * t326 - t231 - (-t244 - t465) * t347 - t404;
t266 = t347 * t295 + t343 * t413;
t555 = -pkin(10) * t266 + t509;
t265 = t343 * t295 - t347 * t413;
t264 = Ifges(5,4) * t265;
t291 = qJD(4) - t294;
t453 = Ifges(6,5) * t265;
t514 = t266 * t563 + t549 * t291 - t264 + t453;
t342 = sin(qJ(6));
t346 = cos(qJ(6));
t174 = t265 * t346 - t266 * t342;
t172 = Ifges(7,4) * t174;
t369 = t265 * t342 + t266 * t346;
t554 = Ifges(7,2) * t369 - t172;
t114 = t343 * t205 + t347 * t219;
t67 = t291 * t350 + t555;
t282 = t291 * qJ(5);
t84 = pkin(10) * t265 + t114;
t75 = t282 + t84;
t20 = -t342 * t75 + t346 * t67;
t21 = t342 * t67 + t346 * t75;
t392 = Ifges(4,6) * t413;
t89 = -pkin(4) * t291 + t509;
t90 = t282 + t114;
t553 = -t323 * mrSges(4,1) - t113 * mrSges(5,1) + t89 * mrSges(6,1) + t20 * mrSges(7,1) + t114 * mrSges(5,2) - t21 * mrSges(7,2) - t90 * mrSges(6,3) + t392 / 0.2e1;
t218 = -pkin(3) * t413 - t248;
t356 = t266 * qJ(5) - t218;
t103 = t265 * pkin(4) - t356;
t457 = Ifges(5,4) * t266;
t148 = -Ifges(5,2) * t265 + Ifges(5,6) * t291 + t457;
t386 = mrSges(6,1) * t343 - mrSges(6,3) * t347;
t388 = mrSges(5,1) * t343 + mrSges(5,2) * t347;
t552 = -t103 * t386 + t148 * t471 - t218 * t388;
t393 = Ifges(4,5) * t413;
t551 = t323 * mrSges(4,2) + t393 / 0.2e1;
t259 = t413 * t310;
t239 = t259 * qJD(1);
t233 = Ifges(7,3) * t239;
t454 = Ifges(7,4) * t369;
t285 = qJD(6) - t291;
t479 = -t285 / 0.2e1;
t487 = -t369 / 0.2e1;
t308 = t344 * t345 - t348 * t349;
t258 = t413 * t308;
t238 = t258 * qJD(1);
t163 = -qJD(4) * t265 - t347 * t238;
t164 = qJD(4) * t266 - t343 * t238;
t46 = qJD(6) * t174 + t163 * t346 + t164 * t342;
t47 = -qJD(6) * t369 - t163 * t342 + t164 * t346;
t528 = Ifges(7,5) * t46 + Ifges(7,6) * t47;
t82 = t265 * t350 + t356;
t550 = (t174 * t20 + t21 * t369) * mrSges(7,3) + (Ifges(7,5) * t174 - Ifges(7,6) * t369) * t479 + (Ifges(7,1) * t174 - t454) * t487 - t82 * (mrSges(7,1) * t369 + mrSges(7,2) * t174) - t233 + t528;
t548 = t549 * t239 + (-Ifges(5,4) + Ifges(6,5)) * t164 + t563 * t163;
t402 = qJD(2) * t496;
t391 = qJD(1) * t402;
t300 = t345 * t391;
t366 = t349 * t391;
t154 = qJD(3) * t249 + t300 * t344 - t348 * t366;
t546 = -qJ(5) * t163 - qJD(5) * t266 + t154;
t368 = t342 * t347 - t343 * t346;
t203 = t368 * t294;
t508 = qJD(4) - qJD(6);
t257 = t508 * t368;
t545 = t203 - t257;
t367 = t342 * t343 + t346 * t347;
t204 = t367 * t294;
t256 = t508 * t367;
t544 = t204 - t256;
t416 = qJD(4) * t347;
t290 = pkin(4) * t417 - qJ(5) * t416 - t343 * qJD(5);
t543 = -pkin(5) * t564 - t290;
t541 = t343 * t549 + t347 * t562;
t451 = Ifges(6,5) * t347;
t455 = Ifges(5,4) * t347;
t538 = t343 * t563 - t451 + t455;
t418 = qJD(2) * t345;
t410 = pkin(2) * t418;
t123 = pkin(3) * t239 + pkin(9) * t238 + qJD(1) * t410;
t153 = qJD(3) * t248 + t348 * t300 + t344 * t366;
t27 = t343 * t123 + t347 * t153 + t205 * t416 - t219 * t417;
t28 = t123 * t347 - t343 * t153 - t205 * t417 - t219 * t416;
t536 = t27 * t347 - t28 * t343;
t17 = t239 * qJ(5) + t291 * qJD(5) + t27;
t23 = -pkin(4) * t239 - t28;
t535 = t17 * t347 + t23 * t343 + t89 * t416;
t375 = Ifges(6,3) * t343 + t451;
t381 = -Ifges(5,2) * t343 + t455;
t263 = Ifges(6,5) * t266;
t145 = Ifges(6,6) * t291 + Ifges(6,3) * t265 + t263;
t428 = t343 * t145;
t480 = t266 / 0.2e1;
t482 = t265 / 0.2e1;
t483 = -t265 / 0.2e1;
t530 = (-t113 * t347 - t114 * t343) * mrSges(5,3) + t375 * t482 + t381 * t483 + t428 / 0.2e1 - t552 + t539 * t480 + t540 * t291 / 0.2e1;
t502 = t46 / 0.2e1;
t501 = t47 / 0.2e1;
t485 = -t239 / 0.2e1;
t305 = t463 * t343;
t240 = t305 * t346 - t306 * t342;
t523 = qJD(6) * t240 + t558 * t342 + t346 * t559;
t241 = t305 * t342 + t306 * t346;
t522 = -qJD(6) * t241 - t342 * t559 + t558 * t346;
t324 = t495 * t343;
t267 = t324 * t346 - t326 * t342;
t521 = qJD(6) * t267 + t342 * t556 + t346 * t557;
t269 = t324 * t342 + t326 * t346;
t520 = -qJD(6) * t269 - t342 * t557 + t346 * t556;
t319 = t346 * qJ(5) + t342 * t350;
t519 = -qJD(6) * t319 - t342 * t555 - t346 * t84;
t318 = -t342 * qJ(5) + t346 * t350;
t518 = qJD(6) * t318 - t342 * t84 + t346 * t555;
t74 = mrSges(5,1) * t164 + mrSges(5,2) * t163;
t517 = -m(5) * t154 - t74;
t223 = t368 * t310;
t430 = t294 * t347;
t421 = pkin(4) * t431 - qJ(5) * t430;
t152 = t249 + t421;
t516 = t152 + t543;
t251 = t312 * t344 - t297;
t156 = t251 + t421;
t409 = t344 * t449;
t515 = t156 - t409 + t543;
t513 = -t152 + t290;
t512 = t348 * t325 + t327 * t344;
t339 = t343 * qJ(5);
t511 = t347 * pkin(4) + t339;
t461 = mrSges(5,3) * t266;
t195 = mrSges(5,1) * t291 - t461;
t196 = -mrSges(6,1) * t291 + mrSges(6,2) * t266;
t422 = -t195 + t196;
t193 = -mrSges(6,2) * t265 + mrSges(6,3) * t291;
t462 = mrSges(5,3) * t265;
t194 = -mrSges(5,2) * t291 - t462;
t423 = t193 + t194;
t92 = -mrSges(6,2) * t164 + mrSges(6,3) * t239;
t93 = mrSges(5,1) * t239 - mrSges(5,3) * t163;
t94 = -t239 * mrSges(6,1) + t163 * mrSges(6,2);
t95 = -mrSges(5,2) * t239 - mrSges(5,3) * t164;
t506 = m(6) * (-t417 * t90 + t535) + m(5) * (-t113 * t416 - t114 * t417 + t536) + (t92 + t95) * t347 + (-t93 + t94) * t343 + (-t343 * t423 + t347 * t422) * qJD(4);
t270 = t325 * t344 - t327 * t348;
t316 = t345 * t402;
t317 = t349 * t402;
t181 = t270 * qJD(3) + t316 * t344 - t348 * t317;
t505 = m(4) / 0.2e1;
t504 = Ifges(7,4) * t502 + Ifges(7,2) * t501 + Ifges(7,6) * t485;
t503 = Ifges(7,1) * t502 + Ifges(7,4) * t501 + Ifges(7,5) * t485;
t71 = Ifges(7,2) * t174 + Ifges(7,6) * t285 + t454;
t500 = -t71 / 0.2e1;
t499 = t71 / 0.2e1;
t72 = Ifges(7,1) * t369 + Ifges(7,5) * t285 + t172;
t498 = -t72 / 0.2e1;
t497 = t72 / 0.2e1;
t494 = pkin(1) * mrSges(3,1);
t493 = pkin(1) * mrSges(3,2);
t492 = t163 / 0.2e1;
t491 = -t164 / 0.2e1;
t490 = t164 / 0.2e1;
t489 = -t174 / 0.2e1;
t488 = t174 / 0.2e1;
t486 = t369 / 0.2e1;
t484 = t239 / 0.2e1;
t481 = -t266 / 0.2e1;
t478 = t285 / 0.2e1;
t477 = -t291 / 0.2e1;
t475 = -t294 / 0.2e1;
t474 = -t295 / 0.2e1;
t472 = -t343 / 0.2e1;
t470 = -t347 / 0.2e1;
t469 = t347 / 0.2e1;
t468 = m(4) * t323;
t466 = pkin(4) * t295;
t464 = pkin(10) * t310;
t460 = Ifges(4,1) * t295;
t459 = Ifges(3,4) * t345;
t289 = Ifges(4,4) * t294;
t458 = Ifges(4,4) * t295;
t450 = Ifges(4,2) * t294;
t447 = t174 * Ifges(7,6);
t446 = t369 * Ifges(7,5);
t442 = t285 * Ifges(7,3);
t441 = t294 * mrSges(4,3);
t440 = t295 * mrSges(4,3);
t439 = Ifges(3,5) * qJD(2);
t438 = Ifges(3,6) * qJD(2);
t436 = qJ(5) * t265;
t435 = qJ(5) * t347;
t434 = qJD(2) * mrSges(3,1);
t433 = qJD(2) * mrSges(3,2);
t432 = t154 * t512;
t426 = t343 * t348;
t425 = t347 * t348;
t424 = -mrSges(4,1) * t413 + mrSges(5,1) * t265 + mrSges(5,2) * t266 + t440;
t247 = pkin(3) * t308 - pkin(9) * t310 + t334;
t167 = t343 * t247 + t347 * t270;
t414 = qJD(5) * t347;
t411 = pkin(3) + t511;
t407 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t406 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t405 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t126 = t308 * qJ(5) + t167;
t333 = -pkin(2) * t348 - pkin(3);
t401 = t439 / 0.2e1;
t400 = -t438 / 0.2e1;
t136 = t213 * t347 - t245;
t143 = t244 * t347 - t231;
t260 = t343 * t270;
t166 = t247 * t347 - t260;
t11 = -pkin(10) * t163 + t239 * t350 - t28;
t12 = pkin(10) * t164 + t17;
t2 = qJD(6) * t20 + t11 * t342 + t12 * t346;
t3 = -qJD(6) * t21 + t11 * t346 - t12 * t342;
t390 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t389 = mrSges(5,1) * t347 - mrSges(5,2) * t343;
t387 = mrSges(6,1) * t347 + mrSges(6,3) * t343;
t380 = Ifges(5,2) * t347 + t456;
t374 = -Ifges(6,3) * t347 + t452;
t373 = pkin(4) * t343 - t435;
t88 = t260 + (-t247 - t464) * t347 + t350 * t308;
t96 = t343 * t464 + t126;
t42 = -t342 * t96 + t346 * t88;
t43 = t342 * t88 + t346 * t96;
t302 = t333 - t511;
t124 = -mrSges(7,2) * t285 + mrSges(7,3) * t174;
t125 = mrSges(7,1) * t285 - mrSges(7,3) * t369;
t370 = t124 * t346 - t125 * t342;
t162 = pkin(3) * t259 + pkin(9) * t258 + t410;
t180 = qJD(3) * t512 + t316 * t348 + t317 * t344;
t49 = t162 * t347 - t343 * t180 - t247 * t417 - t270 * t416;
t365 = t343 * t350 + t435;
t48 = t343 * t162 + t347 * t180 + t247 * t416 - t270 * t417;
t29 = t259 * qJ(5) + t308 * qJD(5) + t48;
t354 = t28 * mrSges(5,1) - t23 * mrSges(6,1) - t27 * mrSges(5,2) + t17 * mrSges(6,3) - t390;
t146 = t266 * Ifges(5,5) - t265 * Ifges(5,6) + t291 * Ifges(5,3);
t147 = t266 * Ifges(6,4) + t291 * Ifges(6,2) + t265 * Ifges(6,6);
t18 = t164 * t350 - t546;
t214 = t392 + t450 + t458;
t215 = t289 + t393 + t460;
t37 = pkin(4) * t164 + t546;
t58 = Ifges(6,5) * t163 + Ifges(6,6) * t239 + Ifges(6,3) * t164;
t59 = Ifges(5,4) * t163 - Ifges(5,2) * t164 + Ifges(5,6) * t239;
t70 = t442 + t446 + t447;
t351 = (-t2 * t367 + t20 * t544 + t21 * t545 + t3 * t368) * mrSges(7,3) + (-mrSges(7,1) * t545 - mrSges(7,2) * t544) * t82 + t538 * t492 + t541 * t484 + (t113 * t430 + t114 * t431 + t536) * mrSges(5,3) + t530 * qJD(4) + (Ifges(4,1) * t474 + t375 * t483 + t381 * t482 + t540 * t477 + t539 * t481 - t551 + t552) * t294 + (-Ifges(4,2) * t475 - Ifges(7,3) * t479 - Ifges(7,5) * t487 - Ifges(7,6) * t489 + Ifges(5,6) * t482 + Ifges(6,6) * t483 + t549 * t481 + (Ifges(5,3) + Ifges(6,2)) * t477 + t553) * t295 + (Ifges(7,1) * t256 - Ifges(7,4) * t257) * t486 + (Ifges(7,4) * t256 - Ifges(7,2) * t257) * t488 + (t214 + t70) * t295 / 0.2e1 - t367 * t504 + (-Ifges(7,5) * t368 - Ifges(7,6) * t367) * t485 + (-Ifges(7,4) * t368 - Ifges(7,2) * t367) * t501 + (-Ifges(7,1) * t368 - Ifges(7,4) * t367) * t502 + t18 * (mrSges(7,1) * t367 - mrSges(7,2) * t368) - t368 * t503 + (Ifges(7,5) * t256 - Ifges(7,6) * t257) * t478 + t514 * (t416 / 0.2e1 - t430 / 0.2e1) + (t147 + t146 - t458) * t474 + (t289 + t428 + t215) * t475 + (Ifges(7,4) * t204 - Ifges(7,2) * t203) * t489 + (-mrSges(4,1) - t389) * t154 + (Ifges(7,1) * t204 - Ifges(7,4) * t203) * t487 + (Ifges(7,5) * t204 - Ifges(7,6) * t203) * t479 + t59 * t469 + t58 * t470 + t248 * t441 - t153 * mrSges(4,2) - t37 * t387 + (-t89 * t430 - t564 * t90 + t535) * mrSges(6,2) + t548 * t471 - Ifges(4,5) * t238 - Ifges(4,6) * t239 + t374 * t490 + t380 * t491 + t256 * t497 + t204 * t498 - t257 * t499 - t203 * t500;
t340 = t347 * pkin(5);
t335 = Ifges(3,4) * t419;
t322 = mrSges(3,3) * t419 - t433;
t321 = -mrSges(3,3) * t420 + t434;
t303 = t340 + t411;
t293 = Ifges(3,1) * t420 + t335 + t439;
t292 = t438 + (Ifges(3,2) * t349 + t459) * qJD(1);
t281 = -t302 + t340;
t276 = t290 + t409;
t273 = -mrSges(4,2) * t413 + t441;
t243 = -mrSges(4,1) * t294 + mrSges(4,2) * t295;
t236 = Ifges(6,2) * t239;
t234 = Ifges(5,3) * t239;
t224 = t367 * t310;
t185 = mrSges(6,1) * t265 - mrSges(6,3) * t266;
t184 = pkin(4) * t266 + t436;
t177 = t310 * t373 - t512;
t161 = Ifges(6,4) * t163;
t160 = Ifges(5,5) * t163;
t159 = Ifges(5,6) * t164;
t158 = Ifges(6,6) * t164;
t138 = t310 * t365 + t512;
t131 = -pkin(4) * t308 - t166;
t121 = t266 * t350 - t436;
t108 = -t143 - t466;
t98 = -t136 - t466;
t81 = -mrSges(7,1) * t174 + mrSges(7,2) * t369;
t73 = mrSges(6,1) * t164 - mrSges(6,3) * t163;
t69 = t256 * t310 + t258 * t368;
t68 = t223 * t508 - t367 * t258;
t57 = -t373 * t258 + (qJD(4) * t511 - t414) * t310 + t181;
t41 = -t365 * t258 + (t414 + (t347 * t350 - t339) * qJD(4)) * t310 - t181;
t40 = -pkin(4) * t259 - t49;
t34 = mrSges(7,2) * t239 + mrSges(7,3) * t47;
t33 = -mrSges(7,1) * t239 - mrSges(7,3) * t46;
t19 = (-t258 * t343 + t310 * t416) * pkin(10) + t29;
t15 = t350 * t259 + (t258 * t347 + t310 * t417) * pkin(10) - t49;
t10 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t5 = -qJD(6) * t43 + t15 * t346 - t19 * t342;
t4 = qJD(6) * t42 + t15 * t342 + t19 * t346;
t1 = [(t238 * t512 - t239 * t270) * mrSges(4,3) + ((Ifges(4,2) + Ifges(7,3) / 0.2e1 + t406) * t239 + t405 * t164 + t407 * t163 + Ifges(4,4) * t238 + t354 + t161 / 0.2e1 + t160 / 0.2e1 - t159 / 0.2e1 + t158 / 0.2e1 + t233 / 0.2e1 + t234 / 0.2e1 + t236 / 0.2e1 - t153 * mrSges(4,3) - t528) * t308 + (-Ifges(4,1) * t238 - Ifges(4,4) * t239 + t58 * t471 + t59 * t472 + t37 * t386 + t375 * t490 + t381 * t491 + (mrSges(4,3) + t388) * t154 + (-t27 * t343 - t28 * t347) * mrSges(5,3) + (-t17 * t343 + t23 * t347) * mrSges(6,2) + (t148 * t470 + t103 * t387 + t218 * t389 + t374 * t483 + t380 * t482 + (t113 * t343 - t114 * t347) * mrSges(5,3) + (-t343 * t89 - t347 * t90) * mrSges(6,2) + t538 * t481 + t541 * t477 + t514 * t472) * qJD(4) + t539 * t492 + t540 * t484 + (qJD(4) * t145 + t548) * t469) * t310 + t334 * (mrSges(4,1) * t239 - mrSges(4,2) * t238) - t512 * t74 - (t530 - t248 * mrSges(4,3) + (-t343 * t90 + t347 * t89) * mrSges(6,2) + t215 / 0.2e1 + t289 / 0.2e1 + t460 / 0.2e1 + t514 * t469 + t551) * t258 + (-t214 / 0.2e1 - t249 * mrSges(4,3) - t447 / 0.2e1 - t446 / 0.2e1 - t442 / 0.2e1 + t147 / 0.2e1 + t146 / 0.2e1 - t458 / 0.2e1 - t450 / 0.2e1 - t70 / 0.2e1 + t405 * t265 + t406 * t291 + t407 * t266 - t553) * t259 + m(5) * (t113 * t49 + t114 * t48 + t166 * t28 + t167 * t27 + t181 * t218 - t432) + m(4) * (t153 * t270 + t180 * t249 - t181 * t248 - t432) + m(7) * (t138 * t18 + t2 * t43 + t20 * t5 + t21 * t4 + t3 * t42 + t41 * t82) + m(6) * (t103 * t57 + t126 * t17 + t131 * t23 + t177 * t37 + t29 * t90 + t40 * t89) + t29 * t193 + t48 * t194 + t49 * t195 + t40 * t196 + t57 * t185 + (-t2 * t223 - t20 * t68 + t21 * t69 - t224 * t3) * mrSges(7,3) + t18 * (mrSges(7,1) * t223 + mrSges(7,2) * t224) + (Ifges(7,5) * t224 - Ifges(7,6) * t223) * t485 + (Ifges(7,4) * t224 - Ifges(7,2) * t223) * t501 + (Ifges(7,1) * t224 - Ifges(7,4) * t223) * t502 + t177 * t73 + t166 * t93 + t167 * t95 + t138 * t10 + t131 * t94 + t4 * t124 + t5 * t125 + t126 * t92 + t41 * t81 + t82 * (-mrSges(7,1) * t69 + mrSges(7,2) * t68) + t42 * t33 + t43 * t34 + (-pkin(7) * t322 - t292 / 0.2e1 + t400 + (-0.2e1 * t494 - 0.3e1 / 0.2e1 * t459 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t349) * qJD(1) + (t243 + 0.2e1 * t468 + qJD(1) * (mrSges(4,1) * t308 + mrSges(4,2) * t310)) * pkin(2)) * t418 + (Ifges(7,5) * t68 + Ifges(7,6) * t69) * t478 + (-pkin(7) * t321 + t293 / 0.2e1 + t401 + (-0.2e1 * t493 + 0.3e1 / 0.2e1 * Ifges(3,4) * t349) * qJD(1)) * t349 * qJD(2) + t180 * t273 + (Ifges(7,1) * t68 + Ifges(7,4) * t69) * t486 + (Ifges(7,4) * t68 + Ifges(7,2) * t69) * t488 + t68 * t497 + t69 * t499 + t224 * t503 - t223 * t504 + t424 * t181; ((t238 * t348 - t239 * t344) * mrSges(4,3) + (t424 * t344 + (t343 * t422 + t347 * t423 + t273) * t348) * qJD(3)) * pkin(2) + (-t156 + t276) * t185 + (t18 * t281 + t2 * t241 + t20 * t522 + t21 * t523 + t240 * t3 + t515 * t82) * m(7) + t523 * t124 + t249 * t440 + t506 * t332 - t517 * t333 + t522 * t125 - m(6) * (t103 * t156 + t89 * t98 + t90 * t97) - m(5) * (t113 * t136 + t114 * t137 + t218 * t251) + m(6) * (t103 * t276 + t302 * t37) + t515 * t81 - t97 * t193 - t137 * t194 - t136 * t195 - t98 * t196 - m(4) * (-t248 * t251 + t249 * t252) + 0.2e1 * ((t153 * t344 - t154 * t348) * t505 + (m(5) * (-t113 * t426 + t114 * t425 + t218 * t344) / 0.2e1 + m(6) * (t425 * t90 + t426 * t89) / 0.2e1 + (-t248 * t344 + t249 * t348) * t505) * qJD(3)) * pkin(2) + ((t401 - t335 / 0.2e1 - t293 / 0.2e1 + qJD(1) * t493 + (t321 - t434) * pkin(7)) * t349 + (t400 + t292 / 0.2e1 + (t494 + t459 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t349) * qJD(1) + (t322 + t433) * pkin(7) + (-t243 - t468) * pkin(2)) * t345) * qJD(1) + t351 + t240 * t33 + t241 * t34 - t252 * t273 + t281 * t10 + t302 * t73 - t424 * t251; t506 * pkin(9) + t520 * t125 - m(5) * (t113 * t143 + t114 * t144 + t218 * t249) + t516 * t81 + t521 * t124 - t106 * t193 - t144 * t194 - t143 * t195 - t108 * t196 + t513 * t185 + (-t424 + t440) * t249 + t351 + t267 * t33 + t269 * t34 - t248 * t273 + t303 * t10 - t411 * t73 + t517 * pkin(3) + (t18 * t303 + t2 * t269 + t20 * t520 + t21 * t521 + t267 * t3 + t516 * t82) * m(7) + (t103 * t513 - t106 * t90 - t108 * t89 - t37 * t411) * m(6); -t550 + (-t549 * t265 - t266 * t562) * t477 + (-t265 * t563 + t145 + t263 - t457) * t481 - t174 * t498 + t554 * t489 + t354 + (-pkin(4) * t23 + qJ(5) * t17 - t103 * t184 - t114 * t89 + t509 * t90) * m(6) + (t265 * t89 + t266 * t90) * mrSges(6,2) + t518 * t124 + (-t121 * t82 + t2 * t319 + t20 * t519 + t21 * t518 + t3 * t318) * m(7) + t519 * t125 + t161 + t160 - t159 + t158 + (-Ifges(5,2) * t266 - t264 + t514) * t482 + t369 * t500 + qJD(5) * t193 - t184 * t185 - t121 * t81 - pkin(4) * t94 + qJ(5) * t92 + t148 * t480 - t103 * (mrSges(6,1) * t266 + mrSges(6,3) * t265) - t218 * (mrSges(5,1) * t266 - mrSges(5,2) * t265) + t234 + t236 + (Ifges(6,3) * t266 - t453) * t483 + t318 * t33 + t319 * t34 + (-t422 + t461) * t114 + (-t423 - t462) * t113; t346 * t33 + t342 * t34 + (t185 - t81) * t266 + t370 * qJD(6) + (-t193 - t370) * t291 + t94 + (t2 * t342 - t266 * t82 + t3 * t346 + t285 * (-t20 * t342 + t21 * t346)) * m(7) + (t103 * t266 - t291 * t90 + t23) * m(6); t71 * t486 - t20 * t124 + t21 * t125 + t390 + (t72 - t554) * t489 + t550;];
tauc  = t1(:);

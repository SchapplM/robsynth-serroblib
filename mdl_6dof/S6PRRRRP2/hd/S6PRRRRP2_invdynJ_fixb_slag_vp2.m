% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:11
% EndTime: 2019-03-09 00:01:50
% DurationCPUTime: 22.46s
% Computational Cost: add. (10429->757), mult. (23266->999), div. (0->0), fcn. (17560->14), ass. (0->354)
t343 = sin(qJ(2));
t338 = sin(pkin(6));
t462 = qJD(1) * t338;
t423 = t343 * t462;
t342 = sin(qJ(3));
t456 = qJD(3) * t342;
t554 = pkin(3) * t456 - t423;
t599 = Ifges(6,1) + Ifges(7,1);
t606 = -Ifges(6,4) + Ifges(7,5);
t579 = Ifges(7,4) + Ifges(6,5);
t578 = Ifges(7,2) + Ifges(6,3);
t598 = Ifges(6,6) - Ifges(7,6);
t341 = sin(qJ(4));
t345 = cos(qJ(4));
t346 = cos(qJ(3));
t289 = t341 * t342 - t345 * t346;
t357 = t289 * qJD(4);
t216 = -qJD(3) * t289 - t357;
t290 = t341 * t346 + t342 * t345;
t358 = t290 * qJD(4);
t217 = qJD(3) * t290 + t358;
t605 = pkin(4) * t217 - pkin(10) * t216 + t554;
t348 = -pkin(9) - pkin(8);
t424 = qJD(3) * t348;
t294 = t342 * t424;
t295 = t346 * t424;
t347 = cos(qJ(2));
t422 = t347 * t462;
t312 = t348 * t342;
t313 = t348 * t346;
t558 = t345 * t312 + t313 * t341;
t561 = qJD(4) * t558 + t289 * t422 + t294 * t345 + t295 * t341;
t449 = qJD(2) * qJD(3);
t296 = qJDD(2) * t346 - t342 * t449;
t297 = qJDD(2) * t342 + t346 * t449;
t144 = -qJD(2) * t357 + t296 * t341 + t297 * t345;
t335 = qJDD(3) + qJDD(4);
t340 = sin(qJ(5));
t344 = cos(qJ(5));
t286 = t290 * qJD(2);
t446 = qJD(3) + qJD(4);
t228 = t340 * t286 - t344 * t446;
t452 = qJD(5) * t228;
t75 = t344 * t144 + t340 * t335 - t452;
t544 = t75 / 0.2e1;
t229 = t344 * t286 + t340 * t446;
t76 = qJD(5) * t229 + t340 * t144 - t344 * t335;
t542 = t76 / 0.2e1;
t145 = -qJD(2) * t358 + t296 * t345 - t297 * t341;
t143 = qJDD(5) - t145;
t540 = t143 / 0.2e1;
t604 = -mrSges(6,1) - mrSges(7,1);
t603 = mrSges(6,2) - mrSges(7,3);
t581 = mrSges(6,3) + mrSges(7,2);
t129 = mrSges(7,1) * t228 - mrSges(7,3) * t229;
t503 = t286 * mrSges(5,3);
t559 = mrSges(5,1) * t446 - mrSges(6,1) * t228 - mrSges(6,2) * t229 - t503;
t602 = t129 - t559;
t285 = t289 * qJD(2);
t276 = qJD(5) + t285;
t300 = qJD(2) * pkin(8) + t423;
t399 = pkin(9) * qJD(2) + t300;
t339 = cos(pkin(6));
t461 = qJD(1) * t339;
t421 = t342 * t461;
t222 = t346 * t399 + t421;
t212 = t345 * t222;
t322 = t346 * t461;
t221 = -t342 * t399 + t322;
t213 = qJD(3) * pkin(3) + t221;
t109 = t341 * t213 + t212;
t103 = pkin(10) * t446 + t109;
t329 = pkin(3) * t346 + pkin(2);
t274 = -qJD(2) * t329 - t422;
t140 = pkin(4) * t285 - pkin(10) * t286 + t274;
t44 = -t103 * t340 + t140 * t344;
t565 = qJD(6) - t44;
t35 = -pkin(5) * t276 + t565;
t241 = t300 * t346 + t421;
t460 = qJD(2) * t338;
t412 = qJD(1) * t460;
t315 = t347 * t412;
t448 = qJDD(1) * t338;
t269 = t343 * t448 + t315;
t257 = qJDD(2) * pkin(8) + t269;
t447 = qJDD(1) * t339;
t121 = -qJD(3) * t241 - t257 * t342 + t346 * t447;
t100 = qJDD(3) * pkin(3) - pkin(9) * t297 + t121;
t120 = qJD(3) * t322 + t346 * t257 - t300 * t456 + t342 * t447;
t104 = pkin(9) * t296 + t120;
t453 = qJD(4) * t345;
t454 = qJD(4) * t341;
t25 = t341 * t100 + t345 * t104 + t213 * t453 - t222 * t454;
t22 = pkin(10) * t335 + t25;
t314 = t343 * t412;
t268 = t347 * t448 - t314;
t256 = -qJDD(2) * pkin(2) - t268;
t201 = -pkin(3) * t296 + t256;
t43 = -pkin(4) * t145 - pkin(10) * t144 + t201;
t450 = qJD(5) * t344;
t451 = qJD(5) * t340;
t6 = -t103 * t451 + t140 * t450 + t344 * t22 + t340 * t43;
t2 = qJ(6) * t143 + qJD(6) * t276 + t6;
t45 = t103 * t344 + t140 * t340;
t7 = -qJD(5) * t45 - t22 * t340 + t344 * t43;
t4 = -pkin(5) * t143 + qJDD(6) - t7;
t389 = t2 * t344 + t340 * t4;
t601 = t35 * t450 + t389;
t600 = mrSges(5,2) - t581;
t597 = t579 * t143 + t599 * t75 + t606 * t76;
t596 = -t228 * t598 + t229 * t579 + t578 * t276;
t224 = Ifges(6,4) * t228;
t506 = t228 * Ifges(7,5);
t576 = t229 * t599 + t579 * t276 - t224 + t506;
t595 = t446 * Ifges(5,5);
t594 = t446 * Ifges(5,6);
t202 = pkin(4) * t289 - pkin(10) * t290 - t329;
t231 = t312 * t341 - t313 * t345;
t562 = t340 * t202 + t344 * t231;
t572 = -qJD(5) * t562 - t340 * t561 + t344 * t605;
t570 = t202 * t450 - t231 * t451 + t340 * t605 + t344 * t561;
t112 = t221 * t341 + t212;
t593 = -pkin(3) * t454 + t112;
t486 = t285 * t344;
t487 = t285 * t340;
t592 = -qJD(6) * t340 + (-t450 - t486) * qJ(6) + (t451 + t487) * pkin(5);
t591 = -t340 * t598 + t344 * t579;
t508 = Ifges(7,5) * t340;
t510 = Ifges(6,4) * t340;
t590 = -t344 * t599 - t508 + t510;
t475 = t338 * t343;
t281 = t339 * t346 - t342 * t475;
t501 = cos(pkin(11));
t401 = t501 * t347;
t337 = sin(pkin(11));
t477 = t337 * t343;
t280 = -t339 * t477 + t401;
t474 = t338 * t346;
t589 = -t280 * t342 + t337 * t474;
t388 = -t7 * t340 + t344 * t6;
t211 = t341 * t222;
t108 = t345 * t213 - t211;
t102 = -pkin(4) * t446 - t108;
t384 = mrSges(7,1) * t340 - mrSges(7,3) * t344;
t385 = mrSges(6,1) * t340 + mrSges(6,2) * t344;
t48 = t228 * pkin(5) - t229 * qJ(6) + t102;
t588 = t102 * t385 + t384 * t48;
t587 = Ifges(7,5) * t544 + Ifges(7,6) * t540 - t75 * Ifges(6,4) / 0.2e1 - t143 * Ifges(6,6) / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t542;
t38 = mrSges(6,1) * t143 - mrSges(6,3) * t75;
t39 = -t143 * mrSges(7,1) + t75 * mrSges(7,2);
t516 = t39 - t38;
t37 = -mrSges(7,2) * t76 + mrSges(7,3) * t143;
t40 = -mrSges(6,2) * t143 - mrSges(6,3) * t76;
t517 = t37 + t40;
t586 = t516 * t340 + t517 * t344;
t505 = t229 * Ifges(6,4);
t95 = -t228 * Ifges(6,2) + t276 * Ifges(6,6) + t505;
t585 = -t95 / 0.2e1;
t584 = -m(6) - m(7);
t583 = t296 / 0.2e1;
t582 = t297 / 0.2e1;
t575 = t346 * Ifges(4,2);
t573 = -pkin(5) * t217 - t572;
t571 = qJ(6) * t217 + qJD(6) * t289 + t570;
t569 = t592 - t593;
t568 = -t109 + t592;
t131 = mrSges(5,1) * t335 - mrSges(5,3) * t144;
t31 = mrSges(6,1) * t76 + mrSges(6,2) * t75;
t567 = t31 - t131;
t402 = t501 * t343;
t476 = t337 * t347;
t278 = t339 * t402 + t476;
t336 = qJ(3) + qJ(4);
t333 = sin(t336);
t334 = cos(t336);
t403 = t338 * t501;
t203 = -t278 * t333 - t334 * t403;
t497 = t203 * t344;
t498 = t203 * t340;
t564 = -pkin(5) * t497 - qJ(6) * t498;
t478 = t337 * t338;
t205 = -t280 * t333 + t334 * t478;
t495 = t205 * t344;
t496 = t205 * t340;
t563 = -pkin(5) * t495 - qJ(6) * t496;
t258 = -t333 * t475 + t334 * t339;
t491 = t258 * t344;
t492 = t258 * t340;
t560 = -pkin(5) * t491 - qJ(6) * t492;
t527 = pkin(3) * t341;
t327 = pkin(10) + t527;
t440 = pkin(3) * t453;
t557 = -t327 * t451 + t344 * t440;
t556 = -t327 * t450 - t340 * t440;
t363 = t216 * t340 + t290 * t450;
t555 = -t44 * t450 - t45 * t451;
t553 = t143 * t578 + t579 * t75 - t598 * t76;
t552 = t120 * t346 - t121 * t342;
t439 = m(4) * pkin(8) + mrSges(4,3);
t550 = -mrSges(5,3) - t439 + mrSges(3,2);
t387 = -mrSges(4,1) * t346 + mrSges(4,2) * t342;
t356 = m(4) * pkin(2) - t387;
t386 = mrSges(5,1) * t334 - mrSges(5,2) * t333;
t549 = t386 + t356 + mrSges(3,1);
t259 = t333 * t339 + t334 * t475;
t548 = -t258 * mrSges(5,1) + t259 * t600 + t491 * t604 + t492 * t603;
t206 = t280 * t334 + t333 * t478;
t547 = -t205 * mrSges(5,1) + t206 * t600 + t495 * t604 + t496 * t603;
t204 = t278 * t334 - t333 * t403;
t546 = -t203 * mrSges(5,1) + t204 * t600 + t497 * t604 + t498 * t603;
t398 = m(7) * pkin(5) - t604;
t393 = -m(7) * qJ(6) + t603;
t545 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t349 = qJD(2) ^ 2;
t543 = -t76 / 0.2e1;
t538 = -t228 / 0.2e1;
t537 = t228 / 0.2e1;
t536 = -t229 / 0.2e1;
t535 = t229 / 0.2e1;
t534 = -t276 / 0.2e1;
t532 = t285 / 0.2e1;
t530 = t286 / 0.2e1;
t528 = t340 / 0.2e1;
t526 = pkin(3) * t345;
t525 = pkin(4) * t334;
t524 = pkin(5) * t286;
t523 = g(3) * t338;
t515 = mrSges(6,3) * t228;
t514 = mrSges(6,3) * t229;
t513 = Ifges(4,4) * t342;
t512 = Ifges(4,4) * t346;
t511 = Ifges(5,4) * t286;
t509 = Ifges(6,4) * t344;
t507 = Ifges(7,5) * t344;
t504 = t285 * mrSges(5,3);
t36 = qJ(6) * t276 + t45;
t502 = t340 * t36;
t493 = t216 * t344;
t277 = -t339 * t401 + t477;
t490 = t277 * t333;
t279 = t339 * t476 + t402;
t489 = t279 * t333;
t480 = t334 * t340;
t479 = t334 * t344;
t473 = t338 * t347;
t471 = t344 * t347;
t200 = pkin(4) * t286 + pkin(10) * t285;
t54 = t344 * t108 + t340 * t200;
t113 = t221 * t345 - t211;
t459 = qJD(2) * t342;
t443 = pkin(3) * t459;
t162 = t200 + t443;
t52 = t344 * t113 + t340 * t162;
t146 = -mrSges(7,2) * t228 + mrSges(7,3) * t276;
t147 = -mrSges(6,2) * t276 - t515;
t469 = t146 + t147;
t148 = mrSges(6,1) * t276 - t514;
t149 = -mrSges(7,1) * t276 + mrSges(7,2) * t229;
t468 = -t148 + t149;
t465 = -t277 * t329 - t278 * t348;
t464 = -t279 * t329 - t280 * t348;
t458 = qJD(2) * t343;
t457 = qJD(2) * t346;
t455 = qJD(3) * t346;
t438 = mrSges(4,3) * t459;
t437 = mrSges(4,3) * t457;
t436 = t36 * t451;
t432 = t333 * t473;
t429 = t338 * t471;
t318 = t340 * t473;
t223 = Ifges(7,5) * t229;
t92 = t276 * Ifges(7,6) + t228 * Ifges(7,3) + t223;
t428 = t92 * t528;
t420 = t338 * t458;
t419 = t347 * t460;
t408 = -t451 / 0.2e1;
t407 = t450 / 0.2e1;
t406 = t203 * pkin(4) + t204 * pkin(10);
t405 = t205 * pkin(4) + pkin(10) * t206;
t404 = t258 * pkin(4) + pkin(10) * t259;
t400 = t449 / 0.2e1;
t391 = t589 * pkin(3);
t381 = t513 + t575;
t380 = -Ifges(6,2) * t340 + t509;
t378 = Ifges(4,5) * t346 - Ifges(4,6) * t342;
t376 = Ifges(7,3) * t340 + t507;
t375 = pkin(5) * t344 + qJ(6) * t340;
t374 = pkin(5) * t340 - qJ(6) * t344;
t372 = -t340 * t45 - t344 * t44;
t53 = -t108 * t340 + t200 * t344;
t51 = -t113 * t340 + t162 * t344;
t110 = t202 * t344 - t231 * t340;
t282 = t339 * t342 + t343 * t474;
t369 = t345 * t281 - t282 * t341;
t168 = t281 * t341 + t282 * t345;
t367 = t281 * pkin(3);
t26 = t100 * t345 - t341 * t104 - t213 * t454 - t222 * t453;
t303 = -pkin(4) - t375;
t136 = t168 * t340 + t429;
t362 = t290 * t451 - t493;
t301 = -qJD(2) * pkin(2) - t422;
t361 = t301 * (mrSges(4,1) * t342 + mrSges(4,2) * t346);
t360 = t342 * (Ifges(4,1) * t346 - t513);
t355 = -t278 * t342 - t346 * t403;
t354 = t391 + t405;
t23 = -pkin(4) * t335 - t26;
t353 = t367 + t404;
t352 = t355 * pkin(3);
t119 = qJD(4) * t231 + t294 * t341 - t345 * t295;
t351 = t352 + t406;
t163 = -Ifges(5,2) * t285 + t511 + t594;
t273 = Ifges(5,4) * t285;
t164 = Ifges(5,1) * t286 - t273 + t595;
t9 = pkin(5) * t76 - qJ(6) * t75 - qJD(6) * t229 + t23;
t350 = (-t229 * t590 + t276 * t591) * qJD(5) / 0.2e1 + t508 * t542 + t510 * t543 - (-Ifges(5,1) * t285 - t511 + t596) * t286 / 0.2e1 + t597 * t528 + (-t23 * mrSges(6,1) - t9 * mrSges(7,1) + Ifges(6,2) * t543 - Ifges(7,3) * t542 + t540 * t598 - t587) * t344 + (t45 * mrSges(6,2) - t36 * mrSges(7,3) + t594 / 0.2e1 + Ifges(6,6) * t537 - t274 * mrSges(5,1) - t44 * mrSges(6,1) + t35 * mrSges(7,1) + Ifges(7,6) * t538 - Ifges(5,2) * t532 + t579 * t536 + t578 * t534) * t286 + (t23 * mrSges(6,2) - t9 * mrSges(7,3) + t579 * t540 + t544 * t599) * t340 + (t428 + t588) * qJD(5) + (t595 / 0.2e1 - t380 * t537 + t274 * mrSges(5,2) - t376 * t538 + t590 * t536 - t591 * t534 + t588) * t285 + (-t380 / 0.2e1 + t376 / 0.2e1) * t452 + (t585 + t92 / 0.2e1) * t487 + (t35 * t486 - t36 * t487 + t601) * mrSges(7,2) + (-t507 + t509) * t544 + t163 * t530 + (-t273 + t164) * t532 - t108 * t504 + (t407 + t486 / 0.2e1) * t576 + Ifges(5,3) * t335 + Ifges(5,5) * t144 + Ifges(5,6) * t145 + (-t44 * t486 - t45 * t487 + t388) * mrSges(6,3) + t95 * t408 - t25 * mrSges(5,2) + t26 * mrSges(5,1);
t330 = Ifges(4,4) * t457;
t328 = -pkin(4) - t526;
t309 = -qJD(3) * mrSges(4,2) + t437;
t308 = qJD(3) * mrSges(4,1) - t438;
t293 = t329 * t473;
t292 = t387 * qJD(2);
t288 = t303 - t526;
t284 = Ifges(4,1) * t459 + Ifges(4,5) * qJD(3) + t330;
t283 = Ifges(4,6) * qJD(3) + qJD(2) * t381;
t272 = t286 * qJ(6);
t271 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t297;
t270 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t296;
t242 = -mrSges(5,2) * t446 - t504;
t240 = -t300 * t342 + t322;
t225 = -mrSges(4,1) * t296 + mrSges(4,2) * t297;
t219 = qJD(3) * t281 + t346 * t419;
t218 = -qJD(3) * t282 - t342 * t419;
t208 = t259 * t340 + t429;
t197 = mrSges(5,1) * t285 + mrSges(5,2) * t286;
t137 = t168 * t344 - t318;
t132 = -mrSges(5,2) * t335 + mrSges(5,3) * t145;
t128 = pkin(5) * t229 + qJ(6) * t228;
t126 = t206 * t340 - t279 * t344;
t124 = t204 * t340 - t277 * t344;
t117 = t290 * t374 - t558;
t88 = -pkin(5) * t289 - t110;
t83 = qJ(6) * t289 + t562;
t59 = -mrSges(5,1) * t145 + mrSges(5,2) * t144;
t57 = qJD(4) * t168 - t345 * t218 + t219 * t341;
t56 = qJD(4) * t369 + t218 * t341 + t219 * t345;
t50 = -t53 - t524;
t49 = t272 + t54;
t47 = -t51 - t524;
t46 = t272 + t52;
t33 = -qJD(5) * t318 + t168 * t450 + t340 * t56 - t344 * t420;
t32 = -qJD(5) * t136 + t340 * t420 + t344 * t56;
t30 = mrSges(7,1) * t76 - mrSges(7,3) * t75;
t29 = t374 * t216 + (qJD(5) * t375 - qJD(6) * t344) * t290 + t119;
t1 = [m(2) * qJDD(1) + t168 * t132 + t218 * t308 + t219 * t309 + t56 * t242 + t282 * t270 + t281 * t271 + t468 * t33 + t469 * t32 + t517 * t137 + t516 * t136 + t602 * t57 - (t30 + t567) * t369 + (-m(2) - m(3) - m(4) - m(5) + t584) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t349 - t225 - t59) * t347 + (-mrSges(3,1) * t349 - mrSges(3,2) * qJDD(2) + (t197 + t292) * qJD(2)) * t343) * t338 + m(5) * (-t108 * t57 + t109 * t56 + t369 * t26 + t168 * t25 + (-t201 * t347 + t274 * t458) * t338) + m(4) * (t120 * t282 + t121 * t281 + t218 * t240 + t219 * t241 + (-t256 * t347 + t301 * t458) * t338) + m(3) * (qJDD(1) * t339 ^ 2 + (t268 * t347 + t269 * t343) * t338) + m(7) * (t136 * t4 + t137 * t2 + t32 * t36 + t33 * t35 - t369 * t9 + t48 * t57) + m(6) * (t102 * t57 - t136 * t7 + t137 * t6 - t23 * t369 + t32 * t45 - t33 * t44); (t361 + t378 * qJD(3) / 0.2e1) * qJD(3) + (t102 * t119 + t110 * t7 - t23 * t558 + t44 * t572 + t45 * t570 + t562 * t6) * m(6) + (-t108 * t119 + t109 * t561 - t201 * t329 + t231 * t25 + t26 * t558 + t274 * t554) * m(5) - t567 * t558 + (t217 * t578 - t362 * t579 - t363 * t598) * t276 / 0.2e1 + t512 * t582 + t381 * t583 + t363 * t585 + (-t108 * t216 - t109 * t217) * mrSges(5,3) + t596 * t217 / 0.2e1 + (t343 * t523 - t269 + t315) * mrSges(3,2) + (-(t301 * t343 + (-t240 * t342 + t241 * t346) * t347) * t462 - pkin(2) * t256) * m(4) - (t343 * t439 + t347 * t356) * t523 + t360 * t400 + (-Ifges(6,4) * t362 - Ifges(6,2) * t363 + Ifges(6,6) * t217) * t538 + t284 * t455 / 0.2e1 - t283 * t456 / 0.2e1 + ((mrSges(7,2) * t4 - mrSges(6,3) * t7 + t597 / 0.2e1) * t344 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t587) * t340 + t576 * t408 + (m(5) * t108 - m(6) * t102 - m(7) * t48 - t602) * t422 + t201 * mrSges(5,2) - t26 * mrSges(5,3) + Ifges(5,1) * t144 + Ifges(5,4) * t145 + Ifges(5,5) * t335 + t23 * t385 + t376 * t542 + t380 * t543 + t384 * t9 + t407 * t92 + t591 * t540 - t590 * t544) * t290 + t446 * (Ifges(5,5) * t216 - Ifges(5,6) * t217) / 0.2e1 + (t579 * t217 - t599 * t362 + t363 * t606) * t535 + (Ifges(5,1) * t216 - Ifges(5,4) * t217) * t530 + (-Ifges(7,5) * t362 + Ifges(7,6) * t217 + Ifges(7,3) * t363) * t537 + (-t347 * t523 + t268 + t314) * mrSges(3,1) + t44 * (mrSges(6,1) * t217 + mrSges(6,3) * t362) + t35 * (-mrSges(7,1) * t217 - mrSges(7,2) * t362) + t48 * (mrSges(7,1) * t363 + mrSges(7,3) * t362) + t102 * (mrSges(6,1) * t363 - mrSges(6,2) * t362) + t36 * (-mrSges(7,2) * t363 + mrSges(7,3) * t217) + t45 * (-mrSges(6,2) * t217 - mrSges(6,3) * t363) + t562 * t40 - t292 * t423 + Ifges(3,3) * qJDD(2) + qJDD(3) * (Ifges(4,5) * t342 + Ifges(4,6) * t346) - t329 * t59 - t285 * (Ifges(5,4) * t216 - Ifges(5,2) * t217) / 0.2e1 + t274 * (mrSges(5,1) * t217 + mrSges(5,2) * t216) + t231 * t132 - pkin(2) * t225 + (-m(5) * t293 - t581 * t432 + t584 * (pkin(10) * t432 - t348 * t475 + t473 * t525 + t293) + t393 * (t318 * t334 - t344 * t475) + (-t386 * t347 - (-m(5) * t348 + mrSges(5,3)) * t343 - t398 * (t334 * t471 + t340 * t343)) * t338) * g(3) - t217 * t163 / 0.2e1 + t216 * t164 / 0.2e1 + t29 * t129 + t117 * t30 + (t553 / 0.2e1 + t201 * mrSges(5,1) - t25 * mrSges(5,3) - Ifges(5,4) * t144 - Ifges(5,2) * t145 - Ifges(5,6) * t335 + Ifges(6,6) * t543 + Ifges(7,6) * t542 + t540 * t578 + t544 * t579 + t545) * t289 + t110 * t38 + t88 * t39 + t83 * t37 + t256 * t387 + (m(4) * ((-t240 * t346 - t241 * t342) * qJD(3) + t552) - t308 * t455 - t309 * t456 - t342 * t271 + t346 * t270) * pkin(8) + (-t240 * t455 - t241 * t456 + t552) * mrSges(4,3) + t554 * t197 - t559 * t119 + t561 * t242 + t216 * t428 + t570 * t147 + t571 * t146 + t572 * t148 + t573 * t149 + (t117 * t9 + t2 * t83 + t29 * t48 + t35 * t573 + t36 * t571 + t4 * t88) * m(7) + t576 * t493 / 0.2e1 + (Ifges(4,4) * t582 + Ifges(4,2) * t583 - t309 * t422 + t512 * t400) * t346 + (Ifges(4,1) * t297 + Ifges(4,4) * t583 + t308 * t422 - t400 * t575) * t342 + (-m(5) * t464 + t581 * t489 + t584 * (-pkin(10) * t489 - t279 * t525 + t464) - t398 * (-t279 * t479 + t280 * t340) + t393 * (-t279 * t480 - t280 * t344) + t550 * t280 + t549 * t279) * g(1) + (-m(5) * t465 + t581 * t490 + t584 * (-pkin(10) * t490 - t277 * t525 + t465) - t398 * (-t277 * t479 + t278 * t340) + t393 * (-t277 * t480 - t278 * t344) + t550 * t278 + t549 * t277) * g(2); (t440 - t113) * t242 + (t438 + t308) * t241 + (t437 - t309) * t240 - (-Ifges(4,2) * t459 + t284 + t330) * t457 / 0.2e1 + ((qJD(5) * t372 + t388) * m(6) + ((t344 * t35 - t502) * qJD(5) + t389) * m(7) + t586) * t327 - t349 * t360 / 0.2e1 + ((t25 * t341 + t26 * t345 + (-t108 * t341 + t109 * t345) * qJD(4)) * pkin(3) + t108 * t112 - t109 * t113 - t274 * t443) * m(5) + (-t102 * t112 - t44 * t51 - t45 * t52 + t23 * t328 + (t102 * t341 + (-t340 * t44 + t344 * t45) * t345) * qJD(4) * pkin(3)) * m(6) + (t288 * t9 + (t340 * t35 + t344 * t36) * t440 - t35 * t47 - t36 * t46 + t569 * t48) * m(7) + t283 * t459 / 0.2e1 - t378 * t449 / 0.2e1 + t131 * t526 + t132 * t527 - t197 * t443 - mrSges(7,2) * t436 + t350 - qJD(2) * t361 + t109 * t503 + Ifges(4,3) * qJDD(3) + t328 * t31 + Ifges(4,6) * t296 + Ifges(4,5) * t297 + t288 * t30 - t120 * mrSges(4,2) + t121 * mrSges(4,1) + (-m(5) * t391 - t589 * mrSges(4,1) - (-t280 * t346 - t342 * t478) * mrSges(4,2) - m(6) * t354 - m(7) * (t354 - t563) + t547) * g(1) + t593 * t559 + t555 * mrSges(6,3) + (-t47 - t556) * t149 + (-t51 + t556) * t148 + (-t52 + t557) * t147 + (-t46 + t557) * t146 + (-m(5) * t367 - m(6) * t353 - m(7) * (t353 - t560) - mrSges(4,1) * t281 + mrSges(4,2) * t282 + t548) * g(3) + (-m(5) * t352 - m(6) * t351 - m(7) * (t351 - t564) - t355 * mrSges(4,1) - (-t278 * t346 + t342 * t403) * mrSges(4,2) + t546) * g(2) + t569 * t129; t546 * g(2) + t547 * g(1) + t548 * g(3) + t568 * t129 + (-mrSges(7,2) * t502 + t372 * mrSges(6,3)) * qJD(5) + t350 + t303 * t30 + (t559 + t503) * t109 - t108 * t242 - t49 * t146 - t54 * t147 - t53 * t148 - t50 * t149 - pkin(4) * t31 + ((-t406 + t564) * g(2) + (-t405 + t563) * g(1) + (-t404 + t560) * g(3) + t303 * t9 - t35 * t50 - t36 * t49 + t568 * t48) * m(7) + (-pkin(4) * t23 - g(1) * t405 - g(2) * t406 - g(3) * t404 - t102 * t109 - t44 * t53 - t45 * t54) * m(6) + ((-t340 * t469 + t344 * t468) * qJD(5) + m(7) * (-t436 + t601) + m(6) * (t388 + t555) + t586) * pkin(10); t545 + (-t228 * t599 + t223 - t505 + t92) * t536 + (-pkin(5) * t4 + qJ(6) * t2 - t128 * t48 - t35 * t45 + t36 * t565) * m(7) + (-t228 * t579 - t229 * t598) * t534 + (-t468 + t514) * t45 + (-t469 - t515) * t44 + (Ifges(7,3) * t229 - t506) * t538 + t95 * t535 + (-Ifges(6,2) * t229 - t224 + t576) * t537 + (t393 * (t204 * t344 + t277 * t340) + t398 * t124) * g(2) + (t393 * (t206 * t344 + t279 * t340) + t398 * t126) * g(1) + (t393 * (t259 * t344 - t318) + t398 * t208) * g(3) + (t228 * t35 + t229 * t36) * mrSges(7,2) - t48 * (mrSges(7,1) * t229 + mrSges(7,3) * t228) - t102 * (mrSges(6,1) * t229 - mrSges(6,2) * t228) + qJD(6) * t146 - t128 * t129 + qJ(6) * t37 - pkin(5) * t39 + t553; t229 * t129 - t276 * t146 + (-g(1) * t126 - g(2) * t124 - g(3) * t208 + t48 * t229 - t36 * t276 + t4) * m(7) + t39;];
tau  = t1;

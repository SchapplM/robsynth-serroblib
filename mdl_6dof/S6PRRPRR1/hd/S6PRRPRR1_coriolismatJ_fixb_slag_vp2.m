% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:51:19
% EndTime: 2019-03-08 21:51:38
% DurationCPUTime: 12.47s
% Computational Cost: add. (24684->603), mult. (53749->829), div. (0->0), fcn. (64282->12), ass. (0->332)
t334 = sin(qJ(6));
t337 = cos(qJ(6));
t333 = sin(pkin(6));
t339 = cos(qJ(2));
t452 = t333 * t339;
t336 = sin(qJ(3));
t338 = cos(qJ(3));
t533 = sin(qJ(2));
t423 = t333 * t533;
t483 = cos(pkin(6));
t362 = t336 * t483 + t338 * t423;
t363 = t336 * t423 - t338 * t483;
t481 = sin(pkin(12));
t482 = cos(pkin(12));
t245 = t362 * t482 - t481 * t363;
t335 = sin(qJ(5));
t534 = cos(qJ(5));
t571 = -t362 * t481 - t482 * t363;
t619 = t245 * t534 + t335 * t571;
t127 = -t334 * t452 + t337 * t619;
t163 = t245 * t335 - t534 * t571;
t490 = t337 * mrSges(7,2);
t493 = t334 * mrSges(7,1);
t375 = t493 / 0.2e1 + t490 / 0.2e1;
t371 = t163 * t375;
t311 = t490 + t493;
t468 = t163 * t311;
t492 = t334 * mrSges(7,3);
t428 = t492 / 0.2e1;
t429 = -t492 / 0.2e1;
t605 = t428 + t429;
t631 = t371 + t468 / 0.2e1 + t605 * t127;
t479 = t127 * t337;
t126 = -t334 * t619 - t337 * t452;
t480 = t126 * t334;
t381 = t479 - t480;
t633 = m(7) * (t619 - t381) * t163;
t635 = t633 * qJD(1);
t637 = t631 * qJD(6) + t635;
t636 = qJD(3) + qJD(5);
t391 = mrSges(7,1) * t337 - mrSges(7,2) * t334;
t587 = -t391 - mrSges(6,1);
t406 = t482 * t338;
t302 = -t336 * t481 + t406;
t405 = t481 * t338;
t303 = -t336 * t482 - t405;
t373 = t335 * t302 - t303 * t534;
t596 = -t373 / 0.2e1;
t632 = 0.2e1 * t596;
t511 = t163 * mrSges(6,2);
t629 = t587 * t619;
t630 = t511 + t629;
t555 = -t163 / 0.2e1;
t628 = pkin(5) * t619;
t527 = -qJ(4) - pkin(8);
t392 = t527 * t481;
t304 = t336 * t392;
t279 = -t527 * t406 + t304;
t529 = t302 * pkin(9);
t239 = t279 + t529;
t393 = t527 * t482;
t378 = t336 * t393;
t277 = t405 * t527 + t378;
t528 = t303 * pkin(9);
t351 = t277 + t528;
t154 = t239 * t335 - t351 * t534;
t626 = t154 * t619;
t278 = t338 * t392 + t378;
t624 = -t277 + t278;
t535 = t337 / 0.2e1;
t538 = -t334 / 0.2e1;
t548 = t373 / 0.2e1;
t324 = Ifges(7,5) * t337;
t517 = Ifges(7,6) * t334;
t582 = t324 - t517;
t269 = t534 * t302 + t303 * t335;
t325 = Ifges(7,4) * t337;
t581 = -Ifges(7,2) * t334 + t325;
t602 = Ifges(7,6) * t373 + t269 * t581;
t521 = Ifges(7,4) * t334;
t317 = Ifges(7,1) * t337 - t521;
t603 = Ifges(7,5) * t373 + t269 * t317;
t617 = -t269 / 0.2e1;
t623 = Ifges(6,4) * t632 + t603 * t535 + t602 * t538 + t582 * t548 + (Ifges(6,2) + Ifges(7,3)) * t617;
t568 = -m(7) / 0.2e1;
t613 = t269 * mrSges(6,2);
t622 = -t613 / 0.2e1;
t611 = t311 * t269;
t621 = t154 * t611;
t276 = t338 * t393 - t304;
t620 = t276 + t279;
t199 = t373 * mrSges(6,1) + t613;
t616 = -t269 / 0.4e1;
t615 = t269 / 0.2e1;
t544 = t269 / 0.4e1;
t614 = -t311 / 0.2e1;
t612 = t269 * mrSges(6,3);
t610 = t333 ^ 2 * t533;
t586 = (-t517 / 0.2e1 + t324 / 0.2e1) * t269;
t186 = t311 * t373;
t588 = t373 * mrSges(6,3);
t584 = t186 + t588;
t422 = t482 * pkin(3);
t322 = t422 + pkin(4);
t424 = pkin(3) * t481;
t292 = t335 * t322 + t424 * t534;
t290 = pkin(10) + t292;
t329 = t334 ^ 2;
t331 = t337 ^ 2;
t438 = t329 + t331;
t400 = t438 * t290;
t607 = t292 - t400;
t606 = t336 ^ 2 + t338 ^ 2;
t291 = t322 * t534 - t335 * t424;
t289 = -pkin(5) - t291;
t604 = t289 * t619 + t381 * t291;
t316 = Ifges(7,1) * t334 + t325;
t440 = t337 * t316;
t314 = Ifges(7,2) * t337 + t521;
t446 = t334 * t314;
t583 = t440 / 0.2e1 - t446 / 0.2e1;
t531 = pkin(5) * t373;
t203 = -pkin(10) * t269 + t531;
t281 = t303 * t452;
t282 = t302 * t452;
t214 = t335 * t281 + t282 * t534;
t204 = -t334 * t214 + t337 * t423;
t205 = t337 * t214 + t334 * t423;
t601 = (t204 * t337 + t205 * t334) * t568 - (m(5) + m(6)) * t423 / 0.2e1;
t597 = t316 / 0.4e1;
t594 = mrSges(7,1) * t373;
t593 = mrSges(7,2) * t373;
t592 = mrSges(7,3) * (t331 / 0.2e1 + t329 / 0.2e1);
t312 = Ifges(7,5) * t334 + Ifges(7,6) * t337;
t460 = t373 * t312;
t296 = t302 * mrSges(5,2);
t580 = t303 * mrSges(5,1) - t199 - t296;
t458 = t373 * t334;
t193 = mrSges(7,2) * t269 - mrSges(7,3) * t458;
t457 = t373 * t337;
t196 = -mrSges(7,1) * t269 - mrSges(7,3) * t457;
t536 = -t337 / 0.2e1;
t579 = t193 * t538 + t196 * t536;
t494 = t331 * mrSges(7,3);
t495 = t329 * mrSges(7,3);
t578 = t495 / 0.2e1 + t494 / 0.2e1;
t81 = t154 * t334 + t203 * t337;
t82 = -t154 * t337 + t203 * t334;
t577 = -t81 * t334 + t82 * t337;
t238 = t278 + t528;
t358 = t276 - t529;
t155 = t238 * t534 + t335 * t358;
t327 = t336 * pkin(3);
t288 = -pkin(4) * t303 + t327;
t157 = t203 + t288;
t69 = -t155 * t334 + t157 * t337;
t70 = t155 * t337 + t157 * t334;
t385 = -t69 * t334 + t70 * t337;
t576 = t336 * mrSges(4,1) + t338 * mrSges(4,2);
t575 = -mrSges(4,1) * t338 + mrSges(4,2) * t336;
t119 = -Ifges(7,5) * t269 + t317 * t373;
t444 = t337 * t119;
t116 = -Ifges(7,6) * t269 + t373 * t581;
t451 = t334 * t116;
t574 = Ifges(6,4) * t269 + t444 / 0.2e1 - t451 / 0.2e1 + Ifges(6,1) * t548;
t573 = t154 * t614 + t582 * t544;
t572 = t339 * t606;
t570 = 0.2e1 * m(7);
t569 = m(6) / 0.2e1;
t567 = m(7) / 0.2e1;
t566 = pkin(5) / 0.2e1;
t565 = pkin(10) / 0.2e1;
t564 = m(5) * pkin(3);
t563 = t69 / 0.2e1;
t562 = -t70 / 0.2e1;
t561 = -t81 / 0.2e1;
t560 = t82 / 0.2e1;
t559 = -t126 / 0.2e1;
t558 = t127 / 0.2e1;
t557 = -t619 / 0.2e1;
t554 = t619 / 0.2e1;
t183 = t391 * t373;
t553 = -t183 / 0.2e1;
t552 = t611 / 0.2e1;
t551 = t186 / 0.2e1;
t192 = -t269 * t492 - t593;
t550 = -t192 / 0.2e1;
t542 = t290 / 0.2e1;
t541 = -t291 / 0.2e1;
t540 = -t391 / 0.2e1;
t539 = t314 / 0.4e1;
t537 = t334 / 0.2e1;
t532 = pkin(5) * t611;
t530 = pkin(5) * t311;
t524 = mrSges(7,3) * t373;
t520 = Ifges(6,5) * t269;
t518 = Ifges(6,6) * t373;
t153 = t238 * t335 - t358 * t534;
t516 = t153 * mrSges(6,1);
t515 = t154 * mrSges(6,2);
t514 = t155 * mrSges(6,2);
t156 = t239 * t534 + t335 * t351;
t513 = t156 * mrSges(6,1);
t509 = t214 * mrSges(6,2);
t501 = t290 * t81;
t500 = t290 * t82;
t323 = -pkin(3) * t338 - pkin(2);
t287 = -pkin(4) * t302 + t323;
t152 = -pkin(5) * t269 - pkin(10) * t373 + t287;
t67 = t152 * t337 - t156 * t334;
t499 = t291 * t67;
t68 = t152 * t334 + t156 * t337;
t498 = t291 * t68;
t497 = t302 * mrSges(5,3);
t496 = t303 * mrSges(5,3);
t489 = t337 * mrSges(7,3);
t478 = t153 * t154;
t477 = t153 * t391;
t213 = -t281 * t534 + t282 * t335;
t476 = t154 * t213;
t475 = t154 * t373;
t473 = t156 * t391;
t470 = t163 * t213;
t469 = t163 * t373;
t465 = t204 * t334;
t464 = t205 * t337;
t421 = t373 * t540;
t437 = t564 / 0.2e1;
t344 = t421 + (t302 * t481 + t303 * t482) * t437 + (t289 * t567 - t291 * t569) * t373 + (t292 * t569 + t400 * t567 + t578) * t269;
t195 = -t269 * t489 + t594;
t349 = t288 * t569 + (t334 * t70 + t337 * t69) * t567 + t192 * t537 + t195 * t535 + t336 * t437;
t23 = t344 - t349 + t580;
t462 = t23 * qJD(2);
t397 = t339 * t610;
t26 = m(7) * (t126 * t204 + t127 * t205 + t470) + m(6) * (t214 * t619 - t397 + t470) + m(5) * (t245 * t282 + t281 * t571 - t397) + m(4) * (t572 * t610 - t397);
t461 = t26 * qJD(1);
t456 = t289 * t183;
t455 = t289 * t311;
t407 = pkin(10) * t438;
t347 = t269 * t592 + t622 + (t269 * t407 - t531) * t567 + t421;
t447 = t334 * t269;
t191 = -mrSges(7,3) * t447 - t593;
t441 = t337 * t269;
t194 = -mrSges(7,3) * t441 + t594;
t350 = (t334 * t82 + t337 * t81) * t568 + t622 + t191 * t538 + t194 * t536;
t30 = mrSges(6,1) * t632 + t347 + t350;
t453 = t30 * qJD(2);
t450 = t334 * t603;
t449 = t334 * t195;
t448 = t334 * t196;
t445 = t337 * t602;
t443 = t337 * t192;
t442 = t337 * t193;
t370 = t375 * t269;
t415 = -t442 / 0.2e1;
t416 = t448 / 0.2e1;
t374 = t415 + t416;
t39 = -t370 + t374;
t439 = t39 * qJD(2);
t434 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t433 = -t524 / 0.2e1;
t427 = t489 / 0.2e1;
t426 = t391 / 0.2e1 + mrSges(6,1) / 0.2e1;
t419 = -t452 / 0.2e1;
t411 = t544 + t616;
t410 = -t316 / 0.2e1 - t581 / 0.2e1;
t409 = -t317 / 0.2e1 + t314 / 0.2e1;
t402 = t291 * t438;
t401 = t438 * t163;
t398 = t317 * t537 + t535 * t581 + t583;
t200 = -mrSges(6,1) * t269 + mrSges(6,2) * t373;
t271 = -mrSges(5,1) * t302 - mrSges(5,2) * t303;
t1 = t153 * t186 + t621 + t68 * t192 + t70 * t193 + t67 * t195 + t69 * t196 + t288 * t200 + t287 * t613 + t323 * t296 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t336 + pkin(3) * t271) * t336 + (-t323 * mrSges(5,1) + t620 * mrSges(5,3) - Ifges(5,4) * t303) * t303 + m(5) * (t276 * t277 + t278 * t279 + t323 * t327) + m(6) * (t155 * t156 + t287 * t288 + t478) + m(7) * (t67 * t69 + t68 * t70 + t478) + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t338 + (Ifges(4,1) - Ifges(4,2)) * t336) * t338 + (Ifges(5,4) * t302 + (-Ifges(5,1) + Ifges(5,2)) * t303 + t624 * mrSges(5,3)) * t302 + (t287 * mrSges(6,1) + Ifges(6,1) * t615 + (t153 - t156) * mrSges(6,3) + t623) * t373 + (t582 * t617 + t574 + (t154 + t155) * mrSges(6,3) - t434 * t373) * t269;
t380 = t153 * t163 + t626;
t386 = t334 * t67 - t337 * t68;
t340 = -m(5) * (t245 * t624 - t327 * t452 + t620 * t571) / 0.2e1 - m(6) * (t155 * t619 - t156 * t163 - t288 * t452 + t380) / 0.2e1 + (t126 * t69 + t127 * t70 + t163 * t386 + t380) * t568 + t195 * t559 + t127 * t550 + t588 * t554 - t163 * t416 - t163 * t415 + t584 * t557 + (t576 - t580) * t452 / 0.2e1 + t611 * t555;
t379 = t464 - t465;
t384 = t204 * t429 + t205 * t427 - t509 / 0.2e1 + (t540 - mrSges(6,1) / 0.2e1) * t213;
t341 = (-t213 * t291 + t214 * t292) * t569 + (t213 * t289 + t290 * t379) * t567 + t281 * mrSges(5,1) / 0.2e1 - t282 * mrSges(5,2) / 0.2e1 + (t281 * t482 + t282 * t481) * t437 + t384 + t576 * t419;
t2 = t340 + t341;
t388 = -t2 * qJD(1) + t1 * qJD(2);
t6 = t82 * t193 + t68 * t191 + t81 * t196 + t67 * t194 + m(7) * (t154 * t156 + t67 * t81 + t68 * t82) + t156 * t186 + t621 + t287 * t199 + t623 * t373 - (t586 + (-Ifges(6,1) / 0.2e1 + t434) * t373 - t574) * t269;
t343 = ((t156 + t386) * t567 + t552 + t374) * t163 + (t126 * t81 + t127 * t82 + t626) * t567 + t126 * t194 / 0.2e1 + t191 * t558 + t619 * t551 + t199 * t419;
t356 = m(7) * (-pkin(5) * t213 + pkin(10) * t379);
t7 = t509 / 0.2e1 + t426 * t213 + (t465 / 0.2e1 - t464 / 0.2e1) * mrSges(7,3) - t356 / 0.2e1 + t343;
t387 = t7 * qJD(1) + t6 * qJD(2);
t359 = t163 * t553 + t193 * t559 + t196 * t558;
t376 = t204 * mrSges(7,1) / 0.2e1 - t205 * mrSges(7,2) / 0.2e1;
t16 = (-t480 / 0.2e1 + t479 / 0.2e1) * t524 + t359 + t376;
t187 = t314 * t373;
t188 = t316 * t373;
t9 = t154 * t183 + t67 * t193 - t68 * t196 + (t386 * mrSges(7,3) + t116 * t536 - t188 * t535 + t312 * t615 + (t119 - t187) * t538) * t373;
t383 = -t16 * qJD(1) + t9 * qJD(2);
t20 = t584 * t373 + (t302 ^ 2 + t303 ^ 2) * mrSges(5,3) + (t442 - t448 + t612) * t269 + m(7) * (-t269 * t386 + t475) + m(6) * (t156 * t269 + t475) + m(5) * (t277 * t303 + t279 * t302);
t342 = (t269 * t381 + t469) * t567 + (t269 * t619 + t469) * t569 + m(5) * (t302 * t245 + t303 * t571) / 0.2e1;
t28 = t342 + t601;
t382 = -qJD(1) * t28 - qJD(2) * t20;
t364 = m(7) * (-t163 * t407 - t628);
t13 = (t163 / 0.2e1 + t555) * mrSges(6,2) + (t163 * t607 + t604) * t567 - t364 / 0.2e1 + t587 * (t554 + t557);
t348 = (t154 * t292 + t156 * t289) * t567 - t513 / 0.2e1 - t473 / 0.2e1 + t289 * t552 + t292 * t551;
t367 = t603 / 0.4e1 - t290 * t194 / 0.2e1 + t196 * t541;
t368 = t602 / 0.4e1 + t191 * t542 + t291 * t193 / 0.2e1;
t4 = t532 / 0.2e1 + (t596 + t548) * Ifges(6,6) + (t615 + t617) * Ifges(6,5) + (t154 / 0.2e1 + t155 / 0.2e1) * mrSges(6,2) + (m(7) * t566 + t426) * t153 + (pkin(10) * t550 - t602 / 0.4e1 - t411 * t316 + (t560 + t562) * mrSges(7,3) + (-pkin(10) * t70 / 0.4e1 + t500 / 0.4e1 + t498 / 0.4e1) * t570 + t368) * t337 + (t195 * t565 - t603 / 0.4e1 + t411 * t314 + (t561 + t563) * mrSges(7,3) + (pkin(10) * t69 / 0.4e1 - t501 / 0.4e1 - t499 / 0.4e1) * t570 + t367) * t334 + t348;
t352 = -t291 * mrSges(6,2) + mrSges(7,3) * t402 + t292 * t587;
t51 = m(7) * (t289 * t292 + t291 * t400) + t352;
t369 = t13 * qJD(1) + t4 * qJD(2) + t51 * qJD(3);
t361 = mrSges(7,1) * t563 + mrSges(7,2) * t562 + Ifges(7,3) * t548;
t10 = -t456 / 0.2e1 + (Ifges(7,5) * t615 - t119 / 0.4e1 + t187 / 0.4e1 + t196 * t542) * t337 + (Ifges(7,6) * t617 + t188 / 0.4e1 + t116 / 0.4e1 + t193 * t542) * t334 + ((t539 - t317 / 0.4e1) * t337 + (t581 / 0.4e1 + t597) * t334 + t290 * t592) * t373 + t361 + t573;
t166 = t398 + t455;
t33 = -t468 / 0.2e1 + t371;
t365 = -t33 * qJD(1) - t10 * qJD(2) + t166 * qJD(3);
t360 = mrSges(7,1) * t561 + mrSges(7,2) * t560 + Ifges(7,3) * t596;
t353 = -t451 / 0.4e1 + t444 / 0.4e1 - t334 * t188 / 0.4e1 - t337 * t187 / 0.4e1 - t573 + t605 * t68 + (-t314 / 0.4e1 + t317 / 0.4e1) * t457 - (t581 + t316) * t458 / 0.4e1;
t345 = t353 + pkin(5) * t553 + (-t373 * t592 + t579) * pkin(10);
t14 = t345 + t360 - t586;
t207 = t409 * t334 + t410 * t337 + t530;
t36 = (t614 + t375) * t163;
t95 = (t566 - t289 / 0.2e1) * t311 + (mrSges(7,2) * t541 + t410) * t337 + (mrSges(7,1) * t541 + t409) * t334;
t357 = t36 * qJD(1) - t14 * qJD(2) + t95 * qJD(3) + t207 * qJD(5);
t96 = -t530 / 0.2e1 + t455 / 0.2e1 - t375 * t291 + t398;
t40 = -t370 - t374;
t31 = t347 - t350;
t27 = t342 - t601;
t25 = t344 + t349;
t17 = t126 * t373 * t428 + t433 * t479 - t359 + t376;
t15 = t345 + Ifges(7,5) * t441 / 0.2e1 - Ifges(7,6) * t447 / 0.2e1 - t360;
t12 = t604 * t567 + t364 / 0.2e1 + t511 / 0.2e1 + t629 + (mrSges(6,2) / 0.2e1 + t607 * t567 - t592 - t578) * t163;
t11 = t353 + t456 / 0.2e1 + t586 + t361 + t433 * t400 + t579 * t290;
t8 = t356 / 0.2e1 + t343 + t384;
t5 = (mrSges(7,3) * t560 + (t498 + t500) * t567 + t269 * t597 + t368) * t337 + ((-t499 - t501) * t567 + mrSges(7,3) * t561 - t269 * t539 + t367) * t334 + (-pkin(5) * t153 + t385 * pkin(10)) * t567 + t348 + t443 * t565 + t70 * t427 + t440 * t544 + t69 * t429 + t446 * t616 - pkin(10) * t449 / 0.2e1 + (t312 / 0.4e1 - Ifges(6,6) / 0.2e1) * t373 + t445 / 0.4e1 + t450 / 0.4e1 - t477 / 0.2e1 + t460 / 0.4e1 + Ifges(6,5) * t615 + t520 / 0.2e1 - t518 / 0.2e1 - t532 / 0.2e1 + t515 / 0.2e1 - t514 / 0.2e1 - t516 / 0.2e1;
t3 = -t340 + t341;
t18 = [t26 * qJD(2) + t633 * t636, t3 * qJD(3) + t27 * qJD(4) + t8 * qJD(5) + t17 * qJD(6) + t461 + (t204 * t196 + m(7) * (t204 * t67 + t205 * t68 + t476) + t205 * t193 + t214 * t612 + m(6) * (t156 * t214 + t476) + m(5) * (t277 * t281 + t279 * t282) + m(4) * (-t533 * pkin(2) + pkin(8) * t572) * t333 + t282 * t497 + t281 * t496 + (m(5) * t323 + m(6) * t287 - mrSges(3,1) + t200 + t271 + t575) * t423 + t584 * t213 + (mrSges(4,3) * t606 - mrSges(3,2)) * t452) * qJD(2), t3 * qJD(2) + ((-t245 * t482 + t481 * t571) * t564 + t363 * mrSges(4,2) - t571 * mrSges(5,2) - t245 * mrSges(5,1) - t362 * mrSges(4,1) + (-m(6) * t291 + m(7) * t289) * t619 - (m(6) * t292 + m(7) * t400 + t494 + t495) * t163 + t630) * qJD(3) + t12 * qJD(5) + t637, qJD(2) * t27, t8 * qJD(2) + t12 * qJD(3) + (m(7) * (-pkin(10) * t401 - t628) - mrSges(7,3) * t401 + t630) * qJD(5) + t637, t17 * qJD(2) + (-mrSges(7,1) * t127 - mrSges(7,2) * t126) * qJD(6) + t636 * t631; -qJD(3) * t2 + qJD(4) * t28 + qJD(5) * t7 - qJD(6) * t16 - t461, qJD(3) * t1 + qJD(4) * t20 + qJD(5) * t6 + qJD(6) * t9 ((m(7) * t153 + t611) * t289 - t291 * t612 + (m(7) * t385 + t443 - t449) * t290 - t518 + t520 + t445 / 0.2e1 + m(6) * (-t153 * t291 + t155 * t292) + t583 * t269 + t385 * mrSges(7,3) + t575 * pkin(8) + (t276 * t482 + t278 * t481) * t564 + t424 * t496 - t514 - t516 + t460 / 0.2e1 + Ifges(4,5) * t338 - Ifges(4,6) * t336 + Ifges(5,5) * t302 + Ifges(5,6) * t303 + t276 * mrSges(5,1) - t278 * mrSges(5,2) - t422 * t497 + t450 / 0.2e1 - t292 * t588 - t477) * qJD(3) + t25 * qJD(4) + t5 * qJD(5) + t11 * qJD(6) + t388, qJD(3) * t25 + qJD(5) * t31 + qJD(6) * t40 - t382, t5 * qJD(3) + t31 * qJD(4) + t15 * qJD(6) + t387 + (-t473 + t603 * t537 + t602 * t535 + t515 - t513 - (-Ifges(6,5) - t583) * t269 + (t312 / 0.2e1 - Ifges(6,6)) * t373 + (-m(7) * t156 - t611) * pkin(5) + (m(7) * t577 + t337 * t191 - t334 * t194) * pkin(10) + t577 * mrSges(7,3)) * qJD(5), t11 * qJD(3) + t40 * qJD(4) + t15 * qJD(5) + (-mrSges(7,1) * t68 - mrSges(7,2) * t67 - t460) * qJD(6) + t383; qJD(2) * t2 + qJD(5) * t13 - qJD(6) * t33 - t635, qJD(4) * t23 + qJD(5) * t4 - qJD(6) * t10 - t388, qJD(5) * t51 + qJD(6) * t166, t462 (m(7) * (-pkin(5) * t292 + pkin(10) * t402) + t352) * qJD(5) + t96 * qJD(6) + t369, t96 * qJD(5) + (-t290 * t391 + t582) * qJD(6) + t365; -qJD(2) * t28, -qJD(3) * t23 - qJD(5) * t30 - qJD(6) * t39 + t382, -t462, 0, -t453, -t311 * qJD(6) - t439; -qJD(2) * t7 - qJD(3) * t13 - qJD(6) * t36 - t635, -qJD(3) * t4 + qJD(4) * t30 + qJD(6) * t14 - t387, -qJD(6) * t95 - t369, t453, -t207 * qJD(6) (-pkin(10) * t391 + t582) * qJD(6) - t357; t16 * qJD(2) + t33 * qJD(3) + t36 * qJD(5), qJD(3) * t10 + qJD(4) * t39 - qJD(5) * t14 - t383, qJD(5) * t95 - t365, t439, t357, 0;];
Cq  = t18;

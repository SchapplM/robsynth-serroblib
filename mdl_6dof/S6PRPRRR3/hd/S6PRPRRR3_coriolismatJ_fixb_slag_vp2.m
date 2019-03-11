% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:57
% EndTime: 2019-03-08 20:32:12
% DurationCPUTime: 8.92s
% Computational Cost: add. (23324->504), mult. (51193->699), div. (0->0), fcn. (61706->12), ass. (0->310)
t310 = sin(pkin(12));
t312 = cos(pkin(12));
t511 = sin(qJ(4));
t514 = cos(qJ(4));
t286 = -t310 * t511 + t312 * t514;
t287 = -t310 * t514 - t312 * t511;
t314 = sin(qJ(5));
t513 = cos(qJ(5));
t355 = t286 * t314 - t287 * t513;
t566 = -t355 / 0.2e1;
t622 = 0.2e1 * t566;
t313 = sin(qJ(6));
t315 = cos(qJ(6));
t311 = sin(pkin(6));
t512 = sin(qJ(2));
t406 = t311 * t512;
t469 = cos(pkin(6));
t344 = t310 * t469 + t312 * t406;
t345 = -t310 * t406 + t312 * t469;
t235 = t344 * t514 + t345 * t511;
t329 = -t344 * t511 + t345 * t514;
t324 = t235 * t513 + t314 * t329;
t316 = cos(qJ(2));
t439 = t311 * t316;
t123 = -t313 * t439 + t315 * t324;
t151 = t235 * t314 - t329 * t513;
t475 = t315 * mrSges(7,2);
t477 = t313 * mrSges(7,1);
t358 = t475 / 0.2e1 + t477 / 0.2e1;
t350 = t151 * t358;
t290 = t475 + t477;
t456 = t151 * t290;
t476 = t313 * mrSges(7,3);
t409 = t476 / 0.2e1;
t410 = -t476 / 0.2e1;
t576 = t409 + t410;
t602 = t350 + t456 / 0.2e1 + t576 * t123;
t467 = t123 * t315;
t122 = -t313 * t324 - t315 * t439;
t468 = t122 * t313;
t363 = t467 - t468;
t615 = m(7) * (t324 - t363) * t151;
t619 = t615 * qJD(1);
t621 = qJD(6) * t602 + t619;
t620 = qJD(4) + qJD(5);
t516 = t315 / 0.2e1;
t518 = -t313 / 0.2e1;
t525 = t355 / 0.2e1;
t304 = Ifges(7,5) * t315;
t495 = Ifges(7,6) * t313;
t548 = t304 - t495;
t263 = t286 * t513 + t287 * t314;
t501 = Ifges(7,4) * t313;
t503 = Ifges(7,1) * t315;
t374 = -t501 + t503;
t571 = Ifges(7,5) * t355 + t263 * t374;
t498 = Ifges(7,2) * t313;
t500 = Ifges(7,4) * t315;
t373 = -t498 + t500;
t572 = Ifges(7,6) * t355 + t263 * t373;
t590 = -t263 / 0.2e1;
t618 = Ifges(6,4) * t622 + t516 * t571 + t518 * t572 + t548 * t525 + (Ifges(6,2) + Ifges(7,3)) * t590;
t505 = pkin(8) + qJ(3);
t386 = t505 * t312;
t387 = t505 * t310;
t267 = -t386 * t511 - t387 * t514;
t234 = t287 * pkin(9) + t267;
t268 = t386 * t514 - t387 * t511;
t337 = pkin(9) * t286 + t268;
t142 = -t234 * t513 + t314 * t337;
t185 = t290 * t355;
t301 = -pkin(3) * t312 - pkin(2);
t273 = -pkin(4) * t286 + t301;
t541 = t234 * t314 + t337 * t513;
t581 = t290 * t263;
t585 = t263 * mrSges(6,2);
t599 = mrSges(6,1) * t355 + t585;
t616 = t142 * t581 + t185 * t541 + t273 * t599;
t291 = Ifges(7,5) * t313 + Ifges(7,6) * t315;
t551 = t355 * t291;
t560 = Ifges(6,6) * t355;
t375 = mrSges(7,1) * t315 - mrSges(7,2) * t313;
t578 = t541 * t375;
t583 = t541 * mrSges(6,1);
t603 = t315 * t572;
t604 = t313 * t571;
t609 = t142 * mrSges(6,2);
t614 = t551 / 0.2e1 - t560 - t578 - t583 + t609 + t603 / 0.2e1 + t604 / 0.2e1;
t294 = Ifges(7,1) * t313 + t500;
t423 = t315 * t294;
t293 = Ifges(7,2) * t315 + t501;
t431 = t313 * t293;
t613 = (t431 / 0.4e1 - t423 / 0.4e1) * t263 - t551 / 0.4e1 + t560 / 0.2e1 + t578 / 0.2e1 + t583 / 0.2e1 - t603 / 0.4e1 - t604 / 0.4e1 - t609 / 0.2e1;
t612 = t581 / 0.2e1;
t611 = -t585 / 0.2e1;
t610 = m(6) * t273;
t608 = mrSges(6,1) + t375;
t607 = t142 * t313;
t606 = t142 * t314;
t605 = t142 * t315;
t466 = t541 * t142;
t553 = t324 * t375;
t558 = t324 * mrSges(6,1);
t586 = t151 * mrSges(6,2);
t601 = -t553 - t558 + t586;
t600 = m(7) * t541 + t581;
t598 = t553 / 0.2e1 + t558 / 0.2e1 - t586 / 0.2e1;
t308 = t313 ^ 2;
t309 = t315 ^ 2;
t420 = t308 + t309;
t580 = t420 * t151;
t596 = -pkin(5) * t324 - pkin(10) * t580;
t509 = pkin(4) * t314;
t302 = pkin(10) + t509;
t417 = t513 * pkin(4);
t303 = -t417 - pkin(5);
t594 = -t302 * t580 + t303 * t324;
t589 = t263 / 0.2e1;
t588 = -t290 / 0.2e1;
t421 = t310 ^ 2 + t312 ^ 2;
t584 = t421 * mrSges(4,3);
t582 = t151 * t541;
t550 = (t304 / 0.2e1 - t495 / 0.2e1) * t263;
t577 = mrSges(6,3) * t355 + t185;
t508 = pkin(5) * t355;
t203 = -pkin(10) * t263 + t508;
t563 = mrSges(7,1) * t355;
t562 = mrSges(7,2) * t355;
t556 = t142 * t324;
t554 = t311 ^ 2 * t316;
t528 = t185 / 0.2e1;
t552 = t324 * t528;
t265 = -t287 * mrSges(5,1) + mrSges(5,2) * t286;
t549 = t265 + t599;
t444 = t355 * t313;
t192 = mrSges(7,2) * t263 - mrSges(7,3) * t444;
t443 = t355 * t315;
t195 = -mrSges(7,1) * t263 - mrSges(7,3) * t443;
t517 = -t315 / 0.2e1;
t547 = t192 * t518 + t195 * t517;
t80 = t203 * t315 + t607;
t81 = t203 * t313 - t605;
t367 = -t313 * t80 + t315 * t81;
t354 = t374 * t355;
t116 = -Ifges(7,5) * t263 + t354;
t428 = t315 * t116;
t353 = t373 * t355;
t113 = -Ifges(7,6) * t263 + t353;
t438 = t313 * t113;
t545 = Ifges(6,4) * t263 + t428 / 0.2e1 - t438 / 0.2e1 + Ifges(6,1) * t525;
t544 = t263 * t548 / 0.4e1 + t142 * t588;
t271 = t287 * t439;
t272 = t286 * t439;
t213 = -t271 * t513 + t272 * t314;
t543 = -t608 * t213 / 0.2e1;
t542 = mrSges(7,3) * t420;
t214 = t271 * t314 + t272 * t513;
t206 = t214 * t315 + t313 * t406;
t453 = t206 * t315;
t205 = -t214 * t313 + t315 * t406;
t454 = t205 * t313;
t485 = t214 * mrSges(6,2);
t540 = t485 / 0.2e1 + (-t453 / 0.2e1 + t454 / 0.2e1) * mrSges(7,3);
t474 = t315 * mrSges(7,3);
t194 = -t263 * t474 + t563;
t362 = t556 + t582;
t145 = -pkin(5) * t263 - pkin(10) * t355 + t273;
t63 = t145 * t315 - t313 * t541;
t64 = t145 * t313 + t315 * t541;
t369 = t313 * t63 - t315 * t64;
t506 = t287 * pkin(4);
t191 = -t263 * t476 - t562;
t527 = t191 / 0.2e1;
t531 = t122 / 0.2e1;
t536 = m(7) / 0.2e1;
t538 = m(6) / 0.2e1;
t154 = t203 - t506;
t68 = t154 * t315 + t607;
t69 = t154 * t313 - t605;
t539 = (t439 * t506 + t362 - t582) * t538 + (t122 * t68 + t123 * t69 + t362) * t536 + t194 * t531 + t123 * t527 + t552 + (-t142 * t538 + (t525 + t566) * mrSges(6,3)) * t324 + (t369 * t536 + t612) * t151;
t537 = -m(7) / 0.2e1;
t535 = m(6) * pkin(4);
t534 = m(7) * pkin(4);
t533 = mrSges(7,1) / 0.2e1;
t532 = -mrSges(7,2) / 0.2e1;
t520 = -t302 / 0.2e1;
t519 = t302 / 0.2e1;
t510 = m(4) * t311;
t507 = pkin(5) * t290;
t504 = mrSges(7,3) * t355;
t499 = Ifges(6,5) * t263;
t482 = t263 * mrSges(6,3);
t479 = t308 * mrSges(7,3);
t478 = t309 * mrSges(7,3);
t473 = t68 * t313;
t472 = t69 * t315;
t464 = t142 * t213;
t463 = t142 * t355;
t458 = t151 * t213;
t457 = t151 * t355;
t352 = t421 * t512;
t379 = t512 * t554;
t23 = m(7) * (t122 * t205 + t123 * t206 + t458) + m(6) * (t214 * t324 - t379 + t458) + m(5) * (t235 * t272 + t271 * t329 - t379) + m(4) * (t352 * t554 - t379);
t450 = t23 * qJD(1);
t382 = t420 * t263;
t396 = t375 * t525;
t419 = t535 / 0.2e1;
t323 = (t302 * t382 + t303 * t355) * t536 - t396 + (t263 * t314 - t355 * t513) * t419 + t589 * t542;
t335 = (t313 * t69 + t315 * t68) * t536 + t313 * t527 + t194 * t516 - m(6) * t506 / 0.2e1;
t24 = t323 - t335 - t549;
t448 = t24 * qJD(2);
t388 = -t309 / 0.2e1 - t308 / 0.2e1;
t376 = mrSges(7,3) * t388;
t333 = -t263 * t376 + t611 + (pkin(10) * t382 - t508) * t536 - t396;
t432 = t313 * t263;
t190 = -mrSges(7,3) * t432 - t562;
t424 = t315 * t263;
t193 = -mrSges(7,3) * t424 + t563;
t336 = (t313 * t81 + t315 * t80) * t537 + t611 + t190 * t518 + t193 * t517;
t29 = mrSges(6,1) * t622 + t333 + t336;
t441 = t29 * qJD(2);
t182 = t375 * t355;
t440 = t303 * t182;
t435 = t313 * t193;
t434 = t313 * t194;
t433 = t313 * t195;
t427 = t315 * t190;
t426 = t315 * t191;
t425 = t315 * t192;
t349 = t358 * t263;
t357 = t433 / 0.2e1 - t425 / 0.2e1;
t38 = -t349 + t357;
t422 = t38 * qJD(2);
t418 = mrSges(7,3) * t472;
t415 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t412 = -t479 / 0.2e1;
t411 = -t478 / 0.2e1;
t408 = t474 / 0.2e1;
t405 = t313 * t513;
t404 = t315 * t513;
t395 = -t439 / 0.2e1;
t392 = t425 / 0.2e1;
t385 = t421 * qJ(3);
t377 = -t405 / 0.2e1;
t200 = -mrSges(6,1) * t263 + mrSges(6,2) * t355;
t1 = t64 * t191 + t69 * t192 + t63 * t194 + t68 * t195 + t301 * t265 + (-Ifges(5,4) * t287 - pkin(4) * t200) * t287 - t506 * t610 + m(7) * (t63 * t68 + t64 * t69 + t466) + (Ifges(5,4) * t286 + (-Ifges(5,1) + Ifges(5,2)) * t287) * t286 + (Ifges(6,1) * t589 + t618) * t355 + (-t415 * t355 + t548 * t590 + t545) * t263 + t616;
t361 = t453 - t454;
t326 = (t213 * t303 + t302 * t361) * t536 + t271 * mrSges(5,1) / 0.2e1 - t272 * mrSges(5,2) / 0.2e1 + (-t213 * t513 + t214 * t314) * t419;
t3 = -t326 + t549 * t395 - (t392 - t433 / 0.2e1) * t151 + t539 + t540 - t543;
t371 = t3 * qJD(1) + t1 * qJD(2);
t6 = t80 * t195 + t63 * t193 + m(7) * (t63 * t80 + t64 * t81 + t466) + t81 * t192 + t64 * t190 + t618 * t355 - (t550 + (-Ifges(6,1) / 0.2e1 + t415) * t355 - t545) * t263 + t616;
t321 = ((t541 + t369) * t536 + t612 + t357) * t151 + (t122 * t80 + t123 * t81 + t556) * t536 + t193 * t531 + t123 * t190 / 0.2e1 + t552 + t599 * t395;
t341 = m(7) * (-pkin(5) * t213 + pkin(10) * t361);
t8 = (t375 / 0.2e1 + mrSges(6,1) / 0.2e1) * t213 - t341 / 0.2e1 + t321 + t540;
t370 = t8 * qJD(1) + t6 * qJD(2);
t368 = t472 - t473;
t366 = t205 * t410 + t206 * t408 - t485 / 0.2e1 + t543;
t334 = (-t467 / 0.2e1 + t468 / 0.2e1) * t504 + t192 * t531 - t123 * t195 / 0.2e1 + t151 * t182 / 0.2e1;
t359 = t205 * t533 + t206 * t532;
t14 = t334 - t359;
t186 = t355 * t293;
t187 = t355 * t294;
t9 = t142 * t182 + t63 * t192 - t64 * t195 + (t369 * mrSges(7,3) + t113 * t517 - t187 * t516 + t291 * t589 + (t116 - t186) * t518) * t355;
t365 = qJD(1) * t14 + qJD(2) * t9;
t20 = t577 * t355 + (t286 ^ 2 + t287 ^ 2) * mrSges(5,3) + t584 + (t425 - t433 + t482) * t263 + m(7) * (-t263 * t369 + t463) + m(6) * (t263 * t541 + t463) + m(5) * (t267 * t287 + t268 * t286) + m(4) * t385;
t319 = m(5) * (t286 * t235 + t287 * t329) / 0.2e1 + (t263 * t324 + t457) * t538 + (t263 * t363 + t457) * t536 + t352 * t510 / 0.2e1;
t330 = (t205 * t315 + t206 * t313) * t537 - (m(4) + m(5) + m(6)) * t406 / 0.2e1;
t26 = t319 + t330;
t364 = -qJD(1) * t26 - qJD(2) * t20;
t356 = t431 / 0.2e1 - t423 / 0.2e1;
t351 = t420 * t513;
t332 = -t608 * t509 + (-mrSges(6,2) + t542) * t417;
t167 = (t302 * t351 + t303 * t314) * t534 + t332;
t322 = ((-t122 * t405 + t123 * t404 + t151 * t314) * pkin(4) + t594) * t536 + t151 * t412 + t151 * t411 - t598;
t327 = t596 * t537 - (t412 + t411) * t151 + t598;
t17 = t322 + t327;
t317 = (t303 * t541 + t367 * t302 + (t404 * t64 - t405 * t63 + t606) * pkin(4)) * t536 + Ifges(6,5) * t589 + t303 * t612 + t435 * t520 + t509 * t528 + t427 * t519 + t80 * t410 + t81 * t408 + pkin(4) * t195 * t377 + t392 * t417 - t613;
t318 = -t499 / 0.2e1 + t68 * t409 - t418 / 0.2e1 + (-t537 * t541 + t612) * pkin(5) + (t368 * t537 + t434 / 0.2e1 - t426 / 0.2e1) * pkin(10) + t613;
t4 = t318 + t317;
t348 = t17 * qJD(1) + t4 * qJD(2) + t167 * qJD(4);
t343 = Ifges(7,3) * t525 + t532 * t69 + t533 * t68;
t10 = -t440 / 0.2e1 + (Ifges(7,5) * t589 - t116 / 0.4e1 + t186 / 0.4e1 + t195 * t519) * t315 + (Ifges(7,6) * t590 + t187 / 0.4e1 + t113 / 0.4e1 + t192 * t519) * t313 + ((t293 / 0.4e1 - t503 / 0.4e1) * t315 - t302 * t376 + (t294 / 0.4e1 + t500 / 0.2e1 - t498 / 0.4e1) * t313) * t355 + t343 + t544;
t278 = t303 * t290;
t331 = t373 * t517 + t374 * t518 + t356;
t204 = -t278 + t331;
t34 = -t456 / 0.2e1 + t350;
t346 = -qJD(1) * t34 - qJD(2) * t10 - qJD(4) * t204;
t342 = Ifges(7,3) * t566 - t80 * mrSges(7,1) / 0.2e1 + t81 * mrSges(7,2) / 0.2e1;
t328 = t278 / 0.2e1 - t507 / 0.2e1 - t331;
t339 = (mrSges(7,1) * t377 + t404 * t532) * pkin(4);
t124 = t339 - t328;
t338 = -t438 / 0.4e1 + t428 / 0.4e1 - t294 * t444 / 0.4e1 - t293 * t443 / 0.4e1 - t544 + t576 * t64 - (t187 + t353) * t313 / 0.4e1 + (-t186 + t354) * t315 / 0.4e1;
t325 = (t388 * t504 + t547) * pkin(10) - pkin(5) * t182 / 0.2e1 + t338;
t13 = t325 + t342 - t550;
t208 = t331 + t507;
t36 = (t588 + t358) * t151;
t340 = t36 * qJD(1) - t13 * qJD(2) + t124 * qJD(4) + t208 * qJD(5);
t125 = t339 + t328;
t39 = -t349 - t357;
t31 = t323 + t335;
t30 = t333 - t336;
t25 = t319 - t330;
t16 = t322 - t327;
t15 = t334 + t359;
t12 = t325 + Ifges(7,5) * t424 / 0.2e1 - Ifges(7,6) * t432 / 0.2e1 - t342;
t11 = t550 + t440 / 0.2e1 + t338 + t343 + t420 * t504 * t520 + t547 * t302;
t7 = t341 / 0.2e1 + t321 + t366;
t5 = -t318 + t317;
t2 = t357 * t151 + (-t265 / 0.2e1 - t599 / 0.2e1) * t439 + t326 + t366 + t539;
t18 = [t23 * qJD(2) + t615 * t620, t25 * qJD(3) + t2 * qJD(4) + t7 * qJD(5) + t15 * qJD(6) + t450 + (m(7) * (t205 * t63 + t206 * t64 + t464) + m(6) * (t214 * t541 + t464) + t206 * t192 + t205 * t195 + t214 * t482 + m(5) * (t267 * t271 + t268 * t272) + (-pkin(2) * t512 + t316 * t385) * t510 + (m(5) * t301 - mrSges(4,1) * t312 - mrSges(5,1) * t286 + mrSges(4,2) * t310 - mrSges(5,2) * t287 - mrSges(3,1) + t200 + t610) * t406 + t577 * t213 + (t271 * t287 + t272 * t286) * mrSges(5,3) + (-mrSges(3,2) + t584) * t439) * qJD(2), qJD(2) * t25, t2 * qJD(2) + (m(7) * t594 - t324 * t513 * t535 - t329 * mrSges(5,2) - t235 * mrSges(5,1) - (t314 * t535 + t478 + t479) * t151 + t601) * qJD(4) + t16 * qJD(5) + t621, t7 * qJD(2) + t16 * qJD(4) + (m(7) * t596 - mrSges(7,3) * t580 + t601) * qJD(5) + t621, t15 * qJD(2) + (-mrSges(7,1) * t123 - mrSges(7,2) * t122) * qJD(6) + t620 * t602; qJD(3) * t26 + qJD(4) * t3 + qJD(5) * t8 + qJD(6) * t14 - t450, qJD(3) * t20 + qJD(4) * t1 + qJD(5) * t6 + qJD(6) * t9, qJD(4) * t31 + qJD(5) * t30 + qJD(6) * t39 - t364, t31 * qJD(3) + (t499 + Ifges(5,5) * t286 + Ifges(5,6) * t287 + (-t513 * t541 - t606) * t535 - t267 * mrSges(5,2) - t268 * mrSges(5,1) + t423 * t589 + t431 * t590 - mrSges(7,3) * t473 + t418 + t600 * t303 + (m(7) * t368 + t426 - t434) * t302 + (-t263 * t417 - t355 * t509) * mrSges(6,3) + t614) * qJD(4) + t5 * qJD(5) + t11 * qJD(6) + t371, t30 * qJD(3) + t5 * qJD(4) + t12 * qJD(6) + t370 + (-(-Ifges(6,5) + t356) * t263 - t600 * pkin(5) + (m(7) * t367 + t427 - t435) * pkin(10) + t367 * mrSges(7,3) + t614) * qJD(5), t39 * qJD(3) + t11 * qJD(4) + t12 * qJD(5) + (-mrSges(7,1) * t64 - mrSges(7,2) * t63 - t551) * qJD(6) + t365; -qJD(2) * t26, -qJD(4) * t24 - qJD(5) * t29 - qJD(6) * t38 + t364, 0, -t448, -t441, -qJD(6) * t290 - t422; -qJD(2) * t3 + qJD(5) * t17 - qJD(6) * t34 - t619, qJD(3) * t24 + qJD(5) * t4 - qJD(6) * t10 - t371, t448, qJD(5) * t167 - qJD(6) * t204 ((-pkin(5) * t314 + pkin(10) * t351) * t534 + t332) * qJD(5) + t125 * qJD(6) + t348, t125 * qJD(5) + (-t302 * t375 + t548) * qJD(6) + t346; -qJD(2) * t8 - qJD(4) * t17 - qJD(6) * t36 - t619, qJD(3) * t29 - qJD(4) * t4 + qJD(6) * t13 - t370, t441, -qJD(6) * t124 - t348, -t208 * qJD(6) (-pkin(10) * t375 + t548) * qJD(6) - t340; -qJD(2) * t14 + qJD(4) * t34 + qJD(5) * t36, qJD(3) * t38 + qJD(4) * t10 - qJD(5) * t13 - t365, t422, qJD(5) * t124 - t346, t340, 0;];
Cq  = t18;

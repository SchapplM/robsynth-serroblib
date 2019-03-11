% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:41:11
% EndTime: 2019-03-08 20:41:24
% DurationCPUTime: 8.52s
% Computational Cost: add. (13270->503), mult. (28355->668), div. (0->0), fcn. (30485->10), ass. (0->330)
t459 = qJD(4) + qJD(5);
t324 = sin(qJ(4));
t327 = cos(qJ(4));
t505 = cos(pkin(6));
t422 = t505 * t327;
t321 = sin(pkin(6));
t548 = cos(qJ(2));
t443 = t321 * t548;
t273 = -t324 * t443 + t422;
t323 = sin(qJ(5));
t353 = -t505 * t324 - t327 * t443;
t547 = cos(qJ(5));
t200 = t273 * t323 - t547 * t353;
t347 = t547 * t273 + t323 * t353;
t322 = sin(qJ(6));
t326 = cos(qJ(6));
t325 = sin(qJ(2));
t481 = t321 * t325;
t143 = t322 * t481 + t326 * t347;
t471 = t326 * t143;
t142 = -t322 * t347 + t326 * t481;
t480 = t322 * t142;
t377 = t347 - t471 + t480;
t625 = t377 * m(7) * t200;
t626 = qJD(1) * t625;
t317 = t322 ^ 2;
t319 = t326 ^ 2;
t462 = t317 + t319;
t577 = mrSges(7,3) * t462;
t552 = t322 / 0.2e1;
t624 = -t326 / 0.2e1;
t510 = t326 * mrSges(7,2);
t515 = t322 * mrSges(7,1);
t381 = t510 / 0.2e1 + t515 / 0.2e1;
t363 = t381 * t200;
t513 = t322 * mrSges(7,3);
t446 = -t513 / 0.2e1;
t479 = t322 * t143;
t298 = t510 + t515;
t498 = t200 * t298;
t564 = mrSges(7,3) / 0.2e1;
t619 = t363 + t143 * t446 + t479 * t564 + t498 / 0.2e1;
t622 = t619 * qJD(6);
t620 = Ifges(7,3) - Ifges(6,1) + Ifges(6,2);
t315 = Ifges(7,4) * t326;
t403 = Ifges(7,2) * t322 - t315;
t534 = Ifges(7,4) * t322;
t404 = Ifges(7,1) * t326 - t534;
t550 = t326 / 0.2e1;
t299 = Ifges(7,2) * t326 + t534;
t300 = Ifges(7,1) * t322 + t315;
t600 = t299 * t552 + t300 * t624;
t340 = -t322 * t404 / 0.2e1 + t403 * t550 + t600;
t289 = t323 * t324 - t547 * t327;
t290 = t323 * t327 + t547 * t324;
t156 = -Ifges(7,6) * t289 + t403 * t290;
t158 = -Ifges(7,5) * t289 - t404 * t290;
t401 = Ifges(7,5) * t322 + Ifges(7,6) * t326;
t372 = t289 * t401;
t368 = t156 * t550 + t158 * t552 - t372 / 0.2e1 + Ifges(6,6) * t289 + (-Ifges(6,5) + t600) * t290;
t511 = t326 * mrSges(7,1);
t514 = t322 * mrSges(7,2);
t297 = -t511 + t514;
t328 = -pkin(2) - pkin(8);
t540 = -pkin(9) + t328;
t296 = t540 * t324;
t423 = t540 * t327;
t354 = t547 * t296 + t323 * t423;
t585 = t354 * t297;
t590 = t354 * mrSges(6,1);
t234 = t296 * t323 - t547 * t423;
t606 = t234 * mrSges(6,2);
t617 = t368 + t585 - t590 + t606;
t586 = t347 * t297;
t591 = t347 * mrSges(6,1);
t607 = t200 * mrSges(6,2);
t616 = t586 - t591 + t607;
t615 = t586 / 0.2e1 - t591 / 0.2e1 + t607 / 0.2e1;
t614 = -t585 / 0.2e1 + t590 / 0.2e1 - t606 / 0.2e1;
t601 = t462 * t200;
t613 = -pkin(5) * t347 - pkin(10) * t601;
t545 = pkin(4) * t323;
t310 = pkin(10) + t545;
t457 = t547 * pkin(4);
t311 = -t457 - pkin(5);
t612 = -t310 * t601 + t311 * t347;
t253 = t289 * t481;
t254 = t290 * t481;
t611 = (t297 / 0.2e1 - mrSges(6,1) / 0.2e1) * t253 - t254 * mrSges(6,2) / 0.2e1;
t610 = t142 / 0.2e1;
t609 = -t200 / 0.2e1;
t559 = -t289 / 0.2e1;
t308 = t324 * pkin(4) + qJ(3);
t608 = m(6) * t308;
t405 = t324 * mrSges(5,1) + t327 * mrSges(5,2);
t421 = -t290 * mrSges(6,1) + t289 * mrSges(6,2);
t355 = -t405 + t421;
t605 = mrSges(4,3) - t355;
t604 = t234 * t322;
t603 = t234 * t323;
t602 = t234 * t326;
t589 = t234 * t347;
t493 = t234 * t354;
t486 = t290 * t322;
t222 = mrSges(7,2) * t289 + mrSges(7,3) * t486;
t468 = t326 * t222;
t485 = t290 * t326;
t224 = -mrSges(7,1) * t289 + mrSges(7,3) * t485;
t475 = t322 * t224;
t582 = -t468 / 0.2e1 + t475 / 0.2e1;
t593 = pkin(5) * t354;
t587 = t311 * t354;
t584 = t377 * t289;
t542 = pkin(5) * t290;
t219 = pkin(10) * t289 + t308 + t542;
t106 = t219 * t326 - t322 * t354;
t107 = t219 * t322 + t326 * t354;
t392 = t106 * t322 - t107 * t326;
t378 = t354 + t392;
t482 = t311 * t298;
t105 = -t340 + t482;
t38 = -t498 / 0.2e1 + t363;
t359 = -t298 / 0.2e1 + t381;
t128 = t359 * t289;
t460 = t128 * qJD(3);
t217 = t289 * t299;
t218 = t289 * t300;
t549 = t326 / 0.4e1;
t551 = t322 / 0.4e1;
t341 = t299 * t549 - t326 * t404 / 0.4e1 + (-t403 + t300) * t551;
t411 = mrSges(7,3) * (t319 / 0.2e1 + t317 / 0.2e1);
t336 = t310 * t411 + t341;
t233 = -pkin(5) * t289 + pkin(10) * t290;
t544 = pkin(4) * t327;
t221 = t233 + t544;
t386 = t221 * t322 - t602;
t387 = t221 * t326 + t604;
t566 = mrSges(7,2) / 0.2e1;
t568 = -mrSges(7,1) / 0.2e1;
t344 = t386 * t566 + t387 * t568;
t509 = t326 * mrSges(7,3);
t225 = t290 * mrSges(7,1) + t289 * t509;
t466 = t326 * t225;
t223 = -t290 * mrSges(7,2) + t289 * t513;
t476 = t322 * t223;
t375 = -t476 / 0.2e1 - t466 / 0.2e1;
t215 = t297 * t289;
t561 = t215 / 0.2e1;
t435 = t311 * t561;
t555 = t298 / 0.2e1;
t438 = t234 * t555;
t314 = Ifges(7,5) * t326;
t451 = t314 / 0.2e1;
t519 = t290 * Ifges(7,5);
t159 = -t404 * t289 + t519;
t469 = t326 * t159;
t518 = t290 * Ifges(7,6);
t157 = t403 * t289 + t518;
t477 = t322 * t157;
t532 = Ifges(7,6) * t322;
t553 = t314 / 0.4e1;
t563 = Ifges(7,3) / 0.2e1;
t9 = t438 + t469 / 0.4e1 + t218 * t551 - t477 / 0.4e1 + t217 * t549 + t435 + t375 * t310 + (-0.3e1 / 0.4e1 * t532 + t553 + t451) * t290 + (t563 + t336) * t289 + t344;
t583 = -t38 * qJD(1) + t9 * qJD(2) + t105 * qJD(4) - t460;
t467 = t326 * t223;
t474 = t322 * t225;
t374 = -t467 / 0.2e1 + t474 / 0.2e1;
t487 = t289 * t323;
t581 = t547 * t290 + t487;
t541 = pkin(5) * t298;
t116 = t340 + t541;
t337 = pkin(10) * t411 + t341;
t376 = t290 * t553 + t438;
t424 = t159 / 0.4e1 + t217 / 0.4e1;
t560 = -t225 / 0.2e1;
t380 = pkin(10) * t560 + t424;
t117 = t233 * t326 + t604;
t118 = t233 * t322 - t602;
t384 = t117 * t568 + t118 * t566;
t425 = -t157 / 0.4e1 + t218 / 0.4e1;
t454 = -pkin(10) * t223 / 0.2e1;
t455 = -pkin(5) * t215 / 0.2e1;
t12 = t455 + (t519 / 0.2e1 + t380) * t326 + (-0.3e1 / 0.4e1 * t518 + t454 + t425) * t322 + (t563 + t337) * t289 + t376 + t384;
t41 = t359 * t200;
t442 = t322 * t547;
t412 = -t442 / 0.2e1;
t441 = t326 * t547;
t567 = -mrSges(7,2) / 0.2e1;
t348 = (mrSges(7,1) * t412 + t441 * t567) * pkin(4);
t79 = (pkin(5) / 0.2e1 - t311 / 0.2e1) * t298 + t348 + t340;
t580 = t41 * qJD(1) - t12 * qJD(2) + t79 * qJD(4) + t116 * qJD(5) + t460;
t402 = t532 - t314;
t256 = t290 * t297;
t579 = -t289 * t577 + t256;
t216 = t290 * t298;
t576 = t216 * t609 + t143 * t222 / 0.2e1 + t224 * t610;
t320 = t327 ^ 2;
t575 = m(5) / 0.2e1;
t574 = -m(6) / 0.2e1;
t573 = m(6) / 0.2e1;
t572 = -m(7) / 0.2e1;
t571 = m(7) / 0.2e1;
t570 = m(6) * pkin(4);
t569 = m(7) * pkin(4);
t565 = -mrSges(6,3) / 0.2e1;
t558 = t289 / 0.2e1;
t556 = t290 / 0.2e1;
t554 = -t310 / 0.2e1;
t417 = t462 * t289;
t546 = m(7) * (-pkin(10) * t417 - t542);
t543 = pkin(5) * t216;
t370 = t289 * t298;
t360 = t370 / 0.2e1;
t408 = t462 * t234;
t427 = t467 / 0.2e1;
t19 = (t474 - t216) * t559 + (t378 * t572 + t427) * t289 + ((t234 - t408) * t572 + t360 + t582) * t290;
t346 = t378 * t571 + t374;
t391 = -t117 * t322 + t118 * t326;
t21 = ((t234 + t391) * t572 + t582) * t290 + (t216 - t346) * t289;
t539 = -t19 * qJD(4) - t21 * qJD(5);
t538 = mrSges(5,1) * t327;
t537 = mrSges(5,2) * t324;
t535 = Ifges(6,4) * t289;
t531 = Ifges(7,3) * t289;
t521 = t289 * mrSges(6,3);
t413 = (0.1e1 - t462) * t290;
t64 = m(7) * t289 * t413;
t463 = t64 * qJD(3);
t507 = -t128 * qJD(6) - t463;
t129 = (t381 + t555) * t289;
t506 = t129 * qJD(6) + t463;
t499 = t200 * t253;
t497 = t347 * t289;
t211 = -t322 * t254 + t326 * t443;
t495 = t211 * t322;
t212 = t326 * t254 + t322 * t443;
t494 = t212 * t326;
t492 = t234 * t253;
t489 = t253 * t289;
t488 = t254 * t290;
t303 = t321 ^ 2 * t325 * t548;
t30 = m(7) * (t142 * t211 + t143 * t212 + t499) + m(6) * (t254 * t347 + t303 + t499) + m(5) * (t303 + (-t320 * t443 + (-t422 + t273) * t324) * t481);
t484 = t30 * qJD(1);
t483 = t311 * t216;
t478 = t322 * t156;
t470 = t326 * t158;
t343 = (t289 * t411 + t375) * t290 + t215 * t558;
t382 = -t514 / 0.2e1 + t511 / 0.2e1;
t34 = t343 - t382;
t464 = t34 * qJD(2);
t461 = t324 ^ 2 + t320;
t23 = ((t200 - t601) * t290 + t584) * t571;
t27 = (t200 * t413 + t584) * t571;
t458 = t23 * qJD(4) + t27 * qJD(5);
t456 = t546 / 0.2e1;
t453 = t310 * t475;
t452 = t310 * t468;
t448 = t521 / 0.2e1;
t445 = t509 / 0.2e1;
t437 = t486 / 0.2e1;
t436 = -t485 / 0.2e1;
t434 = -t481 / 0.2e1;
t433 = t481 / 0.2e1;
t410 = t461 * t481;
t232 = -mrSges(6,1) * t289 - mrSges(6,2) * t290;
t409 = t232 * t433 + t347 * t448 + t576;
t407 = -t106 * t224 - t107 * t222 + t234 * t216 - t308 * t232 - t354 * t521;
t406 = -t537 + t538;
t399 = t211 * t446 + t212 * t445 + t611;
t11 = t106 * t223 - t107 * t225 + t234 * t215 + (t401 * t556 + t157 * t550 + t218 * t624 - t392 * mrSges(7,3) + (t217 + t159) * t552) * t289;
t342 = (t471 / 0.2e1 - t480 / 0.2e1) * t289 * mrSges(7,3) + t223 * t610 + t143 * t560 + t200 * t561;
t383 = t211 * mrSges(7,1) / 0.2e1 + t212 * t567;
t16 = t342 - t383;
t398 = t16 * qJD(1) + t11 * qJD(2);
t397 = t23 * qJD(1) - t19 * qJD(2);
t396 = t27 * qJD(1) - t21 * qJD(2);
t388 = t494 - t495;
t335 = (t388 * t290 + t489) * t571 + (t488 + t489) * t573 + t410 * t575;
t373 = m(7) * (t326 * t142 + t479);
t37 = (t574 - m(5) / 0.2e1) * t481 - t373 / 0.2e1 + t335;
t50 = t476 + t466 + (m(5) + m(4)) * qJ(3) + m(7) * (t106 * t326 + t107 * t322) + t608 + t605;
t395 = qJD(1) * t37 - qJD(2) * t50;
t394 = t23 * qJD(3) + t626;
t393 = t27 * qJD(3) + t626;
t389 = t200 * t354 + t589;
t385 = t451 - t532 / 0.2e1;
t379 = -t518 / 0.4e1 + t425;
t369 = t421 + t579;
t367 = t289 * (-mrSges(6,3) - t298);
t366 = t462 * t547;
t329 = (t481 * t544 + t389 - t589) * t574 + (t387 * t142 + t386 * t143 + t389) * t572 + t406 * t434 + (t448 + t360) * t347 + (-t354 * t574 + t392 * t572 - t374) * t200;
t333 = (t253 * t311 + t388 * t310) * t571 + (-t547 * t253 + t254 * t323) * t570 / 0.2e1 + t434 * t537 + t433 * t538 + t399;
t2 = t232 * t434 + t497 * t565 + t329 + t333 - t576;
t278 = Ifges(6,4) * t290;
t3 = t320 * Ifges(5,4) + t354 * t370 - qJ(3) * t406 - t387 * t225 - t386 * t223 - m(7) * (t106 * t387 + t107 * t386 + t493) + (-Ifges(5,4) * t324 + (Ifges(5,1) - Ifges(5,2)) * t327) * t324 + (t470 / 0.2e1 - t478 / 0.2e1 + t354 * mrSges(6,3) + (Ifges(6,4) - t385) * t289) * t289 + t407 + (t421 - t608) * t544 + (t469 / 0.2e1 - t477 / 0.2e1 - t278 + t620 * t289 + t385 * t290) * t290;
t365 = -t2 * qJD(1) - t3 * qJD(2) - t19 * qJD(3);
t345 = (t565 - t381) * t497 + t409;
t349 = m(7) * (-pkin(5) * t253 + t388 * pkin(10));
t356 = t117 * t142 + t118 * t143 + t589;
t4 = t374 * t200 + (-t494 / 0.2e1 + t495 / 0.2e1) * mrSges(7,3) + (t378 * t200 + t356) * t571 - t349 / 0.2e1 + t345 - t611;
t7 = -t407 - t531 * t556 + m(7) * (t106 * t117 + t107 * t118 + t493) + t159 * t436 + t157 * t437 + t118 * t223 + t117 * t225 + t354 * t367 + (t478 - t535) * t558 + (t402 * t289 + t470 + t535) * t559 + (-Ifges(6,2) * t558 + t402 * t556 + t278 + (-Ifges(6,1) + t620) * t559) * t290;
t364 = t4 * qJD(1) + t7 * qJD(2) - t21 * qJD(3);
t362 = t387 * t513;
t361 = t386 * t509;
t358 = t311 * t290 - t310 * t417;
t352 = Ifges(7,5) * t436 + Ifges(7,6) * t437 - t531 / 0.2e1 + t376;
t332 = ((-t142 * t442 + t143 * t441 + t200 * t323) * pkin(4) + t612) * t571 + t577 * t609 + t615;
t334 = -t564 * t601 + t571 * t613 + t615;
t14 = t332 - t334;
t338 = m(7) * ((t366 * t290 + t487) * pkin(4) + t358);
t53 = t456 - t338 / 0.2e1;
t330 = (t587 + t391 * t310 + (-t106 * t442 + t107 * t441 + t603) * pkin(4)) * t571 - t483 / 0.2e1 + t117 * t446 + t118 * t445 - t453 / 0.2e1 + t452 / 0.2e1 - t370 * t545 / 0.2e1 + pkin(4) * t225 * t412 + t427 * t457 - t614;
t331 = -t572 * t593 - t543 / 0.2e1 + t362 / 0.2e1 - t361 / 0.2e1 + (-t408 * t572 + t582) * pkin(10) + t614;
t8 = t330 + t331;
t339 = (t297 - mrSges(6,1)) * t545 + (-mrSges(6,2) + t577) * t457;
t86 = (t366 * t310 + t311 * t323) * t569 + t339;
t351 = t14 * qJD(1) + t8 * qJD(2) - t53 * qJD(3) + t86 * qJD(4);
t307 = qJ(3) * t443;
t80 = -t541 / 0.2e1 + t482 / 0.2e1 + t348 - t340;
t47 = t338 / 0.2e1 + t456 + t369;
t36 = m(4) * t481 + t373 / 0.2e1 + t335 + (m(6) + m(5)) * t433;
t35 = t343 + t382;
t17 = t342 + t383;
t15 = t332 + t334;
t13 = t455 + (t454 + t379) * t322 + t380 * t326 + t337 * t289 + t352 - t384;
t10 = t435 + (t223 * t554 + t379) * t322 + (t225 * t554 + t424) * t326 + t336 * t289 - t344 + t352;
t6 = -t331 + t330 + t368;
t5 = t356 * t571 + t349 / 0.2e1 + t346 * t200 + t345 + t399;
t1 = t333 - t329 + t409;
t18 = [t30 * qJD(2) + t459 * t625, t36 * qJD(3) + t1 * qJD(4) + t5 * qJD(5) + t17 * qJD(6) + t484 + (-mrSges(6,3) * t488 + t211 * t225 + t212 * t223 + t253 * t367 + ((-mrSges(5,3) * t461 - mrSges(3,1) + mrSges(4,2)) * t325 + (-mrSges(3,2) + t605) * t548) * t321 + 0.2e1 * (t106 * t211 + t107 * t212 + t492) * t571 + 0.2e1 * (t254 * t354 + t308 * t443 + t492) * t573 + 0.2e1 * (t328 * t410 + t307) * t575 + m(4) * (-pkin(2) * t481 + t307)) * qJD(2), qJD(2) * t36 + t458, t1 * qJD(2) + (m(7) * t612 - t347 * t547 * t570 - t353 * mrSges(5,2) - t273 * mrSges(5,1) - (t323 * t570 + t577) * t200 + t616) * qJD(4) + t15 * qJD(5) + t622 + t394, t5 * qJD(2) + t15 * qJD(4) + (m(7) * t613 - mrSges(7,3) * t601 + t616) * qJD(5) + t622 + t393, t17 * qJD(2) + (-mrSges(7,1) * t143 - mrSges(7,2) * t142) * qJD(6) + t459 * t619; -qJD(3) * t37 - qJD(4) * t2 + qJD(5) * t4 + qJD(6) * t16 - t484, qJD(3) * t50 - qJD(4) * t3 + qJD(5) * t7 + qJD(6) * t11, qJD(6) * t35 - t395 + t539 ((-t354 * t547 - t603) * t570 + t361 + m(7) * (-t310 * t408 + t587) + t452 - t453 - t362 - Ifges(5,6) * t327 - Ifges(5,5) * t324 - t483 + t581 * mrSges(6,3) * pkin(4) - t405 * t328 + t617) * qJD(4) + t6 * qJD(5) + t10 * qJD(6) + t365, t6 * qJD(4) + (-m(7) * t593 + t391 * mrSges(7,3) + t543 + (m(7) * t391 + t468 - t475) * pkin(10) + t617) * qJD(5) + t13 * qJD(6) + t364, t35 * qJD(3) + t10 * qJD(4) + t13 * qJD(5) + (-mrSges(7,1) * t107 - mrSges(7,2) * t106 + t372) * qJD(6) + t398; qJD(2) * t37 + t458, qJD(6) * t34 + t395 + t539, t459 * t64 (m(7) * t358 - t570 * t581 + t355 + t579) * qJD(4) + t47 * qJD(5) + t397 + t506, t47 * qJD(4) + (t369 + t546) * qJD(5) + t396 + t506, qJD(6) * t256 + t129 * t459 + t464; qJD(2) * t2 + qJD(5) * t14 - qJD(6) * t38 - t394, qJD(5) * t8 + qJD(6) * t9 - t365, -qJD(5) * t53 - t397 + t507, qJD(5) * t86 + qJD(6) * t105 ((-pkin(5) * t323 + pkin(10) * t366) * t569 + t339) * qJD(5) + t80 * qJD(6) + t351, t80 * qJD(5) + (t297 * t310 - t402) * qJD(6) + t583; -qJD(2) * t4 - qJD(4) * t14 - qJD(6) * t41 - t393, -qJD(4) * t8 + qJD(6) * t12 - t364, qJD(4) * t53 - t396 + t507, -qJD(6) * t79 - t351, -t116 * qJD(6) (pkin(10) * t297 - t402) * qJD(6) - t580; -t16 * qJD(2) + t38 * qJD(4) + t41 * qJD(5), -qJD(3) * t34 - qJD(4) * t9 - qJD(5) * t12 - t398, t128 * t459 - t464, qJD(5) * t79 - t583, t580, 0;];
Cq  = t18;

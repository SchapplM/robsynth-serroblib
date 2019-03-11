% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:16
% EndTime: 2019-03-09 03:56:24
% DurationCPUTime: 5.86s
% Computational Cost: add. (5144->435), mult. (10615->548), div. (0->0), fcn. (7979->12), ass. (0->198)
t504 = sin(pkin(10));
t505 = cos(pkin(10));
t509 = sin(qJ(3));
t513 = cos(qJ(3));
t534 = t504 * t513 + t505 * t509;
t444 = t534 * qJD(1);
t512 = cos(qJ(5));
t428 = t512 * t444;
t575 = qJD(1) * t513;
t576 = qJD(1) * t509;
t447 = -t504 * t576 + t505 * t575;
t508 = sin(qJ(5));
t396 = t447 * t508 + t428;
t511 = cos(qJ(6));
t570 = qJD(6) * t511;
t636 = t396 * t511 + t570;
t517 = qJD(1) ^ 2;
t514 = cos(qJ(1));
t495 = g(2) * t514;
t510 = sin(qJ(1));
t496 = g(1) * t510;
t581 = t496 - t495;
t526 = -qJ(2) * t517 - t581;
t573 = qJD(3) * t513;
t574 = qJD(3) * t509;
t445 = t504 * t574 - t505 * t573;
t446 = t534 * qJD(3);
t533 = t504 * t509 - t505 * t513;
t535 = -t508 * t533 + t512 * t534;
t366 = qJD(5) * t535 - t445 * t508 + t512 * t446;
t407 = -t508 * t534 - t512 * t533;
t497 = qJDD(3) + qJDD(5);
t498 = qJD(3) + qJD(5);
t635 = -t366 * t498 + t407 * t497;
t536 = -t444 * t508 + t512 * t447;
t565 = qJD(1) * qJD(3);
t403 = -qJDD(1) * t534 + t533 * t565;
t562 = qJDD(1) * t513;
t563 = qJDD(1) * t509;
t404 = t504 * t563 - t505 * t562 + t534 * t565;
t572 = qJD(5) * t508;
t352 = -qJD(5) * t428 + t508 * t403 - t512 * t404 - t447 * t572;
t507 = sin(qJ(6));
t555 = t511 * t352 + t507 * t497 + t498 * t570;
t571 = qJD(6) * t507;
t334 = -t536 * t571 + t555;
t333 = t334 * t511;
t381 = t498 * t507 + t511 * t536;
t475 = t511 * t497;
t335 = qJD(6) * t381 + t352 * t507 - t475;
t379 = -t511 * t498 + t507 * t536;
t634 = -t507 * t335 - t636 * t379 + t333;
t332 = t334 * t507;
t353 = qJD(5) * t536 - t512 * t403 - t404 * t508;
t351 = qJDD(6) + t353;
t348 = t507 * t351;
t597 = t396 * t498;
t599 = t536 * t498;
t601 = t381 * t536;
t626 = qJD(6) + t396;
t633 = t497 * MDP(20) + (-t353 + t599) * MDP(19) - t396 ^ 2 * MDP(17) + (t396 * MDP(16) + MDP(17) * t536 - MDP(27) * t626) * t536 + (t352 + t597) * MDP(18) + (t636 * t381 + t332) * MDP(23) + (t636 * t626 + t348 - t601) * MDP(25);
t515 = -pkin(1) - pkin(7);
t467 = qJD(1) * t515 + qJD(2);
t440 = -qJ(4) * t576 + t467 * t509;
t419 = t504 * t440;
t441 = -qJ(4) * t575 + t513 * t467;
t423 = qJD(3) * pkin(3) + t441;
t384 = t505 * t423 - t419;
t608 = pkin(8) * t447;
t368 = qJD(3) * pkin(4) + t384 - t608;
t594 = t505 * t440;
t385 = t504 * t423 + t594;
t609 = pkin(8) * t444;
t369 = t385 - t609;
t343 = t368 * t512 - t369 * t508;
t339 = -pkin(5) * t498 - t343;
t632 = t339 * t396;
t465 = qJDD(1) * t515 + qJDD(2);
t456 = t513 * t465;
t564 = qJD(1) * qJD(4);
t383 = -t513 * t564 - t467 * t574 + qJDD(3) * pkin(3) + t456 + (t509 * t565 - t562) * qJ(4);
t389 = (-qJ(4) * qJD(1) + t467) * t573 + (-qJ(4) * qJDD(1) + t465 - t564) * t509;
t356 = t505 * t383 - t389 * t504;
t342 = qJDD(3) * pkin(4) + pkin(8) * t404 + t356;
t344 = t368 * t508 + t369 * t512;
t357 = t504 * t383 + t505 * t389;
t345 = pkin(8) * t403 + t357;
t614 = qJD(5) * t344 - t512 * t342 + t508 * t345;
t324 = -pkin(5) * t497 + t614;
t490 = qJ(3) + pkin(10) + qJ(5);
t481 = cos(t490);
t560 = t481 * t496;
t629 = t324 + t560;
t480 = sin(t490);
t628 = g(3) * t480 + t481 * t495;
t458 = pkin(3) * t576 + qJD(1) * qJ(2) + qJD(4);
t411 = pkin(4) * t444 + t458;
t472 = g(3) * t481;
t615 = (qJD(5) * t368 + t345) * t512 + t508 * t342 - t369 * t572;
t625 = t396 * t411 + t581 * t480 + t472 - t615;
t622 = pkin(5) * t536;
t602 = t379 * t536;
t621 = t626 * (t626 * pkin(9) + t622);
t482 = pkin(3) * t505 + pkin(4);
t610 = pkin(3) * t504;
t583 = t508 * t482 + t512 * t610;
t499 = qJDD(1) * qJ(2);
t540 = g(1) * t514 + g(2) * t510;
t500 = qJD(1) * qJD(2);
t559 = 0.2e1 * t500;
t620 = 0.2e1 * t499 + t559 - t540;
t340 = pkin(9) * t498 + t344;
t354 = pkin(5) * t396 - pkin(9) * t536 + t411;
t327 = -t340 * t507 + t354 * t511;
t619 = -t327 * t536 + t339 * t571 + t628 * t511;
t328 = t340 * t511 + t354 * t507;
t618 = t328 * t536 + t339 * t570 + t629 * t507;
t617 = -t411 * t536 - t560 - t614 + t628;
t521 = qJD(5) * t407 - t512 * t445 - t446 * t508;
t613 = -t497 * t535 - t498 * t521;
t588 = qJ(4) - t515;
t436 = -qJD(4) * t513 + t574 * t588;
t460 = t588 * t513;
t437 = -qJD(3) * t460 - qJD(4) * t509;
t387 = t505 * t436 - t437 * t504;
t370 = pkin(8) * t446 + t387;
t388 = t504 * t436 + t505 * t437;
t371 = pkin(8) * t445 + t388;
t459 = t588 * t509;
t409 = t459 * t504 - t505 * t460;
t390 = pkin(8) * t533 + t409;
t410 = -t505 * t459 - t504 * t460;
t391 = -pkin(8) * t534 + t410;
t537 = t390 * t512 - t391 * t508;
t329 = qJD(5) * t537 + t370 * t508 + t371 * t512;
t359 = t390 * t508 + t391 * t512;
t493 = t509 * pkin(3);
t589 = qJ(2) + t493;
t424 = pkin(4) * t534 + t589;
t360 = pkin(5) * t535 - pkin(9) * t407 + t424;
t323 = pkin(9) * t497 + t615;
t546 = qJD(6) * t354 + t323;
t612 = t324 * t407 - t339 * t366 - t359 * t351 - (qJD(6) * t360 + t329) * t626 - t535 * t546;
t607 = g(3) * t509;
t606 = pkin(1) * qJDD(1);
t604 = t339 * t407;
t603 = t360 * t351;
t600 = t381 * t507;
t593 = t507 * t510;
t592 = t507 * t514;
t591 = t510 * t511;
t349 = t511 * t351;
t590 = t511 * t514;
t393 = t505 * t441 - t419;
t392 = -t441 * t504 - t594;
t372 = t392 + t609;
t373 = t393 - t608;
t529 = t482 * t512 - t508 * t610;
t585 = -t529 * qJD(5) + t372 * t508 + t373 * t512;
t584 = qJD(5) * t583 + t372 * t512 - t373 * t508;
t582 = t514 * pkin(1) + t510 * qJ(2);
t503 = t513 ^ 2;
t580 = t509 ^ 2 - t503;
t516 = qJD(3) ^ 2;
t579 = -t516 - t517;
t577 = qJD(1) * t458;
t568 = pkin(3) * t573 + qJD(2);
t561 = qJDD(3) * t509;
t554 = t513 * t565;
t415 = pkin(3) * t575 + pkin(4) * t447;
t548 = t626 * t507;
t532 = qJDD(4) + t499 + t500 + (t554 + t563) * pkin(3);
t374 = -pkin(4) * t403 + t532;
t326 = pkin(5) * t353 - pkin(9) * t352 + t374;
t545 = qJD(6) * t340 - t326;
t439 = pkin(9) + t583;
t543 = pkin(9) * t396 + qJD(6) * t439 + t415 + t622;
t542 = qJDD(2) - t606;
t412 = -pkin(4) * t445 + t568;
t538 = -t351 * t439 + t632;
t530 = t349 - (t396 * t507 + t571) * t626;
t528 = -t366 * t511 - t407 * t571;
t527 = 0.2e1 * qJ(2) * t565 + qJDD(3) * t515;
t525 = -pkin(9) * t351 + t343 * t626 + t632;
t522 = -t356 * t533 + t357 * t534 - t384 * t446 - t385 * t445 - t581;
t520 = -t515 * t516 + t620;
t506 = -qJ(4) - pkin(7);
t492 = t514 * qJ(2);
t489 = qJDD(3) * t513;
t438 = -pkin(5) - t529;
t435 = t480 * t590 - t593;
t434 = t480 * t592 + t591;
t433 = t480 * t591 + t592;
t432 = -t480 * t593 + t590;
t336 = pkin(5) * t521 + pkin(9) * t366 + t412;
t330 = qJD(5) * t359 - t370 * t512 + t371 * t508;
t325 = t511 * t326;
t1 = [(-t542 * pkin(1) - g(1) * (-pkin(1) * t510 + t492) - g(2) * t582 + (t559 + t499) * qJ(2)) * MDP(6) + 0.2e1 * (-t509 * t562 + t565 * t580) * MDP(8) + t581 * MDP(2) + (-t329 * t498 + t352 * t424 - t359 * t497 - t366 * t411 + t374 * t407 + t412 * t536 - t481 * t540) * MDP(22) + (t352 * t407 - t366 * t536) * MDP(16) + (-t352 * t535 - t353 * t407 + t366 * t396 - t521 * t536) * MDP(17) + (-(-t379 * t511 - t600) * t366 + (-t332 - t335 * t511 + (t379 * t507 - t381 * t511) * qJD(6)) * t407) * MDP(24) + (-t513 * t516 - t561) * MDP(10) + (qJDD(1) * t503 - 0.2e1 * t509 * t554) * MDP(7) + t540 * MDP(3) + (t509 * t520 + t513 * t527) * MDP(12) + (-t509 * t527 + t513 * t520) * MDP(13) + t613 * MDP(19) + qJDD(1) * MDP(1) + t620 * MDP(5) + (-t330 * t498 + t353 * t424 + t374 * t535 + t396 * t412 + t411 * t521 - t480 * t540 + t497 * t537) * MDP(21) + (-t407 * t348 - t335 * t535 - t521 * t379 + (t366 * t507 - t407 * t570) * t626) * MDP(26) + (t334 * t535 + t349 * t407 + t381 * t521 + t528 * t626) * MDP(25) + (t351 * t535 + t521 * t626) * MDP(27) + (g(1) * t434 - g(2) * t432 - t328 * t521 + t330 * t381 - t537 * t334 + (-(-qJD(6) * t359 + t336) * t626 - t603 + t545 * t535 - qJD(6) * t604) * t507 + t612 * t511) * MDP(29) + (-g(1) * t435 - g(2) * t433 + t325 * t535 + t327 * t521 + t330 * t379 - t537 * t335 + (t336 * t626 + t603 + (-t340 * t535 - t359 * t626 + t604) * qJD(6)) * t511 + t612 * t507) * MDP(28) + (t333 * t407 + t381 * t528) * MDP(23) + (qJDD(2) - t581 - 0.2e1 * t606) * MDP(4) + (t357 * t410 + t385 * t388 + t356 * t409 + t384 * t387 + t532 * t589 + t458 * t568 - g(1) * (t514 * t493 + t492 + (-pkin(1) + t506) * t510) - g(2) * (t493 * t510 - t506 * t514 + t582)) * MDP(15) + (-t387 * t447 - t388 * t444 + t403 * t410 + t404 * t409 - t522) * MDP(14) + t635 * MDP(18) + (-t509 * t516 + t489) * MDP(9); qJDD(1) * MDP(4) - t517 * MDP(5) + (t542 + t526) * MDP(6) + (t509 * t579 + t489) * MDP(12) + (t513 * t579 - t561) * MDP(13) + (t403 * t534 - t404 * t533 + t444 * t445 + t446 * t447) * MDP(14) + (t522 - t577) * MDP(15) + (-qJD(1) * t396 + t635) * MDP(21) + (-qJD(1) * t536 + t613) * MDP(22) + (-t335 * t407 - t348 * t535 + t366 * t379) * MDP(28) + (-t334 * t407 - t349 * t535 + t366 * t381) * MDP(29) + ((-qJD(1) * t511 - t507 * t521 - t535 * t570) * MDP(28) + (qJD(1) * t507 - t511 * t521 + t535 * t571) * MDP(29)) * t626; -MDP(10) * t563 + (g(3) * t513 + (-t465 - t526) * t509) * MDP(13) + (MDP(7) * t509 * t513 - MDP(8) * t580) * t517 + (-t384 * t392 - t385 * t393 + (t607 + t356 * t505 + t357 * t504 + (-t581 - t577) * t513) * pkin(3)) * MDP(15) + qJDD(3) * MDP(11) + (t438 * t335 - t629 * t511 + t538 * t507 + t584 * t379 + (t507 * t585 - t511 * t543) * t626 + t619) * MDP(28) + (t438 * t334 + t538 * t511 - t628 * t507 + t584 * t381 + (t507 * t543 + t511 * t585) * t626 + t618) * MDP(29) + (-t415 * t536 - t497 * t583 + t498 * t585 + t625) * MDP(22) + MDP(9) * t562 + (-t415 * t396 + t497 * t529 - t498 * t584 + t617) * MDP(21) + (t530 + t602) * MDP(26) + (-t600 * t626 + t634) * MDP(24) + (t513 * t526 + t456 + t607) * MDP(12) + ((t385 + t392) * t447 - (t384 - t393) * t444 + (t403 * t504 + t404 * t505) * pkin(3)) * MDP(14) + t633; (-t444 ^ 2 - t447 ^ 2) * MDP(14) + (t384 * t447 + t385 * t444 + t532 - t540) * MDP(15) + (t353 + t599) * MDP(21) + (t352 - t597) * MDP(22) + (t530 - t602) * MDP(28) + (-t511 * t626 ^ 2 - t348 - t601) * MDP(29); (t344 * t498 + t617) * MDP(21) + (t343 * t498 + t625) * MDP(22) + (-t381 * t548 + t634) * MDP(24) + (-t548 * t626 + t349 + t602) * MDP(26) + (-pkin(5) * t335 - t344 * t379 + t525 * t507 + (-t629 - t621) * t511 + t619) * MDP(28) + (-pkin(5) * t334 - t344 * t381 + t525 * t511 + (-t628 + t621) * t507 + t618) * MDP(29) + t633; t381 * t379 * MDP(23) + (-t379 ^ 2 + t381 ^ 2) * MDP(24) + (t379 * t626 + t555) * MDP(25) + (t381 * t626 + t475) * MDP(26) + t351 * MDP(27) + (-g(1) * t432 - g(2) * t434 + t328 * t626 - t339 * t381 + t325) * MDP(28) + (g(1) * t433 - g(2) * t435 + t327 * t626 + t339 * t379) * MDP(29) + ((-t323 + t472) * MDP(29) + (-MDP(26) * t536 - MDP(28) * t340 - MDP(29) * t354) * qJD(6)) * t511 + (-qJD(6) * t536 * MDP(25) + (-qJD(6) * t498 - t352) * MDP(26) + (-t546 + t472) * MDP(28) + t545 * MDP(29)) * t507;];
tau  = t1;

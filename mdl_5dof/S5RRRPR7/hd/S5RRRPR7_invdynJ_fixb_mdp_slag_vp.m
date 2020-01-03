% Calculate vector of inverse dynamics joint torques for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:30
% EndTime: 2019-12-31 21:17:39
% DurationCPUTime: 6.74s
% Computational Cost: add. (4895->459), mult. (11002->604), div. (0->0), fcn. (7999->14), ass. (0->212)
t515 = sin(qJ(3));
t519 = cos(qJ(2));
t613 = cos(qJ(3));
t569 = qJD(1) * t613;
t516 = sin(qJ(2));
t585 = qJD(1) * t516;
t624 = -t515 * t585 + t519 * t569;
t439 = qJD(5) - t624;
t591 = t515 * t519;
t447 = -qJD(1) * t591 - t516 * t569;
t508 = qJD(2) + qJD(3);
t512 = sin(pkin(9));
t513 = cos(pkin(9));
t426 = t447 * t512 + t513 * t508;
t518 = cos(qJ(5));
t424 = t447 * t513 - t508 * t512;
t514 = sin(qJ(5));
t604 = t424 * t514;
t623 = t426 * t518 + t604;
t626 = t623 * t439;
t625 = t624 * t508;
t549 = t424 * t518 - t426 * t514;
t621 = t439 * t549;
t456 = t512 * t518 + t513 * t514;
t444 = t456 * qJD(5);
t620 = -t456 * t624 + t444;
t455 = t512 * t514 - t518 * t513;
t589 = t439 * t455;
t511 = qJ(2) + qJ(3);
t503 = sin(t511);
t520 = cos(qJ(1));
t595 = t503 * t520;
t517 = sin(qJ(1));
t596 = t503 * t517;
t504 = cos(t511);
t610 = g(3) * t504;
t619 = g(1) * t595 + g(2) * t596 - t610;
t614 = pkin(7) + pkin(6);
t555 = g(1) * t517 - g(2) * t520;
t618 = t555 * t503;
t408 = -pkin(3) * t447 - qJ(4) * t624;
t399 = pkin(2) * t585 + t408;
t471 = t614 * t519;
t461 = qJD(1) * t471;
t448 = t515 * t461;
t470 = t614 * t516;
t459 = qJD(1) * t470;
t416 = -t613 * t459 - t448;
t363 = t513 * t399 - t416 * t512;
t568 = qJD(3) * t613;
t481 = pkin(2) * t568 + qJD(4);
t617 = -t481 * t512 - t363;
t364 = t512 * t399 + t513 * t416;
t563 = t481 * t513 - t364;
t608 = qJD(2) * pkin(2);
t452 = -t459 + t608;
t411 = t613 * t452 - t448;
t366 = t512 * t408 + t513 * t411;
t616 = -qJD(4) * t513 + t366;
t449 = t613 * t461;
t415 = -t515 * t459 + t449;
t584 = qJD(3) * t515;
t558 = pkin(2) * t584 - t415;
t587 = t504 * pkin(3) + t503 * qJ(4);
t506 = qJDD(2) + qJDD(3);
t615 = -pkin(3) * t506 + qJDD(4);
t612 = pkin(2) * t519;
t496 = g(3) * t503;
t609 = t513 * pkin(4);
t505 = t513 * pkin(8);
t607 = qJ(4) * t513;
t565 = qJDD(1) * t613;
t577 = qJDD(1) * t519;
t386 = t515 * t577 + t516 * t565 + t625;
t458 = t613 * t516 + t591;
t419 = t508 * t458;
t578 = qJDD(1) * t516;
t554 = t515 * t578 - t519 * t565;
t387 = t419 * qJD(1) + t554;
t500 = pkin(1) + t612;
t579 = qJD(1) * qJD(2);
t567 = t516 * t579;
t440 = pkin(2) * t567 - t500 * qJDD(1);
t338 = pkin(3) * t387 - qJ(4) * t386 + qJD(4) * t447 + t440;
t566 = t519 * t579;
t421 = qJDD(2) * pkin(2) + t614 * (-t566 - t578);
t427 = t614 * (-t567 + t577);
t529 = t515 * t421 + t613 * t427 + t452 * t568 - t461 * t584;
t346 = t506 * qJ(4) + t508 * qJD(4) + t529;
t324 = t512 * t338 + t513 * t346;
t322 = t324 * t513;
t606 = t387 * t512;
t542 = -t515 * t516 + t613 * t519;
t418 = t508 * t542;
t605 = t418 * t512;
t603 = t624 * t512;
t602 = t624 * t513;
t601 = t458 * t512;
t600 = t458 * t513;
t493 = pkin(2) * t515 + qJ(4);
t597 = t493 * t513;
t594 = t504 * t512;
t593 = t504 * t517;
t592 = t504 * t520;
t575 = t516 * t608;
t362 = pkin(3) * t419 - qJ(4) * t418 - qJD(4) * t458 + t575;
t571 = qJD(2) * t614;
t460 = t516 * t571;
t462 = t519 * t571;
t543 = -t613 * t470 - t515 * t471;
t381 = t543 * qJD(3) - t613 * t460 - t515 * t462;
t335 = t512 * t362 + t513 * t381;
t469 = t500 * qJD(1);
t395 = -pkin(3) * t624 + qJ(4) * t447 - t469;
t412 = t515 * t452 + t449;
t400 = qJ(4) * t508 + t412;
t356 = t512 * t395 + t513 * t400;
t398 = -t508 * pkin(3) + qJD(4) - t411;
t588 = t398 - t481;
t410 = -pkin(3) * t542 - qJ(4) * t458 - t500;
t429 = -t515 * t470 + t613 * t471;
t371 = t512 * t410 + t513 * t429;
t509 = t516 ^ 2;
t586 = -t519 ^ 2 + t509;
t582 = qJD(5) * t514;
t581 = qJD(5) * t518;
t580 = -qJD(4) + t398;
t576 = pkin(8) * t603;
t375 = t386 * t512 - t513 * t506;
t376 = t386 * t513 + t506 * t512;
t573 = -t514 * t375 + t518 * t376 + t426 * t581;
t572 = g(1) * t592 + g(2) * t593 + t496;
t323 = t513 * t338 - t346 * t512;
t318 = pkin(4) * t387 - pkin(8) * t376 + t323;
t319 = -pkin(8) * t375 + t324;
t564 = t518 * t318 - t319 * t514;
t334 = t513 * t362 - t381 * t512;
t562 = t518 * t375 + t514 * t376;
t355 = t513 * t395 - t400 * t512;
t365 = t513 * t408 - t411 * t512;
t370 = t513 * t410 - t429 * t512;
t561 = -t613 * t421 + t515 * t427 + t452 * t584 + t461 * t568;
t499 = -t613 * pkin(2) - pkin(3);
t560 = -t447 * pkin(4) - pkin(8) * t602;
t435 = pkin(4) * t603;
t559 = -t435 + t558;
t557 = -pkin(2) * t516 - pkin(3) * t503;
t556 = g(1) * t520 + g(2) * t517;
t553 = t318 * t514 + t319 * t518;
t552 = -t323 * t512 + t322;
t339 = -pkin(4) * t624 + pkin(8) * t424 + t355;
t340 = pkin(8) * t426 + t356;
t320 = t339 * t518 - t340 * t514;
t321 = t339 * t514 + t340 * t518;
t352 = -pkin(4) * t542 - pkin(8) * t600 + t370;
t361 = -pkin(8) * t601 + t371;
t551 = t352 * t518 - t361 * t514;
t550 = t352 * t514 + t361 * t518;
t548 = t500 + t587;
t547 = t355 * t602 + t356 * t603 + t322 - t572;
t546 = t503 * t556;
t545 = t555 * t504;
t544 = -0.2e1 * pkin(1) * t579 - pkin(6) * qJDD(2);
t451 = t505 + t597;
t541 = qJD(5) * t451 + t560 - t617;
t450 = (-pkin(8) - t493) * t512;
t540 = -qJD(5) * t450 - t563 - t576;
t468 = t505 + t607;
t539 = qJD(4) * t512 + qJD(5) * t468 + t365 + t560;
t467 = (-pkin(8) - qJ(4)) * t512;
t538 = -qJD(5) * t467 - t576 + t616;
t328 = t424 * t582 + t573;
t350 = t561 + t615;
t535 = -t561 + t619;
t522 = qJD(2) ^ 2;
t534 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t522 + t555;
t523 = qJD(1) ^ 2;
t533 = pkin(1) * t523 - pkin(6) * qJDD(1) + t556;
t532 = t350 * t458 + t398 * t418 - t556;
t531 = t355 * t447 + (-t350 + t619) * t513;
t530 = -t469 * t447 + t535;
t329 = -t549 * qJD(5) + t562;
t528 = g(3) * t594 - t356 * t447 + (t350 - t546) * t512;
t333 = pkin(4) * t375 + t350;
t369 = -pkin(4) * t426 + t398;
t507 = pkin(9) + qJ(5);
t502 = cos(t507);
t527 = t320 * t447 + t333 * t455 + t620 * t369 + t619 * t502;
t382 = t429 * qJD(3) - t515 * t460 + t613 * t462;
t384 = qJDD(5) + t387;
t526 = (-t328 * t455 - t329 * t456 + t549 * t620 - t589 * t623) * MDP(23) + (t328 * t456 + t549 * t589) * MDP(22) + (t384 * t456 - t589 * t439 - t447 * t549) * MDP(24) + (-t384 * t455 - t439 * t620 + t447 * t623) * MDP(25) + (t386 - t625) * MDP(13) + (-t554 + (-qJD(1) * t458 - t447) * t508) * MDP(14) + (t447 ^ 2 - t624 ^ 2) * MDP(12) + t506 * MDP(15) + (MDP(11) * t624 + t439 * MDP(26)) * t447;
t501 = sin(t507);
t525 = -t321 * t447 + t333 * t456 + (-t546 + t610) * t501 - t589 * t369;
t524 = t469 * t624 - t529 + t572;
t494 = -pkin(3) - t609;
t473 = qJ(4) * t592;
t472 = qJ(4) * t593;
t466 = t499 - t609;
t433 = t501 * t517 + t502 * t592;
t432 = -t501 * t592 + t502 * t517;
t431 = t501 * t520 - t502 * t593;
t430 = t501 * t593 + t502 * t520;
t404 = t455 * t458;
t403 = t456 * t458;
t397 = pkin(4) * t601 - t543;
t383 = t435 + t412;
t354 = pkin(4) * t605 + t382;
t344 = t456 * t418 + t581 * t600 - t582 * t601;
t343 = -t455 * t418 - t444 * t458;
t330 = -pkin(8) * t605 + t335;
t327 = pkin(4) * t419 - t418 * t505 + t334;
t1 = [0.2e1 * (t516 * t577 - t586 * t579) * MDP(5) + (qJDD(1) * t509 + 0.2e1 * t516 * t566) * MDP(4) + t555 * MDP(2) + t556 * MDP(3) + (t544 * t516 + t534 * t519) * MDP(9) + (-t534 * t516 + t544 * t519) * MDP(10) + (-t328 * t403 + t329 * t404 + t343 * t623 + t344 * t549) * MDP(23) + ((t327 * t518 - t330 * t514) * t439 + t551 * t384 - t564 * t542 + t320 * t419 - t354 * t623 + t397 * t329 + t333 * t403 + t369 * t344 - g(1) * t431 - g(2) * t433 + (t321 * t542 - t550 * t439) * qJD(5)) * MDP(27) + (t329 * t542 - t344 * t439 - t384 * t403 + t419 * t623) * MDP(25) + (t386 * t542 - t387 * t458 + t418 * t624 + t419 * t447) * MDP(12) + (-t382 * t508 - t387 * t500 - t419 * t469 - t440 * t542 + t506 * t543 - t575 * t624 + t545) * MDP(16) + (-t323 * t542 - t334 * t624 + t355 * t419 + t370 * t387 - t375 * t543 - t382 * t426 + t532 * t512 + t513 * t545) * MDP(18) + (t324 * t542 + t335 * t624 - t356 * t419 - t371 * t387 - t376 * t543 - t382 * t424 + t532 * t513 - t555 * t594) * MDP(19) + (qJDD(2) * t516 + t519 * t522) * MDP(6) + (qJDD(2) * t519 - t516 * t522) * MDP(7) + (t418 * t508 + t458 * t506) * MDP(13) + (t386 * t458 - t418 * t447) * MDP(11) + (-t381 * t508 - t386 * t500 - t418 * t469 - t429 * t506 + t440 * t458 - t447 * t575 - t618) * MDP(17) + (-t328 * t404 - t343 * t549) * MDP(22) + (t334 * t424 + t335 * t426 - t370 * t376 - t371 * t375 + t618 + (-t323 * t513 - t324 * t512) * t458 + (-t355 * t513 - t356 * t512) * t418) * MDP(20) + (-t419 * t508 + t506 * t542) * MDP(14) + (-t384 * t542 + t419 * t439) * MDP(26) + (-(t327 * t514 + t330 * t518) * t439 - t550 * t384 + t553 * t542 - t321 * t419 - t354 * t549 + t397 * t328 - t333 * t404 + t369 * t343 - g(1) * t430 - g(2) * t432 + (t320 * t542 - t551 * t439) * qJD(5)) * MDP(28) + (-t328 * t542 + t343 * t439 - t384 * t404 - t419 * t549) * MDP(24) + (t323 * t370 + t324 * t371 + t355 * t334 + t356 * t335 - t350 * t543 + t398 * t382 + (-g(1) * t614 - g(2) * t548) * t520 + (g(1) * t548 - g(2) * t614) * t517) * MDP(21) + qJDD(1) * MDP(1); (t350 * t499 - g(1) * (t557 * t520 + t473) - g(2) * (t557 * t517 + t472) - g(3) * (t587 + t612) + t552 * t493 + t558 * t398 + t563 * t356 + t617 * t355) * MDP(21) + (t415 * t508 + (t613 * t506 - t508 * t584 + t585 * t624) * pkin(2) + t530) * MDP(16) + (t416 * t508 + (t447 * t585 - t506 * t515 - t508 * t568) * pkin(2) + t524) * MDP(17) + (-(t450 * t514 + t451 * t518) * t384 + t466 * t328 + (t541 * t514 + t540 * t518) * t439 - t559 * t549 + t525) * MDP(28) + t526 + ((t450 * t518 - t451 * t514) * t384 + t466 * t329 + (t540 * t514 - t541 * t518) * t439 - t559 * t623 + t527) * MDP(27) + MDP(7) * t577 + MDP(6) * t578 + (-t375 * t597 - t363 * t424 + t563 * t426 + (t376 * t493 - t424 * t481 - t323) * t512 + t547) * MDP(20) + (-t387 * t597 + t376 * t499 - t558 * t424 - (t588 * t513 + t364) * t624 + t528) * MDP(19) + (-t493 * t606 + t375 * t499 - t558 * t426 - (t588 * t512 - t363) * t624 + t531) * MDP(18) + qJDD(2) * MDP(8) + (g(3) * t516 + t533 * t519) * MDP(10) + (-g(3) * t519 + t533 * t516) * MDP(9) + (-t516 * t519 * MDP(4) + t586 * MDP(5)) * t523; (-t350 * pkin(3) - t356 * t366 - t355 * t365 - t398 * t412 - g(1) * (-pkin(3) * t595 + t473) - g(2) * (-pkin(3) * t596 + t472) - g(3) * t587 + (-t355 * t512 + t356 * t513) * qJD(4) + t552 * qJ(4)) * MDP(21) + (t411 * t508 + t524) * MDP(17) + (-(t467 * t514 + t468 * t518) * t384 + t494 * t328 + t383 * t549 + (t539 * t514 + t538 * t518) * t439 + t525) * MDP(28) + (t412 * t508 + t530) * MDP(16) + t526 + ((t467 * t518 - t468 * t514) * t384 + t494 * t329 + t383 * t623 + (t538 * t514 - t539 * t518) * t439 + t527) * MDP(27) + (-t375 * t607 - t365 * t424 - t616 * t426 + (qJ(4) * t376 - qJD(4) * t424 - t323) * t512 + t547) * MDP(20) + (-t387 * t607 - pkin(3) * t376 + t412 * t424 - (t580 * t513 + t366) * t624 + t528) * MDP(19) + (-qJ(4) * t606 - pkin(3) * t375 + t412 * t426 - (t580 * t512 - t365) * t624 + t531) * MDP(18); (t424 * t624 + t375) * MDP(18) + (-t426 * t624 + t376) * MDP(19) + (-t424 ^ 2 - t426 ^ 2) * MDP(20) + (-t355 * t424 - t356 * t426 - t535 + t615) * MDP(21) + (t329 - t621) * MDP(27) + (t328 + t626) * MDP(28); t549 * t623 * MDP(22) + (t549 ^ 2 - t623 ^ 2) * MDP(23) + (t573 - t626) * MDP(24) + (-t562 - t621) * MDP(25) + t384 * MDP(26) + (-g(1) * t432 + g(2) * t430 + t321 * t439 + t369 * t549 + t501 * t496 + t564) * MDP(27) + (g(1) * t433 - g(2) * t431 + t320 * t439 - t369 * t623 + t502 * t496 - t553) * MDP(28) + (MDP(24) * t604 + t549 * MDP(25) - t321 * MDP(27) - t320 * MDP(28)) * qJD(5);];
tau = t1;

% Calculate vector of inverse dynamics joint torques for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:31
% EndTime: 2022-01-20 12:08:39
% DurationCPUTime: 3.36s
% Computational Cost: add. (3626->363), mult. (5481->469), div. (0->0), fcn. (3946->16), ass. (0->199)
t495 = qJ(1) + qJ(2);
t484 = cos(t495);
t499 = sin(qJ(2));
t611 = pkin(1) * t499;
t478 = qJD(2) * t611;
t504 = cos(qJ(2));
t610 = pkin(1) * t504;
t576 = qJD(1) * t478 - qJDD(1) * t610;
t620 = -g(2) * t484 - t576;
t491 = qJD(1) + qJD(2);
t502 = cos(qJ(4));
t503 = cos(qJ(3));
t581 = t502 * t503;
t557 = t491 * t581;
t497 = sin(qJ(4));
t498 = sin(qJ(3));
t589 = t497 * t498;
t558 = t491 * t589;
t402 = -t557 + t558;
t501 = cos(qJ(5));
t423 = t497 * t503 + t498 * t502;
t404 = t423 * t491;
t496 = sin(qJ(5));
t602 = t404 * t496;
t363 = -t402 * t501 - t602;
t562 = -qJD(4) - qJD(5);
t480 = qJD(3) - t562;
t619 = t363 * t480;
t530 = t402 * t496 - t404 * t501;
t618 = t480 * t530;
t489 = qJDD(1) + qJDD(2);
t609 = pkin(2) * t489;
t617 = -t609 - t620;
t422 = -t581 + t589;
t561 = qJD(1) * t611;
t438 = pkin(7) * t491 + t561;
t551 = pkin(8) * t491 + t438;
t393 = t551 * t503;
t386 = t502 * t393;
t392 = t551 * t498;
t387 = qJD(3) * pkin(3) - t392;
t531 = -t387 * t497 - t386;
t606 = pkin(9) * t402;
t340 = -t531 - t606;
t474 = -pkin(3) * t503 - pkin(2);
t573 = qJD(1) * t504;
t560 = pkin(1) * t573;
t405 = t474 * t491 - t560;
t366 = pkin(4) * t402 + t405;
t494 = qJ(3) + qJ(4);
t485 = qJ(5) + t494;
t469 = sin(t485);
t566 = qJD(5) * t496;
t470 = cos(t485);
t595 = t484 * t470;
t482 = sin(t495);
t600 = t470 * t482;
t522 = g(1) * t595 + g(2) * t600 + g(3) * t469 + t340 * t566 - t366 * t363;
t490 = qJD(3) + qJD(4);
t569 = qJD(3) * t503;
t552 = t491 * t569;
t346 = qJD(4) * t557 + t423 * t489 - t490 * t558 + t502 * t552;
t488 = qJDD(3) + qJDD(4);
t563 = qJDD(1) * t499;
t572 = qJD(2) * t504;
t407 = pkin(7) * t489 + (qJD(1) * t572 + t563) * pkin(1);
t593 = t489 * t498;
t350 = -t438 * t569 + qJDD(3) * pkin(3) - t407 * t498 + (-t552 - t593) * pkin(8);
t570 = qJD(3) * t498;
t553 = t491 * t570;
t592 = t489 * t503;
t351 = -t438 * t570 + t407 * t503 + (-t553 + t592) * pkin(8);
t513 = qJD(4) * t531 + t350 * t502 - t497 * t351;
t314 = pkin(4) * t488 - pkin(9) * t346 + t513;
t377 = t490 * t423;
t537 = t422 * t489;
t347 = t377 * t491 + t537;
t568 = qJD(4) * t497;
t613 = (qJD(4) * t387 + t351) * t502 + t497 * t350 - t393 * t568;
t315 = -pkin(9) * t347 + t613;
t598 = t482 * t469;
t601 = t469 * t484;
t515 = g(1) * t601 + g(2) * t598 - g(3) * t470 + t314 * t501 - t496 * t315 + t366 * t530;
t479 = qJDD(5) + t488;
t616 = t479 * MDP(25) + t363 * t530 * MDP(21) + (-t363 ^ 2 + t530 ^ 2) * MDP(22);
t612 = -pkin(7) - pkin(8);
t554 = qJD(3) * t612;
t429 = t498 * t554;
t430 = t503 * t554;
t455 = t612 * t498;
t486 = t503 * pkin(8);
t456 = pkin(7) * t503 + t486;
t577 = t455 * t497 + t456 * t502;
t615 = -qJD(4) * t577 + t423 * t560 - t429 * t497 + t430 * t502;
t567 = qJD(4) * t502;
t614 = -t422 * t560 - t429 * t502 - t430 * t497 - t455 * t567 + t456 * t568;
t399 = t404 * pkin(9);
t384 = t497 * t393;
t546 = t387 * t502 - t384;
t339 = -t399 + t546;
t471 = pkin(7) + t611;
t603 = -pkin(8) - t471;
t419 = t603 * t498;
t420 = t471 * t503 + t486;
t579 = t419 * t497 + t420 * t502;
t575 = g(1) * t484 + g(2) * t482;
t548 = t346 * t496 + t347 * t501;
t320 = -qJD(5) * t530 + t548;
t608 = pkin(2) * t491;
t376 = t490 * t422;
t607 = pkin(9) * t376;
t605 = pkin(9) * t423;
t465 = g(1) * t482;
t481 = sin(t494);
t599 = t481 * t484;
t597 = t482 * t481;
t483 = cos(t494);
t596 = t482 * t483;
t594 = t484 * t483;
t591 = t491 * t498;
t590 = t497 * t479;
t588 = t497 * t501;
t585 = t498 * t503;
t337 = pkin(4) * t490 + t339;
t584 = t501 * t337;
t583 = t501 * t340;
t582 = t501 * t479;
t580 = -t392 * t502 - t384;
t439 = -t560 - t608;
t578 = t439 * t570 + t465 * t503;
t492 = t498 ^ 2;
t574 = -t503 ^ 2 + t492;
t571 = qJD(3) * t491;
t565 = qJD(5) * t501;
t559 = pkin(1) * t572;
t477 = pkin(3) * t570;
t556 = t346 * t501 - t347 * t496 - t402 * t565;
t555 = t439 * t569 + t498 * t617;
t369 = pkin(4) * t377 + t477;
t549 = qJD(3) * t603;
t545 = t392 * t497 - t386;
t544 = t419 * t502 - t420 * t497;
t543 = t455 * t502 - t456 * t497;
t541 = -qJD(5) * t337 - t315;
t540 = t491 * t561;
t538 = t369 - t561;
t367 = t543 - t605;
t374 = t377 * pkin(9);
t536 = -qJD(5) * t367 + t374 + t614;
t418 = t422 * pkin(9);
t368 = -t418 + t577;
t535 = qJD(5) * t368 - t607 - t615;
t534 = -t496 * t337 - t583;
t356 = t544 - t605;
t357 = -t418 + t579;
t533 = t356 * t501 - t357 * t496;
t532 = t356 * t496 + t357 * t501;
t372 = t422 * t501 + t423 * t496;
t373 = -t422 * t496 + t423 * t501;
t473 = -pkin(2) - t610;
t506 = qJD(3) ^ 2;
t529 = t471 * t506 + t473 * t489;
t401 = pkin(4) * t422 + t474;
t528 = t465 + t620;
t330 = -qJD(5) * t372 - t376 * t501 - t377 * t496;
t370 = pkin(3) * t553 + t474 * t489 + t576;
t332 = pkin(4) * t347 + t370;
t527 = -g(1) * t598 + g(2) * t601 + t330 * t366 + t332 * t373;
t526 = -g(1) * t597 + g(2) * t599 + t370 * t423 - t376 * t405;
t331 = qJD(5) * t373 - t376 * t496 + t377 * t501;
t525 = g(1) * t600 - g(2) * t595 + t331 * t366 + t332 * t372;
t524 = g(1) * t596 - g(2) * t594 + t370 * t422 + t377 * t405;
t523 = -qJDD(3) * t471 + t473 * t571;
t388 = t498 * t549 + t503 * t559;
t389 = -t498 * t559 + t503 * t549;
t521 = t388 * t502 + t389 * t497 + t419 * t567 - t420 * t568;
t319 = -t404 * t566 + t556;
t519 = -t477 + t561;
t517 = -t439 * t491 - t407 + t575;
t516 = pkin(7) * t506 - t540 - t609;
t514 = -pkin(7) * qJDD(3) + (t560 - t608) * qJD(3);
t512 = -qJD(4) * t579 - t388 * t497 + t389 * t502;
t510 = (-t319 * t372 - t320 * t373 + t330 * t363 + t331 * t530) * MDP(22) + (t319 * t373 - t330 * t530) * MDP(21) + (-t346 * t422 - t347 * t423 + t376 * t402 - t377 * t404) * MDP(15) + (t330 * t480 + t373 * t479) * MDP(23) + (-t331 * t480 - t372 * t479) * MDP(24) + (t346 * t423 - t376 * t404) * MDP(14) + (-t376 * t490 + t423 * t488) * MDP(16) + (-t377 * t490 - t422 * t488) * MDP(17) + 0.2e1 * (t489 * t585 - t571 * t574) * MDP(8) + (t489 * t492 + 0.2e1 * t498 * t552) * MDP(7) + (qJDD(3) * t503 - t498 * t506) * MDP(10) + (qJDD(3) * t498 + t503 * t506) * MDP(9) + t489 * MDP(4);
t509 = g(1) * t594 + g(2) * t596 + g(3) * t481 + t405 * t402 - t613;
t508 = t404 * t402 * MDP(14) + (t319 - t619) * MDP(23) + (-t320 - t618) * MDP(24) + (t402 * t490 + t346) * MDP(16) - t537 * MDP(17) + (-t402 ^ 2 + t404 ^ 2) * MDP(15) + t488 * MDP(18) + t616;
t507 = g(1) * t599 + g(2) * t597 - g(3) * t483 - t405 * t404 + t513;
t505 = cos(qJ(1));
t500 = sin(qJ(1));
t472 = pkin(3) * t502 + pkin(4);
t448 = t474 - t610;
t431 = t478 + t477;
t391 = t401 - t610;
t378 = pkin(3) * t591 + pkin(4) * t404;
t365 = t369 + t478;
t342 = -t399 + t580;
t341 = t545 + t606;
t329 = t512 + t607;
t328 = -t374 + t521;
t1 = [(-t365 * t363 + t391 * t320 + (-qJD(5) * t532 - t328 * t496 + t329 * t501) * t480 + t533 * t479 + t525) * MDP(26) + (-t365 * t530 + t391 * t319 - (qJD(5) * t533 + t328 * t501 + t329 * t496) * t480 - t532 * t479 + t527) * MDP(27) + (t448 * t347 + t431 * t402 + t488 * t544 + t490 * t512 + t524) * MDP(19) + (t489 * t504 * MDP(5) + (-qJDD(1) - t489) * MDP(6) * t499 + ((-MDP(12) * t503 + MDP(13) * t498 - MDP(5)) * t499 * t491 + ((-qJD(1) - t491) * MDP(6) + (-MDP(12) * t498 - MDP(13) * t503) * qJD(3)) * t504) * qJD(2)) * pkin(1) + ((-t529 - t617) * MDP(12) + t523 * MDP(13)) * t503 + (t523 * MDP(12) + (t529 - t465) * MDP(13)) * t498 + t528 * MDP(5) + t575 * MDP(6) + t510 + qJDD(1) * MDP(1) + (t448 * t346 + t431 * t404 - t488 * t579 - t490 * t521 + t526) * MDP(20) + (g(1) * t500 - g(2) * t505) * MDP(2) + (g(1) * t505 + g(2) * t500) * MDP(3) + t578 * MDP(12) + t555 * MDP(13); (t514 * t498 + (-t516 - t617) * t503 + t578) * MDP(12) + (t514 * t503 + (t516 - t465) * t498 + t555) * MDP(13) + ((-t563 + (-qJD(2) + t491) * t573) * pkin(1) + t575) * MDP(6) + (t474 * t346 - t519 * t404 - t577 * t488 + t490 * t614 + t526) * MDP(20) + (t474 * t347 - t519 * t402 + t543 * t488 + t490 * t615 + t524) * MDP(19) + (t401 * t319 - (t367 * t496 + t368 * t501) * t479 + (t496 * t535 + t501 * t536) * t480 - t538 * t530 + t527) * MDP(27) + (t401 * t320 + (t367 * t501 - t368 * t496) * t479 + (t496 * t536 - t501 * t535) * t480 - t538 * t363 + t525) * MDP(26) + t510 + (t528 + t540) * MDP(5); (t580 * t490 + (-t404 * t591 - t488 * t497 - t490 * t567) * pkin(3) + t509) * MDP(20) + (-g(3) * t503 + t498 * t517) * MDP(12) + (g(3) * t498 + t503 * t517) * MDP(13) + (-t545 * t490 + (-t402 * t591 + t488 * t502 - t490 * t568) * pkin(3) + t507) * MDP(19) + t508 + (t378 * t530 + (-t472 * t479 - t314 + (-pkin(3) * t497 * t562 + t341) * t480) * t496 + (-pkin(3) * t590 + (-pkin(3) * t567 - qJD(5) * t472 + t342) * t480 + t541) * t501 + t522) * MDP(27) + MDP(10) * t592 + MDP(9) * t593 + qJDD(3) * MDP(11) + (t472 * t582 + t378 * t363 - (t341 * t501 - t342 * t496) * t480 + (-t496 * t590 + (-t496 * t502 - t588) * t480 * qJD(4)) * pkin(3) + ((-pkin(3) * t588 - t472 * t496) * t480 + t534) * qJD(5) + t515) * MDP(26) + (-MDP(7) * t585 + MDP(8) * t574) * t491 ^ 2; (t490 * t546 + t509) * MDP(20) + (-t490 * t531 + t507) * MDP(19) + (-(-t339 * t496 - t583) * t480 + t534 * qJD(5) + (t363 * t404 - t480 * t566 + t582) * pkin(4) + t515) * MDP(26) + t508 + ((-t340 * t480 - t314) * t496 + (t339 * t480 + t541) * t501 + (t404 * t530 - t496 * t479 - t480 * t565) * pkin(4) + t522) * MDP(27); (t556 - t619) * MDP(23) + (-t548 - t618) * MDP(24) + (-t480 * t534 + t515) * MDP(26) + (-t501 * t315 - t496 * t314 + (-t340 * t496 + t584) * t480 + t522) * MDP(27) + (-MDP(23) * t602 + MDP(24) * t530 + t534 * MDP(26) - MDP(27) * t584) * qJD(5) + t616;];
tau = t1;

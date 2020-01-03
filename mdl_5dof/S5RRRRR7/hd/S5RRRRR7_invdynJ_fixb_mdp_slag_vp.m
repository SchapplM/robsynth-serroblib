% Calculate vector of inverse dynamics joint torques for
% S5RRRRR7
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
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:39
% EndTime: 2019-12-31 22:22:49
% DurationCPUTime: 6.18s
% Computational Cost: add. (5640->416), mult. (13195->549), div. (0->0), fcn. (10007->14), ass. (0->198)
t496 = sin(qJ(3));
t497 = sin(qJ(2));
t578 = qJD(1) * t497;
t559 = t496 * t578;
t501 = cos(qJ(3));
t502 = cos(qJ(2));
t577 = qJD(1) * t502;
t560 = t501 * t577;
t423 = -t559 + t560;
t424 = -t496 * t577 - t501 * t578;
t495 = sin(qJ(4));
t500 = cos(qJ(4));
t393 = t500 * t423 + t424 * t495;
t499 = cos(qJ(5));
t571 = qJD(5) * t499;
t631 = -t393 * t499 + t571;
t530 = t423 * t495 - t500 * t424;
t490 = qJD(2) + qJD(3);
t566 = qJDD(1) * t502;
t568 = qJD(1) * qJD(2);
t557 = t502 * t568;
t567 = qJDD(1) * t497;
t622 = t557 + t567;
t378 = qJD(3) * t560 - t490 * t559 + t496 * t566 + t501 * t622;
t434 = t496 * t502 + t497 * t501;
t403 = t490 * t434;
t534 = t496 * t567 - t501 * t566;
t379 = t403 * qJD(1) + t534;
t573 = qJD(4) * t500;
t574 = qJD(4) * t495;
t343 = t500 * t378 - t495 * t379 + t423 * t573 + t424 * t574;
t489 = qJDD(2) + qJDD(3);
t484 = qJDD(4) + t489;
t485 = qJD(4) + t490;
t494 = sin(qJ(5));
t563 = t499 * t343 + t494 * t484 + t485 * t571;
t572 = qJD(5) * t494;
t331 = -t530 * t572 + t563;
t329 = t331 * t494;
t330 = t331 * t499;
t382 = t485 * t494 + t499 * t530;
t461 = t499 * t484;
t332 = t382 * qJD(5) + t343 * t494 - t461;
t344 = t530 * qJD(4) + t378 * t495 + t500 * t379;
t342 = qJDD(5) + t344;
t338 = t494 * t342;
t339 = t499 * t342;
t380 = -t499 * t485 + t494 * t530;
t621 = qJD(5) - t393;
t628 = t621 * t494;
t630 = t484 * MDP(22) - t344 * MDP(21) - t393 ^ 2 * MDP(19) + (-t393 * t485 + t343) * MDP(20) + (-MDP(18) * t393 + MDP(19) * t530 + MDP(21) * t485 - MDP(29) * t621) * t530 + (t631 * t382 + t329) * MDP(25) + (-t382 * t530 + t631 * t621 + t338) * MDP(27) + (-t494 * t332 - t631 * t380 - t382 * t628 + t330) * MDP(26) + (t380 * t530 - t621 * t628 + t339) * MDP(28);
t493 = qJ(2) + qJ(3);
t488 = qJ(4) + t493;
t476 = sin(t488);
t498 = sin(qJ(1));
t503 = cos(qJ(1));
t536 = g(1) * t503 + g(2) * t498;
t629 = t536 * t476;
t419 = t424 * pkin(8);
t607 = pkin(6) + pkin(7);
t452 = t607 * t502;
t442 = qJD(1) * t452;
t425 = t496 * t442;
t451 = t607 * t497;
t440 = qJD(1) * t451;
t601 = qJD(2) * pkin(2);
t431 = -t440 + t601;
t550 = t501 * t431 - t425;
t376 = t419 + t550;
t369 = pkin(3) * t490 + t376;
t429 = t501 * t442;
t529 = -t431 * t496 - t429;
t605 = pkin(8) * t423;
t377 = -t529 + t605;
t597 = t377 * t495;
t349 = t500 * t369 - t597;
t347 = -pkin(4) * t485 - t349;
t600 = t347 * t393;
t404 = qJDD(2) * pkin(2) - t607 * t622;
t558 = t497 * t568;
t405 = t607 * (-t558 + t566);
t517 = t529 * qJD(3) + t501 * t404 - t496 * t405;
t336 = pkin(3) * t489 - pkin(8) * t378 + t517;
t576 = qJD(3) * t496;
t611 = (qJD(3) * t431 + t405) * t501 + t496 * t404 - t442 * t576;
t341 = -pkin(8) * t379 + t611;
t588 = t500 * t377;
t350 = t495 * t369 + t588;
t610 = t350 * qJD(4) - t500 * t336 + t495 * t341;
t317 = -pkin(4) * t484 + t610;
t477 = cos(t488);
t602 = g(3) * t477;
t624 = t317 + t602;
t481 = -pkin(2) * t502 - pkin(1);
t450 = t481 * qJD(1);
t406 = -pkin(3) * t423 + t450;
t470 = g(3) * t476;
t612 = (qJD(4) * t369 + t341) * t500 + t495 * t336 - t377 * t574;
t511 = -t406 * t393 + t477 * t536 + t470 - t612;
t618 = pkin(4) * t530;
t617 = (pkin(9) * t621 + t618) * t621;
t582 = -t496 * t451 + t501 * t452;
t348 = pkin(9) * t485 + t350;
t353 = -pkin(4) * t393 - pkin(9) * t530 + t406;
t326 = -t348 * t494 + t353 * t499;
t616 = -t326 * t530 + t347 * t572 + t499 * t629;
t327 = t348 * t499 + t353 * t494;
t615 = t327 * t530 + t347 * t571 + t494 * t624;
t508 = -t406 * t530 - t602 - t610 + t629;
t561 = qJD(2) * t607;
t441 = t497 * t561;
t443 = t502 * t561;
t575 = qJD(3) * t501;
t522 = -t501 * t441 - t496 * t443 - t451 * t575 - t452 * t576;
t360 = -pkin(8) * t403 + t522;
t433 = t496 * t497 - t501 * t502;
t402 = t490 * t433;
t516 = -t582 * qJD(3) + t441 * t496 - t501 * t443;
t361 = pkin(8) * t402 + t516;
t548 = -t501 * t451 - t452 * t496;
t387 = -pkin(8) * t434 + t548;
t388 = -pkin(8) * t433 + t582;
t531 = t387 * t500 - t388 * t495;
t323 = t531 * qJD(4) + t360 * t500 + t361 * t495;
t400 = t500 * t433 + t434 * t495;
t357 = -t400 * qJD(4) - t402 * t500 - t403 * t495;
t401 = -t433 * t495 + t434 * t500;
t409 = pkin(3) * t433 + t481;
t362 = pkin(4) * t400 - pkin(9) * t401 + t409;
t364 = t387 * t495 + t388 * t500;
t316 = pkin(9) * t484 + t612;
t543 = qJD(5) * t353 + t316;
t609 = t317 * t401 - t364 * t342 + t347 * t357 - (qJD(5) * t362 + t323) * t621 - t543 * t400;
t606 = pkin(3) * t424;
t599 = t347 * t401;
t598 = t362 * t342;
t594 = t494 * t498;
t593 = t494 * t503;
t592 = t495 * t496;
t591 = t496 * t500;
t590 = t498 * t499;
t589 = t499 * t503;
t549 = t440 * t496 - t429;
t384 = t549 - t605;
t583 = -t501 * t440 - t425;
t385 = t419 + t583;
t480 = t501 * pkin(2) + pkin(3);
t585 = t384 * t495 + t385 * t500 - t480 * t573 - (-t496 * t574 + (t500 * t501 - t592) * qJD(3)) * pkin(2);
t584 = t384 * t500 - t385 * t495 + t480 * t574 + (t496 * t573 + (t495 * t501 + t591) * qJD(3)) * pkin(2);
t581 = pkin(2) * t591 + t495 * t480;
t491 = t497 ^ 2;
t580 = -t502 ^ 2 + t491;
t483 = t497 * t601;
t395 = pkin(3) * t403 + t483;
t420 = pkin(2) * t558 + t481 * qJDD(1);
t368 = pkin(3) * t379 + t420;
t321 = pkin(4) * t344 - pkin(9) * t343 + t368;
t542 = qJD(5) * t348 - t321;
t359 = -pkin(9) * t393 - t606 + t618;
t418 = pkin(9) + t581;
t482 = pkin(2) * t578;
t540 = qJD(5) * t418 + t359 + t482;
t478 = pkin(3) * t495 + pkin(9);
t539 = qJD(5) * t478 + t359;
t351 = t376 * t495 + t588;
t538 = pkin(3) * t574 - t351;
t352 = t376 * t500 - t597;
t537 = -pkin(3) * t573 + t352;
t535 = g(1) * t498 - g(2) * t503;
t533 = -t342 * t418 - t600;
t532 = -t342 * t478 - t600;
t528 = -pkin(2) * t592 + t480 * t500;
t526 = -0.2e1 * pkin(1) * t568 - pkin(6) * qJDD(2);
t525 = t357 * t499 - t401 * t572;
t521 = -pkin(9) * t342 + t349 * t621 - t600;
t504 = qJD(2) ^ 2;
t519 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t504 + t535;
t505 = qJD(1) ^ 2;
t518 = pkin(1) * t505 - pkin(6) * qJDD(1) + t536;
t514 = -t499 * t624 + t616;
t512 = -t494 * t629 + t615;
t486 = sin(t493);
t487 = cos(t493);
t510 = g(3) * t486 - t450 * t423 + t487 * t536 - t611;
t507 = -g(3) * t487 + t450 * t424 + t486 * t536 + t517;
t506 = t424 * t423 * MDP(11) + (-t423 * t490 + t378) * MDP(13) + (-t534 + (-qJD(1) * t434 - t424) * t490) * MDP(14) + (-t423 ^ 2 + t424 ^ 2) * MDP(12) + t489 * MDP(15) + t630;
t479 = -pkin(3) * t500 - pkin(4);
t417 = -pkin(4) - t528;
t415 = t477 * t589 + t594;
t414 = -t477 * t593 + t590;
t413 = -t477 * t590 + t593;
t412 = t477 * t594 + t589;
t407 = t482 - t606;
t358 = t401 * qJD(4) - t402 * t495 + t500 * t403;
t325 = pkin(4) * t358 - pkin(9) * t357 + t395;
t324 = t364 * qJD(4) + t360 * t495 - t361 * t500;
t318 = t499 * t321;
t1 = [0.2e1 * (t497 * t566 - t580 * t568) * MDP(5) + (t481 * t379 + t450 * t403 + t420 * t433 - t423 * t483 + t535 * t487 + t548 * t489 + t516 * t490) * MDP(16) + (qJDD(1) * t491 + 0.2e1 * t497 * t557) * MDP(4) + (-t323 * t485 + t343 * t409 + t357 * t406 - t364 * t484 + t368 * t401 + t395 * t530 - t535 * t476) * MDP(24) + (t343 * t401 + t357 * t530) * MDP(18) + t536 * MDP(3) + t535 * MDP(2) + (t526 * t497 + t519 * t502) * MDP(9) + (-t519 * t497 + t526 * t502) * MDP(10) + qJDD(1) * MDP(1) + (t342 * t400 + t358 * t621) * MDP(29) + (-g(1) * t413 - g(2) * t415 + t318 * t400 + t324 * t380 + t326 * t358 - t531 * t332 + (t325 * t621 + t598 + (-t348 * t400 - t364 * t621 + t599) * qJD(5)) * t499 + t609 * t494) * MDP(30) + (-g(1) * t412 - g(2) * t414 + t324 * t382 - t327 * t358 - t531 * t331 + (-(-qJD(5) * t364 + t325) * t621 - t598 + t542 * t400 - qJD(5) * t599) * t494 + t609 * t499) * MDP(31) + (-t401 * t338 - t332 * t400 - t358 * t380 + (-t357 * t494 - t401 * t571) * t621) * MDP(28) + (t331 * t400 + t401 * t339 + t358 * t382 + t525 * t621) * MDP(27) + (-t343 * t400 - t344 * t401 + t357 * t393 - t358 * t530) * MDP(19) + (-t324 * t485 + t344 * t409 + t358 * t406 + t368 * t400 - t393 * t395 + t535 * t477 + t484 * t531) * MDP(23) + (t481 * t378 - t450 * t402 + t420 * t434 - t424 * t483 - t535 * t486 - t582 * t489 - t522 * t490) * MDP(17) + (-t402 * t490 + t434 * t489) * MDP(13) + (-t378 * t433 - t379 * t434 - t402 * t423 + t403 * t424) * MDP(12) + (t378 * t434 + t402 * t424) * MDP(11) + (-t403 * t490 - t433 * t489) * MDP(14) + (-t358 * t485 - t400 * t484) * MDP(21) + (t357 * t485 + t401 * t484) * MDP(20) + (qJDD(2) * t497 + t502 * t504) * MDP(6) + (qJDD(2) * t502 - t497 * t504) * MDP(7) + (t401 * t330 + t525 * t382) * MDP(25) + ((-t380 * t499 - t382 * t494) * t357 + (-t329 - t332 * t499 + (t380 * t494 - t382 * t499) * qJD(5)) * t401) * MDP(26); t506 + qJDD(2) * MDP(8) + (g(3) * t497 + t518 * t502) * MDP(10) + (t417 * t331 + t533 * t499 + t584 * t382 + (t540 * t494 + t585 * t499) * t621 + t512) * MDP(31) + (t417 * t332 + t533 * t494 + t584 * t380 + (t585 * t494 - t540 * t499) * t621 + t514) * MDP(30) + (t393 * t407 + t528 * t484 - t584 * t485 + t508) * MDP(23) + (-t407 * t530 - t581 * t484 + t585 * t485 + t511) * MDP(24) + (-g(3) * t502 + t518 * t497) * MDP(9) + MDP(7) * t566 + MDP(6) * t567 + (t583 * t490 + (t424 * t578 - t496 * t489 - t490 * t575) * pkin(2) + t510) * MDP(17) + (-t549 * t490 + (t423 * t578 + t501 * t489 - t490 * t576) * pkin(2) + t507) * MDP(16) + (-t497 * t502 * MDP(4) + t580 * MDP(5)) * t505; t506 + (t479 * t332 + t532 * t494 + t538 * t380 + (t537 * t494 - t539 * t499) * t621 + t514) * MDP(30) + (t479 * t331 + t532 * t499 + t538 * t382 + (t539 * t494 + t537 * t499) * t621 + t512) * MDP(31) + (t352 * t485 + (t424 * t530 - t484 * t495 - t485 * t573) * pkin(3) + t511) * MDP(24) + (-t529 * t490 + t507) * MDP(16) + (t550 * t490 + t510) * MDP(17) + (t351 * t485 + (-t393 * t424 + t484 * t500 - t485 * t574) * pkin(3) + t508) * MDP(23); (t350 * t485 + t508) * MDP(23) + (t349 * t485 + t511) * MDP(24) + (-pkin(4) * t332 - t350 * t380 + t521 * t494 + (-t624 - t617) * t499 + t616) * MDP(30) + (-pkin(4) * t331 - t350 * t382 + t521 * t499 + (-t629 + t617) * t494 + t615) * MDP(31) + t630; t382 * t380 * MDP(25) + (-t380 ^ 2 + t382 ^ 2) * MDP(26) + (t380 * t621 + t563) * MDP(27) + (t382 * t621 + t461) * MDP(28) + t342 * MDP(29) + (-g(1) * t414 + g(2) * t412 + t327 * t621 - t347 * t382 + t318) * MDP(30) + (g(1) * t415 - g(2) * t413 + t326 * t621 + t347 * t380) * MDP(31) + ((-t316 + t470) * MDP(31) + (-MDP(28) * t530 - MDP(30) * t348 - MDP(31) * t353) * qJD(5)) * t499 + (-qJD(5) * t530 * MDP(27) + (-qJD(5) * t485 - t343) * MDP(28) + (-t543 + t470) * MDP(30) + t542 * MDP(31)) * t494;];
tau = t1;

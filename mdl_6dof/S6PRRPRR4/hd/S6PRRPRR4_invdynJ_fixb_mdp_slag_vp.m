% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:43
% EndTime: 2019-03-08 22:14:53
% DurationCPUTime: 8.69s
% Computational Cost: add. (3272->491), mult. (7318->656), div. (0->0), fcn. (5782->12), ass. (0->228)
t508 = sin(pkin(6));
t515 = sin(qJ(2));
t605 = qJD(1) * qJD(2);
t589 = t515 * t605;
t462 = t508 * t589;
t519 = cos(qJ(2));
t636 = t508 * t519;
t471 = qJDD(1) * t636;
t507 = qJDD(2) * pkin(2);
t423 = t462 - t471 - t507;
t514 = sin(qJ(3));
t495 = t514 * qJD(4);
t518 = cos(qJ(3));
t601 = qJDD(2) * t518;
t604 = qJD(2) * qJD(3);
t585 = t518 * t604;
t602 = qJDD(2) * t514;
t681 = t585 + t602;
t541 = pkin(3) * t601 + qJ(4) * t681 + qJD(2) * t495 - t423;
t586 = t514 * t604;
t377 = pkin(3) * t586 - t541;
t510 = cos(pkin(6));
t649 = sin(pkin(11));
t582 = t649 * t515;
t509 = cos(pkin(11));
t634 = t509 * t519;
t431 = t510 * t634 - t582;
t581 = t649 * t519;
t635 = t509 * t515;
t433 = -t510 * t581 - t635;
t564 = g(1) * t433 + g(2) * t431;
t598 = g(3) * t636;
t539 = -t564 - t598;
t690 = -t377 + t539;
t500 = qJDD(3) - qJDD(5);
t520 = -pkin(3) - pkin(4);
t588 = t519 * t605;
t647 = qJDD(2) * pkin(8);
t424 = t647 + (qJDD(1) * t515 + t588) * t508;
t638 = t508 * t515;
t593 = qJD(1) * t638;
t453 = qJD(2) * pkin(8) + t593;
t613 = qJD(3) * t514;
t587 = qJD(1) * t613;
t603 = qJDD(1) * t510;
t612 = qJD(3) * t518;
t566 = t514 * t424 + t453 * t612 + t510 * t587 - t518 * t603;
t554 = -qJDD(4) - t566;
t358 = -pkin(9) * t681 + t520 * qJDD(3) - t554;
t502 = qJDD(3) * qJ(4);
t503 = qJD(3) * qJD(4);
t619 = qJD(1) * t510;
t475 = t518 * t619;
t595 = qJD(3) * t475 + t518 * t424 + t514 * t603;
t369 = -t453 * t613 + t502 + t503 + t595;
t680 = t586 - t601;
t359 = pkin(9) * t680 + t369;
t594 = t520 * qJD(3);
t617 = qJD(2) * t514;
t441 = t514 * t453;
t419 = -t441 + t475;
t678 = qJD(4) - t419;
t679 = -pkin(9) * t617 + t678;
t393 = t594 + t679;
t420 = t518 * t453 + t514 * t619;
t614 = qJD(2) * t518;
t405 = -pkin(9) * t614 + t420;
t504 = qJD(3) * qJ(4);
t395 = t405 + t504;
t513 = sin(qJ(5));
t517 = cos(qJ(5));
t610 = qJD(5) * t517;
t611 = qJD(5) * t513;
t551 = -t517 * t358 + t513 * t359 + t393 * t611 + t395 * t610;
t342 = pkin(5) * t500 + t551;
t432 = t510 * t635 + t581;
t637 = t508 * t518;
t398 = t432 * t514 + t509 * t637;
t399 = -t509 * t508 * t514 + t432 * t518;
t434 = -t510 * t582 + t634;
t583 = t508 * t649;
t400 = t434 * t514 - t518 * t583;
t401 = t434 * t518 + t514 * t583;
t436 = -t510 * t518 + t514 * t638;
t437 = t510 * t514 + t515 * t637;
t557 = t436 * t517 - t437 * t513;
t544 = g(3) * t557 + g(1) * (t400 * t517 - t401 * t513) + g(2) * (t398 * t517 - t399 * t513);
t689 = t342 + t544;
t474 = qJD(1) * t636;
t511 = qJD(2) * pkin(2);
t454 = -t474 - t511;
t421 = -pkin(3) * t614 - qJ(4) * t617 + t454;
t406 = pkin(4) * t614 - t421;
t591 = t513 * t614;
t440 = t517 * t617 - t591;
t688 = -t406 * t440 - t544 - t551;
t542 = t513 * t680 + t517 * t681;
t444 = t513 * t514 + t517 * t518;
t547 = t444 * qJD(5);
t382 = -qJD(2) * t547 + t542;
t501 = qJD(3) - qJD(5);
t512 = sin(qJ(6));
t516 = cos(qJ(6));
t608 = qJD(6) * t516;
t596 = t516 * t382 - t512 * t500 - t501 * t608;
t609 = qJD(6) * t512;
t349 = -t440 * t609 + t596;
t413 = t440 * t516 - t501 * t512;
t577 = t382 * t512 + t516 * t500;
t350 = qJD(6) * t413 + t577;
t540 = t513 * t612 + t514 * t610;
t383 = qJD(2) * t540 - qJD(5) * t591 + qJDD(2) * t444 - t517 * t586;
t639 = t440 * t512;
t411 = t516 * t501 + t639;
t378 = qJDD(6) + t383;
t632 = t516 * t378;
t633 = t512 * t378;
t643 = t349 * t512;
t666 = t444 * qJD(2);
t677 = qJD(6) + t666;
t669 = t677 * t413;
t670 = t677 ^ 2;
t686 = t411 * t677;
t687 = -(t512 * (t350 + t669) + (-t349 + t686) * t516) * MDP(24) + (t516 * t669 + t643) * MDP(23) + (-t413 * t440 + t516 * t670 + t633) * MDP(25) - (-t411 * t440 + t512 * t670 - t632) * MDP(26) - (t440 * t501 + t383) * MDP(19) - (-t440 ^ 2 + t666 ^ 2) * MDP(17) - t500 * MDP(20) + (MDP(16) * t666 - MDP(27) * t677) * t440;
t565 = t514 * t594;
t623 = qJ(4) * t612 + t495;
t667 = t565 + t623 + t593;
t676 = MDP(10) + MDP(12);
t599 = MDP(13) * qJDD(2);
t360 = t393 * t517 - t395 * t513;
t356 = pkin(5) * t501 - t360;
t671 = t356 * t677;
t622 = t517 * qJ(4) + t513 * t520;
t410 = t504 + t420;
t407 = -qJD(3) * pkin(3) + t678;
t552 = t518 * pkin(3) + t514 * qJ(4) + pkin(2);
t394 = pkin(5) * t440 + pkin(10) * t666;
t664 = t677 * (pkin(10) * qJD(6) + t394) + t689;
t488 = qJ(4) * t614;
t428 = t520 * t617 + t488;
t448 = -pkin(10) + t622;
t663 = t677 * (qJD(6) * t448 - t394 + t428) - t689;
t548 = t513 * t358 + t517 * t359 + t393 * t610 - t395 * t611;
t341 = -pkin(10) * t500 + t548;
t373 = pkin(5) * t666 - pkin(10) * t440 + t406;
t443 = t518 * pkin(4) + t552;
t555 = t513 * t518 - t514 * t517;
t386 = pkin(5) * t444 + pkin(10) * t555 + t443;
t396 = qJD(3) * t444 - t547;
t657 = pkin(8) - pkin(9);
t460 = t657 * t514;
t461 = t657 * t518;
t415 = t460 * t513 + t461 * t517;
t563 = g(1) * t434 + g(2) * t432;
t418 = t444 * t636;
t449 = t657 * t613;
t450 = qJD(3) * t461;
t556 = t460 * t517 - t461 * t513;
t629 = -qJD(1) * t418 + qJD(5) * t556 - t449 * t517 + t450 * t513;
t662 = -(qJD(6) * t373 + t341) * t444 - t342 * t555 + t356 * t396 + (-qJD(6) * t386 - t629) * t677 + g(3) * t638 - t415 * t378 + t563;
t660 = -g(1) * t401 - g(2) * t399 - g(3) * t437 - (t419 + t441) * qJD(3) + t595;
t364 = t398 * t513 + t399 * t517;
t367 = t400 * t513 + t401 * t517;
t391 = t436 * t513 + t437 * t517;
t659 = g(1) * t367 + g(2) * t364 + g(3) * t391 + t406 * t666 - t548;
t521 = qJD(3) ^ 2;
t656 = pkin(8) * t521;
t650 = g(3) * t519;
t648 = pkin(8) * qJDD(3);
t646 = qJDD(3) * pkin(3);
t361 = t393 * t513 + t395 * t517;
t357 = -pkin(10) * t501 + t361;
t560 = t357 * t512 - t373 * t516;
t645 = t560 * t440;
t346 = t357 * t516 + t373 * t512;
t644 = t346 * t440;
t642 = t356 * t555;
t641 = t666 * t501;
t522 = qJD(2) ^ 2;
t631 = t518 * t522;
t630 = qJDD(1) - g(3);
t628 = qJD(5) * t415 - t449 * t513 - t450 * t517 - t555 * t474;
t561 = -qJ(4) * t513 + t517 * t520;
t627 = qJD(5) * t561 - t405 * t513 + t517 * t679;
t626 = t622 * qJD(5) + t405 * t517 + t513 * t679;
t625 = t518 * t462 + t587 * t636;
t505 = t514 ^ 2;
t506 = t518 ^ 2;
t621 = t505 - t506;
t620 = t505 + t506;
t616 = qJD(2) * t515;
t615 = qJD(2) * t516;
t607 = qJD(6) * t519;
t600 = qJDD(2) * t519;
t597 = t514 * t631;
t592 = qJD(2) * t636;
t584 = t519 * t604;
t580 = t454 - t511;
t571 = -qJD(2) * t552 + t421;
t568 = t501 ^ 2;
t567 = t517 * t501;
t397 = -t517 * t613 - t518 * t611 + t540;
t562 = pkin(5) * t397 - pkin(10) * t396 + t667;
t553 = -t515 * t522 + t600;
t550 = t396 * t516 + t555 * t609;
t543 = g(3) * t418 + t444 * t564;
t534 = -pkin(10) * t378 + t360 * t677 + t671;
t532 = t386 * t378 - t543;
t529 = t423 - t507 + t564 + t656;
t528 = g(1) * t400 + g(2) * t398 + g(3) * t436 - t566;
t527 = -t448 * t378 - t627 * t677 - t671;
t526 = qJD(3) * t420 + t528;
t368 = pkin(4) * t601 + qJD(2) * t565 + t541;
t525 = qJDD(2) * t552 - t656 + t690;
t372 = -t554 - t646;
t524 = t369 * t518 + t372 * t514 + (t407 * t518 - t410 * t514) * qJD(3) - t563;
t447 = pkin(5) - t561;
t446 = pkin(3) * t617 - t488;
t429 = pkin(3) * t613 - t623;
t403 = qJD(3) * t437 + t514 * t592;
t402 = -qJD(3) * t436 + t518 * t592;
t380 = t391 * t516 + t512 * t636;
t379 = -t391 * t512 + t516 * t636;
t348 = qJD(5) * t557 + t402 * t517 + t403 * t513;
t347 = qJD(5) * t391 + t402 * t513 - t403 * t517;
t344 = pkin(5) * t383 - pkin(10) * t382 + t368;
t343 = t516 * t344;
t1 = [t630 * MDP(1) + (t369 * t437 + t372 * t436 + t402 * t410 + t403 * t407 - g(3)) * MDP(15) + (t347 * t501 - t500 * t557) * MDP(21) + (t348 * t501 + t391 * t500) * MDP(22) + ((-t348 * t512 - t391 * t608) * t677 + t379 * t378 + t347 * t411 - t557 * t350) * MDP(28) + (-(t348 * t516 - t391 * t609) * t677 - t380 * t378 + t347 * t413 - t557 * t349) * MDP(29) + (-MDP(11) + MDP(14)) * (qJD(3) * t402 + qJDD(3) * t437 + (t514 * t553 + t518 * t584) * t508) + (t436 * t514 + t437 * t518) * t599 + (t402 * t518 + t403 * t514 + (t436 * t518 - t437 * t514) * qJD(3)) * MDP(13) * qJD(2) + (t553 * MDP(3) + (-qJDD(2) * t515 - t519 * t522) * MDP(4) + (-t377 * t519 + t421 * t616) * MDP(15) + (t383 * t519 - t616 * t666) * MDP(21) + (t382 * t519 - t440 * t616) * MDP(22) + ((-t512 * t607 - t515 * t615) * MDP(28) - (-t512 * t616 + t516 * t607) * MDP(29)) * t677 + t676 * (-t514 * t584 - t515 * t631)) * t508 + t676 * (-qJD(3) * t403 - qJDD(3) * t436 + t600 * t637); (t421 * t429 + t524 * pkin(8) + ((-pkin(8) * g(3) - qJD(1) * t421) * t515 + (-t407 * t514 - t410 * t518) * qJD(1) * t519) * t508 + t690 * t552) * MDP(15) + (t443 * t382 + t406 * t396 + t415 * t500 + t667 * t440 + t629 * t501 + (t508 * t650 - t368 + t564) * t555) * MDP(22) + (t368 * t444 + t383 * t443 + t397 * t406 - t500 * t556 + t628 * t501 + t666 * t667 - t543) * MDP(21) + (-t349 * t516 * t555 + t413 * t550) * MDP(23) + ((-t411 * t516 - t413 * t512) * t396 - (-t643 - t350 * t516 + (t411 * t512 - t413 * t516) * qJD(6)) * t555) * MDP(24) + (-t382 * t555 + t396 * t440) * MDP(16) + (-t396 * t501 + t500 * t555) * MDP(18) + (-t382 * t444 + t383 * t555 - t396 * t666 - t397 * t440) * MDP(17) + (t471 + t539) * MDP(3) + ((qJD(3) * t571 - t648) * t514 + (-qJD(2) * t429 + t525) * t518 + t625) * MDP(12) + ((qJD(3) * t580 - t648) * t514 + (-t529 - t598) * t518 + t625) * MDP(10) + ((t648 + (-t571 - t474) * qJD(3)) * t518 + ((-t429 + t593) * qJD(2) + t525) * t514) * MDP(14) + ((-t648 + (t580 + t474) * qJD(3)) * t518 + ((-t589 + t650) * t508 + t529) * t514) * MDP(11) + (t620 * t647 + (-g(3) * t515 - t588 * t620) * t508 + t524) * MDP(13) + (-t346 * t397 - t556 * t349 + t628 * t413 + (-(-qJD(6) * t357 + t344) * t444 + qJD(6) * t642 + (qJD(6) * t415 - t562) * t677 - t532) * t512 + t662 * t516) * MDP(29) + (t343 * t444 - t560 * t397 - t556 * t350 + t628 * t411 + (t562 * t677 + (-t357 * t444 - t415 * t677 - t642) * qJD(6) + t532) * t516 + t662 * t512) * MDP(28) + (t555 * t633 - t350 * t444 - t397 * t411 + (-t396 * t512 + t555 * t608) * t677) * MDP(26) + (t349 * t444 + t397 * t413 + t550 * t677 - t555 * t632) * MDP(25) + (t378 * t444 + t397 * t677) * MDP(27) + (-t630 * t638 + t563) * MDP(4) + qJDD(2) * MDP(2) + 0.2e1 * (t514 * t601 - t604 * t621) * MDP(6) + (qJDD(2) * t505 + 0.2e1 * t514 * t585) * MDP(5) + (t397 * t501 + t444 * t500) * MDP(19) + (qJDD(3) * t514 + t518 * t521) * MDP(7) + (qJDD(3) * t518 - t514 * t521) * MDP(8); (t447 * t349 + t626 * t413 + t512 * t663 + t527 * t516 - t644) * MDP(29) + (t447 * t350 + t626 * t411 + t527 * t512 - t516 * t663 - t645) * MDP(28) + (-t428 * t440 + t500 * t622 + t501 * t627 - t659) * MDP(22) + (0.2e1 * t502 + 0.2e1 * t503 + (t421 * t518 + t446 * t514) * qJD(2) + t660) * MDP(14) + (-t454 * t614 - t660) * MDP(11) + (t369 * qJ(4) - t372 * pkin(3) - t421 * t446 - t407 * t420 - g(1) * (-pkin(3) * t400 + qJ(4) * t401) - g(2) * (-pkin(3) * t398 + qJ(4) * t399) - g(3) * (-pkin(3) * t436 + qJ(4) * t437) + t678 * t410) * MDP(15) + (-pkin(3) * t514 + qJ(4) * t518) * t599 + (qJD(5) * t666 - t542 + t641) * MDP(18) + MDP(8) * t601 + MDP(7) * t602 + (0.2e1 * t646 - qJDD(4) + (-t421 * t514 + t446 * t518) * qJD(2) + t526) * MDP(12) + (-t428 * t666 - t500 * t561 + t501 * t626 - t688) * MDP(21) + qJDD(3) * MDP(9) + (-t454 * t617 + t526) * MDP(10) - MDP(5) * t597 + t621 * t522 * MDP(6) - t687; (-qJDD(3) - t597) * MDP(12) + t514 * t599 + (-t505 * t522 - t521) * MDP(14) + (-qJD(3) * t410 + t421 * t617 + qJDD(4) - t528 - t646) * MDP(15) + (-t500 * t517 - t513 * t568 - t617 * t666) * MDP(21) + (-t440 * t617 + t500 * t513 - t517 * t568) * MDP(22) + (-t517 * t350 + (t512 * t567 - t514 * t615) * t677 + (-t411 * t501 - t608 * t677 - t633) * t513) * MDP(28) + (-t517 * t349 + (t512 * t617 + t516 * t567) * t677 + (-t413 * t501 + t609 * t677 - t632) * t513) * MDP(29); (t382 - t641) * MDP(18) + (-t361 * t501 + t688) * MDP(21) + (-t360 * t501 + t659) * MDP(22) + (-pkin(5) * t350 - t361 * t411 + t534 * t512 - t516 * t664 + t645) * MDP(28) + (-pkin(5) * t349 - t361 * t413 + t512 * t664 + t534 * t516 + t644) * MDP(29) + t687; t413 * t411 * MDP(23) + (-t411 ^ 2 + t413 ^ 2) * MDP(24) + (t596 + t686) * MDP(25) + (-t577 + t669) * MDP(26) + t378 * MDP(27) + (-t512 * t341 + t343 + t346 * t677 - t356 * t413 - g(1) * (-t367 * t512 + t433 * t516) - g(2) * (-t364 * t512 + t431 * t516) - g(3) * t379) * MDP(28) + (-t516 * t341 - t512 * t344 - t560 * t677 + t356 * t411 - g(1) * (-t367 * t516 - t433 * t512) - g(2) * (-t364 * t516 - t431 * t512) + g(3) * t380) * MDP(29) + (-MDP(25) * t639 - MDP(26) * t413 - MDP(28) * t346 + MDP(29) * t560) * qJD(6);];
tau  = t1;

% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:34
% EndTime: 2019-03-08 21:49:40
% DurationCPUTime: 6.88s
% Computational Cost: add. (3527->543), mult. (7639->685), div. (0->0), fcn. (5526->10), ass. (0->238)
t506 = sin(pkin(6));
t510 = sin(qJ(2));
t513 = cos(qJ(2));
t603 = qJD(1) * qJD(2);
t582 = t513 * t603;
t507 = cos(pkin(6));
t615 = qJD(3) * t507;
t652 = qJDD(2) * pkin(8);
t683 = t652 + (qJDD(1) * t510 + t582) * t506 + qJD(1) * t615;
t509 = sin(qJ(3));
t602 = qJD(2) * qJD(3);
t580 = t509 * t602;
t512 = cos(qJ(3));
t598 = qJDD(2) * t512;
t681 = t580 - t598;
t622 = qJD(1) * t506;
t589 = t510 * t622;
t465 = qJD(2) * pkin(8) + t589;
t621 = qJD(1) * t507;
t625 = -t509 * t465 + t512 * t621;
t680 = qJD(4) - t625;
t596 = MDP(12) * qJDD(2);
t665 = pkin(4) + pkin(8);
t578 = -qJ(4) * t509 - pkin(2);
t508 = sin(qJ(5));
t511 = cos(qJ(5));
t617 = qJD(2) * t512;
t584 = t508 * t617;
t386 = -qJD(5) * t584 + qJDD(3) * t508 + (qJD(3) * qJD(5) - t681) * t511;
t612 = qJD(3) * t511;
t458 = -t584 + t612;
t618 = qJD(2) * t509;
t491 = qJD(5) + t618;
t644 = t458 * t491;
t682 = -t386 + t644;
t410 = -qJD(3) * pkin(3) + t680;
t679 = t491 ^ 2;
t613 = qJD(3) * t509;
t495 = pkin(3) * t613;
t655 = qJ(4) * t512;
t561 = pkin(9) * t509 - t655;
t610 = qJD(4) * t509;
t530 = qJD(3) * t561 - t610;
t421 = t495 + t530;
t635 = t509 * t513;
t425 = (t508 * t635 + t510 * t511) * t506;
t514 = -pkin(3) - pkin(9);
t453 = t512 * t514 + t578;
t611 = qJD(3) * t512;
t463 = t665 * t611;
t476 = t665 * t509;
t607 = qJD(5) * t511;
t608 = qJD(5) * t508;
t678 = qJD(1) * t425 - t511 * t421 + t453 * t608 - t508 * t463 - t476 * t607;
t677 = t511 * t453 + t508 * t476;
t419 = t512 * t465 + t509 * t621;
t502 = qJD(3) * qJ(4);
t411 = -t502 - t419;
t614 = qJD(3) * t508;
t456 = t511 * t617 + t614;
t609 = qJD(5) * t456;
t385 = -t511 * qJDD(3) - t508 * t681 + t609;
t579 = t512 * t602;
t599 = qJDD(2) * t509;
t543 = t579 + t599;
t455 = qJDD(5) + t543;
t441 = t511 * t455;
t676 = -t491 * t608 + t441;
t634 = t511 * t513;
t675 = t508 * t510 - t509 * t634;
t605 = pkin(4) * t618 + t680;
t674 = -MDP(10) + MDP(13);
t673 = MDP(11) - MDP(14);
t672 = MDP(21) + MDP(23);
t597 = -MDP(22) + MDP(25);
t588 = t513 * t622;
t468 = -pkin(3) * t512 + t578;
t620 = qJD(2) * t468;
t420 = -t588 + t620;
t671 = t420 * t618 + qJDD(4);
t670 = t465 * t611 + t683 * t509;
t658 = cos(pkin(10));
t576 = t658 * t510;
t505 = sin(pkin(10));
t642 = t505 * t513;
t443 = t507 * t576 + t642;
t577 = t506 * t658;
t399 = t443 * t512 - t509 * t577;
t575 = t658 * t513;
t643 = t505 * t510;
t445 = -t507 * t643 + t575;
t401 = t505 * t506 * t509 + t445 * t512;
t640 = t506 * t512;
t449 = t507 * t509 + t510 * t640;
t536 = g(1) * t401 + g(2) * t399 + g(3) * t449;
t583 = t510 * t603;
t639 = t506 * t513;
t560 = -qJDD(1) * t639 + t506 * t583;
t442 = -t507 * t575 + t643;
t444 = t507 * t642 + t576;
t568 = g(1) * t444 + g(2) * t442;
t515 = qJD(3) ^ 2;
t663 = pkin(8) * t515;
t669 = 0.2e1 * qJDD(2) * pkin(2) + t506 * (-g(3) * t513 + t583) - t560 + t568 - t663;
t544 = -qJ(4) * t611 - t610;
t548 = pkin(3) * t580 + t560;
t600 = qJDD(2) * t468;
t371 = qJD(2) * t544 + t548 + t600;
t440 = t495 + t544;
t533 = -g(3) * t639 + t568;
t668 = qJD(2) * (-t440 + t589) - t371 + t533 - t600 - t663;
t667 = -qJD(5) * t677 - t508 * t421 + t511 * t463 + t675 * t622;
t666 = t458 ^ 2;
t664 = pkin(5) * t455;
t659 = qJD(2) * pkin(2);
t657 = pkin(8) * qJDD(3);
t654 = qJ(6) * t455;
t651 = qJDD(3) * pkin(3);
t387 = qJD(3) * t514 + t605;
t408 = qJD(2) * t453 - t588;
t360 = t387 * t508 + t408 * t511;
t650 = t360 * t491;
t649 = t385 * t511;
t648 = t442 * t512;
t647 = t444 * t512;
t646 = t456 * t458;
t645 = t456 * t491;
t641 = t506 * t510;
t638 = t508 * t509;
t636 = t509 * t511;
t633 = t512 * t513;
t632 = t514 * t455;
t631 = qJDD(1) - g(3);
t630 = qJ(6) * t611 + qJD(6) * t509 - t678;
t629 = -pkin(5) * t611 - t667;
t407 = pkin(4) * t617 + t419;
t496 = pkin(3) * t618;
t435 = qJD(2) * t561 + t496;
t628 = t508 * t407 + t511 * t435;
t563 = pkin(5) * t511 + qJ(6) * t508;
t550 = -pkin(4) - t563;
t627 = qJD(5) * t563 - qJD(6) * t511 - t550 * t618 + t680;
t477 = t665 * t512;
t503 = t509 ^ 2;
t504 = t512 ^ 2;
t624 = t503 - t504;
t623 = t503 + t504;
t619 = qJD(2) * t506;
t616 = qJD(3) * t419;
t606 = qJD(5) * t514;
t359 = t387 * t511 - t408 * t508;
t604 = qJD(6) - t359;
t601 = qJDD(1) * t507;
t595 = t509 * t641;
t516 = qJD(2) ^ 2;
t593 = t509 * t512 * t516;
t535 = -t512 * t601 + t670;
t532 = qJDD(4) + t535;
t355 = pkin(4) * t543 + qJDD(3) * t514 + t532;
t363 = qJD(2) * t530 + qJDD(2) * t453 + t548;
t592 = -t508 * t355 - t511 * t363 - t387 * t607;
t591 = -pkin(3) * t648 + t442 * t578;
t590 = -pkin(3) * t647 + t444 * t578;
t587 = t510 * t619;
t586 = t513 * t619;
t573 = t511 * t355 - t508 * t363 - t387 * t608 - t408 * t607;
t572 = t465 * t613 - t509 * t601 - t683 * t512;
t393 = t502 + t407;
t571 = pkin(2) * t639 + pkin(8) * t641 + (pkin(3) * t633 + qJ(4) * t635) * t506;
t570 = t456 * t588;
t569 = t458 * t588;
t567 = g(1) * t445 + g(2) * t443;
t541 = -t408 * t608 - t592;
t346 = qJD(6) * t491 + t541 + t654;
t352 = -pkin(5) * t491 + t604;
t565 = -t352 * t618 - t346;
t347 = qJDD(6) - t573 - t664;
t353 = qJ(6) * t491 + t360;
t564 = -t353 * t618 + t347;
t562 = -pkin(5) * t508 + qJ(6) * t511;
t558 = t352 * t508 + t353 * t511;
t555 = -t453 * t508 + t476 * t511;
t467 = qJ(4) - t562;
t500 = qJDD(3) * qJ(4);
t501 = qJD(3) * qJD(4);
t364 = -t500 - t501 + t572;
t448 = -t507 * t512 + t595;
t405 = -t448 * t508 + t506 * t634;
t404 = t448 * t511 + t508 * t639;
t545 = -t455 * t508 - t491 * t607;
t378 = t442 * t636 + t443 * t508;
t380 = t444 * t636 + t445 * t508;
t424 = t675 * t506;
t539 = -g(1) * t380 - g(2) * t378 - g(3) * t424;
t379 = -t442 * t638 + t443 * t511;
t381 = -t444 * t638 + t445 * t511;
t538 = -g(1) * t381 - g(2) * t379 - g(3) * t425;
t398 = t443 * t509 + t512 * t577;
t400 = t445 * t509 - t505 * t640;
t537 = g(1) * t400 + g(2) * t398 + g(3) * t448;
t366 = pkin(5) * t456 - qJ(6) * t458 + t393;
t531 = t366 * t491 + t632;
t527 = t491 * t606 + t536;
t372 = -t398 * t511 + t442 * t508;
t374 = -t400 * t511 + t444 * t508;
t526 = g(1) * t374 + g(2) * t372 - g(3) * t404 + t573;
t356 = -pkin(4) * t681 - t364;
t348 = pkin(5) * t386 + qJ(6) * t385 - qJD(6) * t458 + t356;
t525 = -t348 + t527;
t466 = -t588 - t659;
t524 = -t657 + (t466 + t588 - t659) * qJD(3);
t523 = t657 + (-t420 - t588 - t620) * qJD(3);
t373 = t398 * t508 + t442 * t511;
t375 = t400 * t508 + t444 * t511;
t522 = -g(1) * t375 - g(2) * t373 + g(3) * t405 + t541;
t521 = -qJD(3) * t625 - t536 - t572;
t520 = -t535 + t537;
t518 = t366 * t458 + qJDD(6) - t526;
t365 = t532 - t651;
t517 = -t364 * t512 + t365 * t509 + (t410 * t512 + t411 * t509) * qJD(3) - t567;
t462 = t665 * t613;
t461 = -qJ(4) * t617 + t496;
t439 = t448 * pkin(3);
t417 = t512 * t563 + t477;
t403 = qJD(3) * t449 + t509 * t586;
t402 = -qJD(3) * t595 + (t586 + t615) * t512;
t397 = pkin(5) * t458 + qJ(6) * t456;
t396 = t400 * pkin(3);
t395 = t398 * pkin(3);
t389 = -pkin(5) * t509 - t555;
t388 = qJ(6) * t509 + t677;
t376 = (qJD(5) * t562 + qJD(6) * t508) * t512 + (-pkin(8) + t550) * t613;
t369 = t645 - t385;
t368 = -pkin(5) * t617 - t407 * t511 + t435 * t508;
t367 = qJ(6) * t617 + t628;
t358 = qJD(5) * t404 + t403 * t508 + t511 * t587;
t357 = -qJD(5) * t405 - t403 * t511 + t508 * t587;
t1 = [t631 * MDP(1) + (-t364 * t449 + t365 * t448 - t402 * t411 + t403 * t410 - g(3)) * MDP(15) + (t357 * t458 - t358 * t456 + t385 * t404 + t386 * t405) * MDP(24) + (-t346 * t405 - t347 * t404 + t348 * t449 + t352 * t357 + t353 * t358 + t366 * t402 - g(3)) * MDP(26) + (t448 * t509 + t449 * t512) * t596 + (t402 * t512 + t403 * t509 + (t448 * t512 - t449 * t509) * qJD(3)) * MDP(12) * qJD(2) + ((qJD(2) * t420 * MDP(15) - qJDD(2) * MDP(4) + (t509 * t673 + t512 * t674 - MDP(3)) * t516) * t510 + (-t371 * MDP(15) + qJDD(2) * MDP(3) - t516 * MDP(4) - t673 * t543 + t674 * t681) * t513) * t506 - t673 * (qJD(3) * t402 + qJDD(3) * t449) + t674 * (qJD(3) * t403 + qJDD(3) * t448) + t672 * (-t357 * t491 + t449 * t386 + t402 * t456 + t404 * t455) + t597 * (t358 * t491 + t385 * t449 - t402 * t458 - t405 * t455); qJDD(2) * MDP(2) + (t631 * t639 + t568) * MDP(3) + (-t631 * t641 + t567) * MDP(4) + (qJDD(2) * t503 + 0.2e1 * t509 * t579) * MDP(5) + 0.2e1 * (t509 * t598 - t602 * t624) * MDP(6) + (qJDD(3) * t509 + t512 * t515) * MDP(7) + (qJDD(3) * t512 - t509 * t515) * MDP(8) + (t524 * t509 + t512 * t669) * MDP(10) + (-t509 * t669 + t524 * t512) * MDP(11) + (t623 * t652 + (-g(3) * t510 - t582 * t623) * t506 + t517) * MDP(12) + (t523 * t509 - t512 * t668) * MDP(13) + (t509 * t668 + t523 * t512) * MDP(14) + (t371 * t468 + t420 * t440 - g(1) * t590 - g(2) * t591 - g(3) * t571 + (-t420 * t510 + (-t410 * t509 + t411 * t512) * t513) * t622 + t517 * pkin(8)) * MDP(15) + (t385 * t508 * t512 + (t508 * t613 - t512 * t607) * t458) * MDP(16) + ((-t456 * t508 + t458 * t511) * t613 + (t649 + t386 * t508 + (t456 * t511 + t458 * t508) * qJD(5)) * t512) * MDP(17) + ((t491 * t614 - t385) * t509 + (qJD(3) * t458 + t545) * t512) * MDP(18) + ((t491 * t612 - t386) * t509 + (-qJD(3) * t456 - t676) * t512) * MDP(19) + (t455 * t509 + t491 * t611) * MDP(20) + (t555 * t455 - t462 * t456 + t477 * t386 + (-t393 * t612 + t573) * t509 + t667 * t491 + (qJD(3) * t359 + t356 * t511 - t393 * t608 - t570) * t512 + t538) * MDP(21) + (-t677 * t455 - t462 * t458 - t477 * t385 + ((qJD(3) * t393 + qJD(5) * t408) * t508 + t592) * t509 + t678 * t491 + (-qJD(3) * t360 - t356 * t508 - t393 * t607 - t569) * t512 - t539) * MDP(22) + (t376 * t456 + t386 * t417 - t389 * t455 + (-t366 * t612 - t347) * t509 - t629 * t491 + (-qJD(3) * t352 + t348 * t511 - t366 * t608 - t570) * t512 + t538) * MDP(23) + (-t385 * t389 - t386 * t388 + t629 * t458 - t630 * t456 + t558 * t613 + (-t346 * t511 - t347 * t508 + (-t352 * t511 + t353 * t508) * qJD(5) + t533) * t512) * MDP(24) + (-t376 * t458 + t385 * t417 + t388 * t455 + (-t366 * t614 + t346) * t509 + t630 * t491 + (qJD(3) * t353 + t348 * t508 + t366 * t607 + t569) * t512 + t539) * MDP(25) + (t346 * t388 + t348 * t417 + t347 * t389 - g(1) * (pkin(5) * t381 - pkin(9) * t647 + qJ(6) * t380 + t445 * t665 + t590) - g(2) * (pkin(5) * t379 - pkin(9) * t648 + qJ(6) * t378 + t443 * t665 + t591) - g(3) * (pkin(5) * t425 + qJ(6) * t424 + (pkin(4) * t510 + pkin(9) * t633) * t506 + t571) + (-t512 * t588 + t376) * t366 + t630 * t353 + t629 * t352) * MDP(26); -MDP(5) * t593 + t624 * t516 * MDP(6) + MDP(7) * t599 + MDP(8) * t598 + qJDD(3) * MDP(9) + (-t466 * t618 + t520 + t616) * MDP(10) + (-t466 * t617 - t521) * MDP(11) + (-pkin(3) * t509 + t655) * t596 + (-0.2e1 * t651 - t616 + (-qJD(2) * t461 - t601) * t512 - t537 + t670 + t671) * MDP(13) + (0.2e1 * t500 + 0.2e1 * t501 + (t420 * t512 + t461 * t509) * qJD(2) + t521) * MDP(14) + (-t364 * qJ(4) - t365 * pkin(3) - t420 * t461 - t410 * t419 - g(1) * (qJ(4) * t401 - t396) - g(2) * (qJ(4) * t399 - t395) - g(3) * (qJ(4) * t449 - t439) - t680 * t411) * MDP(15) + (-t508 * t644 - t649) * MDP(16) + ((-t386 - t644) * t511 + (t385 + t645) * t508) * MDP(17) + ((-t458 * t512 - t491 * t638) * qJD(2) + t676) * MDP(18) + ((t456 * t512 - t491 * t636) * qJD(2) + t545) * MDP(19) - t491 * MDP(20) * t617 + (-t359 * t617 + qJ(4) * t386 + t605 * t456 + (t632 + (t393 - t407) * t491) * t511 + (t356 + (t435 - t606) * t491 - t536) * t508) * MDP(21) + (-qJ(4) * t385 + t628 * t491 + t360 * t617 + t605 * t458 + (-t393 * t491 - t632) * t508 + (t356 - t527) * t511) * MDP(22) + (t352 * t617 + t368 * t491 + t386 * t467 + t456 * t627 - t508 * t525 + t511 * t531) * MDP(23) + (t367 * t456 - t368 * t458 + (t385 * t514 + (-t456 * t514 - t353) * qJD(5) + t564) * t511 + (-t386 * t514 + (t458 * t514 - t352) * qJD(5) + t565) * t508 + t537) * MDP(24) + (-t353 * t617 - t367 * t491 + t385 * t467 - t458 * t627 + t508 * t531 + t511 * t525) * MDP(25) + (-t353 * t367 - t352 * t368 - g(1) * (-pkin(9) * t400 - t396) - g(2) * (-pkin(9) * t398 - t395) - g(3) * (-pkin(9) * t448 - t439) + t627 * t366 + (qJD(5) * t558 + t346 * t508 - t347 * t511) * t514 + (t348 - t536) * t467) * MDP(26); t509 * t596 + (qJDD(3) + t593) * MDP(13) + (-t503 * t516 - t515) * MDP(14) + (-t520 - t651 + t671) * MDP(15) + t441 * MDP(21) - t537 * MDP(26) + (t411 * MDP(15) - t366 * MDP(26) - t456 * t672 + t597 * t458) * qJD(3) + (t455 * MDP(23) + (-t456 * t618 + t385 - t609) * MDP(24) + (qJD(5) * t353 - t564) * MDP(26) + t597 * t679) * t511 + (t682 * MDP(24) + (qJD(5) * t352 - t565) * MDP(26) + t597 * t455 - t672 * t679) * t508; MDP(16) * t646 + (-t456 ^ 2 + t666) * MDP(17) + t369 * MDP(18) + t682 * MDP(19) + t455 * MDP(20) + (-t393 * t458 + t526 + t650) * MDP(21) + (t359 * t491 + t393 * t456 - t522) * MDP(22) + (-t397 * t456 - t518 + t650 + 0.2e1 * t664) * MDP(23) + (pkin(5) * t385 - qJ(6) * t386 + (t353 - t360) * t458 + (t352 - t604) * t456) * MDP(24) + (0.2e1 * t654 - t366 * t456 + t397 * t458 + (0.2e1 * qJD(6) - t359) * t491 + t522) * MDP(25) + (t346 * qJ(6) - t347 * pkin(5) - t366 * t397 - t352 * t360 - g(1) * (-pkin(5) * t374 + qJ(6) * t375) - g(2) * (-pkin(5) * t372 + qJ(6) * t373) - g(3) * (pkin(5) * t404 - qJ(6) * t405) + t604 * t353) * MDP(26); (-t455 + t646) * MDP(23) + t369 * MDP(24) + (-t666 - t679) * MDP(25) + (-t353 * t491 + t518 - t664) * MDP(26);];
tau  = t1;

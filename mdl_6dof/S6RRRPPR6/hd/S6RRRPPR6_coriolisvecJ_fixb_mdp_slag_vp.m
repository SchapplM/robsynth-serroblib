% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:54:11
% EndTime: 2019-03-09 15:54:26
% DurationCPUTime: 9.54s
% Computational Cost: add. (8768->558), mult. (23472->743), div. (0->0), fcn. (18331->10), ass. (0->234)
t558 = sin(qJ(3));
t554 = sin(pkin(6));
t562 = cos(qJ(2));
t638 = qJD(1) * t562;
t618 = t554 * t638;
t599 = t558 * t618;
t634 = qJD(3) * t558;
t688 = t599 - t634;
t559 = sin(qJ(2));
t561 = cos(qJ(3));
t640 = qJD(1) * t554;
t619 = t559 * t640;
t556 = cos(pkin(6));
t639 = qJD(1) * t556;
t624 = pkin(1) * t639;
t506 = -pkin(8) * t619 + t562 * t624;
t580 = t554 * (pkin(2) * t559 - pkin(9) * t562);
t507 = qJD(1) * t580;
t606 = -t506 * t558 + t561 * t507;
t669 = -qJ(4) - pkin(9);
t611 = qJD(3) * t669;
t652 = t561 * t562;
t687 = -(pkin(3) * t559 - qJ(4) * t652) * t640 - t606 - qJD(4) * t558 + t561 * t611;
t643 = t561 * t506 + t558 * t507;
t686 = -qJ(4) * t599 - qJD(4) * t561 - t558 * t611 + t643;
t532 = -qJD(3) + t618;
t553 = sin(pkin(11));
t555 = cos(pkin(11));
t633 = qJD(3) * t561;
t654 = t555 * t561;
t642 = t688 * t553 + t555 * t633 - t618 * t654;
t541 = qJD(2) + t639;
t600 = t558 * t619;
t489 = -t561 * t541 + t600;
t491 = t541 * t558 + t561 * t619;
t450 = t555 * t489 + t491 * t553;
t685 = pkin(5) * t450;
t560 = cos(qJ(6));
t557 = sin(qJ(6));
t661 = t532 * t557;
t424 = -t560 * t450 - t661;
t583 = -t489 * t553 + t555 * t491;
t674 = qJD(6) + t583;
t684 = t424 * t674;
t426 = t450 * t557 - t532 * t560;
t683 = t426 * t674;
t682 = t450 * t532;
t646 = t553 * t687 - t686 * t555;
t519 = t553 * t561 + t555 * t558;
t645 = t532 * t519;
t537 = t559 * t624;
t509 = pkin(8) * t618 + t537;
t681 = -t688 * pkin(3) - t509;
t680 = t583 ^ 2;
t550 = t554 ^ 2;
t625 = qJD(1) * qJD(2);
t679 = -0.2e1 * t550 * t625;
t678 = pkin(5) * t583;
t677 = MDP(5) * (t559 ^ 2 - t562 ^ 2);
t477 = -pkin(2) * t541 - t506;
t448 = pkin(3) * t489 + qJD(4) + t477;
t564 = -qJ(5) * t583 + t448;
t388 = pkin(4) * t450 + t564;
t676 = t388 * t583;
t605 = t674 * t557;
t651 = qJ(5) * t619 - t646;
t647 = t686 * t553 + t555 * t687;
t656 = t554 * t562;
t671 = pkin(1) * t559;
t503 = pkin(8) * t656 + (pkin(9) + t671) * t556;
t504 = (-pkin(2) * t562 - pkin(9) * t559 - pkin(1)) * t554;
t644 = t561 * t503 + t558 * t504;
t675 = t642 * qJ(5) + qJD(5) * t519 - t681;
t478 = pkin(9) * t541 + t509;
t485 = qJD(1) * t504;
t440 = t478 * t561 + t485 * t558;
t417 = -qJ(4) * t489 + t440;
t410 = t553 * t417;
t439 = -t478 * t558 + t561 * t485;
t416 = -qJ(4) * t491 + t439;
t382 = t416 * t555 - t410;
t629 = -qJD(5) + t382;
t657 = t554 * t559;
t516 = -t556 * t561 + t558 * t657;
t635 = qJD(2) * t562;
t616 = t554 * t635;
t470 = -qJD(3) * t516 + t561 * t616;
t517 = t556 * t558 + t561 * t657;
t508 = qJD(2) * t580;
t542 = pkin(8) * t657;
t670 = pkin(1) * t562;
t510 = (t556 * t670 - t542) * qJD(2);
t567 = -qJD(3) * t644 + t561 * t508 - t510 * t558;
t637 = qJD(2) * t559;
t617 = t554 * t637;
t392 = pkin(3) * t617 - qJ(4) * t470 - qJD(4) * t517 + t567;
t615 = t558 * t635;
t469 = qJD(3) * t517 + t554 * t615;
t576 = -t503 * t634 + t504 * t633 + t558 * t508 + t561 * t510;
t397 = -qJ(4) * t469 - qJD(4) * t516 + t576;
t360 = t553 * t392 + t555 * t397;
t357 = -t554 * (qJ(5) * t637 - qJD(5) * t562) - t360;
t408 = -pkin(3) * t532 + t416;
t655 = t555 * t417;
t376 = t553 * t408 + t655;
t371 = qJ(5) * t532 - t376;
t364 = -t371 - t685;
t381 = t416 * t553 + t655;
t612 = t554 * t625;
t597 = t562 * t612;
t464 = -qJD(3) * t600 + t541 * t633 + t561 * t597;
t465 = (t559 * t633 + t615) * t640 + t541 * t634;
t419 = t464 * t555 - t465 * t553;
t548 = -pkin(3) * t555 - pkin(4);
t544 = -pkin(10) + t548;
t673 = t544 * t419 + (t364 - t381 + t685) * t674;
t672 = pkin(4) + pkin(10);
t668 = t381 * t583;
t418 = t464 * t553 + t555 * t465;
t598 = t559 * t612;
t631 = qJD(6) * t560;
t622 = t557 * t418 + t450 * t631 + t560 * t598;
t632 = qJD(6) * t557;
t384 = t532 * t632 + t622;
t667 = t384 * t560;
t666 = t424 * t450;
t665 = t426 * t450;
t664 = t489 * t532;
t663 = t491 * t532;
t518 = t553 * t558 - t654;
t662 = t518 * t557;
t660 = t532 * t558;
t659 = t532 * t561;
t563 = qJD(1) ^ 2;
t658 = t550 * t563;
t653 = t557 * t419;
t415 = t560 * t419;
t498 = qJD(1) * t508;
t499 = qJD(1) * t510;
t566 = -qJD(3) * t440 + t561 * t498 - t558 * t499;
t377 = pkin(3) * t598 - qJ(4) * t464 - qJD(4) * t491 + t566;
t577 = -t478 * t634 + t485 * t633 + t558 * t498 + t561 * t499;
t383 = -qJ(4) * t465 - qJD(4) * t489 + t577;
t355 = t553 * t377 + t555 * t383;
t650 = -pkin(4) * t619 + t647;
t649 = pkin(4) * t645 + t675;
t607 = -t503 * t558 + t561 * t504;
t422 = -pkin(3) * t656 - qJ(4) * t517 + t607;
t432 = -qJ(4) * t516 + t644;
t396 = t553 * t422 + t555 * t432;
t648 = pkin(5) * t645 - t651;
t500 = pkin(8) * t597 + qJD(2) * t537;
t511 = t556 * pkin(1) * t637 + pkin(8) * t616;
t636 = qJD(2) * t561;
t630 = qJD(2) - t541;
t628 = -t629 + t678;
t623 = t560 * t656;
t620 = -pkin(3) * t561 - pkin(2);
t354 = t377 * t555 - t553 * t383;
t602 = t672 * t657;
t590 = qJD(2) * t602;
t349 = pkin(5) * t419 - qJD(1) * t590 - t354;
t436 = pkin(3) * t465 + t500;
t569 = -qJ(5) * t419 - qJD(5) * t583 + t436;
t356 = t418 * t672 + t569;
t610 = t560 * t349 - t356 * t557;
t359 = t392 * t555 - t553 * t397;
t375 = t408 * t555 - t410;
t395 = t422 * t555 - t553 * t432;
t609 = t557 * t619 - t560 * t645;
t608 = t557 * t645 + t560 * t619;
t533 = t669 * t558;
t534 = t669 * t561;
t475 = -t555 * t533 - t534 * t553;
t604 = t674 * t560;
t601 = t550 * t562 * t559 * MDP(4);
t352 = -qJ(5) * t598 + qJD(5) * t532 - t355;
t596 = pkin(3) * t469 + t511;
t595 = pkin(1) * t679;
t391 = pkin(4) * t656 - t395;
t594 = qJD(5) - t375;
t582 = -qJ(5) * t519 + t620;
t454 = t518 * t672 + t582;
t593 = -pkin(5) * t642 - qJD(1) * t602 + qJD(6) * t454 + t647;
t457 = pkin(5) * t519 + t475;
t592 = -qJD(6) * t457 + t645 * t672 + t675;
t591 = pkin(3) * t491 + qJ(5) * t450;
t589 = t349 * t557 + t356 * t560;
t363 = t532 * t672 + t594 + t678;
t368 = t450 * t672 + t564;
t345 = t363 * t560 - t368 * t557;
t346 = t363 * t557 + t368 * t560;
t462 = -t516 * t553 + t517 * t555;
t369 = pkin(5) * t462 + pkin(10) * t656 + t391;
t461 = t555 * t516 + t517 * t553;
t502 = t542 + (-pkin(2) - t670) * t556;
t570 = pkin(3) * t516 + t502;
t565 = -qJ(5) * t462 + t570;
t387 = t461 * t672 + t565;
t587 = t369 * t560 - t387 * t557;
t586 = t369 * t557 + t387 * t560;
t476 = t533 * t553 - t534 * t555;
t585 = -t418 * t476 + t419 * t475;
t390 = qJ(5) * t656 - t396;
t445 = t461 * t560 + t557 * t656;
t348 = -pkin(5) * t418 - t352;
t575 = t348 + (-qJD(6) * t544 + t583 * t672 + t591) * t674;
t574 = -t560 * t418 + t557 * t598;
t573 = t518 * t632 - t609;
t572 = t518 * t631 - t608;
t353 = -pkin(4) * t598 - t354;
t431 = -t469 * t553 + t470 * t555;
t568 = -qJ(5) * t431 - qJD(5) * t462 + t596;
t362 = pkin(4) * t418 + t569;
t545 = pkin(3) * t553 + qJ(5);
t466 = pkin(4) * t518 + t582;
t458 = -pkin(5) * t518 + t476;
t446 = t461 * t557 - t623;
t430 = t555 * t469 + t470 * t553;
t405 = pkin(4) * t461 + t565;
t402 = pkin(4) * t583 + t591;
t399 = qJD(6) * t445 + t430 * t557 + t560 * t617;
t398 = -t560 * t430 - qJD(6) * t623 + (qJD(6) * t461 + t617) * t557;
t385 = qJD(6) * t426 + t574;
t372 = -pkin(5) * t461 - t390;
t370 = pkin(4) * t532 + t594;
t365 = pkin(4) * t430 + t568;
t361 = t430 * t672 + t568;
t358 = -pkin(4) * t617 - t359;
t351 = -pkin(5) * t430 - t357;
t350 = pkin(5) * t431 - t359 - t590;
t344 = -qJD(6) * t346 + t610;
t343 = qJD(6) * t345 + t589;
t1 = [(-(qJD(6) * t587 + t350 * t557 + t361 * t560) * t674 - t586 * t419 - t343 * t462 - t346 * t431 + t351 * t426 + t372 * t384 + t348 * t446 + t364 * t399) * MDP(30) + ((-qJD(6) * t586 + t350 * t560 - t361 * t557) * t674 + t587 * t419 + t344 * t462 + t345 * t431 + t351 * t424 + t372 * t385 - t348 * t445 + t364 * t398) * MDP(29) + (-t385 * t462 - t398 * t674 + t419 * t445 - t424 * t431) * MDP(27) + (t384 * t462 + t399 * t674 + t419 * t446 + t426 * t431) * MDP(26) + (t419 * t462 + t431 * t674) * MDP(28) + t677 * t679 + (-t499 * t556 - t510 * t541 + t562 * t595) * MDP(10) + (-t500 * t556 - t511 * t541 + t559 * t595) * MDP(9) + (t354 * t395 + t355 * t396 + t375 * t359 + t376 * t360 + t436 * t570 + t448 * t596) * MDP(19) + (-t532 * t554 - t550 * t638) * MDP(15) * t637 + (-t464 * t516 - t465 * t517 - t469 * t491 - t470 * t489) * MDP(12) + (t464 * t517 + t470 * t491) * MDP(11) + (t352 * t461 + t353 * t462 + t357 * t450 + t358 * t583 + t370 * t431 + t371 * t430 + t390 * t418 + t391 * t419) * MDP(20) + (-t354 * t462 - t355 * t461 - t359 * t583 - t360 * t450 - t375 * t431 - t376 * t430 - t395 * t419 - t396 * t418) * MDP(18) + (t357 * t532 - t362 * t462 - t365 * t583 - t388 * t431 - t405 * t419 + (t352 * t562 + (-qJD(1) * t390 - t371) * t637) * t554) * MDP(22) + (MDP(6) * t616 - MDP(7) * t617) * (t541 + t639) + (t384 * t445 - t385 * t446 - t398 * t426 - t399 * t424) * MDP(25) + (t384 * t446 + t399 * t426) * MDP(24) + (t352 * t390 + t353 * t391 + t357 * t371 + t358 * t370 + t362 * t405 + t365 * t388) * MDP(23) + 0.2e1 * t601 * t625 + (-t470 * t532 + (-t464 * t562 + (qJD(1) * t517 + t491) * t637) * t554) * MDP(13) + (t469 * t532 + (t465 * t562 + (-qJD(1) * t516 - t489) * t637) * t554) * MDP(14) + (-t567 * t532 + t511 * t489 + t502 * t465 + t500 * t516 + t477 * t469 + (-t566 * t562 + (qJD(1) * t607 + t439) * t637) * t554) * MDP(16) + (-t358 * t532 - t362 * t461 - t365 * t450 - t388 * t430 - t405 * t418 + (-t353 * t562 + (qJD(1) * t391 + t370) * t637) * t554) * MDP(21) + (t576 * t532 + t511 * t491 + t502 * t464 + t500 * t517 + t477 * t470 + (t577 * t562 + (-qJD(1) * t644 - t440) * t637) * t554) * MDP(17); (t419 * t519 + t642 * t674) * MDP(28) + ((-t454 * t557 + t457 * t560) * t419 + t344 * t519 + t458 * t385 - t348 * t560 * t518 + (t557 * t592 - t560 * t593) * t674 + t648 * t424 + t642 * t345 + t573 * t364) * MDP(29) + (-t385 * t519 + t415 * t518 - t424 * t642 - t573 * t674) * MDP(27) + (t384 * t519 + t426 * t642 + t518 * t653 + t572 * t674) * MDP(26) + (-(t454 * t560 + t457 * t557) * t419 - t343 * t519 + t458 * t384 + t348 * t662 + (t557 * t593 + t560 * t592) * t674 + t648 * t426 - t642 * t346 + t572 * t364) * MDP(30) + (-t354 * t475 + t355 * t476 + t647 * t375 + t646 * t376 + t436 * t620 + t448 * t681) * MDP(19) + (-t362 * t518 + t645 * t388 - t418 * t466 + t649 * t450 + t650 * t532) * MDP(21) + (t532 * MDP(15) - t630 * MDP(7) + (qJD(2) * t475 - t370) * MDP(21) + (qJD(2) * t476 + t371) * MDP(22)) * t619 + t658 * t677 + t630 * MDP(6) * t618 - t563 * t601 + (-t362 * t519 - t642 * t388 - t419 * t466 + t651 * t532 + t583 * t649) * MDP(22) + (-t354 * t519 - t355 * t518 - t375 * t642 + t376 * t645 - t450 * t646 - t583 * t647 + t585) * MDP(18) + (t352 * t518 + t353 * t519 + t370 * t642 - t371 * t645 + t450 * t651 - t583 * t650 + t585) * MDP(20) + (-t352 * t476 + t353 * t475 + t362 * t466 - t370 * t650 + t371 * t651 - t388 * t649) * MDP(23) + (-t532 * t633 + (t532 * t652 + (qJD(2) * t558 - t491) * t559) * t640) * MDP(13) + (t464 * t558 - t491 * t659) * MDP(11) + (-pkin(2) * t465 - t500 * t561 + t606 * t532 - t509 * t489 + (pkin(9) * t659 + t477 * t558) * qJD(3) + (-t439 * t559 + (-pkin(9) * t637 - t477 * t562) * t558) * t640) * MDP(16) + (t532 * t634 + (-t562 * t660 + (t489 + t636) * t559) * t640) * MDP(14) + (-pkin(2) * t464 + t500 * t558 - t643 * t532 - t509 * t491 + (-pkin(9) * t660 + t477 * t561) * qJD(3) + (-t477 * t652 + (-pkin(9) * t636 + t440) * t559) * t640) * MDP(17) + (t384 * t662 + t426 * t572) * MDP(24) + ((t464 + t664) * t561 + (-t465 + t663) * t558) * MDP(12) + (t609 * t426 + t608 * t424 + (t667 - t385 * t557 + (-t424 * t560 - t426 * t557) * qJD(6)) * t518) * MDP(25) + (pkin(8) * t598 + t506 * t541 + (-t556 * t625 + t658) * t670) * MDP(10) + (t509 * t541 + t658 * t671 - t500) * MDP(9); t491 * t489 * MDP(11) + (-t489 ^ 2 + t491 ^ 2) * MDP(12) + (t464 - t664) * MDP(13) + (-t465 - t663) * MDP(14) + MDP(15) * t598 + (-t440 * t532 - t477 * t491 + t566) * MDP(16) + (-t439 * t532 + t477 * t489 - t577) * MDP(17) + (t376 * t583 - t668 + (-t418 * t553 - t419 * t555) * pkin(3) + (-t375 + t382) * t450) * MDP(18) + (t375 * t381 - t376 * t382 + (t354 * t555 + t355 * t553 - t448 * t491) * pkin(3)) * MDP(19) + (-t371 * t583 - t418 * t545 + t419 * t548 - t668 + (t370 + t629) * t450) * MDP(20) + (t381 * t532 + t676 + t402 * t450 + (-pkin(4) + t548) * t598 - t354) * MDP(21) + (-t388 * t450 + t402 * t583 + t532 * t629 + t545 * t598 - t352) * MDP(22) + (-t352 * t545 + t353 * t548 - t370 * t381 + t371 * t629 - t388 * t402) * MDP(23) + (-t426 * t605 + t667) * MDP(24) + ((-t385 - t683) * t560 + (-t384 + t684) * t557) * MDP(25) + (-t605 * t674 + t415 + t665) * MDP(26) + (-t604 * t674 - t653 - t666) * MDP(27) + t674 * t450 * MDP(28) + (t345 * t450 + t545 * t385 + t628 * t424 + t575 * t557 + t560 * t673) * MDP(29) + (-t346 * t450 + t545 * t384 + t628 * t426 - t557 * t673 + t575 * t560) * MDP(30); (t375 * t583 + t376 * t450 + t436) * MDP(19) + (t532 * t583 - t418) * MDP(21) + (-t419 - t682) * MDP(22) + (-t370 * t583 - t371 * t450 + t362) * MDP(23) + (-t653 + t666) * MDP(29) + (-t415 + t665) * MDP(30) + (-MDP(29) * t604 + MDP(30) * t605) * t674 + (MDP(18) + MDP(20)) * (-t450 ^ 2 - t680); (t419 - t682) * MDP(20) + (-t450 * t583 + t598) * MDP(21) + (-t532 ^ 2 - t680) * MDP(22) + (-t371 * t532 + t353 + t676) * MDP(23) + (t424 * t532 + t415) * MDP(29) + (t426 * t532 - t653) * MDP(30) - (MDP(29) * t557 + MDP(30) * t560) * t674 ^ 2; t426 * t424 * MDP(24) + (-t424 ^ 2 + t426 ^ 2) * MDP(25) + (t622 + t684) * MDP(26) + (-t574 + t683) * MDP(27) + t419 * MDP(28) + (t346 * t674 - t364 * t426 + t610) * MDP(29) + (t345 * t674 + t364 * t424 - t589) * MDP(30) + (MDP(26) * t661 - MDP(27) * t426 - MDP(29) * t346 - MDP(30) * t345) * qJD(6);];
tauc  = t1;

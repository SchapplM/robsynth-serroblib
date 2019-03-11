% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:57
% EndTime: 2019-03-08 22:28:13
% DurationCPUTime: 10.50s
% Computational Cost: add. (6685->525), mult. (18443->765), div. (0->0), fcn. (15653->14), ass. (0->223)
t542 = cos(pkin(7));
t550 = cos(qJ(3));
t551 = cos(qJ(2));
t649 = t550 * t551;
t546 = sin(qJ(3));
t547 = sin(qJ(2));
t653 = t546 * t547;
t568 = -t542 * t653 + t649;
t628 = qJD(3) * t550;
t608 = t542 * t628;
t540 = sin(pkin(6));
t636 = qJD(1) * t540;
t681 = -pkin(2) * t608 + t568 * t636;
t539 = sin(pkin(7));
t583 = pkin(3) * t546 - qJ(4) * t550;
t556 = qJD(3) * t583 - qJD(4) * t546;
t616 = t547 * t636;
t680 = (t556 - t616) * t539;
t629 = qJD(3) * t546;
t611 = t539 * t629;
t679 = pkin(9) * t611 - qJD(4) * t542 + t681;
t538 = sin(pkin(13));
t541 = cos(pkin(13));
t645 = t679 * t538 + t680 * t541;
t644 = t680 * t538 - t679 * t541;
t564 = (-pkin(10) * t541 * t550 + pkin(4) * t546) * t539;
t559 = qJD(3) * t564;
t678 = t559 + t645;
t610 = t539 * t628;
t593 = t538 * t610;
t677 = -pkin(10) * t593 + t644;
t633 = qJD(2) * t542;
t530 = qJD(3) + t633;
t634 = qJD(2) * t539;
t613 = t546 * t634;
t479 = t530 * t541 - t538 * t613;
t480 = t530 * t538 + t541 * t613;
t545 = sin(qJ(5));
t549 = cos(qJ(5));
t433 = -t549 * t479 + t480 * t545;
t431 = qJD(6) + t433;
t676 = t431 ^ 2;
t632 = qJD(2) * t550;
t612 = t539 * t632;
t525 = -qJD(5) + t612;
t675 = t433 * t525;
t651 = t547 * t550;
t652 = t546 * t551;
t570 = t542 * t651 + t652;
t609 = t542 * t629;
t640 = -pkin(2) * t609 - pkin(9) * t610 + t570 * t636;
t650 = t549 * t541;
t505 = t538 * t545 - t650;
t655 = t539 * t550;
t560 = t505 * t655;
t639 = qJD(2) * t560 - t505 * qJD(5);
t506 = t538 * t549 + t541 * t545;
t561 = t506 * t655;
t638 = -qJD(2) * t561 + t506 * qJD(5);
t615 = t551 * t636;
t543 = cos(pkin(6));
t635 = qJD(1) * t543;
t617 = t539 * t635;
t674 = qJD(2) * t615 + qJD(3) * t617;
t664 = qJD(2) * pkin(2);
t517 = t615 + t664;
t673 = t517 * t542 + t617;
t484 = pkin(9) * t655 + (pkin(2) * t546 + qJ(4)) * t542;
t584 = -pkin(3) * t550 - qJ(4) * t546;
t485 = (-pkin(2) + t584) * t539;
t438 = -t484 * t538 + t541 * t485;
t656 = t539 * t546;
t498 = t538 * t542 + t541 * t656;
t407 = -pkin(4) * t655 - pkin(10) * t498 + t438;
t439 = t541 * t484 + t538 * t485;
t496 = t538 * t656 - t541 * t542;
t418 = -pkin(10) * t496 + t439;
t646 = t545 * t407 + t549 * t418;
t672 = qJD(5) * t646 + t677 * t545 - t549 * t678;
t626 = qJD(5) * t549;
t627 = qJD(5) * t545;
t671 = t407 * t626 - t418 * t627 + t545 * t678 + t677 * t549;
t504 = pkin(9) * t634 + t616;
t494 = t546 * t504;
t440 = t550 * t673 - t494;
t491 = t583 * t634;
t408 = -t440 * t538 + t541 * t491;
t387 = qJD(2) * t564 + t408;
t409 = t541 * t440 + t538 * t491;
t595 = t538 * t612;
t396 = -pkin(10) * t595 + t409;
t665 = pkin(10) + qJ(4);
t523 = t665 * t538;
t524 = t665 * t541;
t574 = -t523 * t549 - t524 * t545;
t670 = -qJD(4) * t505 + qJD(5) * t574 - t545 * t387 - t549 * t396;
t465 = -t523 * t545 + t524 * t549;
t669 = qJD(4) * t506 + qJD(5) * t465 + t387 * t549 - t396 * t545;
t642 = pkin(4) * t593 - t640;
t575 = t479 * t545 + t480 * t549;
t668 = qJD(5) * t575;
t441 = t550 * t504 + t546 * t673;
t426 = qJ(4) * t530 + t441;
t529 = t542 * t635;
t453 = t529 + (qJD(2) * t584 - t517) * t539;
t384 = -t426 * t538 + t541 * t453;
t375 = -pkin(4) * t612 - pkin(10) * t480 + t384;
t385 = t541 * t426 + t538 * t453;
t378 = pkin(10) * t479 + t385;
t356 = t375 * t545 + t378 * t549;
t596 = t542 * t616;
t573 = qJD(2) * t596;
t618 = t517 * t608 + t550 * t674;
t395 = qJD(4) * t530 + (-qJD(3) * t504 - t573) * t546 + t618;
t450 = (t556 + t616) * t634;
t376 = -t395 * t538 + t541 * t450;
t367 = qJD(2) * t559 + t376;
t377 = t541 * t395 + t538 * t450;
t621 = qJD(2) * qJD(3);
t607 = t539 * t621;
t589 = t550 * t607;
t572 = t538 * t589;
t371 = -pkin(10) * t572 + t377;
t554 = -qJD(5) * t356 + t549 * t367 - t545 * t371;
t590 = t546 * t607;
t346 = -pkin(5) * t590 - t554;
t667 = t431 * (pkin(5) * t575 + pkin(11) * t431) + t346;
t535 = t539 ^ 2;
t666 = (-t546 * t550 * MDP(5) + (t546 ^ 2 - t550 ^ 2) * MDP(6)) * t535;
t641 = t479 * t626 + t589 * t650;
t397 = (-qJD(5) * t480 - t572) * t545 + t641;
t544 = sin(qJ(6));
t548 = cos(qJ(6));
t624 = qJD(6) * t548;
t619 = t548 * t397 - t525 * t624 + t544 * t590;
t625 = qJD(6) * t544;
t363 = -t575 * t625 + t619;
t662 = t363 * t544;
t659 = t575 * t544;
t414 = t548 * t525 + t659;
t661 = t414 * t431;
t416 = -t525 * t544 + t548 * t575;
t660 = t416 * t431;
t658 = t506 * t548;
t558 = qJD(3) * t561;
t398 = qJD(2) * t558 + t668;
t654 = t544 * t398;
t394 = t548 * t398;
t648 = -pkin(5) * t611 + t672;
t643 = pkin(5) * t613 + t669;
t631 = qJD(3) * t538;
t630 = qJD(3) * t541;
t623 = qJD(3) - t530;
t620 = t544 * t655;
t534 = -pkin(4) * t541 - pkin(3);
t614 = t535 * t632;
t566 = t545 * t367 + t549 * t371 + t375 * t626 - t378 * t627;
t345 = pkin(11) * t590 + t566;
t402 = t504 * t628 + t517 * t609 + t546 * t674 + t550 * t573;
t390 = pkin(4) * t572 + t402;
t358 = pkin(5) * t398 - pkin(11) * t397 + t390;
t605 = -t345 * t544 + t548 * t358;
t603 = t397 * t544 - t548 * t590;
t602 = -t544 * t639 - t548 * t613;
t601 = t544 * t613 - t548 * t639;
t600 = t548 * t431;
t422 = pkin(4) * t595 + t441;
t594 = t540 * t547 * t634;
t454 = t549 * t496 + t498 * t545;
t455 = -t496 * t545 + t498 * t549;
t487 = pkin(9) * t656 + (-pkin(2) * t550 - pkin(3)) * t542;
t456 = pkin(4) * t496 + t487;
t379 = pkin(5) * t454 - pkin(11) * t455 + t456;
t588 = -pkin(11) * t611 - qJD(6) * t379 - t671;
t458 = pkin(5) * t505 - pkin(11) * t506 + t534;
t587 = pkin(11) * t613 - qJD(6) * t458 - t670;
t369 = -pkin(11) * t655 + t646;
t411 = -qJD(3) * t560 - qJD(5) * t454;
t412 = qJD(5) * t455 + t558;
t586 = -pkin(5) * t412 + pkin(11) * t411 + qJD(6) * t369 - t642;
t585 = -pkin(5) * t638 + pkin(11) * t639 + qJD(6) * t465 + t422;
t582 = t345 * t548 + t358 * t544;
t352 = -pkin(11) * t525 + t356;
t423 = -pkin(3) * t530 + qJD(4) - t440;
t399 = -pkin(4) * t479 + t423;
t361 = pkin(5) * t433 - pkin(11) * t575 + t399;
t348 = t352 * t548 + t361 * t544;
t581 = t352 * t544 - t361 * t548;
t355 = t375 * t549 - t378 * t545;
t569 = t542 * t652 + t651;
t460 = t540 * t569 + t543 * t656;
t497 = -t539 * t540 * t551 + t542 * t543;
t427 = -t460 * t538 + t497 * t541;
t428 = t460 * t541 + t497 * t538;
t381 = t427 * t545 + t428 * t549;
t571 = t542 * t649 - t653;
t459 = -t540 * t571 - t543 * t655;
t580 = t381 * t548 + t459 * t544;
t579 = -t381 * t544 + t459 * t548;
t577 = t407 * t549 - t418 * t545;
t576 = t427 * t549 - t428 * t545;
t429 = t455 * t544 + t548 * t655;
t563 = t506 * t624 - t602;
t562 = -t506 * t625 - t601;
t351 = pkin(5) * t525 - t355;
t557 = -pkin(11) * t398 + (t351 + t355) * t431;
t555 = -qJ(4) * t629 + (-pkin(3) * qJD(3) + qJD(4) - t423) * t550;
t552 = qJD(2) ^ 2;
t474 = -t517 * t539 + t529;
t430 = t455 * t548 - t620;
t421 = t543 * t610 + (qJD(2) * t568 + qJD(3) * t571) * t540;
t420 = t543 * t611 + (qJD(2) * t570 + qJD(3) * t569) * t540;
t405 = t421 * t541 + t538 * t594;
t404 = -t421 * t538 + t541 * t594;
t374 = -qJD(6) * t620 + t411 * t544 + t455 * t624 - t548 * t611;
t373 = -qJD(6) * t429 + t411 * t548 + t544 * t611;
t368 = pkin(5) * t655 - t577;
t364 = qJD(6) * t416 + t603;
t354 = qJD(5) * t381 - t404 * t549 + t405 * t545;
t353 = qJD(5) * t576 + t404 * t545 + t405 * t549;
t344 = -qJD(6) * t348 + t605;
t343 = -qJD(6) * t581 + t582;
t1 = [-t420 * t479 * MDP(12) + t420 * t480 * MDP(13) + (-t404 * t480 + t405 * t479) * MDP(14) + (t376 * t427 + t377 * t428 + t384 * t404 + t385 * t405 + t402 * t459 + t420 * t423) * MDP(15) + (t354 * t525 + t398 * t459 + t420 * t433) * MDP(21) + (t353 * t525 + t397 * t459 + t420 * t575) * MDP(22) + ((-qJD(6) * t580 - t353 * t544 + t420 * t548) * t431 + t579 * t398 + t354 * t414 - t576 * t364) * MDP(28) + (-(qJD(6) * t579 + t353 * t548 + t420 * t544) * t431 - t580 * t398 + t354 * t416 - t576 * t363) * MDP(29) + (-MDP(10) * t420 - MDP(11) * t421) * t530 + (-MDP(4) * t551 + (-MDP(3) + (-MDP(10) * t550 + MDP(11) * t546) * t535) * t547) * t552 * t540 + ((-MDP(12) * t404 + MDP(13) * t405) * t550 + ((t497 * MDP(11) + (-t427 * t541 - t428 * t538) * MDP(14) + (t538 * MDP(12) + t541 * MDP(13)) * t459) * t550 + (MDP(10) * t497 + MDP(12) * t427 - MDP(13) * t428 + MDP(21) * t576 - MDP(22) * t381) * t546) * qJD(3)) * t634; (-t402 * t542 + t640 * t530 + (t474 * t539 - t535 * t664) * t629) * MDP(10) + (-(-t546 * t573 + t618) * t542 + t681 * t530 + (-pkin(2) * t614 + t542 * t494 + (pkin(9) * t530 * t546 + t474 * t550) * t539) * qJD(3)) * MDP(11) + (t402 * t496 + t640 * t479 + ((qJD(2) * t438 + t384) * t629 + (t423 * t631 - t376 + (t487 * t631 - t645) * qJD(2)) * t550) * t539) * MDP(12) + (t402 * t498 - t640 * t480 + ((-qJD(2) * t439 - t385) * t629 + (t423 * t630 + t377 + (t487 * t630 + t644) * qJD(2)) * t550) * t539) * MDP(13) + (-t376 * t498 - t377 * t496 - t645 * t480 + t644 * t479 + (-t384 * t541 - t385 * t538 + (-t438 * t541 - t439 * t538) * qJD(2)) * t610) * MDP(14) + (t376 * t438 + t377 * t439 + t384 * t645 + t385 * t644 + t402 * t487 - t423 * t640) * MDP(15) + (t397 * t455 + t411 * t575) * MDP(16) + (-t397 * t454 - t398 * t455 - t411 * t433 - t412 * t575) * MDP(17) + (-t411 * t525 + (-t397 * t550 + (qJD(2) * t455 + t575) * t629) * t539) * MDP(18) + (t412 * t525 + (t398 * t550 + (-qJD(2) * t454 - t433) * t629) * t539) * MDP(19) + (-t525 * t539 - t614) * MDP(20) * t629 + (t390 * t454 + t456 * t398 + t399 * t412 + t672 * t525 + t642 * t433 + (-t554 * t550 + (qJD(2) * t577 + t355) * t629) * t539) * MDP(21) + (t390 * t455 + t456 * t397 + t399 * t411 + t671 * t525 + t642 * t575 + (t566 * t550 + (-qJD(2) * t646 - t356) * t629) * t539) * MDP(22) + (t363 * t430 + t373 * t416) * MDP(23) + (-t363 * t429 - t364 * t430 - t373 * t414 - t374 * t416) * MDP(24) + (t363 * t454 + t373 * t431 + t398 * t430 + t412 * t416) * MDP(25) + (-t364 * t454 - t374 * t431 - t398 * t429 - t412 * t414) * MDP(26) + (t398 * t454 + t412 * t431) * MDP(27) + ((-t369 * t544 + t379 * t548) * t398 + t344 * t454 - t581 * t412 + t368 * t364 + t346 * t429 + t351 * t374 + (t544 * t588 - t548 * t586) * t431 + t648 * t414) * MDP(28) + (-(t369 * t548 + t379 * t544) * t398 - t343 * t454 - t348 * t412 + t368 * t363 + t346 * t430 + t351 * t373 + (t544 * t586 + t548 * t588) * t431 + t648 * t416) * MDP(29) - 0.2e1 * t666 * t621 + (MDP(7) * t610 - MDP(8) * t611) * (t530 + t633); t623 * MDP(7) * t612 + (t441 * t530 - t402) * MDP(10) + (t504 * t629 + t440 * t530 + (-t474 * t655 + t546 * t596) * qJD(2) - t618) * MDP(11) + (-t402 * t541 + t441 * t479 + (-t384 * t546 + t408 * t550 + t538 * t555) * t634) * MDP(12) + (t402 * t538 - t441 * t480 + (t385 * t546 - t409 * t550 + t541 * t555) * t634) * MDP(13) + (t408 * t480 - t409 * t479 + (qJD(4) * t479 + t384 * t612 + t377) * t541 + (qJD(4) * t480 + t385 * t612 - t376) * t538) * MDP(14) + (-pkin(3) * t402 - t384 * t408 - t385 * t409 - t423 * t441 + (-t384 * t538 + t385 * t541) * qJD(4) + (-t376 * t538 + t377 * t541) * qJ(4)) * MDP(15) + (t397 * t506 + t575 * t639) * MDP(16) + (-t397 * t505 - t398 * t506 - t433 * t639 - t575 * t638) * MDP(17) + (t390 * t505 + t534 * t398 + t638 * t399 - t422 * t433) * MDP(21) + (t390 * t506 + t534 * t397 + t639 * t399 - t422 * t575) * MDP(22) + (t363 * t658 + t416 * t562) * MDP(23) + (t602 * t416 + t601 * t414 + (-t662 - t364 * t548 + (t414 * t544 - t416 * t548) * qJD(6)) * t506) * MDP(24) + (t363 * t505 + t394 * t506 + t416 * t638 + t431 * t562) * MDP(25) + (-t364 * t505 - t414 * t638 - t431 * t563 - t506 * t654) * MDP(26) + (t398 * t505 + t431 * t638) * MDP(27) + ((t458 * t548 - t465 * t544) * t398 + t344 * t505 - t574 * t364 + t346 * t544 * t506 + (t544 * t587 - t548 * t585) * t431 + t643 * t414 - t638 * t581 + t563 * t351) * MDP(28) + (-(t458 * t544 + t465 * t548) * t398 - t343 * t505 - t574 * t363 + t346 * t658 + (t544 * t585 + t548 * t587) * t431 + t643 * t416 - t638 * t348 + t562 * t351) * MDP(29) + (-t639 * MDP(18) + t638 * MDP(19) + MDP(21) * t669 + MDP(22) * t670) * t525 + (-t623 * MDP(8) - t474 * MDP(10) + (qJD(3) * t506 - t575) * MDP(18) + (-qJD(3) * t505 + t433) * MDP(19) + (qJD(3) * t574 - t355) * MDP(21) + (-qJD(3) * t465 + t356) * MDP(22) + t525 * MDP(20)) * t613 + t666 * t552; (-t479 ^ 2 - t480 ^ 2) * MDP(14) + (t384 * t480 - t385 * t479 + t402) * MDP(15) + (t479 * t627 + t480 * t626 - t525 * t575) * MDP(21) + (-t480 * t627 + t641 + t675) * MDP(22) + (-t414 * t575 + t394) * MDP(28) + (-t416 * t575 - t654) * MDP(29) + (-t480 * MDP(12) - t479 * MDP(13) + ((MDP(21) * t545 + MDP(13)) * t541 + (MDP(21) * t549 - MDP(22) * t545 + MDP(12)) * t538) * qJD(3)) * t612 - (MDP(28) * t544 + MDP(29) * t548) * t676; -t433 ^ 2 * MDP(17) + (t397 - t675) * MDP(18) + (-t506 * t589 - t668) * MDP(19) + MDP(20) * t590 + (-t356 * t525 + t554) * MDP(21) + (-t355 * t525 + t399 * t433 - t566) * MDP(22) + (t416 * t600 + t662) * MDP(23) + ((t363 - t661) * t548 + (-t364 - t660) * t544) * MDP(24) + (t431 * t600 + t654) * MDP(25) + (-t676 * t544 + t394) * MDP(26) + (-pkin(5) * t364 - t356 * t414 + t557 * t544 - t548 * t667) * MDP(28) + (-pkin(5) * t363 - t356 * t416 + t544 * t667 + t557 * t548) * MDP(29) + (MDP(16) * t433 + MDP(17) * t575 - t525 * MDP(19) - t399 * MDP(21) - t416 * MDP(25) + t414 * MDP(26) - t431 * MDP(27) + MDP(28) * t581 + t348 * MDP(29)) * t575; t416 * t414 * MDP(23) + (-t414 ^ 2 + t416 ^ 2) * MDP(24) + (t619 + t661) * MDP(25) + (-t603 + t660) * MDP(26) + t398 * MDP(27) + (t348 * t431 - t351 * t416 + t605) * MDP(28) + (t351 * t414 - t431 * t581 - t582) * MDP(29) + (-MDP(25) * t659 - MDP(26) * t416 - MDP(28) * t348 + MDP(29) * t581) * qJD(6);];
tauc  = t1;

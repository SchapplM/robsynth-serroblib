% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:58
% EndTime: 2019-03-09 08:20:08
% DurationCPUTime: 9.02s
% Computational Cost: add. (3179->582), mult. (6703->712), div. (0->0), fcn. (4187->8), ass. (0->243)
t536 = sin(pkin(9));
t537 = cos(pkin(9));
t636 = qJD(2) * t537;
t543 = cos(qJ(2));
t638 = qJD(1) * t543;
t457 = -t536 * t638 + t636;
t624 = qJDD(1) * t543;
t509 = pkin(7) * t624;
t531 = qJDD(2) * qJ(3);
t532 = qJD(2) * qJD(3);
t540 = sin(qJ(2));
t626 = qJD(1) * qJD(2);
t610 = t540 * t626;
t429 = pkin(7) * t610 - t509 - t531 - t532;
t571 = pkin(3) * t624 + qJDD(4) - t429;
t553 = pkin(3) * t610 - t571;
t602 = qJDD(2) * t537 - t536 * t624;
t422 = t536 * t610 + t602;
t665 = qJ(5) * t422;
t548 = qJD(5) * t457 + t553 + t665;
t421 = qJDD(2) * t536 + (-t610 + t624) * t537;
t675 = pkin(4) * t421;
t359 = -t548 + t675;
t541 = sin(qJ(1));
t544 = cos(qJ(1));
t594 = g(1) * t544 + g(2) * t541;
t575 = t594 * t543;
t670 = g(3) * t540;
t684 = t575 + t670;
t685 = t359 - t684;
t677 = pkin(3) + pkin(7);
t609 = t543 * t626;
t625 = qJDD(1) * t540;
t565 = t609 + t625;
t639 = qJD(1) * t540;
t511 = pkin(7) * t639;
t629 = pkin(3) * t639 + qJD(3) + t511;
t622 = MDP(17) + MDP(20);
t682 = t537 * qJD(5) - t629;
t539 = sin(qJ(6));
t542 = cos(qJ(6));
t651 = t542 * t536;
t656 = t537 * t543;
t681 = t539 * t656 - t543 * t651;
t680 = MDP(15) + MDP(19);
t679 = MDP(16) - MDP(21);
t664 = qJ(5) * t536;
t676 = pkin(4) + pkin(5);
t574 = t537 * t676 + t664;
t652 = t541 * t543;
t488 = qJ(3) * t652;
t650 = t543 * t544;
t491 = qJ(3) * t650;
t520 = t540 * qJ(3);
t526 = t543 * pkin(2);
t642 = t526 + t520;
t615 = t543 * qJ(4) + t642;
t669 = pkin(2) + qJ(4);
t678 = t594 * t669 * t540 - g(1) * t491 - g(2) * t488 - g(3) * t615;
t613 = t537 * t638;
t637 = qJD(2) * t536;
t455 = t613 + t637;
t404 = t455 * t539 + t457 * t542;
t603 = -t542 * t421 + t539 * t422;
t361 = qJD(6) * t404 + t603;
t454 = t457 ^ 2;
t674 = g(1) * t541;
t671 = g(2) * t544;
t668 = pkin(8) - t669;
t667 = pkin(7) * qJDD(2);
t666 = qJ(3) * t543;
t663 = qJ(5) * t537;
t662 = qJDD(2) * pkin(2);
t659 = t457 * t539;
t402 = -t542 * t455 + t659;
t495 = -qJD(6) + t639;
t661 = t402 * t495;
t660 = t404 * t495;
t534 = t540 ^ 2;
t546 = qJD(1) ^ 2;
t658 = t534 * t546;
t413 = t536 * t421;
t657 = t536 * t543;
t655 = t540 * t541;
t654 = t540 * t544;
t653 = t540 * t546;
t494 = pkin(2) * t610;
t587 = qJ(4) * t540 - t666;
t633 = qJD(3) * t540;
t549 = qJD(2) * t587 - qJD(4) * t543 - t633;
t606 = -pkin(1) - t520;
t561 = -t543 * t669 + t606;
t371 = qJD(1) * t549 + qJDD(1) * t561 + t494;
t493 = pkin(7) * t609;
t508 = pkin(7) * t625;
t607 = qJDD(3) + t493 + t508;
t395 = pkin(3) * t565 - qJD(2) * qJD(4) - qJDD(2) * t669 + t607;
t358 = t537 * t371 + t536 * t395;
t635 = qJD(2) * t540;
t515 = pkin(2) * t635;
t412 = t515 + t549;
t634 = qJD(2) * t543;
t467 = t677 * t634;
t380 = t537 * t412 + t536 * t467;
t427 = t561 * qJD(1);
t431 = -qJD(2) * t669 + t629;
t377 = t537 * t427 + t536 * t431;
t516 = pkin(2) * t639;
t437 = qJD(1) * t587 + t516;
t513 = pkin(7) * t638;
t466 = pkin(3) * t638 + t513;
t394 = t537 * t437 + t536 * t466;
t452 = -pkin(1) - t615;
t478 = t677 * t540;
t406 = t537 * t452 + t536 * t478;
t580 = -t537 * t539 + t651;
t614 = t537 * t639;
t648 = t580 * qJD(6) + t539 * t614 - t639 * t651;
t579 = t539 * t536 + t542 * t537;
t443 = t579 * qJD(6);
t567 = t540 * t579;
t647 = -qJD(1) * t567 + t443;
t608 = t669 * t625;
t646 = (-t609 * t669 - t608) * t537;
t645 = -t574 * t639 + t682;
t590 = pkin(4) * t537 + t664;
t644 = t590 * t639 - t682;
t479 = t677 * t543;
t527 = t544 * pkin(7);
t641 = t544 * pkin(3) + t527;
t535 = t543 ^ 2;
t640 = t534 - t535;
t533 = qJD(2) * qJ(3);
t632 = qJD(4) * t455;
t631 = qJD(6) * t539;
t630 = qJD(6) * t542;
t518 = t540 * qJD(5);
t441 = qJD(4) + t533 + t466;
t563 = qJ(5) * t457 - t441;
t381 = pkin(4) * t455 - t563;
t628 = -qJD(4) + t381;
t627 = qJD(4) - t441;
t621 = t676 * t540;
t619 = t543 * t653;
t617 = t539 * t421 + t542 * t422 + t455 * t630;
t373 = qJ(5) * t639 + t377;
t383 = qJ(5) * t638 + t394;
t398 = t540 * qJ(5) + t406;
t616 = g(1) * t654 + g(2) * t655 - g(3) * t543;
t612 = t669 * t634;
t611 = qJD(5) * t657;
t605 = -qJD(2) * pkin(2) + qJD(3);
t357 = -t536 * t371 + t395 * t537;
t588 = qJDD(5) - t357;
t356 = -pkin(4) * t565 + t588;
t352 = -pkin(5) * t565 - pkin(8) * t422 + t356;
t355 = qJ(5) * t565 + qJD(1) * t518 + t358;
t353 = pkin(8) * t421 + t355;
t604 = t542 * t352 - t539 * t353;
t379 = -t536 * t412 + t467 * t537;
t376 = -t536 * t427 + t431 * t537;
t393 = -t536 * t437 + t466 * t537;
t405 = -t536 * t452 + t478 * t537;
t369 = qJ(5) * t634 + t380 + t518;
t601 = t544 * pkin(1) + pkin(2) * t650 + t541 * pkin(7) + qJ(3) * t654;
t600 = t413 * t669 + t616;
t599 = -t508 + t616;
t598 = t536 * t608;
t545 = qJD(2) ^ 2;
t597 = pkin(7) * t545 + t671;
t445 = t536 * t541 - t537 * t654;
t447 = t536 * t544 + t537 * t655;
t596 = -g(1) * t447 - g(2) * t445;
t446 = t536 * t654 + t537 * t541;
t448 = -t536 * t655 + t537 * t544;
t595 = -g(1) * t448 - g(2) * t446;
t593 = t373 * t639 - t356;
t592 = t377 * t639 + t357;
t591 = qJD(5) - t376;
t589 = pkin(4) * t536 - t663;
t586 = t539 * t352 + t542 * t353;
t362 = -pkin(8) * t457 - qJD(1) * t621 + t591;
t363 = pkin(8) * t455 + t373;
t350 = t362 * t542 - t363 * t539;
t351 = t362 * t539 + t363 * t542;
t378 = pkin(8) * t657 - t405 - t621;
t384 = pkin(8) * t656 + t398;
t585 = t378 * t542 - t384 * t539;
t584 = t378 * t539 + t384 * t542;
t583 = t447 * t542 - t448 * t539;
t582 = t447 * t539 + t448 * t542;
t472 = t511 + t605;
t477 = -t513 - t533;
t581 = t472 * t543 + t477 * t540;
t578 = qJD(4) * t457 + t422 * t669;
t577 = t606 - t526;
t576 = t541 * pkin(3) + qJ(4) * t650 + t601;
t573 = -0.2e1 * pkin(1) * t626 - t667;
t451 = t577 * qJD(1);
t572 = t451 * t639 + qJDD(3) - t599;
t570 = -qJ(3) * t634 - t633;
t469 = t668 * t536;
t562 = -pkin(8) * t536 * t540 - t543 * t676;
t569 = -qJD(1) * t562 + qJD(4) * t537 - qJD(6) * t469 + t393;
t470 = t668 * t537;
t568 = pkin(8) * t614 - qJD(4) * t536 - qJD(6) * t470 - t383;
t566 = t580 * t540;
t564 = t457 * t631 - t617;
t560 = 0.2e1 * qJDD(1) * pkin(1) - t597;
t474 = -pkin(1) - t642;
t559 = t667 + (-qJD(1) * t474 - t451) * qJD(2);
t556 = t561 * t674;
t554 = -t553 - t684;
t399 = qJD(1) * t570 + qJDD(1) * t577 + t494;
t439 = t515 + t570;
t552 = qJD(1) * t439 + qJDD(1) * t474 + t399 + t597;
t436 = t607 - t662;
t550 = qJD(2) * t581 - t429 * t543 + t436 * t540;
t547 = (-pkin(3) * t626 - g(3)) * t540 - t575 + t571;
t498 = g(1) * t652;
t471 = qJ(3) + t589;
t465 = t677 * t635;
t463 = -qJ(3) * t638 + t516;
t459 = -qJDD(6) + t565;
t444 = -t536 * t676 - qJ(3) + t663;
t432 = t579 * t543;
t419 = t543 * t590 + t479;
t409 = -t543 * t574 - t479;
t400 = -pkin(4) * t540 - t405;
t397 = t611 + (-t590 - t677) * t635;
t392 = t445 * t539 + t446 * t542;
t391 = t445 * t542 - t446 * t539;
t387 = qJD(2) * t566 + t443 * t543;
t386 = qJD(2) * t567 + qJD(6) * t681;
t385 = -pkin(4) * t638 - t393;
t382 = -t611 + (t574 + t677) * t635;
t374 = -pkin(4) * t634 - t379;
t372 = -pkin(4) * t639 + t591;
t366 = -t455 * t676 + t563;
t365 = -pkin(8) * t537 * t635 + t369;
t364 = qJD(2) * t562 - t379;
t354 = -t421 * t676 + t548;
t1 = [(-t387 * t495 - t404 * t634 - t459 * t681 + t540 * t564) * MDP(25) + ((t364 * t539 + t365 * t542) * t495 + t584 * t459 + t586 * t540 + t351 * t634 + t382 * t404 - t409 * t564 + t354 * t681 + t366 * t387 - g(1) * t583 - g(2) * t391 + (t350 * t540 + t495 * t585) * qJD(6)) * MDP(29) + (t387 * t404 - t564 * t681) * MDP(23) + (-t361 * t681 - t386 * t404 - t387 * t402 - t432 * t564) * MDP(24) + (t540 * t559 + t543 * t552 - t498) * MDP(12) + (-g(1) * t641 - g(2) * t576 + t357 * t405 + t358 * t406 + t376 * t379 + t377 * t380 - t441 * t465 - t479 * t553 - t556) * MDP(18) + (t421 * t479 - t455 * t465 + (-t553 * t537 + (qJD(1) * t405 + t376) * qJD(2)) * t543 + (qJD(1) * t379 + qJDD(1) * t405 - t441 * t636 + t357) * t540 + t595) * MDP(15) + (t422 * t479 - t457 * t465 + (t553 * t536 + (-qJD(1) * t406 - t377) * qJD(2)) * t543 + (-qJD(1) * t380 - qJDD(1) * t406 + t441 * t637 - t358) * t540 - t596) * MDP(16) + qJDD(1) * MDP(1) + (pkin(7) * t550 - g(1) * t527 - g(2) * t601 + t399 * t474 + t451 * t439 - t577 * t674) * MDP(14) + (t573 * t543 + (-t560 - t674) * t540) * MDP(10) + (t559 * t543 + (-t552 + t674) * t540) * MDP(13) + (-t671 + t674) * MDP(2) + (-t379 * t457 - t380 * t455 - t405 * t422 - t406 * t421 + t498 + (-t376 * t536 + t377 * t537) * t635 + (t357 * t536 - t358 * t537 - t671) * t543) * MDP(17) + (-t369 * t455 + t374 * t457 - t398 * t421 + t400 * t422 + t498 + (t372 * t536 + t373 * t537) * t635 + (-t355 * t537 - t356 * t536 - t671) * t543) * MDP(20) + 0.2e1 * (t540 * t624 - t626 * t640) * MDP(5) + (t355 * t398 + t373 * t369 + t359 * t419 + t381 * t397 + t356 * t400 + t372 * t374 - g(1) * (pkin(4) * t448 + qJ(5) * t447 + t641) - g(2) * (pkin(4) * t446 + qJ(5) * t445 + t576) - t556) * MDP(22) + (-(t364 * t542 - t365 * t539) * t495 - t585 * t459 - t604 * t540 - t350 * t634 + t382 * t402 + t409 * t361 - t354 * t432 + t366 * t386 - g(1) * t582 - g(2) * t392 + (t351 * t540 + t495 * t584) * qJD(6)) * MDP(28) + (t459 * t540 + t495 * t634) * MDP(27) + (t361 * t540 + t386 * t495 + t402 * t634 - t432 * t459) * MDP(26) + (t397 * t455 + t419 * t421 + (t359 * t537 + (-qJD(1) * t400 - t372) * qJD(2)) * t543 + (-qJD(1) * t374 - qJDD(1) * t400 - t381 * t636 - t356) * t540 + t595) * MDP(19) + (-t397 * t457 - t419 * t422 + (t359 * t536 + (qJD(1) * t398 + t373) * qJD(2)) * t543 + (qJD(1) * t369 + qJDD(1) * t398 - t381 * t637 + t355) * t540 + t596) * MDP(21) + (qJDD(1) * t534 + 0.2e1 * t540 * t609) * MDP(4) + (t540 * t573 + t543 * t560 + t498) * MDP(9) + ((t534 + t535) * qJDD(1) * pkin(7) + t550 - t594) * MDP(11) + t594 * MDP(3) + (qJDD(2) * t540 + t543 * t545) * MDP(6) + (qJDD(2) * t543 - t540 * t545) * MDP(7); MDP(7) * t624 + MDP(6) * t625 + (t421 * t471 + t644 * t455 + t685 * t536 + (t372 * t543 + (t537 * t628 + t385) * t540) * qJD(1) + t646) * MDP(19) + (-t598 - t422 * t471 - t644 * t457 - t685 * t537 + (-t373 * t543 - t383 * t540 + (t540 * t628 - t612) * t536) * qJD(1)) * MDP(21) + (t404 * t648 - t564 * t579) * MDP(23) + (t404 * t638 - t459 * t579 - t495 * t648) * MDP(25) + (-(-t469 * t539 - t470 * t542) * t459 + t444 * t361 - t354 * t580 - g(3) * t566 + (t539 * t568 - t542 * t569) * t495 + t645 * t402 + t647 * t366 + (t350 * qJD(1) - t580 * t594) * t543) * MDP(28) + (-t361 * t579 - t402 * t648 - t404 * t647 - t564 * t580) * MDP(24) + (-t402 * t638 - t459 * t580 + t495 * t647) * MDP(26) + (-t553 * qJ(3) - t377 * t394 - t376 * t393 - (t357 * t537 + t358 * t536) * t669 + t629 * t441 + (-t376 * t537 - t377 * t536) * qJD(4) + t678) * MDP(18) + (t598 + qJ(3) * t422 + t629 * t457 + t554 * t537 + (t377 * t543 + t394 * t540 + (t540 * t627 + t612) * t536) * qJD(1)) * MDP(16) + t640 * MDP(5) * t546 - t495 * MDP(27) * t638 + qJDD(2) * MDP(8) + (t670 - t509 + (pkin(1) * t546 + t594) * t543) * MDP(10) + (-t463 * t638 + t572 - 0.2e1 * t662) * MDP(12) + ((-pkin(2) * t540 + t666) * qJDD(1) + ((-t477 - t533) * t540 + (-t472 + t605) * t543) * qJD(1)) * MDP(11) + (-t429 * qJ(3) - t477 * qJD(3) - t436 * pkin(2) - t451 * t463 - g(1) * (-pkin(2) * t654 + t491) - g(2) * (-pkin(2) * t655 + t488) - g(3) * t642 - t581 * qJD(1) * pkin(7)) * MDP(14) + (pkin(1) * t653 + t599) * MDP(9) + (t359 * t471 - t373 * t383 - t372 * t385 - (t355 * t536 - t356 * t537) * t669 + t644 * t381 + (t372 * t537 - t373 * t536) * qJD(4) - t684 * t589 + t678) * MDP(22) + (qJ(3) * t421 + t629 * t455 + t554 * t536 + (-t376 * t543 + (-t537 * t627 - t393) * t540) * qJD(1) + t646) * MDP(15) + (t383 * t455 - t385 * t457 + (-t372 * t639 - t355 + t632) * t536 + (t578 - t593) * t537 + t600) * MDP(20) + (t393 * t457 + t394 * t455 + (t376 * t639 - t358 + t632) * t536 + (t578 - t592) * t537 + t600) * MDP(17) - MDP(4) * t619 + (t509 + 0.2e1 * t531 + 0.2e1 * t532 + (qJD(1) * t463 - g(3)) * t540 + (qJD(1) * t451 - t594) * t543) * MDP(13) + ((t469 * t542 - t470 * t539) * t459 - t444 * t564 + (t539 * t569 + t542 * t568) * t495 + t645 * t404 + t648 * t366 - t351 * t638 + (t354 + t684) * t579) * MDP(29); MDP(11) * t625 + (qJDD(2) + t619) * MDP(12) + (-t545 - t658) * MDP(13) + (t493 + t572 - t662) * MDP(14) + (MDP(28) * t579 + MDP(29) * t580) * t459 + (t477 * MDP(14) - t441 * MDP(18) - t381 * MDP(22) + t402 * MDP(28) + t404 * MDP(29) - t457 * t679) * qJD(2) + (t592 * MDP(18) + t593 * MDP(22) - t622 * t422 - t658 * t679) * t537 + (MDP(28) * t648 - MDP(29) * t647) * t495 + t680 * (t537 * t625 - t536 * t658 + (-t455 + t613) * qJD(2)) + (t358 * MDP(18) + t355 * MDP(22) - t679 * t625 + (-t679 * t634 + (-t376 * MDP(18) + MDP(22) * t372 + t457 * t622) * t540) * qJD(1)) * t536 + t622 * (-t455 * t614 - t413) + (-MDP(18) - MDP(22)) * t616; (t376 * t457 + t377 * t455 + t547) * MDP(18) + (t675 - t665 + t373 * t455 + (-qJD(5) - t372) * t457 + t547) * MDP(22) + (-t361 + t660) * MDP(28) + (t564 - t661) * MDP(29) + t680 * (t457 * t639 + t421) + t679 * ((-t455 + t637) * t639 + t602) + t622 * (-t455 ^ 2 - t454); (t455 * t457 - t625) * MDP(19) + t602 * MDP(20) + (-t454 - t658) * MDP(21) + (-pkin(4) * t625 - g(1) * t445 + g(2) * t447 - g(3) * t656 + t381 * t457 + t588) * MDP(22) + (-t402 * t457 - t542 * t459 + t495 * t631) * MDP(28) + (-t404 * t457 + t539 * t459 + t495 * t630) * MDP(29) + ((-MDP(22) * pkin(4) - MDP(19)) * t634 + ((t455 + t637) * MDP(20) - t373 * MDP(22) + (-MDP(28) * t539 - MDP(29) * t542) * t495) * t540) * qJD(1); t404 * t402 * MDP(23) + (-t402 ^ 2 + t404 ^ 2) * MDP(24) + (t617 - t661) * MDP(25) + (-t603 - t660) * MDP(26) - t459 * MDP(27) + (-g(1) * t391 + g(2) * t583 - g(3) * t432 - t351 * t495 - t366 * t404 + t604) * MDP(28) + (g(1) * t392 - g(2) * t582 + g(3) * t681 - t350 * t495 + t366 * t402 - t586) * MDP(29) + (-MDP(25) * t659 - MDP(26) * t404 - MDP(28) * t351 - MDP(29) * t350) * qJD(6);];
tau  = t1;

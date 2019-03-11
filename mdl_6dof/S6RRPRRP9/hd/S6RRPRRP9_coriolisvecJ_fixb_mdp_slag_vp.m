% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:34:17
% EndTime: 2019-03-09 12:34:33
% DurationCPUTime: 10.68s
% Computational Cost: add. (10663->561), mult. (28160->767), div. (0->0), fcn. (22717->10), ass. (0->236)
t603 = cos(pkin(11));
t609 = cos(qJ(4));
t700 = t609 * t603;
t601 = sin(pkin(11));
t606 = sin(qJ(4));
t706 = t601 * t606;
t565 = -t700 + t706;
t602 = sin(pkin(6));
t610 = cos(qJ(2));
t703 = t602 * t610;
t619 = t565 * t703;
t531 = qJD(1) * t619;
t561 = t565 * qJD(4);
t722 = t561 - t531;
t566 = t601 * t609 + t603 * t606;
t620 = t566 * t703;
t684 = -qJD(1) * t620 + t566 * qJD(4);
t607 = sin(qJ(2));
t637 = pkin(2) * t607 - qJ(3) * t610;
t682 = qJD(1) * t602;
t553 = t637 * t682;
t663 = t607 * t682;
t604 = cos(pkin(6));
t681 = qJD(1) * t604;
t668 = pkin(1) * t681;
t554 = -pkin(8) * t663 + t610 * t668;
t499 = t603 * t553 - t554 * t601;
t702 = t603 * t610;
t623 = t602 * (pkin(3) * t607 - pkin(9) * t702);
t477 = qJD(1) * t623 + t499;
t500 = t601 * t553 + t603 * t554;
t680 = qJD(1) * t610;
t591 = t602 * t680;
t642 = t601 * t591;
t484 = -pkin(9) * t642 + t500;
t718 = pkin(9) + qJ(3);
t578 = t718 * t601;
t579 = t718 * t603;
t635 = -t578 * t609 - t579 * t606;
t736 = -qJD(3) * t565 + qJD(4) * t635 - t606 * t477 - t609 * t484;
t605 = sin(qJ(5));
t608 = cos(qJ(5));
t592 = qJD(2) + t681;
t542 = t592 * t603 - t601 * t663;
t543 = t592 * t601 + t603 * t663;
t636 = t542 * t606 + t543 * t609;
t671 = t591 - qJD(4);
t470 = -t605 * t671 + t608 * t636;
t490 = -t609 * t542 + t543 * t606;
t727 = qJD(5) + t490;
t710 = t727 * t605;
t735 = t470 * t710;
t729 = pkin(10) * t663 - t736;
t555 = pkin(8) * t591 + t607 * t668;
t523 = pkin(3) * t642 + t555;
t734 = t684 * pkin(4) + t722 * pkin(10) - t523;
t532 = qJ(3) * t592 + t555;
t550 = (-pkin(2) * t610 - qJ(3) * t607 - pkin(1)) * t602;
t537 = qJD(1) * t550;
t478 = -t532 * t601 + t603 * t537;
t436 = -pkin(3) * t591 - pkin(9) * t543 + t478;
t535 = (qJD(2) * t637 - qJD(3) * t607) * t602;
t519 = qJD(1) * t535;
t670 = qJD(1) * qJD(2);
t658 = t602 * t670;
t640 = t607 * t658;
t667 = pkin(1) * qJD(2) * t604;
t645 = qJD(1) * t667;
t627 = -pkin(8) * t640 + t610 * t645;
t520 = qJD(3) * t592 + t627;
t473 = t603 * t519 - t520 * t601;
t618 = qJD(2) * t623;
t437 = qJD(1) * t618 + t473;
t479 = t603 * t532 + t601 * t537;
t450 = pkin(9) * t542 + t479;
t474 = t601 * t519 + t603 * t520;
t639 = t610 * t658;
t634 = t601 * t639;
t451 = -pkin(9) * t634 + t474;
t677 = qJD(4) * t609;
t678 = qJD(4) * t606;
t626 = -t436 * t677 - t606 * t437 + t450 * t678 - t609 * t451;
t384 = pkin(10) * t640 - t626;
t405 = t606 * t436 + t609 * t450;
t402 = -pkin(10) * t671 + t405;
t525 = -pkin(2) * t592 + qJD(3) - t554;
t488 = -pkin(3) * t542 + t525;
t414 = pkin(4) * t490 - pkin(10) * t636 + t488;
t387 = t402 * t608 + t414 * t605;
t617 = qJD(2) * t620;
t616 = qJD(1) * t617;
t720 = qJD(4) * t636;
t449 = t616 + t720;
t548 = pkin(8) * t639 + t607 * t645;
t515 = pkin(3) * t634 + t548;
t688 = t542 * t677 + t639 * t700;
t613 = (-qJD(4) * t543 - t634) * t606 + t688;
t400 = t449 * pkin(4) - pkin(10) * t613 + t515;
t376 = -qJD(5) * t387 - t384 * t605 + t608 * t400;
t650 = t608 * t671;
t676 = qJD(5) * t605;
t410 = qJD(5) * t650 - t605 * t640 - t608 * t613 + t636 * t676;
t371 = pkin(5) * t449 + qJ(6) * t410 - qJD(6) * t470 + t376;
t468 = t605 * t636 + t650;
t380 = -qJ(6) * t468 + t387;
t733 = t380 * t727 + t371;
t675 = qJD(5) * t608;
t375 = t608 * t384 + t605 * t400 - t402 * t676 + t414 * t675;
t612 = t605 * t613 - t608 * t640;
t411 = qJD(5) * t470 + t612;
t372 = -qJ(6) * t411 - qJD(6) * t468 + t375;
t386 = -t402 * t605 + t608 * t414;
t379 = -qJ(6) * t470 + t386;
t378 = pkin(5) * t727 + t379;
t732 = -t378 * t727 + t372;
t731 = t468 * t490;
t730 = t490 * t671;
t717 = -qJ(6) - pkin(10);
t728 = -qJ(6) * t490 + qJD(5) * t717;
t598 = t602 ^ 2;
t726 = -0.2e1 * t598 * t670;
t600 = t610 ^ 2;
t725 = MDP(5) * (t607 ^ 2 - t600);
t549 = pkin(8) * t703 + (pkin(1) * t607 + qJ(3)) * t604;
t495 = -t549 * t601 + t603 * t550;
t704 = t602 * t607;
t560 = t601 * t604 + t603 * t704;
t458 = -pkin(3) * t703 - pkin(9) * t560 + t495;
t496 = t603 * t549 + t601 * t550;
t559 = t601 * t704 - t604 * t603;
t476 = -pkin(9) * t559 + t496;
t692 = t606 * t458 + t609 * t476;
t417 = -pkin(10) * t703 + t692;
t503 = t609 * t559 + t560 * t606;
t504 = -t559 * t606 + t560 * t609;
t552 = pkin(8) * t704 + (-pkin(1) * t610 - pkin(2)) * t604;
t510 = pkin(3) * t559 + t552;
t426 = pkin(4) * t503 - pkin(10) * t504 + t510;
t694 = t608 * t417 + t605 * t426;
t529 = -t578 * t606 + t579 * t609;
t724 = qJD(3) * t566 + qJD(4) * t529 + t477 * t609;
t723 = t734 * t608;
t597 = -pkin(3) * t603 - pkin(2);
t512 = pkin(4) * t565 - pkin(10) * t566 + t597;
t721 = t512 * t675 + t605 * t734 - t608 * t729;
t719 = t470 ^ 2;
t715 = t405 * t607;
t714 = t410 * t605;
t712 = t478 * t607;
t711 = t479 * t607;
t709 = t566 * t605;
t708 = t566 * t608;
t611 = qJD(1) ^ 2;
t707 = t598 * t611;
t705 = t601 * t610;
t701 = t605 * t449;
t443 = t608 * t449;
t517 = t608 * t529;
t699 = t378 - t379;
t502 = -t531 * t608 + t605 * t663;
t633 = qJ(6) * t561 - qJD(6) * t566;
t698 = qJ(6) * t502 - t517 * qJD(5) + t633 * t608 + t723 + ((qJ(6) * t566 - t512) * qJD(5) + t729) * t605 + t684 * pkin(5);
t501 = -t531 * t605 - t608 * t663;
t660 = t566 * t675;
t697 = (-qJD(5) * t529 + t633) * t605 + t721 + (t501 - t660) * qJ(6);
t404 = t609 * t436 - t606 * t450;
t429 = pkin(4) * t636 + pkin(10) * t490;
t696 = t608 * t404 + t605 * t429;
t695 = -t605 * t411 - t468 * t675;
t480 = t606 * t484;
t632 = pkin(4) * t663 - t480;
t690 = t632 + t724;
t689 = t605 * t512 + t517;
t679 = qJD(2) * t607;
t662 = t602 * t679;
t631 = -pkin(8) * t662 + t610 * t667;
t541 = qJD(3) * t604 + t631;
t483 = t601 * t535 + t603 * t541;
t687 = qJD(6) * t608 + t605 * t728 - t696;
t428 = t608 * t429;
t686 = -pkin(5) * t636 - t428 + t728 * t608 + (-qJD(6) + t404) * t605;
t661 = qJD(2) * t703;
t556 = pkin(8) * t661 + t607 * t667;
t674 = t607 * MDP(19);
t669 = pkin(1) * t707;
t666 = t606 * t705;
t665 = t605 * t703;
t641 = t601 * t661;
t524 = pkin(3) * t641 + t556;
t656 = -t417 * t605 + t608 * t426;
t655 = t458 * t609 - t606 * t476;
t654 = t561 * t605 + t501;
t653 = t561 * t608 + t502;
t652 = t608 * t512 - t529 * t605;
t482 = t603 * t535 - t541 * t601;
t651 = t727 ^ 2;
t648 = t727 * t608;
t647 = t671 * t602;
t644 = t598 * t607 * t610 * MDP(4);
t643 = -t436 * t678 + t609 * t437 - t450 * t677 - t606 * t451;
t638 = pkin(1) * t726;
t416 = pkin(4) * t703 - t655;
t455 = t482 + t618;
t472 = -pkin(9) * t641 + t483;
t630 = t455 * t609 - t458 * t678 - t606 * t472 - t476 * t677;
t629 = -t543 * t678 + t688;
t485 = t504 * t605 + t608 * t703;
t401 = pkin(4) * t671 - t404;
t628 = -pkin(10) * t449 + t401 * t727;
t625 = t606 * t455 + t458 * t677 + t609 * t472 - t476 * t678;
t391 = pkin(10) * t662 + t625;
t463 = -qJD(2) * t619 - qJD(4) * t503;
t464 = qJD(4) * t504 + t617;
t409 = pkin(4) * t464 - pkin(10) * t463 + t524;
t624 = t608 * t391 + t605 * t409 - t417 * t676 + t426 * t675;
t622 = -t654 + t660;
t621 = -t566 * t676 - t653;
t385 = -pkin(4) * t640 - t643;
t392 = -pkin(4) * t662 - t630;
t615 = -qJ(3) * t679 + (-pkin(2) * qJD(2) + qJD(3) - t525) * t610;
t614 = -qJD(5) * t694 - t391 * t605 + t608 * t409;
t377 = pkin(5) * t411 + t385;
t584 = t717 * t608;
t583 = t717 * t605;
t486 = t504 * t608 - t665;
t465 = t468 ^ 2;
t440 = -qJ(6) * t709 + t689;
t430 = pkin(5) * t565 - qJ(6) * t708 + t652;
t419 = -qJD(5) * t665 + t463 * t605 + t504 * t675 - t608 * t662;
t418 = qJD(5) * t485 - t608 * t463 - t605 * t662;
t395 = t468 * pkin(5) + qJD(6) + t401;
t389 = -qJ(6) * t485 + t694;
t382 = pkin(5) * t503 - qJ(6) * t486 + t656;
t374 = -qJ(6) * t419 - qJD(6) * t485 + t624;
t373 = pkin(5) * t464 + qJ(6) * t418 - qJD(6) * t486 + t614;
t1 = [(-t598 * t680 - t647) * qJD(2) * t674 + (t372 * t389 + t380 * t374 + t371 * t382 + t378 * t373 + t377 * (pkin(5) * t485 + t416) + t395 * (pkin(5) * t419 + t392)) * MDP(30) + (-t592 * t631 - t604 * t627 + t610 * t638) * MDP(10) + (-t548 * t604 - t556 * t592 + t607 * t638) * MDP(9) + (-t463 * t671 + (-t629 * t610 + (t636 * t607 + (t600 * t602 * t706 + t504 * t607) * qJD(1)) * qJD(2)) * t602) * MDP(17) + (t625 * t671 + t524 * t636 + t510 * t629 + t515 * t504 + t488 * t463 + (-t626 * t610 + (-t715 + (-t510 * t666 - t607 * t692) * qJD(1)) * qJD(2)) * t602) * MDP(21) + (MDP(6) * t661 - MDP(7) * t662) * (t592 + t681) + (t463 * t636 + t504 * t613) * MDP(15) + (-t504 * t449 - t463 * t490 - t464 * t636 - t503 * t613) * MDP(16) + t725 * t726 + (t449 * t503 + t464 * t727) * MDP(26) + (-t411 * t503 - t419 * t727 - t449 * t485 - t464 * t468) * MDP(25) + (-t410 * t503 - t418 * t727 + t449 * t486 + t464 * t470) * MDP(24) + (-t375 * t503 + t385 * t486 - t387 * t464 + t392 * t470 - t401 * t418 - t416 * t410 - t449 * t694 - t624 * t727) * MDP(28) + (t376 * t503 + t385 * t485 + t386 * t464 + t392 * t468 + t401 * t419 + t416 * t411 + t449 * t656 + t614 * t727) * MDP(27) + (t543 * t556 + t548 * t560 + ((qJD(1) * t483 + t474) * t610 + (t525 * t702 - t711 + (-t496 * t607 + t552 * t702) * qJD(1)) * qJD(2)) * t602) * MDP(12) + (-t542 * t556 + t548 * t559 + ((-qJD(1) * t482 - t473) * t610 + (t525 * t705 + t712 + (t495 * t607 + t552 * t705) * qJD(1)) * qJD(2)) * t602) * MDP(11) + (t464 * t671 + (t449 * t610 + (-qJD(1) * t503 - t490) * t679) * t602) * MDP(18) + (-t630 * t671 + t524 * t490 + t510 * t449 + t515 * t503 + t488 * t464 + (-t643 * t610 + (qJD(1) * t655 + t404) * t679) * t602) * MDP(20) + (t473 * t495 + t474 * t496 + t478 * t482 + t479 * t483 + t525 * t556 + t548 * t552) * MDP(14) + 0.2e1 * t644 * t670 + (-t473 * t560 - t474 * t559 - t482 * t543 + t483 * t542 + (-t478 * t603 - t479 * t601 + (-t495 * t603 - t496 * t601) * qJD(1)) * t661) * MDP(13) + (-t410 * t486 - t418 * t470) * MDP(22) + (t410 * t485 - t411 * t486 + t418 * t468 - t419 * t470) * MDP(23) + (-t371 * t486 - t372 * t485 - t373 * t470 - t374 * t468 + t378 * t418 - t380 * t419 + t382 * t410 - t389 * t411) * MDP(29); (-pkin(2) * t548 - t478 * t499 - t479 * t500 - t525 * t555 + (-t478 * t601 + t479 * t603) * qJD(3) + (-t473 * t601 + t474 * t603) * qJ(3)) * MDP(14) + (t372 * t440 + t371 * t430 + t377 * (pkin(5) * t709 - t635) + ((qJD(3) * t603 - qJD(4) * t578) * t606 + (qJD(3) * t601 + qJD(4) * t579 + t477) * t609 + t622 * pkin(5) + t632) * t395 + t697 * t380 + t698 * t378) * MDP(30) + t707 * t725 + ((-qJD(2) * t565 + t490) * t663 + t684 * t671) * MDP(18) + ((qJD(2) * t566 - t636) * t663 + t722 * t671) * MDP(17) + (t652 * t449 + t376 * t565 - t635 * t411 + t385 * t709 + (-t529 * t675 + (-qJD(5) * t512 + t729) * t605 + t723) * t727 + t690 * t468 + t684 * t386 + t622 * t401) * MDP(27) + (-t689 * t449 - t375 * t565 + t635 * t410 + t385 * t708 + (t529 * t676 - t721) * t727 + t690 * t470 - t684 * t387 + t621 * t401) * MDP(28) + (-t411 * t565 - t468 * t684 - t566 * t701 - t622 * t727) * MDP(25) + (-t410 * t565 + t443 * t566 + t470 * t684 + t621 * t727) * MDP(24) + (t449 * t565 + t684 * t727) * MDP(26) + (t566 * t613 - t636 * t722) * MDP(15) + (-t566 * t449 + t490 * t722 - t565 * t613 - t636 * t684) * MDP(16) + (MDP(6) * t591 - MDP(7) * t663) * (qJD(2) - t592) + (t597 * t449 + t515 * t565 - t523 * t490 + t684 * t488 + (qJD(2) * t635 - t404) * t663 + (-t480 + t724) * t671) * MDP(20) + (t597 * t629 + t515 * t566 - t523 * t636 - t722 * t488 + (t715 + (-t529 * t607 - t597 * t666) * qJD(2)) * t682 + t736 * t671) * MDP(21) + (t654 * t470 + t653 * t468 + (t714 - t411 * t608 + (t468 * t605 - t470 * t608) * qJD(5)) * t566) * MDP(23) + (-t543 * t555 + t548 * t601 + (-t500 * t610 + t603 * t615 + t711) * t682) * MDP(12) + (t542 * t555 - t548 * t603 + (t499 * t610 + t601 * t615 - t712) * t682) * MDP(11) + (-t410 * t708 + t470 * t621) * MDP(22) + (t410 * t430 - t411 * t440 - t698 * t470 - t697 * t468 + t654 * t380 + t653 * t378 + (-t371 * t608 - t372 * t605 + (t378 * t605 - t380 * t608) * qJD(5)) * t566) * MDP(29) + (t555 * t592 + t607 * t669 - t548) * MDP(9) + (t554 * t592 + t610 * t669 - t627) * MDP(10) + (t499 * t543 - t500 * t542 + (qJD(3) * t542 + t478 * t591 + t474) * t603 + (qJD(3) * t543 + t479 * t591 - t473) * t601) * MDP(13) + qJD(1) * t647 * t674 - t611 * t644; (-t542 ^ 2 - t543 ^ 2) * MDP(13) + (t478 * t543 - t479 * t542 + t548) * MDP(14) + (t616 + 0.2e1 * t720) * MDP(20) + (t613 + t730) * MDP(21) + (-t468 * t636 - t605 * t651 + t443) * MDP(27) + (-t470 * t636 - t608 * t651 - t701) * MDP(28) + ((t410 - t731) * t608 + t735 + t695) * MDP(29) + (-t395 * t636 + t605 * t732 + t608 * t733) * MDP(30) + ((qJD(2) * t601 - t543) * MDP(11) + (qJD(2) * t603 - t542) * MDP(12) - t636 * MDP(20)) * t591; -t490 ^ 2 * MDP(16) + (t613 - t730) * MDP(17) + (-t566 * t639 - t720) * MDP(18) + MDP(19) * t640 + (-t405 * t671 + t643) * MDP(20) + (-t404 * t671 + t488 * t490 + t626) * MDP(21) + (t470 * t648 - t714) * MDP(22) + ((-t410 - t731) * t608 - t735 + t695) * MDP(23) + (t648 * t727 + t701) * MDP(24) + (-t710 * t727 + t443) * MDP(25) + (-pkin(4) * t411 - t385 * t608 - t405 * t468 + (-pkin(10) * t675 - t428) * t727 + (t404 * t727 + t628) * t605) * MDP(27) + (pkin(4) * t410 + t385 * t605 - t405 * t470 + (pkin(10) * t676 + t696) * t727 + t628 * t608) * MDP(28) + (t410 * t583 + t411 * t584 - t687 * t468 - t686 * t470 - t605 * t733 + t608 * t732) * MDP(29) + (-t372 * t584 + t371 * t583 + t377 * (-pkin(5) * t608 - pkin(4)) + (pkin(5) * t710 - t405) * t395 + t687 * t380 + t686 * t378) * MDP(30) + (t490 * MDP(15) + MDP(16) * t636 - t671 * MDP(18) - t488 * MDP(20) - t470 * MDP(24) + t468 * MDP(25) - MDP(26) * t727 - t386 * MDP(27) + t387 * MDP(28)) * t636; t470 * t468 * MDP(22) + (-t465 + t719) * MDP(23) + (t468 * t727 - t410) * MDP(24) + (-t612 + (-qJD(5) + t727) * t470) * MDP(25) + t449 * MDP(26) + (t387 * t727 - t401 * t470 + t376) * MDP(27) + (t386 * t727 + t401 * t468 - t375) * MDP(28) + (pkin(5) * t410 - t468 * t699) * MDP(29) + (t699 * t380 + (-t395 * t470 + t371) * pkin(5)) * MDP(30); (-t465 - t719) * MDP(29) + (t378 * t470 + t380 * t468 + t377) * MDP(30);];
tauc  = t1;

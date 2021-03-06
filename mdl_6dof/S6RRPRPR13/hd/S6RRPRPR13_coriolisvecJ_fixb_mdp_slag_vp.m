% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:51
% EndTime: 2019-03-09 11:30:09
% DurationCPUTime: 10.88s
% Computational Cost: add. (7675->626), mult. (19843->857), div. (0->0), fcn. (14864->10), ass. (0->251)
t621 = cos(qJ(2));
t613 = sin(pkin(6));
t618 = sin(qJ(2));
t695 = qJD(1) * t618;
t675 = t613 * t695;
t615 = cos(pkin(6));
t696 = qJD(1) * t615;
t682 = pkin(1) * t696;
t699 = -pkin(8) * t675 + t621 * t682;
t685 = qJD(3) - t699;
t599 = qJD(2) + t696;
t620 = cos(qJ(4));
t617 = sin(qJ(4));
t697 = qJD(1) * t613;
t674 = t621 * t697;
t653 = t617 * t674;
t540 = t599 * t620 - t653;
t585 = qJD(4) + t675;
t612 = sin(pkin(11));
t614 = cos(pkin(11));
t477 = t540 * t612 - t585 * t614;
t619 = cos(qJ(6));
t712 = t619 * t477;
t479 = t540 * t614 + t585 * t612;
t616 = sin(qJ(6));
t727 = t479 * t616;
t425 = t712 + t727;
t538 = t599 * t617 + t620 * t674;
t535 = qJD(6) + t538;
t743 = t425 * t535;
t647 = pkin(4) * t620 + qJ(5) * t617;
t742 = (-pkin(3) - t647) * t675 - qJD(4) * t647 + qJD(5) * t620 - t685;
t736 = -t612 * t616 + t614 * t619;
t735 = t736 * qJD(6);
t640 = t477 * t616 - t479 * t619;
t741 = t535 * t640;
t529 = t612 * t617 * t675 - t614 * t674;
t690 = qJD(4) * t617;
t740 = -t612 * t690 - t529;
t610 = t618 ^ 2;
t739 = MDP(5) * (-t621 ^ 2 + t610);
t598 = pkin(2) * t675;
t646 = pkin(9) * t618 - qJ(3) * t621;
t522 = t646 * t697 + t598;
t550 = pkin(8) * t674 + t618 * t682;
t524 = pkin(3) * t674 + t550;
t704 = t522 * t620 + t524 * t617;
t437 = qJ(5) * t674 + t704;
t738 = -t437 * t612 + t614 * t742;
t622 = -pkin(2) - pkin(9);
t688 = qJD(4) * t622;
t669 = t620 * t688;
t708 = (t437 - t669) * t614 + t742 * t612;
t733 = pkin(1) * t618;
t604 = t615 * t733;
t716 = t613 * t621;
t737 = pkin(8) * t716 + t604;
t715 = t614 * t616;
t562 = t612 * t619 + t715;
t626 = t562 * qJD(6);
t686 = pkin(3) * t675 + t685;
t734 = pkin(3) + pkin(8);
t732 = pkin(10) * t614;
t731 = pkin(10) * t620;
t730 = pkin(10) + qJ(5);
t683 = qJD(1) * qJD(2);
t667 = t613 * t683;
t588 = t621 * t667;
t464 = t599 * t622 + t686;
t652 = t618 * t667;
t580 = pkin(2) * t652;
t692 = qJD(3) * t618;
t624 = (qJD(2) * t646 - t692) * t613;
t474 = qJD(1) * t624 + t580;
t666 = -qJ(3) * t618 - pkin(1);
t520 = (t621 * t622 + t666) * t613;
t498 = qJD(1) * t520;
t681 = pkin(1) * qJD(2) * t615;
t657 = qJD(1) * t681;
t541 = pkin(8) * t588 + t618 * t657;
t504 = pkin(3) * t588 + t541;
t689 = qJD(4) * t620;
t655 = t464 * t690 + t474 * t617 + t498 * t689 - t504 * t620;
t395 = -pkin(4) * t588 + t655;
t729 = t395 * t612;
t728 = t395 * t614;
t484 = -qJD(4) * t653 + t599 * t689 - t620 * t652;
t726 = t484 * t736;
t725 = t484 * t562;
t724 = t484 * t612;
t723 = t484 * t614;
t722 = t538 * t585;
t721 = t538 * t612;
t720 = t540 * t585;
t632 = t585 * t620;
t609 = t613 ^ 2;
t719 = t609 * qJD(1) ^ 2;
t717 = t613 * t618;
t714 = t617 * t618;
t713 = t617 * t622;
t710 = t620 * t622;
t628 = -t464 * t689 - t474 * t620 + t498 * t690 - t504 * t617;
t392 = qJ(5) * t588 + qJD(5) * t585 - t628;
t700 = -pkin(8) * t652 + t621 * t657;
t510 = -qJD(3) * t599 - t700;
t475 = -pkin(3) * t652 - t510;
t483 = -qJD(4) * t538 + t617 * t652;
t411 = pkin(4) * t484 - qJ(5) * t483 - qJD(5) * t540 + t475;
t379 = t392 * t614 + t411 * t612;
t694 = qJD(2) * t618;
t673 = t613 * t694;
t593 = pkin(2) * t673;
t492 = t593 + t624;
t601 = pkin(8) * t717;
t676 = -pkin(1) * t621 - pkin(2);
t502 = pkin(3) * t717 + t601 + (-pkin(9) + t676) * t615;
t525 = (t716 * t734 + t604) * qJD(2);
t627 = t492 * t620 + t502 * t689 - t520 * t690 + t525 * t617;
t693 = qJD(2) * t621;
t403 = (qJ(5) * t693 + qJD(5) * t618) * t613 + t627;
t597 = t621 * t681;
t607 = t615 * qJD(3);
t501 = -t673 * t734 + t597 + t607;
t555 = t615 * t617 + t620 * t716;
t515 = -qJD(4) * t555 + t617 * t673;
t678 = t617 * t716;
t516 = -qJD(4) * t678 + t615 * t689 - t620 * t673;
t556 = t615 * t620 - t678;
t418 = pkin(4) * t516 - qJ(5) * t515 - qJD(5) * t556 + t501;
t382 = t403 * t614 + t418 * t612;
t709 = t612 * t669 + t738;
t431 = t464 * t617 + t498 * t620;
t422 = qJ(5) * t585 + t431;
t589 = t599 * qJ(3);
t489 = t589 + t524;
t432 = pkin(4) * t538 - qJ(5) * t540 + t489;
t394 = t422 * t614 + t432 * t612;
t430 = t464 * t620 - t498 * t617;
t462 = pkin(4) * t540 + qJ(5) * t538;
t410 = t430 * t614 + t462 * t612;
t705 = t502 * t617 + t520 * t620;
t434 = qJ(5) * t717 + t705;
t542 = -qJ(3) * t615 - t737;
t519 = pkin(3) * t716 - t542;
t443 = pkin(4) * t555 - qJ(5) * t556 + t519;
t406 = t434 * t614 + t443 * t612;
t530 = (t612 * t621 + t614 * t714) * t697;
t707 = t529 * t616 - t530 * t619 - t620 * t626 - t690 * t736;
t706 = -t530 * t616 + t619 * t740 + t620 * t735 - t690 * t715;
t660 = -t522 * t617 + t524 * t620;
t438 = -pkin(4) * t674 - t660;
t664 = pkin(5) * t612 - t622;
t703 = -pkin(5) * t529 - t664 * t690 - t438;
t702 = t538 * t736 + t735;
t701 = t538 * t562 + t626;
t568 = pkin(4) * t617 - qJ(5) * t620 + qJ(3);
t533 = t568 * t612 + t614 * t713;
t691 = qJD(4) * t535;
t421 = -pkin(4) * t585 + qJD(5) - t430;
t684 = -qJD(5) + t421;
t680 = t618 * t632;
t679 = t621 * t719;
t453 = t483 * t612 - t588 * t614;
t454 = t483 * t614 + t588 * t612;
t677 = -qJD(6) * t712 - t453 * t616 + t454 * t619;
t672 = t613 * t693;
t671 = t585 * t689;
t668 = t609 * t683;
t665 = -t612 * t622 + pkin(5);
t378 = -t392 * t612 + t411 * t614;
t375 = pkin(5) * t484 - pkin(10) * t454 + t378;
t376 = -pkin(10) * t453 + t379;
t663 = t375 * t619 - t376 * t616;
t381 = -t403 * t612 + t418 * t614;
t393 = -t422 * t612 + t432 * t614;
t409 = -t430 * t612 + t462 * t614;
t405 = -t434 * t612 + t443 * t614;
t662 = t453 * t619 + t454 * t616;
t661 = t502 * t620 - t520 * t617;
t656 = t618 * t679;
t654 = t620 * t675;
t651 = -0.2e1 * pkin(1) * t668;
t650 = t550 * t599 - t541;
t503 = -t612 * t731 + t533;
t649 = -pkin(5) * t654 - pkin(10) * t530 + qJD(6) * t503 - (t617 * t732 + t620 * t665) * qJD(4) + t738;
t558 = t614 * t568;
t490 = -t614 * t731 + t617 * t665 + t558;
t648 = pkin(10) * t740 - qJD(6) * t490 + t708;
t645 = t375 * t616 + t376 * t619;
t644 = -t378 * t612 + t379 * t614;
t383 = pkin(5) * t538 - pkin(10) * t479 + t393;
t385 = -pkin(10) * t477 + t394;
t372 = t383 * t619 - t385 * t616;
t373 = t383 * t616 + t385 * t619;
t514 = t556 * t614 + t612 * t717;
t388 = pkin(5) * t555 - pkin(10) * t514 + t405;
t513 = t556 * t612 - t614 * t717;
t396 = -pkin(10) * t513 + t406;
t643 = t388 * t619 - t396 * t616;
t642 = t388 * t616 + t396 * t619;
t641 = -t393 * t612 + t394 * t614;
t439 = t513 * t619 + t514 * t616;
t440 = -t513 * t616 + t514 * t619;
t551 = t737 * qJD(2);
t639 = t541 * t615 + t551 * t599;
t638 = -pkin(8) * t673 + t597;
t637 = -t492 * t617 - t502 * t690 - t520 * t689 + t525 * t620;
t636 = -t599 * t674 + t588;
t635 = t489 * t618 + t622 * t693;
t634 = t585 * t617;
t578 = t730 * t614;
t631 = pkin(5) * t540 + qJD(5) * t612 + qJD(6) * t578 + t538 * t732 + t409;
t577 = t730 * t612;
t630 = pkin(10) * t721 - qJD(5) * t614 + qJD(6) * t577 + t410;
t386 = -qJD(6) * t727 + t677;
t436 = -pkin(4) * t717 - t661;
t543 = (-pkin(2) * t621 + t666) * t613;
t625 = (-qJ(3) * t693 - t692) * t613;
t404 = -pkin(4) * t672 - t637;
t387 = -qJD(6) * t640 + t662;
t606 = -pkin(5) * t614 - pkin(4);
t566 = t620 * t588;
t559 = t664 * t620;
t548 = -qJ(3) * t674 + t598;
t546 = t736 * t620;
t545 = t562 * t620;
t544 = t615 * t676 + t601;
t534 = -t607 - t638;
t532 = -t612 * t713 + t558;
t531 = qJD(1) * t543;
t528 = t599 * t612 - t614 * t654;
t527 = t599 * t614 + t612 * t654;
t526 = t593 + t625;
t521 = -t589 - t550;
t517 = -pkin(2) * t599 + t685;
t507 = qJD(1) * t625 + t580;
t499 = t531 * t675;
t466 = t515 * t614 + t612 * t672;
t465 = t515 * t612 - t614 * t672;
t419 = pkin(5) * t513 + t436;
t417 = -pkin(5) * t721 + t431;
t412 = pkin(5) * t477 + t421;
t399 = qJD(6) * t440 + t465 * t619 + t466 * t616;
t398 = -qJD(6) * t439 - t465 * t616 + t466 * t619;
t390 = pkin(5) * t465 + t404;
t384 = pkin(5) * t453 + t395;
t380 = -pkin(10) * t465 + t382;
t377 = pkin(5) * t516 - pkin(10) * t466 + t381;
t371 = -qJD(6) * t373 + t663;
t370 = qJD(6) * t372 + t645;
t1 = [((-qJD(6) * t642 + t377 * t619 - t380 * t616) * t535 + t643 * t484 + t371 * t555 + t372 * t516 + t390 * t425 + t419 * t387 + t384 * t439 + t412 * t399) * MDP(31) + (-(qJD(6) * t643 + t377 * t616 + t380 * t619) * t535 - t642 * t484 - t370 * t555 - t373 * t516 - t390 * t640 + t419 * t386 + t384 * t440 + t412 * t398) * MDP(32) + (t386 * t555 + t398 * t535 + t440 * t484 - t516 * t640) * MDP(28) + (t386 * t440 - t398 * t640) * MDP(26) + (-t386 * t439 - t387 * t440 - t398 * t425 + t399 * t640) * MDP(27) + (MDP(6) * t672 - MDP(7) * t673) * (t599 + t696) + (t618 * t651 - t639) * MDP(9) + (-t599 * t638 - t615 * t700 + t621 * t651) * MDP(10) + (t585 * t613 + t609 * t695) * MDP(19) * t693 + (t515 * t585 + (t483 * t618 + (qJD(1) * t556 + t540) * t693) * t613) * MDP(17) + (-t516 * t585 + (-t484 * t618 + (-qJD(1) * t555 - t538) * t693) * t613) * MDP(18) + (t637 * t585 + t501 * t538 + t519 * t484 + t475 * t555 + t489 * t516 + (-t655 * t618 + (qJD(1) * t661 + t430) * t693) * t613) * MDP(20) + (-t510 * t615 - t534 * t599 + (-t531 * t693 - t507 * t618 + (-t526 * t618 - t543 * t693) * qJD(1)) * t613) * MDP(13) + ((-t531 * t694 + t507 * t621 + (t526 * t621 - t543 * t694) * qJD(1)) * t613 + t639) * MDP(12) + (-t510 * t621 + t541 * t618 + (t517 * t621 + t521 * t618) * qJD(2) + (-t534 * t621 + t551 * t618 + (t542 * t618 + t544 * t621) * qJD(2)) * qJD(1)) * t613 * MDP(11) + (-t627 * t585 + t501 * t540 + t519 * t483 + t475 * t556 + t489 * t515 + (t628 * t618 + (-qJD(1) * t705 - t431) * t693) * t613) * MDP(21) + (-t483 * t555 - t484 * t556 - t515 * t538 - t516 * t540) * MDP(16) + (t483 * t556 + t515 * t540) * MDP(15) + (t507 * t543 + t510 * t542 + t517 * t551 + t521 * t534 + t526 * t531 + t541 * t544) * MDP(14) + (t484 * t555 + t516 * t535) * MDP(30) + (-t387 * t555 - t399 * t535 - t425 * t516 - t439 * t484) * MDP(29) + (t378 * t555 + t381 * t538 + t393 * t516 + t395 * t513 + t404 * t477 + t405 * t484 + t421 * t465 + t436 * t453) * MDP(22) + (-t379 * t555 - t382 * t538 - t394 * t516 + t395 * t514 + t404 * t479 - t406 * t484 + t421 * t466 + t436 * t454) * MDP(23) + (-t378 * t514 - t379 * t513 - t381 * t479 - t382 * t477 - t393 * t466 - t394 * t465 - t405 * t454 - t406 * t453) * MDP(24) + 0.2e1 * (MDP(4) * t618 * t621 - t739) * t668 + (t378 * t405 + t379 * t406 + t381 * t393 + t382 * t394 + t395 * t436 + t404 * t421) * MDP(25); ((-t484 - t720) * t620 + (-t483 + t722) * t617) * MDP(16) + (pkin(1) * t679 + t599 * t699 - t700) * MDP(10) + (-t387 * t617 - t425 * t632 - t484 * t545 - t535 * t706) * MDP(29) + (-pkin(2) * t541 - qJ(3) * t510 - t517 * t550 - t521 * t685 - t531 * t548) * MDP(14) + (t483 * t620 - t540 * t634) * MDP(15) + (qJ(3) * t484 + t475 * t617 - t660 * t585 + t686 * t538 + (t489 * t620 - t585 * t713) * qJD(4) + (-t430 * t621 + t620 * t635) * t697) * MDP(20) + (-t585 * t690 + t566 + (-t540 * t621 - t585 * t714) * t697) * MDP(17) + (-t548 * t674 + t499 - t650) * MDP(12) + (-qJD(2) + t599) * MDP(7) * t675 + (t719 * t733 + t650) * MDP(9) - MDP(4) * t656 + (t386 * t546 - t640 * t707) * MDP(26) + (t386 * t617 + t484 * t546 + t535 * t707 - t632 * t640) * MDP(28) + (-t386 * t545 - t387 * t546 - t425 * t707 + t640 * t706) * MDP(27) + (-(t490 * t616 + t503 * t619) * t484 - t370 * t617 + t559 * t386 + t384 * t546 + (t616 * t649 + t619 * t648) * t535 - t703 * t640 + t707 * t412 - t373 * t632) * MDP(32) + (-t671 + (-t680 + (-qJD(2) * t617 + t538) * t621) * t697) * MDP(18) + (t685 * t599 + (t531 * t621 + t548 * t618) * t697 - t510) * MDP(13) + (t484 * t617 + t535 * t632) * MDP(30) - t585 * MDP(19) * t674 + t636 * MDP(6) + (-t421 * t530 - t438 * t479 - t484 * t533 + t708 * t538 + (-t379 + (-t421 * t614 + t479 * t622) * qJD(4)) * t617 + (-t394 * t585 - t454 * t622 + t728) * t620) * MDP(23) + (-t421 * t529 - t438 * t477 + t484 * t532 - t709 * t538 + (t378 + (-t421 * t612 + t477 * t622) * qJD(4)) * t617 + (t393 * t585 - t453 * t622 + t729) * t620) * MDP(22) + ((t490 * t619 - t503 * t616) * t484 + t371 * t617 + t559 * t387 + t384 * t545 + (t616 * t648 - t619 * t649) * t535 + t703 * t425 + t706 * t412 + t372 * t632) * MDP(31) + (t393 * t530 + t394 * t529 - t453 * t533 - t454 * t532 + (-t378 * t614 - t379 * t612) * t620 + t709 * t479 + t708 * t477 + (t393 * t614 + t394 * t612) * t690) * MDP(24) + (qJ(3) * t483 + t475 * t620 + t704 * t585 + t686 * t540 + (-t489 * t617 - t585 * t710) * qJD(4) + (t431 * t621 - t617 * t635) * t697) * MDP(21) + (-t395 * t710 + t378 * t532 + t379 * t533 + (t617 * t688 - t438) * t421 - t708 * t394 - t709 * t393) * MDP(25) + ((-qJ(3) * qJD(2) - t521 - t550) * t618 + (-pkin(2) * qJD(2) - t517 + t685) * t621) * MDP(11) * t697 + t719 * t739; t636 * MDP(11) + MDP(12) * t656 + (-t599 ^ 2 - t610 * t719) * MDP(13) + (t521 * t599 + t499 + t541) * MDP(14) + (-t538 * t599 - t585 * t634 + t566) * MDP(20) + (-t671 - t540 * t599 + (-t617 * t693 - t680) * t697) * MDP(21) + (-t453 * t620 + (-t612 * t689 - t527) * t538 + (t477 * t585 - t724) * t617) * MDP(22) + (-t454 * t620 + (-t614 * t689 + t528) * t538 + (t479 * t585 - t723) * t617) * MDP(23) + (t477 * t528 + t479 * t527 + (-t453 * t614 + t454 * t612) * t617 + (-t477 * t614 + t479 * t612) * t689) * MDP(24) + (-t393 * t527 - t394 * t528 + (qJD(4) * t641 - t395) * t620 + (t421 * t585 + t644) * t617) * MDP(25) + (-(t527 * t619 - t528 * t616) * t535 + (-t562 * t691 - t387) * t620 + (t425 * t585 - t535 * t735 - t725) * t617) * MDP(31) + ((t527 * t616 + t528 * t619) * t535 + (-t691 * t736 - t386) * t620 + (t535 * t626 - t585 * t640 - t726) * t617) * MDP(32); -t538 ^ 2 * MDP(16) + (t483 + t722) * MDP(17) + (-t484 + t720) * MDP(18) + MDP(19) * t588 + (t431 * t585 - t655) * MDP(20) + (t430 * t585 + t489 * t538 + t628) * MDP(21) + (-qJ(5) * t724 - pkin(4) * t453 - t728 - t431 * t477 + (t612 * t684 - t409) * t538) * MDP(22) + (-qJ(5) * t723 - pkin(4) * t454 + t729 - t431 * t479 + (t614 * t684 + t410) * t538) * MDP(23) + (t409 * t479 + t410 * t477 + (-qJ(5) * t453 - qJD(5) * t477 - t393 * t538 + t379) * t614 + (qJ(5) * t454 + qJD(5) * t479 - t394 * t538 - t378) * t612) * MDP(24) + (-pkin(4) * t395 + qJ(5) * t644 + qJD(5) * t641 - t393 * t409 - t394 * t410 - t421 * t431) * MDP(25) + (t386 * t562 - t640 * t702) * MDP(26) + (t386 * t736 - t387 * t562 - t425 * t702 + t640 * t701) * MDP(27) + (t535 * t702 + t725) * MDP(28) + (-t535 * t701 + t726) * MDP(29) + ((-t577 * t619 - t578 * t616) * t484 + t606 * t387 - t384 * t736 - t417 * t425 + (t616 * t630 - t619 * t631) * t535 + t701 * t412) * MDP(31) + (-(-t577 * t616 + t578 * t619) * t484 + t606 * t386 + t384 * t562 + t417 * t640 + (t616 * t631 + t619 * t630) * t535 + t702 * t412) * MDP(32) + (MDP(15) * t538 + MDP(16) * t540 - MDP(20) * t489 - MDP(22) * t393 + MDP(23) * t394 + MDP(28) * t640 + MDP(29) * t425 - MDP(30) * t535 - MDP(31) * t372 + MDP(32) * t373) * t540; (t479 * t538 + t453) * MDP(22) + (-t477 * t538 + t454) * MDP(23) + (-t477 ^ 2 - t479 ^ 2) * MDP(24) + (t393 * t479 + t394 * t477 + t395) * MDP(25) + (t387 - t741) * MDP(31) + (t386 - t743) * MDP(32); -t640 * t425 * MDP(26) + (-t425 ^ 2 + t640 ^ 2) * MDP(27) + (t677 + t743) * MDP(28) + (-t662 - t741) * MDP(29) + t484 * MDP(30) + (t373 * t535 + t412 * t640 + t663) * MDP(31) + (t372 * t535 + t412 * t425 - t645) * MDP(32) + (-MDP(28) * t727 + MDP(29) * t640 - MDP(31) * t373 - MDP(32) * t372) * qJD(6);];
tauc  = t1;

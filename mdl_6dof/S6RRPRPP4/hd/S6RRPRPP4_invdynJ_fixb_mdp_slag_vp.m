% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:59
% EndTime: 2019-03-09 10:02:08
% DurationCPUTime: 8.55s
% Computational Cost: add. (6291->607), mult. (12853->733), div. (0->0), fcn. (7911->10), ass. (0->281)
t611 = sin(qJ(2));
t697 = qJD(1) * qJD(2);
t675 = t611 * t697;
t614 = cos(qJ(2));
t695 = qJDD(1) * t614;
t781 = t675 - t695;
t762 = pkin(3) + pkin(7);
t613 = cos(qJ(4));
t610 = sin(qJ(4));
t711 = qJD(2) * t610;
t714 = qJD(1) * t614;
t523 = t613 * t714 + t711;
t676 = t610 * t714;
t709 = qJD(2) * t613;
t525 = -t676 + t709;
t607 = sin(pkin(9));
t608 = cos(pkin(9));
t470 = t608 * t523 + t525 * t607;
t715 = qJD(1) * t611;
t566 = qJD(4) + t715;
t780 = t470 * t566;
t704 = qJD(4) * t613;
t677 = t608 * t704;
t682 = t613 * t715;
t705 = qJD(4) * t610;
t743 = t607 * t610;
t722 = -t607 * t705 + t608 * t682 - t715 * t743 + t677;
t650 = t607 * t613 + t608 * t610;
t721 = t607 * t704 + t608 * t705 + t650 * t715;
t583 = pkin(7) * t715;
t779 = qJD(3) + t583;
t772 = MDP(22) + MDP(25);
t652 = -t523 * t607 + t608 * t525;
t778 = t652 ^ 2;
t612 = sin(qJ(1));
t615 = cos(qJ(1));
t657 = g(1) * t615 + g(2) * t612;
t777 = t657 * t614;
t728 = t610 * pkin(4) + qJ(3);
t588 = pkin(2) * t715;
t752 = qJ(3) * t614;
t655 = pkin(8) * t611 - t752;
t498 = qJD(1) * t655 + t588;
t584 = pkin(7) * t714;
t531 = pkin(3) * t714 + t584;
t514 = t613 * t531;
t742 = t610 * t611;
t446 = -t498 * t610 + t514 + (pkin(4) * t614 - qJ(5) * t742) * qJD(1);
t720 = t613 * t498 + t610 * t531;
t453 = qJ(5) * t682 + t720;
t616 = -pkin(2) - pkin(8);
t727 = qJ(5) - t616;
t669 = t727 * t613;
t494 = -qJD(4) * t669 - qJD(5) * t610;
t632 = -qJD(5) * t613 + t705 * t727;
t725 = (-t446 + t632) * t608 + (t453 - t494) * t607;
t592 = t611 * qJ(3);
t597 = t614 * pkin(2);
t717 = t597 + t592;
t536 = -pkin(1) - t717;
t518 = -pkin(8) * t614 + t536;
t547 = t762 * t611;
t719 = t613 * t518 + t610 * t547;
t609 = -qJ(5) - pkin(8);
t741 = t610 * t612;
t694 = pkin(4) * t741;
t738 = t611 * t612;
t776 = t609 * t738 + t614 * t694;
t739 = t610 * t615;
t569 = pkin(4) * t739;
t736 = t611 * t615;
t775 = t614 * t569 + t609 * t736;
t703 = qJD(4) * t614;
t678 = t610 * t703;
t774 = t611 * t709 + t678;
t674 = t614 * t697;
t696 = qJDD(1) * t611;
t638 = t674 + t696;
t520 = qJDD(4) + t638;
t501 = t613 * t520;
t773 = -t566 * t705 + t501;
t465 = -qJD(4) * t676 + qJDD(2) * t610 + (qJD(2) * qJD(4) - t781) * t613;
t770 = pkin(4) * t465 + qJDD(5);
t464 = -qJD(4) * t523 + t613 * qJDD(2) + t610 * t781;
t431 = t464 * t607 + t608 * t465;
t432 = t464 * t608 - t465 * t607;
t455 = t608 * t494 + t607 * t632;
t533 = t727 * t610;
t476 = -t533 * t607 + t608 * t669;
t477 = -t608 * t533 - t607 * t669;
t769 = -t477 * t431 + t432 * t476 - t455 * t470;
t604 = qJD(2) * qJ(3);
t506 = t604 + t531;
t475 = pkin(4) * t523 + qJD(5) + t506;
t423 = pkin(5) * t470 - qJ(6) * t652 + t475;
t601 = qJ(4) + pkin(9);
t589 = sin(t601);
t590 = cos(t601);
t489 = t589 * t612 - t590 * t736;
t491 = t589 * t615 + t590 * t738;
t600 = g(3) * t614;
t671 = -pkin(1) - t592;
t636 = t614 * t616 + t671;
t486 = t636 * qJD(1);
t700 = pkin(3) * t715 + t779;
t493 = qJD(2) * t616 + t700;
t451 = t486 * t613 + t493 * t610;
t565 = pkin(2) * t675;
t707 = qJD(3) * t611;
t625 = qJD(2) * t655 - t707;
t449 = qJD(1) * t625 + qJDD(1) * t636 + t565;
t564 = pkin(7) * t674;
t580 = pkin(7) * t696;
t673 = qJDD(3) + t564 + t580;
t468 = t638 * pkin(3) + qJDD(2) * t616 + t673;
t666 = -t449 * t610 + t613 * t468;
t623 = -qJD(4) * t451 + t666;
t404 = pkin(4) * t520 - qJ(5) * t464 - qJD(5) * t525 + t623;
t687 = -t613 * t449 - t610 * t468 - t493 * t704;
t407 = -qJ(5) * t465 - qJD(5) * t523 - t486 * t705 - t687;
t395 = t608 * t404 - t607 * t407;
t672 = -qJDD(6) + t395;
t768 = g(1) * t489 - g(2) * t491 - t423 * t652 + t590 * t600 + t672;
t767 = t506 * t566 + t616 * t520;
t754 = g(3) * t611;
t766 = t754 + t777;
t579 = pkin(4) * t613 + pkin(3);
t765 = pkin(4) * t704 + t579 * t715 + t779;
t396 = t607 * t404 + t608 * t407;
t688 = t520 * qJ(6) + t396;
t393 = qJD(6) * t566 + t688;
t759 = pkin(5) * t520;
t394 = -t672 - t759;
t450 = -t486 * t610 + t613 * t493;
t440 = -qJ(5) * t525 + t450;
t436 = pkin(4) * t566 + t440;
t441 = -qJ(5) * t523 + t451;
t748 = t441 * t607;
t410 = t436 * t608 - t748;
t408 = -pkin(5) * t566 + qJD(6) - t410;
t438 = t608 * t441;
t411 = t607 * t436 + t438;
t409 = qJ(6) * t566 + t411;
t649 = -t608 * t613 + t743;
t683 = -g(1) * t736 - g(2) * t738 + t600;
t764 = t393 * t650 + t394 * t649 + t408 * t721 + t409 * t722 + t683;
t763 = -t395 * t649 + t396 * t650 - t410 * t721 + t411 * t722 + t683;
t761 = pkin(2) * t611;
t758 = g(1) * t612;
t755 = g(2) * t615;
t753 = pkin(7) * qJDD(2);
t751 = qJDD(2) * pkin(2);
t413 = t440 * t607 + t438;
t750 = t413 * t652;
t747 = t464 * t613;
t746 = t520 * t610;
t745 = t523 * t566;
t744 = t525 * t566;
t740 = t610 * t614;
t737 = t611 * t613;
t618 = qJD(1) ^ 2;
t735 = t611 * t618;
t734 = t612 * t613;
t733 = t612 * t614;
t732 = t613 * t614;
t731 = t613 * t615;
t730 = t614 * t615;
t710 = qJD(2) * t611;
t587 = pkin(2) * t710;
t480 = t587 + t625;
t708 = qJD(2) * t614;
t532 = t762 * t708;
t515 = t613 * t532;
t667 = qJ(5) * t614 - t518;
t701 = qJD(5) * t614;
t419 = pkin(4) * t708 + t515 + t667 * t704 + (-qJ(5) * t710 - qJD(4) * t547 - t480 + t701) * t610;
t637 = t613 * t480 - t518 * t705 + t610 * t532 + t547 * t704;
t424 = qJ(5) * t774 - t613 * t701 + t637;
t401 = t607 * t419 + t608 * t424;
t726 = pkin(5) * t722 + qJ(6) * t721 + qJD(6) * t649 + t765;
t421 = t607 * t446 + t608 * t453;
t724 = pkin(5) * t714 - t725;
t416 = qJ(6) * t714 + t421;
t723 = t455 - t416;
t527 = t613 * t547;
t456 = pkin(4) * t611 + t610 * t667 + t527;
t461 = -qJ(5) * t732 + t719;
t429 = t607 * t456 + t608 * t461;
t691 = t611 * t734;
t718 = pkin(4) * t691 + t569;
t548 = t762 * t614;
t605 = t611 ^ 2;
t606 = t614 ^ 2;
t716 = t605 - t606;
t713 = qJD(2) * t523;
t712 = qJD(2) * t525;
t706 = qJD(4) * t486;
t702 = qJD(4) * t616;
t414 = t440 * t608 - t748;
t698 = qJD(6) - t414;
t693 = g(3) * t732;
t692 = t610 * t736;
t690 = t611 * t731;
t689 = t614 * t735;
t598 = t615 * pkin(7);
t686 = t615 * t579 + t609 * t733 + t598;
t581 = pkin(7) * t695;
t602 = qJDD(2) * qJ(3);
t603 = qJD(2) * qJD(3);
t685 = t581 + t602 + t603;
t568 = pkin(4) * t732;
t684 = t568 + t548;
t681 = t610 * t710;
t668 = -qJD(2) * pkin(2) + qJD(3);
t665 = pkin(4) * t678;
t664 = t615 * pkin(1) + pkin(2) * t730 + t612 * pkin(7) + qJ(3) * t736;
t663 = -t580 - t683;
t662 = pkin(3) * t695 + t685;
t530 = t762 * t710;
t561 = qJ(3) * t733;
t661 = -pkin(2) * t738 + t561;
t563 = qJ(3) * t730;
t660 = -pkin(2) * t736 + t563;
t572 = g(1) * t733;
t659 = -g(2) * t730 + t572;
t617 = qJD(2) ^ 2;
t658 = pkin(7) * t617 + t755;
t656 = pkin(5) * t589 - qJ(6) * t590;
t400 = t419 * t608 - t424 * t607;
t428 = t456 * t608 - t461 * t607;
t535 = t583 + t668;
t541 = -t584 - t604;
t651 = t535 * t614 + t541 * t611;
t648 = qJD(1) * t530;
t647 = qJD(2) * (-pkin(7) - t579);
t646 = t566 * t610;
t645 = t671 - t597;
t644 = pkin(4) * t742 - t609 * t614 + t717;
t642 = -0.2e1 * pkin(1) * t697 - t753;
t507 = t645 * qJD(1);
t641 = t507 * t715 + qJDD(3) - t663;
t640 = -t566 * t704 - t746;
t639 = -qJ(3) * t708 - t707;
t496 = t650 * t614;
t635 = 0.2e1 * qJDD(1) * pkin(1) - t658;
t633 = t662 + t770;
t629 = t753 + (-qJD(1) * t536 - t507) * qJD(2);
t626 = pkin(4) * t692 + t612 * t579 - t609 * t730 + t664;
t474 = -t648 + t662;
t624 = t474 - t766;
t466 = qJD(1) * t639 + qJDD(1) * t645 + t565;
t500 = t587 + t639;
t622 = qJD(1) * t500 + qJDD(1) * t536 + t466 + t658;
t487 = pkin(7) * t675 - t685;
t497 = t673 - t751;
t620 = qJD(2) * t651 - t487 * t614 + t497 * t611;
t619 = pkin(5) * t431 - qJ(6) * t432 - qJD(6) * t652 + t633;
t578 = -pkin(4) * t608 - pkin(5);
t574 = pkin(4) * t607 + qJ(6);
t545 = pkin(4) * t690;
t528 = -qJ(3) * t714 + t588;
t511 = -t610 * t738 + t731;
t510 = t691 + t739;
t509 = t692 + t734;
t508 = t690 - t741;
t495 = t649 * t614;
t492 = -t589 * t738 + t590 * t615;
t490 = t589 * t736 + t590 * t612;
t467 = pkin(5) * t650 + qJ(6) * t649 + t728;
t458 = -t607 * t774 - t608 * t681 + t614 * t677;
t457 = qJD(4) * t496 - t649 * t710;
t443 = -pkin(5) * t495 + qJ(6) * t496 + t684;
t435 = t474 + t770;
t430 = pkin(4) * t525 + pkin(5) * t652 + qJ(6) * t470;
t427 = -pkin(5) * t611 - t428;
t426 = qJ(6) * t611 + t429;
t412 = -pkin(5) * t457 + qJ(6) * t458 + qJD(6) * t496 + t611 * t647 - t665;
t399 = -pkin(5) * t708 - t400;
t398 = qJ(6) * t708 + qJD(6) * t611 + t401;
t397 = -t648 + t619;
t1 = [(-t637 * t566 - t719 * t520 - t530 * t525 + t548 * t464 + g(1) * t510 - g(2) * t508 + ((qJD(2) * t506 + t706) * t610 + t687) * t611 + (-qJD(2) * t451 - t474 * t610 - t506 * t704) * t614) * MDP(21) + (t395 * t496 + t396 * t495 - t400 * t652 - t401 * t470 + t410 * t458 + t411 * t457 - t428 * t432 - t429 * t431 + t659) * MDP(22) + (t393 * t495 - t394 * t496 - t398 * t470 + t399 * t652 - t408 * t458 + t409 * t457 - t426 * t431 + t427 * t432 + t659) * MDP(25) + (-g(1) * t491 - g(2) * t489 + t393 * t611 + t397 * t496 + t398 * t566 + t409 * t708 - t412 * t652 + t423 * t458 + t426 * t520 - t432 * t443) * MDP(26) + ((-t480 * t610 + t515) * t566 + (-t518 * t610 + t527) * t520 + t666 * t611 - t530 * t523 + t548 * t465 + t474 * t732 - g(1) * t511 - g(2) * t509 + (t450 * t614 - t506 * t737) * qJD(2) + (-t451 * t611 - t506 * t740 - t566 * t719) * qJD(4)) * MDP(20) + ((t566 * t709 - t465) * t611 + (-t713 - t773) * t614) * MDP(18) + (t393 * t426 + t409 * t398 + t397 * t443 + t423 * t412 + t394 * t427 + t408 * t399 - g(1) * (pkin(5) * t492 + qJ(6) * t491 + t686) - g(2) * (pkin(5) * t490 + qJ(6) * t489 + t626) - (-t611 * t728 - pkin(1) - t597) * t758) * MDP(27) + (t396 * t429 + t411 * t401 + t395 * t428 + t410 * t400 + t435 * t684 - t475 * t665 - g(1) * (-pkin(1) * t612 - pkin(2) * t733 + t686) - g(2) * t626 + (t475 * t647 + t728 * t758) * t611) * MDP(23) + (pkin(7) * t620 - g(1) * t598 - g(2) * t664 + t466 * t536 + t507 * t500 - t645 * t758) * MDP(14) + (t642 * t614 + (-t635 - t758) * t611) * MDP(10) + (-t755 + t758) * MDP(2) + (t629 * t614 + (-t622 + t758) * t611) * MDP(13) + ((-t523 * t610 + t525 * t613) * t710 + (-t747 + t465 * t610 + (t523 * t613 + t525 * t610) * qJD(4)) * t614) * MDP(16) + (-t464 * t740 + (-t613 * t703 + t681) * t525) * MDP(15) + 0.2e1 * (t611 * t695 - t697 * t716) * MDP(5) + ((t605 + t606) * qJDD(1) * pkin(7) + t620 - t657) * MDP(11) + t657 * MDP(3) + (t611 * t629 + t614 * t622 - t572) * MDP(12) + (t520 * t611 + t566 * t708) * MDP(19) + ((t566 * t711 + t464) * t611 + (t640 + t712) * t614) * MDP(17) + (qJDD(1) * t605 + 0.2e1 * t611 * t674) * MDP(4) + (-g(1) * t492 - g(2) * t490 - t394 * t611 - t397 * t495 - t399 * t566 - t408 * t708 + t412 * t470 - t423 * t457 - t427 * t520 + t431 * t443) * MDP(24) + (t611 * t642 + t614 * t635 + t572) * MDP(9) + qJDD(1) * MDP(1) + (qJDD(2) * t611 + t614 * t617) * MDP(6) + (qJDD(2) * t614 - t611 * t617) * MDP(7); t716 * MDP(5) * t618 - MDP(4) * t689 - t566 * MDP(19) * t714 + MDP(6) * t696 + MDP(7) * t695 + qJDD(2) * MDP(8) + (pkin(1) * t735 + t663) * MDP(9) + (t754 - t581 + (pkin(1) * t618 + t657) * t614) * MDP(10) + ((t752 - t761) * qJDD(1) + ((-t541 - t604) * t611 + (-t535 + t668) * t614) * qJD(1)) * MDP(11) + (-t528 * t714 + t641 - 0.2e1 * t751) * MDP(12) + (t581 + 0.2e1 * t602 + 0.2e1 * t603 + (qJD(1) * t528 - g(3)) * t611 + (qJD(1) * t507 - t657) * t614) * MDP(13) + (-pkin(7) * qJD(1) * t651 - t497 * pkin(2) - g(1) * t660 - g(2) * t661 - g(3) * t717 - t487 * qJ(3) - t541 * qJD(3) - t507 * t528) * MDP(14) + (-t525 * t646 + t747) * MDP(15) + ((-t465 - t744) * t613 + (-t464 + t745) * t610) * MDP(16) + ((-t525 * t614 - t566 * t742) * qJD(1) + t773) * MDP(17) + ((t523 * t614 - t566 * t737) * qJD(1) + t640) * MDP(18) + (-t450 * t714 + qJ(3) * t465 - t514 * t566 + t700 * t523 + t767 * t613 + ((t498 - t702) * t566 + t624) * t610) * MDP(20) + (qJ(3) * t464 + t720 * t566 + t451 * t714 + t700 * t525 - t767 * t610 + (-t566 * t702 + t624) * t613) * MDP(21) + (t421 * t470 - t652 * t725 - t763 + t769) * MDP(22) + (t396 * t477 - t395 * t476 + t435 * t728 - g(1) * (t660 + t775) - g(2) * (t661 + t776) - g(3) * t644 + t765 * t475 + (t455 - t421) * t411 + t725 * t410) * MDP(23) + (t397 * t650 + t408 * t714 + t423 * t722 + t431 * t467 + t470 * t726 - t476 * t520 - t566 * t724 - t589 * t766) * MDP(24) + (t416 * t470 + t652 * t724 - t764 + t769) * MDP(25) + (t397 * t649 - t409 * t714 + t721 * t423 - t432 * t467 + t477 * t520 + t723 * t566 + t590 * t766 - t726 * t652) * MDP(26) + (t393 * t477 + t397 * t467 + t394 * t476 - g(1) * (t563 + t775) - g(2) * (t561 + t776) - g(3) * (t611 * t656 + t644) + t726 * t423 + t723 * t409 + t724 * t408 + t657 * (-t614 * t656 + t761)) * MDP(27); MDP(11) * t696 + (qJDD(2) + t689) * MDP(12) + (-t605 * t618 - t617) * MDP(13) + (qJD(2) * t541 + t564 + t641 - t751) * MDP(14) + (t501 - t713) * MDP(20) + (-t712 - t746) * MDP(21) + (-qJD(2) * t475 + t763) * MDP(23) + (-qJD(2) * t470 - t520 * t649) * MDP(24) + t650 * t520 * MDP(26) + (-qJD(2) * t423 + t764) * MDP(27) + (qJD(2) * MDP(26) + t721 * t772) * t652 + (-MDP(21) * t566 * t613 - MDP(20) * t646 - MDP(24) * t721 + MDP(26) * t722) * t566 + t772 * (-t650 * t431 + t432 * t649 - t470 * t722); t525 * t523 * MDP(15) + (-t523 ^ 2 + t525 ^ 2) * MDP(16) + (t464 + t745) * MDP(17) + (-t465 + t744) * MDP(18) + t520 * MDP(19) + (-g(1) * t508 - g(2) * t510 + t451 * t566 - t506 * t525 + t623 + t693) * MDP(20) + (g(1) * t509 - g(2) * t511 + t450 * t566 + t506 * t523 + (t706 - t600) * t610 + t687) * MDP(21) + (t411 * t652 - t750 + (-t431 * t607 - t432 * t608) * pkin(4) + (-t410 + t414) * t470) * MDP(22) + (-t411 * t414 + t410 * t413 - g(1) * t545 - g(2) * t718 + (g(1) * t741 + t395 * t608 + t396 * t607 - t475 * t525 + t693) * pkin(4)) * MDP(23) + (t413 * t566 - t430 * t470 + (pkin(5) - t578) * t520 + t768) * MDP(24) + (t409 * t652 - t431 * t574 + t432 * t578 - t750 + (t408 - t698) * t470) * MDP(25) + (t589 * t600 - g(1) * t490 + g(2) * t492 - t423 * t470 + t430 * t652 + t520 * t574 + (0.2e1 * qJD(6) - t414) * t566 + t688) * MDP(26) + (t393 * t574 + t394 * t578 - t423 * t430 - t408 * t413 - g(1) * (-pkin(5) * t489 + qJ(6) * t490 + t545 - t694) - g(2) * (pkin(5) * t491 - qJ(6) * t492 + t718) - g(3) * (-t568 + (-pkin(5) * t590 - qJ(6) * t589) * t614) + t698 * t409) * MDP(27); (t410 * t652 + t411 * t470 + t633) * MDP(23) + (t566 * t652 + t431) * MDP(24) + (-t432 + t780) * MDP(26) + (-t408 * t652 + t409 * t470 + t619) * MDP(27) + t772 * (-t470 ^ 2 - t778) + (t777 + (t697 * t762 + g(3)) * t611) * (-MDP(23) - MDP(27)); (t470 * t652 - t520) * MDP(24) + (t432 + t780) * MDP(25) + (-t566 ^ 2 - t778) * MDP(26) + (-t409 * t566 - t759 - t768) * MDP(27);];
tau  = t1;

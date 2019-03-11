% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:34
% EndTime: 2019-03-09 10:06:45
% DurationCPUTime: 10.26s
% Computational Cost: add. (4284->612), mult. (8516->702), div. (0->0), fcn. (4952->6), ass. (0->263)
t668 = MDP(23) - MDP(28);
t577 = cos(qJ(4));
t574 = sin(qJ(4));
t687 = qJD(2) * t574;
t578 = cos(qJ(2));
t690 = qJD(1) * t578;
t489 = t577 * t690 + t687;
t575 = sin(qJ(2));
t671 = qJD(1) * qJD(2);
t651 = t575 * t671;
t669 = qJDD(1) * t578;
t753 = -t651 + t669;
t425 = qJD(4) * t489 - t577 * qJDD(2) + t574 * t753;
t691 = qJD(1) * t575;
t530 = qJD(4) + t691;
t718 = t489 * t530;
t755 = t425 - t718;
t760 = t755 * t668;
t759 = t425 + t718;
t652 = t574 * t690;
t426 = -qJD(4) * t652 + qJDD(2) * t574 + (qJD(2) * qJD(4) + t753) * t577;
t685 = qJD(2) * t577;
t491 = -t652 + t685;
t717 = t491 * t530;
t758 = -t426 - t717;
t757 = -t426 + t717;
t746 = MDP(22) + MDP(26);
t756 = MDP(24) + MDP(27);
t738 = pkin(3) + pkin(7);
t739 = pkin(2) + pkin(8);
t704 = qJ(6) - t739;
t650 = t578 * t671;
t670 = qJDD(1) * t575;
t606 = t650 + t670;
t488 = qJDD(4) + t606;
t754 = t488 * qJ(5) + t530 * qJD(5);
t548 = pkin(7) * t691;
t675 = pkin(3) * t691 + qJD(3) + t548;
t752 = t426 * qJ(6) + t489 * qJD(6);
t681 = qJD(4) * t577;
t751 = -t574 * t488 - t530 * t681;
t750 = 0.2e1 * t754;
t467 = t577 * t488;
t682 = qJD(4) * t574;
t499 = t530 * t682;
t749 = -t499 + t467;
t747 = -qJD(5) * t577 + t675;
t576 = sin(qJ(1));
t579 = cos(qJ(1));
t634 = g(1) * t579 + g(2) * t576;
t557 = t575 * qJ(3);
t648 = -pkin(1) - t557;
t603 = -t578 * t739 + t648;
t449 = t603 * qJD(1);
t452 = -qJD(2) * t739 + t675;
t411 = -t574 * t449 + t577 * t452;
t673 = qJD(5) - t411;
t479 = t488 * pkin(4);
t745 = t479 - qJDD(5);
t567 = g(3) * t575;
t707 = t578 * t579;
t710 = t576 * t578;
t744 = -g(1) * t707 - g(2) * t710 - t567;
t550 = pkin(7) * t690;
t497 = pkin(3) * t690 + t550;
t571 = qJD(2) * qJ(3);
t469 = t571 + t497;
t616 = qJ(5) * t491 - t469;
t737 = pkin(4) + pkin(5);
t401 = -t489 * t737 + qJD(6) + t616;
t708 = t577 * t579;
t471 = t574 * t576 - t575 * t708;
t711 = t576 * t577;
t473 = t574 * t579 + t575 * t711;
t529 = pkin(2) * t651;
t730 = qJ(3) * t578;
t626 = pkin(8) * t575 - t730;
t683 = qJD(3) * t575;
t594 = qJD(2) * t626 - t683;
t410 = qJD(1) * t594 + qJDD(1) * t603 + t529;
t527 = pkin(7) * t650;
t545 = pkin(7) * t670;
t649 = qJDD(3) + t527 + t545;
t430 = pkin(3) * t606 - qJDD(2) * t739 + t649;
t640 = -t574 * t410 + t577 * t430 - t449 * t681 - t452 * t682;
t709 = t577 * t578;
t596 = g(1) * t471 - g(2) * t473 + g(3) * t709 + t640;
t592 = t596 + t745;
t725 = qJ(6) * t425;
t743 = (qJD(6) + t401) * t491 + t592 - t725;
t726 = qJ(5) * t577;
t742 = t574 * t737 - t726;
t727 = qJ(5) * t574;
t610 = t577 * t737 + t727;
t741 = -0.2e1 * pkin(1);
t740 = t489 ^ 2;
t487 = t491 ^ 2;
t521 = t530 ^ 2;
t736 = pkin(5) * t488;
t582 = qJD(2) ^ 2;
t735 = pkin(7) * t582;
t734 = g(1) * t576;
t568 = g(3) * t578;
t563 = t578 * pkin(2);
t731 = pkin(7) * qJDD(2);
t729 = qJ(5) * t426;
t728 = qJ(5) * t489;
t724 = qJDD(2) * pkin(2);
t546 = pkin(7) * t669;
t569 = qJDD(2) * qJ(3);
t570 = qJD(2) * qJD(3);
t450 = pkin(7) * t651 - t546 - t569 - t570;
t431 = pkin(3) * t753 - t450;
t387 = t426 * pkin(4) + t425 * qJ(5) - t491 * qJD(5) + t431;
t723 = t387 * t577;
t412 = t577 * t449 + t574 * t452;
t722 = t412 * t530;
t721 = t425 * t577;
t720 = t431 * t577;
t719 = t489 * t491;
t716 = t574 * t575;
t715 = t574 * t578;
t714 = t575 * t576;
t713 = t575 * t577;
t712 = t575 * t579;
t583 = qJD(1) ^ 2;
t706 = t578 * t583;
t705 = t739 * t488;
t703 = -t530 * t610 - t747;
t553 = pkin(2) * t691;
t458 = qJD(1) * t626 + t553;
t702 = t577 * t458 + t574 * t497;
t642 = qJD(4) * t704;
t447 = t574 * t458;
t644 = -t497 * t577 + t447;
t676 = qJD(6) * t577;
t701 = t574 * t642 - t676 - (-qJ(6) * t716 - t578 * t737) * qJD(1) - t644;
t416 = qJ(5) * t690 + t702;
t656 = t577 * t691;
t677 = qJD(6) * t574;
t700 = qJ(6) * t656 + t577 * t642 - t416 + t677;
t629 = pkin(4) * t577 + t727;
t699 = t530 * t629 + t747;
t693 = t563 + t557;
t659 = -t578 * pkin(8) - t693;
t482 = -pkin(1) + t659;
t512 = t738 * t575;
t698 = t577 * t482 + t574 * t512;
t523 = qJ(3) * t710;
t666 = pkin(4) * t715;
t697 = t576 * t666 + t523;
t525 = qJ(3) * t707;
t696 = t579 * t666 + t525;
t695 = -g(1) * t710 + g(2) * t707;
t513 = t738 * t578;
t572 = t575 ^ 2;
t573 = t578 ^ 2;
t692 = t572 - t573;
t689 = qJD(2) * t489;
t688 = qJD(2) * t491;
t686 = qJD(2) * t575;
t684 = qJD(2) * t578;
t680 = qJD(4) * t578;
t679 = qJD(4) * t739;
t678 = qJD(5) * t574;
t399 = qJ(6) * t491 + t411;
t674 = qJD(5) - t399;
t665 = t575 * t706;
t664 = t576 * t709;
t663 = t577 * t707;
t662 = -t577 * t410 - t574 * t430 - t452 * t681;
t428 = t575 * qJ(5) + t698;
t661 = t530 * t656 - t751;
t660 = g(1) * t663 + g(2) * t664 + g(3) * t713;
t658 = -g(1) * t712 - g(2) * t714 + t568;
t657 = t574 * t691;
t655 = t574 * t686;
t654 = t401 * t681;
t653 = t577 * t679;
t647 = -qJD(2) * pkin(2) + qJD(3);
t472 = t574 * t712 + t711;
t646 = -t471 * pkin(4) + qJ(5) * t472;
t474 = -t574 * t714 + t708;
t645 = t473 * pkin(4) - qJ(5) * t474;
t643 = -t574 * t482 + t512 * t577;
t639 = t579 * pkin(1) + pkin(2) * t707 + t576 * pkin(7) + qJ(3) * t712;
t638 = -t545 - t658;
t637 = t546 + t744;
t636 = -g(1) * t473 - g(2) * t471;
t635 = -g(1) * t474 - g(2) * t472;
t633 = -g(2) * t579 + t734;
t400 = qJ(6) * t489 + t412;
t605 = -t449 * t682 - t662;
t385 = t605 + t754;
t402 = -pkin(4) * t530 + t673;
t632 = -t402 * t691 - t385;
t386 = -t640 - t745;
t518 = t530 * qJ(5);
t405 = t518 + t412;
t631 = -t405 * t691 + t386;
t628 = -pkin(4) * t574 + t726;
t627 = pkin(5) * t574 - t726;
t624 = t402 * t574 + t405 * t577;
t503 = t548 + t647;
t506 = -t550 - t571;
t623 = t503 * t578 + t506 * t575;
t622 = g(3) * (pkin(4) * t716 - t659);
t619 = qJDD(2) * t575 + t578 * t582;
t618 = t648 - t563;
t564 = t579 * pkin(7);
t617 = t579 * pkin(3) + t474 * pkin(4) + qJ(5) * t473 + t564;
t552 = pkin(2) * t686;
t444 = t552 + t594;
t498 = t738 * t684;
t615 = -t574 * t444 - t482 * t681 + t498 * t577 - t512 * t682;
t470 = t618 * qJD(1);
t609 = t470 * t691 + qJDD(3) - t638;
t384 = -pkin(5) * t426 + qJDD(6) - t387;
t608 = t384 * t577 - t401 * t682;
t607 = -qJ(3) * t684 - t683;
t604 = t577 * t444 - t482 * t682 + t574 * t498 + t512 * t681;
t602 = t610 * t578;
t383 = t385 + t752;
t394 = -t530 * t737 + t674;
t397 = t400 + t518;
t601 = t383 * t574 + t658 + (t656 + t681) * t397 + (t657 + t682) * t394;
t505 = -pkin(1) - t693;
t600 = t731 + (-qJD(1) * t505 - t470) * qJD(2);
t414 = pkin(4) * t489 - t616;
t599 = t414 * t530 - t705;
t391 = qJ(5) * t684 + t575 * qJD(5) + t604;
t427 = qJD(1) * t607 + qJDD(1) * t618 + t529;
t461 = t552 + t607;
t598 = qJD(1) * t461 + qJDD(1) * t505 + t427 + t735;
t597 = -t578 * t634 - t567;
t595 = t576 * pkin(3) + t472 * pkin(4) + pkin(8) * t707 + qJ(5) * t471 + t639;
t591 = g(1) * t472 - g(2) * t474 - g(3) * t715 - t605;
t457 = t649 - t724;
t590 = qJD(2) * t623 - t450 * t578 + t457 * t575;
t589 = t530 * t679 + t597;
t586 = t414 * t491 - t592;
t584 = t411 * t530 + t591;
t504 = qJ(3) - t628;
t502 = t704 * t577;
t501 = t704 * t574;
t496 = t738 * t686;
t494 = -qJ(3) * t690 + t553;
t478 = -qJ(3) - t742;
t443 = t578 * t629 + t513;
t433 = pkin(4) * t491 + t728;
t432 = -t602 - t513;
t429 = -pkin(4) * t575 - t643;
t418 = -pkin(4) * t690 + t644;
t417 = qJ(6) * t709 + t428;
t415 = -t491 * t737 - t728;
t413 = qJ(6) * t715 - t575 * t737 - t643;
t406 = (qJD(4) * t628 + t678) * t578 + (-t629 - t738) * t686;
t398 = (qJD(4) * t742 - t678) * t578 + (t610 + t738) * t686;
t393 = -pkin(4) * t684 - t615;
t389 = t578 * t676 + (-t574 * t680 - t575 * t685) * qJ(6) + t391;
t388 = -qJ(6) * t655 + (qJ(6) * t681 - qJD(2) * t737 + t677) * t578 - t615;
t381 = -qJD(6) * t491 + t386 + t725 - t736;
t1 = [(t383 * t417 + t397 * t389 + t381 * t413 + t394 * t388 + t384 * t432 + t401 * t398 - g(1) * (pkin(5) * t474 + t617) - g(2) * (pkin(5) * t472 - qJ(6) * t707 + t595) - (t578 * t704 + t648) * t734) * MDP(29) + ((t572 + t573) * qJDD(1) * pkin(7) + t590 - t634) * MDP(11) + t634 * MDP(3) + (t600 * t578 + (-t598 + t633) * t575) * MDP(13) + t633 * MDP(2) + t619 * MDP(6) + ((t530 * t685 - t426) * t575 + (-t689 - t749) * t578) * MDP(18) + ((t530 * t687 - t425) * t575 + (t688 + t751) * t578) * MDP(17) + (t615 * t530 + t643 * t488 - t496 * t489 + t513 * t426 + (-t469 * t685 + t640) * t575 + (qJD(2) * t411 - t469 * t682 + t720) * t578 + t635) * MDP(20) + (-t604 * t530 - t698 * t488 - t496 * t491 - t513 * t425 + ((qJD(2) * t469 + qJD(4) * t449) * t574 + t662) * t575 + (-qJD(2) * t412 - t431 * t574 - t469 * t681) * t578 - t636) * MDP(21) + (qJDD(2) * t578 - t575 * t582) * MDP(7) + (t425 * t715 + (-t577 * t680 + t655) * t491) * MDP(15) + (t575 * t600 + t578 * t598 + t695) * MDP(12) + ((-t489 * t574 + t491 * t577) * t686 + (t721 + t426 * t574 + (t489 * t577 + t491 * t574) * qJD(4)) * t578) * MDP(16) + (qJDD(1) * t572 + 0.2e1 * t575 * t650) * MDP(4) + (t488 * t575 + t530 * t684) * MDP(19) + (0.2e1 * pkin(1) * t753 - t619 * pkin(7) - t695) * MDP(9) + qJDD(1) * MDP(1) + (-t388 * t530 - t398 * t489 - t413 * t488 - t426 * t432 + (t401 * t685 - t381) * t575 + (-qJD(2) * t394 - t608) * t578 + t635) * MDP(26) + (t391 * t530 - t406 * t491 + t425 * t443 + t428 * t488 + (-t414 * t687 + t385) * t575 + (qJD(2) * t405 + t387 * t574 + t414 * t681) * t578 + t636) * MDP(24) + (t389 * t530 + t398 * t491 + t417 * t488 - t425 * t432 + (t401 * t687 + t383) * t575 + (qJD(2) * t397 - t384 * t574 - t654) * t578 + t636) * MDP(27) + (-t393 * t530 + t406 * t489 + t426 * t443 - t429 * t488 + (-t414 * t685 - t386) * t575 + (-qJD(2) * t402 - t414 * t682 + t723) * t578 + t635) * MDP(22) + (pkin(7) * t590 - g(1) * t564 - g(2) * t639 + t427 * t505 + t470 * t461 - t618 * t734) * MDP(14) + (-g(1) * t617 - g(2) * t595 + t385 * t428 + t386 * t429 + t387 * t443 + t405 * t391 + t402 * t393 + t414 * t406 - t603 * t734) * MDP(25) + ((t671 * t741 - t731) * t578 + (qJDD(1) * t741 - t633 + t735) * t575) * MDP(10) + 0.2e1 * (t575 * t669 - t671 * t692) * MDP(5) + (-t388 * t491 + t389 * t489 + t413 * t425 + t417 * t426 + (-t394 * t574 - t397 * t577) * t686 + (t381 * t574 + t383 * t577 + (t394 * t577 - t397 * t574) * qJD(4)) * t578 + t695) * MDP(28) + (-t391 * t489 + t393 * t491 - t425 * t429 - t426 * t428 + t624 * t686 + (-t385 * t577 - t386 * t574 + (-t402 * t577 + t405 * t574) * qJD(4)) * t578 - t695) * MDP(23); (t383 * t501 - t381 * t502 + t384 * t478 - g(1) * t696 - g(2) * t697 - t622 + (g(3) * qJ(6) - t627 * t634) * t578 + t703 * t401 + t700 * t397 + t701 * t394 + (-g(3) * t627 - t634 * t704) * t575) * MDP(29) + (t387 * t504 - t405 * t416 - t402 * t418 - g(1) * (-qJ(5) * t663 + t696) - g(2) * (-qJ(5) * t664 + t697) - t622 + (g(3) * t726 + t634 * t739) * t575 + t699 * t414 - (qJD(4) * t624 + t385 * t574 - t386 * t577) * t739) * MDP(25) + (t416 * t489 - t418 * t491 + (-t425 * t739 + (t489 * t739 - t405) * qJD(4) + t631) * t577 + (t426 * t739 + (-t491 * t739 - t402) * qJD(4) + t632) * t574 - t658) * MDP(23) + (-t574 * t717 - t721) * MDP(15) + (t412 * t690 - qJ(3) * t425 + t720 + (t653 + t702) * t530 + t675 * t491 + (-t469 * t530 + t705) * t574 - t660) * MDP(21) + (-t405 * t690 - t723 + t425 * t504 + (-t416 - t653) * t530 - t699 * t491 + t599 * t574 + t660) * MDP(24) + (t759 * t574 + t758 * t577) * MDP(16) + (pkin(1) * t575 * t583 + t638) * MDP(9) + (0.2e1 * t569 + 0.2e1 * t570 + (t470 * t578 + t494 * t575) * qJD(1) + t637) * MDP(13) + ((-t491 * t578 - t530 * t716) * qJD(1) + t749) * MDP(17) + (-t654 - t426 * t478 + t488 * t502 - t701 * t530 - t703 * t489 + (t394 * t578 - t401 * t713) * qJD(1) + (-t384 + t597) * t574) * MDP(26) + (-t381 * t577 - t425 * t502 + t426 * t501 + t489 * t700 - t491 * t701 + t601) * MDP(28) + (-t411 * t690 + qJ(3) * t426 + t675 * t489 - t705 * t577 + (t431 + t589) * t574 + (t447 + (t469 - t497) * t577) * t530) * MDP(20) + (-t450 * qJ(3) - t506 * qJD(3) - t457 * pkin(2) - t470 * t494 - g(1) * (-pkin(2) * t712 + t525) - g(2) * (-pkin(2) * t714 + t523) - g(3) * t693 - t623 * qJD(1) * pkin(7)) * MDP(14) + (-t425 * t478 + t488 * t501 + t700 * t530 + t703 * t491 + (-t397 * t578 - t401 * t716) * qJD(1) + t608 + t660) * MDP(27) + MDP(7) * t669 + MDP(6) * t670 + (pkin(1) * t706 - t637) * MDP(10) + (t402 * t690 + t418 * t530 + t426 * t504 + t699 * t489 + t599 * t577 + (t387 + t589) * t574) * MDP(22) - t530 * MDP(19) * t690 + t692 * MDP(5) * t583 - MDP(4) * t665 + qJDD(2) * MDP(8) + (-t494 * t690 + t609 - 0.2e1 * t724) * MDP(12) + ((-pkin(2) * t575 + t730) * qJDD(1) + ((-t506 - t571) * t575 + (-t503 + t647) * t578) * qJD(1)) * MDP(11) + (t489 * t690 - t661) * MDP(18); MDP(11) * t670 + (qJDD(2) + t665) * MDP(12) + (-t572 * t583 - t582) * MDP(13) + (qJD(2) * t506 + t527 + t609 - t724) * MDP(14) + (t467 - t689) * MDP(20) - MDP(21) * t688 + (-qJD(2) * t414 + t658) * MDP(25) + (qJD(2) * t401 + t601) * MDP(29) + ((qJD(4) * t405 - t631) * MDP(25) - t381 * MDP(29) - MDP(21) * t521 + t746 * t488 + t760) * t577 + (-t488 * MDP(21) + (qJD(4) * t402 - t632) * MDP(25) - MDP(20) * t521 + t757 * t668) * t574 + t746 * (-t530 * t657 - t499 - t689) + t756 * (t661 + t688); MDP(15) * t719 + (t487 - t740) * MDP(16) - t755 * MDP(17) + t757 * MDP(18) + t488 * MDP(19) + (-t469 * t491 + t596 + t722) * MDP(20) + (t469 * t489 + t584) * MDP(21) + (-t433 * t489 + t479 - t586 + t722) * MDP(22) + (pkin(4) * t425 - t729 + (t405 - t412) * t491 + (t402 - t673) * t489) * MDP(23) + (-t414 * t489 + t433 * t491 - t584 + t750) * MDP(24) + (-t386 * pkin(4) - g(1) * t646 - g(2) * t645 + t385 * qJ(5) - t402 * t412 + t405 * t673 - t414 * t433 + t568 * t629) * MDP(25) + (t400 * t530 + t415 * t489 + (pkin(5) + t737) * t488 + t743) * MDP(26) + (-t399 * t530 + t401 * t489 - t415 * t491 - t591 + t750 + t752) * MDP(27) + (t729 - t425 * t737 + (-t397 + t400) * t491 + (-t394 + t674) * t489) * MDP(28) + (t383 * qJ(5) - t381 * t737 - t394 * t400 - t401 * t415 - g(1) * (-pkin(5) * t471 + t646) - g(2) * (pkin(5) * t473 + t645) + t674 * t397 + g(3) * t602) * MDP(29); (-t405 * t530 + t586) * MDP(25) + (-t397 * t530 - t736 - t743) * MDP(29) + t746 * (-t488 + t719) + t756 * (-t521 - t487) - t760; t758 * MDP(26) - t759 * MDP(27) + (-t487 - t740) * MDP(28) + (t394 * t491 - t397 * t489 + t384 - t744) * MDP(29);];
tau  = t1;

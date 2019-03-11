% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:48:00
% EndTime: 2019-03-09 10:48:12
% DurationCPUTime: 9.86s
% Computational Cost: add. (5872->529), mult. (12928->668), div. (0->0), fcn. (9185->12), ass. (0->248)
t617 = cos(qJ(4));
t618 = cos(qJ(2));
t613 = sin(qJ(4));
t614 = sin(qJ(2));
t725 = t613 * t614;
t525 = t617 * t618 + t725;
t509 = t525 * qJD(1);
t706 = qJD(1) * t618;
t707 = qJD(1) * t614;
t511 = -t613 * t706 + t617 * t707;
t608 = sin(pkin(10));
t609 = cos(pkin(10));
t655 = -t509 * t608 + t609 * t511;
t639 = t525 * qJD(4);
t695 = qJD(1) * qJD(2);
t683 = t614 * t695;
t693 = qJDD(1) * t618;
t763 = t683 - t693;
t682 = t618 * t695;
t694 = qJDD(1) * t614;
t764 = t682 + t694;
t443 = -qJD(1) * t639 + t613 * t763 + t617 * t764;
t702 = qJD(4) * t617;
t703 = qJD(4) * t613;
t704 = qJD(2) * t618;
t777 = t613 * t704 + t614 * t702 - t618 * t703;
t444 = qJD(1) * t777 + qJDD(1) * t525 - t617 * t683;
t423 = t443 * t609 - t444 * t608;
t599 = qJDD(2) - qJDD(4);
t600 = qJD(2) - qJD(4);
t612 = sin(qJ(6));
t616 = cos(qJ(6));
t698 = qJD(6) * t616;
t689 = t616 * t423 - t612 * t599 - t600 * t698;
t699 = qJD(6) * t612;
t409 = -t655 * t699 + t689;
t447 = -t600 * t612 + t616 * t655;
t561 = t616 * t599;
t410 = qJD(6) * t447 + t423 * t612 + t561;
t733 = t447 * t655;
t737 = t409 * t612;
t422 = -t443 * t608 - t609 * t444;
t421 = qJDD(6) - t422;
t726 = t612 * t421;
t459 = t609 * t509 + t511 * t608;
t775 = qJD(6) + t459;
t780 = t775 * t616;
t762 = -t775 * t780 - t726;
t769 = t447 * t775;
t445 = t616 * t600 + t612 * t655;
t770 = t445 * t775;
t788 = -((t410 + t769) * t612 - (t409 - t770) * t616) * MDP(25) + (t447 * t780 + t737) * MDP(24) + (-t733 - t762) * MDP(26) + (-t509 * t600 + t443) * MDP(17) - t655 * MDP(28) * t775 - (t511 * t600 + t444) * MDP(18) + t511 * t509 * MDP(15) - (t509 ^ 2 - t511 ^ 2) * MDP(16) - t599 * MDP(19);
t582 = pkin(7) * t707;
t782 = -pkin(8) * t707 + qJD(3) + t582;
t583 = pkin(7) * t706;
t534 = -pkin(8) * t706 + t583;
t620 = -pkin(2) - pkin(3);
t678 = -qJ(3) * t613 + t617 * t620;
t779 = qJD(4) * t678 - t613 * t534 + t617 * t782;
t536 = qJ(3) * t617 + t613 * t620;
t778 = qJD(4) * t536 + t617 * t534 + t613 * t782;
t577 = pkin(7) * t694;
t681 = pkin(7) * t682 + qJDD(3) + t577;
t466 = -pkin(8) * t764 + t620 * qJDD(2) + t681;
t578 = pkin(7) * t693;
t602 = qJDD(2) * qJ(3);
t603 = qJD(2) * qJD(3);
t488 = -pkin(7) * t683 + t578 + t602 + t603;
t468 = pkin(8) * t763 + t488;
t686 = t620 * qJD(2);
t489 = t686 + t782;
t604 = qJD(2) * qJ(3);
t512 = t534 + t604;
t656 = t489 * t613 + t512 * t617;
t776 = -t656 * qJD(4) + t617 * t466 - t613 * t468;
t747 = pkin(4) * t511;
t774 = pkin(5) * t655 + pkin(9) * t459 + t747;
t676 = t617 * t489 - t512 * t613;
t739 = qJ(5) * t511;
t439 = t676 - t739;
t436 = -pkin(4) * t600 + t439;
t740 = qJ(5) * t509;
t440 = t656 - t740;
t437 = t609 * t440;
t414 = t608 * t436 + t437;
t412 = -pkin(9) * t600 + t414;
t513 = -qJD(1) * pkin(1) - pkin(2) * t706 - qJ(3) * t707;
t487 = pkin(3) * t706 - t513;
t455 = pkin(4) * t509 + qJD(5) + t487;
t417 = pkin(5) * t459 - pkin(9) * t655 + t455;
t400 = -t412 * t612 + t417 * t616;
t772 = t400 * t655;
t401 = t412 * t616 + t417 * t612;
t771 = t401 * t655;
t734 = t445 * t655;
t767 = -t739 + t779;
t766 = t740 - t778;
t615 = sin(qJ(1));
t652 = t613 * t618 - t614 * t617;
t496 = t652 * t615;
t619 = cos(qJ(1));
t720 = t618 * t619;
t691 = t613 * t720;
t723 = t614 * t619;
t498 = -t617 * t723 + t691;
t744 = g(3) * t525;
t765 = -g(1) * t498 - g(2) * t496 + t487 * t511 - t744 - t776;
t711 = t618 * pkin(2) + t614 * qJ(3);
t759 = -pkin(1) - t711;
t748 = pkin(7) - pkin(8);
t542 = t748 * t614;
t543 = t748 * t618;
t714 = t613 * t542 + t617 * t543;
t598 = g(1) * t619;
t756 = g(2) * t615 + t598;
t569 = pkin(4) * t608 + pkin(9);
t405 = -pkin(4) * t599 - qJ(5) * t443 - qJD(5) * t511 + t776;
t642 = -t613 * t466 - t617 * t468 - t489 * t702 + t512 * t703;
t408 = -qJ(5) * t444 - qJD(5) * t509 - t642;
t396 = t405 * t609 - t408 * t608;
t394 = pkin(5) * t599 - t396;
t601 = qJ(4) + pkin(10);
t586 = sin(t601);
t587 = cos(t601);
t653 = t614 * t586 + t618 * t587;
t721 = t615 * t618;
t724 = t614 * t615;
t635 = g(1) * (t586 * t720 - t587 * t723) + g(2) * (t586 * t721 - t587 * t724) + g(3) * t653 - t394;
t753 = t775 * (qJD(6) * t569 + t774) - t635;
t531 = -pkin(4) + t678;
t476 = t608 * t531 + t609 * t536;
t474 = -pkin(9) + t476;
t574 = qJ(3) * t706;
t687 = qJD(1) * t620;
t495 = t614 * t687 + t574;
t752 = t775 * (qJD(6) * t474 + t495 - t774) + t635;
t742 = pkin(7) * qJDD(2);
t750 = qJD(2) * (qJD(1) * t759 + t513) - t742;
t705 = qJD(2) * t614;
t477 = -t617 * t705 + t777;
t533 = t748 * t705;
t535 = qJD(2) * t543;
t641 = -t617 * t533 + t613 * t535 + t542 * t702 - t543 * t703;
t425 = -qJ(5) * t477 - qJD(5) * t525 + t641;
t478 = qJD(2) * t525 - t639;
t629 = -qJD(4) * t714 + t533 * t613 + t617 * t535;
t623 = -qJ(5) * t478 + qJD(5) * t652 + t629;
t403 = t609 * t425 + t608 * t623;
t735 = t440 * t608;
t413 = t436 * t609 - t735;
t411 = pkin(5) * t600 - t413;
t469 = t609 * t525 - t608 * t652;
t470 = -t525 * t608 - t609 * t652;
t521 = t618 * pkin(3) - t759;
t648 = pkin(4) * t525 + t521;
t429 = pkin(5) * t469 - pkin(9) * t470 + t648;
t454 = -qJ(5) * t525 + t714;
t673 = t617 * t542 - t543 * t613;
t645 = qJ(5) * t652 + t673;
t431 = t609 * t454 + t608 * t645;
t434 = -t477 * t608 + t478 * t609;
t397 = t608 * t405 + t609 * t408;
t669 = -pkin(9) * t599 + qJD(6) * t417 + t397;
t749 = t394 * t470 + t411 * t434 - t431 * t421 - (qJD(6) * t429 + t403) * t775 - t469 * t669 + t598;
t746 = pkin(4) * t613;
t597 = g(1) * t615;
t745 = g(2) * t619;
t607 = qJDD(1) * pkin(1);
t738 = qJDD(2) * pkin(2);
t736 = t411 * t470;
t732 = t775 * t612;
t576 = pkin(4) * t617 + pkin(3);
t728 = t576 * t618;
t622 = qJD(1) ^ 2;
t722 = t614 * t622;
t419 = t616 * t421;
t719 = t608 * t767 - t609 * t766;
t718 = t608 * t766 + t609 * t767;
t524 = t608 * t617 + t609 * t613;
t717 = t600 * t524;
t523 = t608 * t613 - t609 * t617;
t716 = t600 * t523;
t588 = t614 * qJD(3);
t712 = qJ(3) * t704 + t588;
t605 = t614 ^ 2;
t606 = t618 ^ 2;
t709 = t605 - t606;
t701 = qJD(6) * t412;
t700 = qJD(6) * t655;
t692 = pkin(4) * t725;
t690 = t618 * t722;
t688 = -g(1) * t723 - g(2) * t724 + g(3) * t618;
t680 = t597 - t745;
t679 = -qJD(2) * pkin(2) + qJD(3);
t665 = t600 ^ 2;
t664 = t619 * pkin(1) + pkin(2) * t720 + t615 * pkin(7) + qJ(3) * t723;
t663 = -t577 - t688;
t662 = t614 * t686;
t621 = qJD(2) ^ 2;
t661 = pkin(7) * t621 + t745;
t483 = t653 * t615;
t659 = g(1) * t483 + t429 * t421;
t658 = pkin(2) * t693 + qJ(3) * t764 + qJD(1) * t588 + t607;
t657 = -t701 - t745;
t475 = t531 * t609 - t536 * t608;
t539 = t582 + t679;
t541 = t583 + t604;
t654 = t539 * t618 - t541 * t614;
t651 = qJD(6) * t524 + t707;
t649 = -t459 * t732 - t699 * t775 + t419;
t494 = t681 - t738;
t647 = -0.2e1 * pkin(1) * t695 - t742;
t644 = t434 * t616 - t470 * t699;
t480 = t662 + t712;
t637 = -t661 + 0.2e1 * t607;
t634 = g(2) * t483 + g(3) * (-t586 * t618 + t587 * t614) - t669;
t633 = pkin(4) * t477 + t480;
t416 = t439 * t609 - t735;
t632 = -t569 * t421 + (t411 + t416) * t775;
t630 = -t474 * t421 + (-t411 - t718) * t775;
t465 = pkin(2) * t683 - t658;
t502 = pkin(2) * t705 - t712;
t628 = -qJD(1) * t502 - qJDD(1) * t759 - t465 - t661;
t448 = pkin(3) * t693 + qJD(1) * t662 + t658;
t627 = qJD(2) * t654 + t488 * t618 + t494 * t614;
t625 = pkin(4) * t444 + qJDD(5) + t448;
t497 = t525 * t615;
t499 = t525 * t619;
t624 = g(1) * t499 + g(2) * t497 - g(3) * t652 + t487 * t509 + t642;
t611 = -qJ(5) - pkin(8);
t593 = t619 * pkin(7);
t570 = -pkin(4) * t609 - pkin(5);
t565 = g(1) * t721;
t559 = qJ(3) * t720;
t557 = qJ(3) * t721;
t530 = pkin(2) * t707 - t574;
t485 = t653 * t619;
t473 = pkin(5) - t475;
t472 = t485 * t616 - t612 * t615;
t471 = -t485 * t612 - t615 * t616;
t433 = t609 * t477 + t478 * t608;
t430 = t454 * t608 - t609 * t645;
t415 = t439 * t608 + t437;
t406 = pkin(5) * t433 - pkin(9) * t434 + t633;
t402 = t425 * t608 - t609 * t623;
t399 = -pkin(5) * t422 - pkin(9) * t423 + t625;
t398 = t616 * t399;
t1 = [t680 * MDP(2) + (qJDD(1) * t605 + 0.2e1 * t614 * t682) * MDP(4) + (g(1) * t497 - g(2) * t499 + t521 * t444 + t448 * t525 + t487 * t477 + t480 * t509 - t599 * t673 - t600 * t629) * MDP(20) + (t409 * t470 * t616 + t447 * t644) * MDP(24) + (t614 * t647 + t618 * t637 + t565) * MDP(9) + (pkin(7) * t627 - g(1) * t593 - g(2) * t664 + t513 * t502 + (t465 - t597) * t759) * MDP(14) + (-g(1) * t496 + g(2) * t498 + t521 * t443 - t448 * t652 + t487 * t478 + t480 * t511 + t599 * t714 + t600 * t641) * MDP(21) + (-t478 * t600 + t599 * t652) * MDP(17) + (-t443 * t525 + t444 * t652 - t477 * t511 - t478 * t509) * MDP(16) + (-t443 * t652 + t478 * t511) * MDP(15) + (t397 * t431 + t414 * t403 - t396 * t430 - t413 * t402 + t625 * t648 + t455 * t633 - g(1) * (t611 * t619 + t593) - g(2) * (t576 * t720 + t619 * t692 + t664) + (-g(1) * (t759 - t692 - t728) - g(2) * t611) * t615) * MDP(23) + ((t605 + t606) * qJDD(1) * pkin(7) + t627 - t756) * MDP(12) + t756 * MDP(3) + 0.2e1 * (t614 * t693 - t695 * t709) * MDP(5) + (t614 * t750 + t628 * t618 + t565) * MDP(11) + (-t750 * t618 + (t628 + t597) * t614) * MDP(13) + (t647 * t618 + (-t637 - t597) * t614) * MDP(10) + (t409 * t469 + t419 * t470 + t433 * t447 + t644 * t775) * MDP(26) + (t421 * t469 + t433 * t775) * MDP(28) + (-g(2) * t471 - t401 * t433 + t402 * t447 + t430 * t409 + (-(-qJD(6) * t431 + t406) * t775 - (t399 - t701) * t469 - qJD(6) * t736 - t659) * t612 + t749 * t616) * MDP(30) + (-g(2) * t472 + t398 * t469 + t400 * t433 + t402 * t445 + t430 * t410 + (t406 * t775 + (-t412 * t469 - t431 * t775 + t736) * qJD(6) + t659) * t616 + t749 * t612) * MDP(29) + (-t470 * t726 - t410 * t469 - t433 * t445 + (-t434 * t612 - t470 * t698) * t775) * MDP(27) + (-t396 * t470 - t397 * t469 + t402 * t655 - t403 * t459 - t413 * t434 - t414 * t433 + t422 * t431 + t423 * t430 + t756) * MDP(22) + ((-t445 * t616 - t447 * t612) * t434 + (-t737 - t410 * t616 + (t445 * t612 - t447 * t616) * qJD(6)) * t470) * MDP(25) + (qJDD(2) * t614 + t618 * t621) * MDP(6) + (qJDD(2) * t618 - t614 * t621) * MDP(7) + (t477 * t600 + t525 * t599) * MDP(18) + qJDD(1) * MDP(1); -MDP(4) * t690 + (-t495 * t509 - t678 * t599 + t600 * t778 + t765) * MDP(20) + (-t495 * t511 + t536 * t599 + t600 * t779 - t624) * MDP(21) + (t578 + 0.2e1 * t602 + 0.2e1 * t603 + (qJD(1) * t530 - g(3)) * t614 + (qJD(1) * t513 - t756) * t618) * MDP(13) + (g(3) * t614 - t578 + (pkin(1) * t622 + t756) * t618) * MDP(10) + (t422 * t476 - t423 * t475 + (-t414 + t719) * t655 + (t413 - t718) * t459) * MDP(22) + ((-pkin(2) * t614 + qJ(3) * t618) * qJDD(1) + ((t541 - t604) * t614 + (-t539 + t679) * t618) * qJD(1)) * MDP(12) + (t473 * t409 + t719 * t447 + t612 * t752 + t630 * t616 - t771) * MDP(30) + (t473 * t410 + t719 * t445 + t630 * t612 - t616 * t752 + t772) * MDP(29) + (t397 * t476 + t396 * t475 - t455 * (t574 - t747) - g(1) * (pkin(4) * t691 + t559) - g(2) * (t721 * t746 + t557) - g(3) * (t711 + t728) + t718 * t414 - t719 * t413 + (-g(3) * t746 - t455 * t687 + t756 * (pkin(2) + t576)) * t614) * MDP(23) + MDP(7) * t693 + (t732 * t775 - t419 - t734) * MDP(27) + MDP(6) * t694 + (0.2e1 * t738 - qJDD(3) + (-t513 * t614 + t530 * t618) * qJD(1) + t663) * MDP(11) + t709 * MDP(5) * t622 + (pkin(1) * t722 + t663) * MDP(9) + (t488 * qJ(3) + t541 * qJD(3) - t494 * pkin(2) - t513 * t530 - g(1) * (-pkin(2) * t723 + t559) - g(2) * (-pkin(2) * t724 + t557) - g(3) * t711 - t654 * qJD(1) * pkin(7)) * MDP(14) + qJDD(2) * MDP(8) - t788; (-qJDD(2) - t690) * MDP(11) + MDP(12) * t694 + (-t605 * t622 - t621) * MDP(13) + (-qJD(2) * t541 + t513 * t707 + t494 + t688) * MDP(14) + (-t509 * t707 - t599 * t617 - t613 * t665) * MDP(20) + (-t511 * t707 + t613 * t599 - t617 * t665) * MDP(21) + (t422 * t524 + t423 * t523 - t459 * t716 - t655 * t717) * MDP(22) + (-t396 * t523 + t397 * t524 + t413 * t717 + t414 * t716 - t455 * t707 + t688) * MDP(23) + (-t524 * t726 + t523 * t410 - t717 * t445 + (-t612 * t716 - t616 * t651) * t775) * MDP(29) + (-t524 * t419 + t523 * t409 - t717 * t447 + (t612 * t651 - t616 * t716) * t775) * MDP(30); (-t600 * t656 - t765) * MDP(20) + (-t600 * t676 + t624) * MDP(21) + ((t422 * t608 - t423 * t609) * pkin(4) - (t413 - t416) * t459 + (t414 - t415) * t655) * MDP(22) + (t413 * t415 - t414 * t416 + (t396 * t609 + t397 * t608 - t455 * t511 + t652 * t756 + t744) * pkin(4)) * MDP(23) + (t649 + t734) * MDP(27) + (t570 * t410 - t415 * t445 + t632 * t612 - t616 * t753 - t772) * MDP(29) + (t570 * t409 - t415 * t447 + t612 * t753 + t632 * t616 + t771) * MDP(30) + t788; (-t459 ^ 2 - t655 ^ 2) * MDP(22) + (t413 * t655 + t414 * t459 + t625 + t680) * MDP(23) + (t649 - t734) * MDP(29) + (-t733 + t762) * MDP(30); t447 * t445 * MDP(24) + (-t445 ^ 2 + t447 ^ 2) * MDP(25) + (t689 + t770) * MDP(26) + (-t561 + t769) * MDP(27) + t421 * MDP(28) + (-g(1) * t471 + t401 * t775 - t411 * t447 + t398) * MDP(29) + (g(1) * t472 + t400 * t775 + t411 * t445) * MDP(30) + (-MDP(27) * t700 + MDP(29) * t657 + MDP(30) * t634) * t616 + (-MDP(26) * t700 + (qJD(6) * t600 - t423) * MDP(27) + t634 * MDP(29) + (-t399 - t657) * MDP(30)) * t612;];
tau  = t1;

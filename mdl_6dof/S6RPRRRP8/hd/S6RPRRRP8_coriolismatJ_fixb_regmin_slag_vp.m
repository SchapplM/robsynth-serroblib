% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x31]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRP8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:37
% EndTime: 2019-03-09 06:25:00
% DurationCPUTime: 13.68s
% Computational Cost: add. (12106->670), mult. (21467->796), div. (0->0), fcn. (22619->6), ass. (0->507)
t472 = cos(qJ(5));
t473 = cos(qJ(3));
t764 = cos(qJ(4));
t456 = t764 * t473;
t470 = sin(qJ(4));
t471 = sin(qJ(3));
t689 = t470 * t471;
t426 = -t456 + t689;
t425 = t426 ^ 2;
t427 = t470 * t473 + t471 * t764;
t787 = t427 ^ 2;
t548 = t787 / 0.2e1 + t425 / 0.2e1;
t828 = t472 * t548;
t469 = sin(qJ(5));
t827 = t548 * t469;
t762 = pkin(5) * t426;
t474 = -pkin(1) - pkin(7);
t805 = -pkin(8) + t474;
t430 = t805 * t471;
t431 = t805 * t473;
t564 = t470 * t430 - t764 * t431;
t818 = t564 * t469;
t826 = t762 / 0.2e1 - t818 / 0.2e1;
t262 = t426 * t469;
t401 = pkin(5) * t262;
t688 = t472 * qJ(6);
t544 = t426 * t688 - t401;
t170 = t544 + t564;
t726 = t170 * t469;
t158 = t726 / 0.2e1;
t435 = pkin(5) * t469 - t688;
t701 = t435 * t469;
t582 = -t701 / 0.2e1;
t785 = pkin(5) / 0.2e1;
t825 = t158 + (t582 + t785) * t426 + t826;
t623 = qJD(3) + qJD(4);
t817 = t564 * t472;
t824 = -t817 / 0.2e1;
t268 = t818 / 0.2e1;
t419 = t426 * qJ(6);
t758 = t427 * pkin(9);
t761 = t426 * pkin(4);
t314 = t758 - t761;
t763 = pkin(3) * t473;
t280 = t314 + t763;
t237 = t469 * t280;
t680 = t237 / 0.2e1 + t824;
t822 = t680 - t419;
t693 = t469 * qJ(6);
t756 = t472 * pkin(5);
t537 = -t693 - t756;
t434 = -pkin(4) + t537;
t621 = t764 * pkin(3);
t418 = -t621 + t434;
t387 = t418 * t469;
t424 = t434 * t469;
t671 = t387 / 0.2e1 + t424 / 0.2e1;
t700 = t435 * t472;
t821 = t700 - t671;
t820 = t762 - t818;
t648 = qJD(5) * t472;
t595 = t427 * t648;
t504 = (0.1e1 / 0.2e1 + t548) * t472;
t793 = qJD(1) * t504;
t815 = t262 * t623 - t595 - t793;
t797 = t623 * t427;
t538 = t426 * t797;
t687 = t472 * t280;
t135 = -t687 + t820;
t735 = t135 * t469;
t672 = -t817 + t237;
t133 = -t419 + t672;
t736 = t133 * t472;
t511 = t736 / 0.2e1 + t735 / 0.2e1;
t257 = t469 * t427;
t260 = t472 * t427;
t541 = t430 * t764 + t470 * t431;
t487 = -pkin(5) * t257 + qJ(6) * t260 + t541;
t725 = t487 * t434;
t814 = -t725 / 0.2e1 - t511 * pkin(9);
t454 = pkin(3) * t471 + qJ(2);
t760 = t426 * pkin(9);
t545 = t427 * pkin(4) + t760;
t497 = t454 + t545;
t799 = t541 * t469;
t154 = -t472 * t497 + t799;
t759 = t427 * pkin(5);
t128 = t154 - t759;
t801 = t487 * t469;
t813 = (t128 - t801) * t426;
t812 = (t154 - t799) * t426;
t798 = t541 * t472;
t155 = t469 * t497 + t798;
t811 = (t155 - t798) * t426;
t465 = -t469 / 0.2e1;
t708 = t427 * qJ(6);
t127 = t155 + t708;
t738 = t127 * t472;
t779 = t487 / 0.2e1;
t810 = (t779 - t738 / 0.2e1 + t128 * t465) * t426;
t264 = t426 * t472;
t796 = t623 * t469;
t56 = t264 * qJD(5) + t427 * t796;
t809 = t623 * t564;
t800 = t487 * t472;
t808 = t170 * t260 + (-t127 + t800) * t426;
t807 = t623 * t541;
t804 = t170 * t487;
t757 = t470 * pkin(3);
t457 = pkin(9) + t757;
t784 = -pkin(9) / 0.2e1;
t614 = t457 / 0.2e1 + t784;
t803 = t426 * t614;
t773 = t418 / 0.2e1;
t521 = t621 / 0.2e1 + t773;
t771 = -t434 / 0.2e1;
t802 = t427 * (t771 + t521);
t765 = t472 / 0.2e1;
t549 = t764 * t765;
t609 = t469 * t764;
t551 = -t609 / 0.2e1;
t483 = (pkin(5) * t551 + qJ(6) * t549) * pkin(3);
t770 = t434 / 0.2e1;
t572 = t770 + t773;
t168 = -t435 * t572 + t483;
t255 = t537 * t426;
t575 = -t154 / 0.2e1 + t128 / 0.2e1;
t576 = -t127 / 0.2e1 + t155 / 0.2e1;
t485 = t469 * t576 + t472 * t575;
t769 = t435 / 0.2e1;
t590 = t170 * t769;
t477 = pkin(9) * t485 + t255 * t770 + t590;
t285 = t469 * t314;
t673 = t285 - t817;
t140 = -t419 + t673;
t686 = t472 * t314;
t141 = -t686 + t820;
t781 = t141 / 0.2e1;
t516 = -t140 * qJ(6) / 0.2e1 + pkin(5) * t781;
t6 = t477 + t516;
t579 = -t688 / 0.2e1;
t774 = t401 / 0.2e1;
t189 = t774 + (t579 - t435 / 0.2e1) * t426;
t644 = t189 * qJD(2);
t702 = t435 * t434;
t795 = -t6 * qJD(1) + t168 * qJD(3) - qJD(4) * t702 + t644;
t476 = t255 * t773 + t457 * t485 + t590;
t782 = -t133 / 0.2e1;
t517 = qJ(6) * t782 + t135 * t785;
t4 = t476 + t517;
t712 = t418 * t435;
t794 = -t4 * qJD(1) - qJD(3) * t712 + t644;
t226 = t787 - t425;
t659 = qJD(2) * t504;
t634 = t427 * qJD(1);
t791 = -t634 - qJD(5);
t467 = t469 ^ 2;
t468 = t472 ^ 2;
t767 = t468 / 0.2e1;
t256 = (-t467 / 0.2e1 + t767) * t426;
t690 = t469 * t472;
t600 = qJD(1) * t690;
t362 = t425 * t600;
t790 = t256 * t623 + t362;
t649 = qJD(5) * t469;
t596 = t427 * t649;
t789 = -t264 * t623 - t596;
t447 = t468 - t467;
t554 = t426 * t600;
t228 = t447 * t623 + 0.2e1 * t554;
t601 = t472 * t634;
t788 = -qJD(5) * t504 - t601;
t786 = pkin(4) / 0.2e1;
t732 = t141 * t469;
t733 = t140 * t472;
t510 = t733 / 0.2e1 + t732 / 0.2e1;
t780 = t170 / 0.2e1;
t10 = (t780 + t510) * t427 + t810;
t8 = (t780 + t511) * t427 + t810;
t783 = t8 * qJD(3) + t10 * qJD(4);
t778 = -t255 / 0.2e1;
t777 = -t280 / 0.2e1;
t776 = t817 / 0.2e1;
t775 = -t314 / 0.2e1;
t772 = -t426 / 0.2e1;
t458 = -t621 - pkin(4);
t768 = -t458 / 0.2e1;
t464 = t469 / 0.2e1;
t766 = -t472 / 0.2e1;
t753 = pkin(3) * qJD(4);
t752 = t8 * qJD(1);
t421 = t456 / 0.2e1 - t689 / 0.2e1;
t393 = t468 * t772;
t694 = t467 * t426;
t669 = t694 / 0.2e1 + t393;
t210 = -t421 + t669;
t751 = t210 * qJD(6);
t211 = t421 + t669;
t750 = t211 * qJD(6);
t515 = t693 / 0.2e1 + t756 / 0.2e1;
t737 = t128 * t472;
t739 = t127 * t469;
t24 = t154 * t766 - t739 / 0.2e1 + t155 * t464 + t737 / 0.2e1 + t515 * t427;
t749 = t24 * qJD(5);
t23 = (t759 / 0.2e1 - t575) * t472 + (t708 / 0.2e1 - t576) * t469;
t748 = -t23 * qJD(5) + t260 * qJD(6);
t729 = t155 * t427;
t47 = -t729 + (-t170 * t472 - t255 * t469) * t426;
t747 = qJD(1) * t47;
t730 = t154 * t427;
t48 = -t730 + (t255 * t472 - t726) * t426;
t746 = qJD(1) * t48;
t740 = t127 * t427;
t61 = t170 * t264 + t740;
t745 = qJD(1) * t61;
t62 = -t737 + t739;
t744 = qJD(1) * t62;
t75 = t262 * t564 + t730;
t743 = qJD(1) * t75;
t76 = -t264 * t564 - t729;
t742 = qJD(1) * t76;
t741 = t10 * qJD(1);
t480 = t485 * t427 + t255 * t426 / 0.2e1;
t14 = t480 - t515;
t734 = t14 * qJD(1);
t15 = -t127 * t154 + t128 * t155 + t170 * t255;
t731 = t15 * qJD(1);
t99 = t128 * t260;
t16 = t135 * t264 + t99 + (-t133 * t426 - t740) * t469;
t728 = t16 * qJD(1);
t17 = t141 * t264 + t99 + (-t140 * t426 - t740) * t469;
t727 = t17 * qJD(1);
t19 = t155 * t264 + (-t738 + (-t128 + t154) * t469) * t426;
t720 = t19 * qJD(1);
t20 = t23 * qJD(1);
t25 = t133 * t427 + t808;
t719 = t25 * qJD(1);
t26 = (-t135 - t726) * t427 + t813;
t718 = t26 * qJD(1);
t27 = t140 * t427 + t808;
t717 = t27 * qJD(1);
t28 = (-t141 - t726) * t427 + t813;
t716 = t28 * qJD(1);
t711 = t418 * t472;
t710 = t426 * t457;
t709 = t426 * t470;
t707 = t427 * t418;
t706 = t427 * t434;
t705 = t427 * t458;
t43 = t687 * t427 + t812;
t704 = t43 * qJD(1);
t703 = t434 * t472;
t44 = (-t672 - t817) * t427 + t811;
t699 = t44 * qJD(1);
t45 = t314 * t260 + t812;
t698 = t45 * qJD(1);
t697 = t457 * t427;
t696 = t458 * t426;
t46 = (-t673 - t817) * t427 + t811;
t695 = t46 * qJD(1);
t678 = t687 / 0.2e1 + t268;
t676 = t285 / 0.2e1 + t824;
t675 = -t285 / 0.2e1 + t776;
t267 = t686 / 0.2e1;
t674 = t267 + t268;
t652 = qJD(4) * t472;
t655 = qJD(3) * t472;
t668 = (t652 + t655) * t469;
t617 = t470 * t753;
t445 = t469 * t617;
t462 = t467 * qJD(6);
t667 = t462 - t445;
t666 = -t467 - t468;
t665 = qJD(1) * qJ(2);
t213 = t226 * t469;
t663 = qJD(1) * t213;
t214 = t226 * t472;
t662 = qJD(1) * t214;
t235 = -t426 * t454 + t427 * t763;
t661 = qJD(1) * t235;
t236 = -t426 * t763 - t427 * t454;
t660 = qJD(1) * t236;
t206 = t766 + t828;
t658 = qJD(2) * t206;
t657 = qJD(2) * t427;
t656 = qJD(3) * t469;
t654 = qJD(4) * t454;
t653 = qJD(4) * t469;
t651 = qJD(5) * t154;
t650 = qJD(5) * t427;
t647 = qJD(6) * t469;
t134 = (0.1e1 + t666) * t427 * t426;
t646 = t134 * qJD(2);
t173 = t465 - t827;
t645 = t173 * qJD(1);
t203 = t464 + t827;
t643 = t203 * qJD(1);
t641 = t226 * qJD(1);
t640 = t257 * qJD(1);
t242 = t260 * qJD(1);
t283 = t447 * t425;
t638 = t283 * qJD(1);
t284 = t666 * t426;
t637 = t284 * qJD(1);
t636 = t421 * qJD(1);
t635 = t426 * qJD(1);
t415 = t427 * qJD(6);
t446 = t471 ^ 2 - t473 ^ 2;
t633 = t446 * qJD(1);
t632 = t471 * qJD(1);
t631 = t471 * qJD(3);
t630 = t472 * qJD(6);
t629 = t473 * qJD(1);
t628 = t473 * qJD(3);
t620 = qJD(3) * t757;
t619 = pkin(9) * t649;
t618 = pkin(9) * t648;
t616 = -t762 / 0.2e1;
t615 = t786 + t768;
t587 = t426 * t769;
t190 = t426 * t579 + t587 + t774;
t613 = t190 * qJD(5) - t262 * qJD(6) + t646;
t159 = -t726 / 0.2e1;
t577 = t260 / 0.2e1;
t584 = t264 / 0.2e1;
t612 = t418 * t584 + t457 * t577 + t159;
t611 = pkin(9) * t577 + t434 * t584 + t159;
t610 = -t189 * qJD(5) - t646;
t608 = t764 * t467;
t607 = t764 * t468;
t606 = qJ(2) * t632;
t605 = qJ(2) * t629;
t604 = t454 * t635;
t603 = t454 * t634;
t602 = t469 * t634;
t599 = t469 * t657;
t598 = t457 * t649;
t597 = t457 * t648;
t324 = t426 * t634;
t453 = t469 * t648;
t594 = t426 * t647;
t450 = t469 * t630;
t593 = t471 * t629;
t586 = t710 / 0.2e1;
t585 = -t264 / 0.2e1;
t583 = -t703 / 0.2e1;
t581 = -t696 / 0.2e1;
t578 = -t260 / 0.2e1;
t574 = -t487 / 0.2e1 + t779;
t571 = t767 + t467 / 0.2e1;
t570 = t764 * qJD(3);
t569 = t764 * qJD(4);
t501 = (-t757 / 0.2e1 + t614) * t426;
t34 = t574 * t472 + (t501 - t802) * t469;
t444 = t472 * t620;
t568 = -qJD(1) * t34 + t444;
t563 = t264 * t757;
t357 = t563 / 0.2e1;
t36 = t357 + t574 * t469 + (t802 - t803) * t472;
t443 = t469 * t620;
t567 = -qJD(1) * t36 + t443;
t356 = -t563 / 0.2e1;
t555 = -t621 / 0.2e1;
t522 = t555 + t768;
t494 = (-pkin(4) / 0.2e1 + t522) * t427;
t60 = t356 + (t494 + t803) * t472;
t566 = -qJD(1) * t60 - t443;
t57 = (t494 + t501) * t469;
t565 = -qJD(1) * t57 + t444;
t561 = t623 * t472;
t558 = t472 * t617;
t557 = pkin(9) * t262 / 0.2e1;
t550 = t609 / 0.2e1;
t547 = t778 - t758 / 0.2e1;
t543 = t571 * t457;
t542 = t778 - t697 / 0.2e1;
t166 = qJD(1) * t211 + t668;
t200 = qJD(1) * t256 - t668;
t539 = t469 * t561;
t536 = -t706 + t760;
t11 = t127 * t133 + t128 * t135 + t804;
t535 = t11 * qJD(1) + t8 * qJD(2);
t12 = t127 * t140 + t128 * t141 + t804;
t534 = t12 * qJD(1) + t10 * qJD(2);
t533 = t735 + t736;
t532 = t732 + t733;
t531 = t707 - t710;
t530 = -t705 + t710;
t293 = t701 + t711;
t513 = (t587 - t170 / 0.2e1) * t472;
t478 = t513 + (t418 * t772 + t542) * t469;
t40 = t478 - t822;
t529 = -qJD(1) * t40 + qJD(3) * t293;
t294 = -t387 + t700;
t291 = t418 * t585;
t31 = t291 + (t777 + t542) * t472 + t825;
t528 = -qJD(1) * t31 + qJD(3) * t294;
t30 = (t140 / 0.2e1 + t782) * t472 + (t781 - t135 / 0.2e1) * t469;
t499 = t608 + t607;
t420 = t499 * pkin(3);
t527 = -qJD(1) * t30 - qJD(3) * t420;
t526 = t426 * t791;
t300 = t425 * t468 + t787;
t525 = -qJD(1) * t300 - t650;
t524 = qJD(5) * t435 - t647;
t523 = t450 - t558;
t519 = t426 * t539;
t518 = t796 * t264;
t50 = t616 + t612 + t678;
t509 = -qJD(1) * t50 + t418 * t656;
t340 = t457 * t578;
t67 = t340 + (t581 + t777) * t472;
t508 = -qJD(1) * t67 - t458 * t656;
t482 = (t697 / 0.2e1 + t696 / 0.2e1) * t469 + t776;
t65 = t482 + t680;
t507 = -qJD(1) * t65 - t458 * t655;
t506 = t472 * t526;
t505 = t262 * qJD(5) - t427 * t561;
t219 = -qJD(5) * t421 + t324;
t502 = t426 * t582 + t158 + 0.2e1 * t616;
t500 = -0.2e1 * t519;
t475 = t510 * t457 + (t127 * t549 + t128 * t550 + t470 * t780) * pkin(3) + t487 * t773;
t2 = t475 + t814;
t215 = (t418 * t470 + t457 * t499) * pkin(3);
t493 = t607 / 0.2e1 + t608 / 0.2e1;
t54 = (pkin(3) * t493 + t771 + t773) * t427 + (t757 / 0.2e1 - t543 + t571 * pkin(9)) * t426;
t498 = t2 * qJD(1) + t54 * qJD(2) + t215 * qJD(3);
t440 = pkin(3) * t551;
t194 = t440 + t821;
t326 = -t424 + t700;
t322 = t426 * t583;
t37 = t322 + (t775 + t547) * t472 + t825;
t496 = -qJD(1) * t37 + qJD(3) * t194 + qJD(4) * t326;
t441 = pkin(3) * t549;
t195 = t472 * t572 + t441 + t701;
t325 = t701 + t703;
t479 = t513 + (t426 * t771 + t547) * t469;
t42 = t419 + t479 + t675;
t495 = -qJD(1) * t42 + qJD(3) * t195 + qJD(4) * t325;
t347 = t469 * t615 + t440;
t378 = pkin(9) * t578;
t71 = t378 + (t761 / 0.2e1 + t775) * t472;
t491 = pkin(4) * t653 - qJD(1) * t71 + qJD(3) * t347;
t442 = t472 * t555;
t348 = t472 * t615 + t442;
t484 = (t758 / 0.2e1 - t761 / 0.2e1) * t469 + t776;
t69 = t484 + t676;
t490 = pkin(4) * t652 - qJD(1) * t69 + qJD(3) * t348;
t439 = pkin(3) * t550;
t224 = t439 + t671;
t52 = t616 + t611 + t674;
t488 = -qJD(1) * t52 + qJD(3) * t224 + t434 * t653;
t486 = qJD(5) * t537 + t630;
t481 = (-t764 * t427 / 0.2e1 - t709 / 0.2e1) * pkin(3) + t586;
t437 = t447 * qJD(5);
t402 = t420 * qJD(4);
t363 = t426 * t450;
t361 = t791 * qJ(6);
t350 = pkin(4) * t766 + t458 * t765 + t442;
t349 = pkin(4) * t465 + t458 * t464 + t440;
t298 = t623 * t426;
t295 = t467 * t623 - t554;
t281 = t623 * t421;
t238 = t256 * qJD(5);
t225 = t439 - t671;
t223 = t242 + t648;
t222 = -t640 - t649;
t204 = t464 - t827;
t197 = t440 - t821;
t196 = -t701 - t711 / 0.2e1 + t583 + t441;
t193 = -0.2e1 * t469 * t506;
t176 = t765 - t828;
t175 = t465 + t827;
t169 = t702 / 0.2e1 + t712 / 0.2e1 + t483;
t165 = -t324 * t468 - t238;
t156 = t623 * t284;
t153 = t155 * qJD(5);
t131 = qJD(5) * t260 + t662;
t130 = -qJD(5) * t257 - t663;
t92 = -t260 * t623 - t324 * t469;
t91 = -t238 + (t468 * t635 - t539) * t427;
t88 = -t426 * t796 - t662;
t87 = -t426 * t561 + t663;
t78 = t469 * t526;
t77 = 0.2e1 * (qJD(5) - t634) * t426 * t690 - t447 * t797;
t72 = pkin(4) * t584 + t267 + 0.2e1 * t268 + t378;
t70 = t484 + t675;
t68 = t472 * t581 + t268 + t340 + t678;
t66 = t482 - t680;
t59 = t356 + t799 + pkin(9) * t584 + pkin(4) * t577 + (t427 * t522 + t586) * t472;
t58 = -t798 + t557 + t257 * t786 + (-t705 / 0.2e1 + t481) * t469;
t53 = t707 / 0.2e1 + pkin(9) * t393 + t706 / 0.2e1 + t694 * t784 - t426 * t543 + (t709 / 0.2e1 + t493 * t427) * pkin(3);
t51 = -t686 / 0.2e1 + t611 + t826;
t49 = -t687 / 0.2e1 + t612 + t826;
t41 = -t419 + t479 + t676;
t39 = t478 + t822;
t38 = t472 * t547 + t322 + t502 + t674;
t35 = t357 - t801 + pkin(9) * t585 + t434 * t577 + (-t710 / 0.2e1 + t521 * t427) * t472;
t33 = -t800 + t257 * t771 + t557 + (-t707 / 0.2e1 + t481) * t469;
t32 = t472 * t542 + t291 + t502 + t678;
t29 = t510 + t511;
t18 = t486 - t20;
t13 = t480 + t515;
t5 = t477 - t516;
t3 = t476 - t517;
t1 = t475 - t814;
t7 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t471 * t628, t446 * qJD(3), 0, 0, 0, qJ(2) * t628 + qJD(2) * t471, -qJ(2) * t631 + qJD(2) * t473, t538, t623 * t226, 0, 0, 0, qJD(3) * t235 - t426 * t654 + t657, -qJD(2) * t426 + qJD(3) * t236 - t427 * t654, -t425 * t453 + t468 * t538, -t283 * qJD(5) + t427 * t500, -t214 * t623 + t426 * t596, t213 * t623 + t426 * t595, -t538, qJD(3) * t43 + qJD(4) * t45 + qJD(5) * t76 + t472 * t657, qJD(3) * t44 + qJD(4) * t46 + qJD(5) * t75 - t599, qJD(3) * t26 + qJD(4) * t28 + qJD(5) * t47 + (-t425 * t647 + t657) * t472, -qJD(2) * t284 - qJD(3) * t16 - qJD(4) * t17 - qJD(5) * t19 + t427 * t594, qJD(3) * t25 + qJD(4) * t27 + qJD(5) * t48 + qJD(6) * t300 + t599, qJD(2) * t62 + qJD(3) * t11 + qJD(4) * t12 + qJD(5) * t15 + qJD(6) * t61; 0, 0, 0, 0, qJD(1), t665, 0, 0, 0, 0, 0, t632, t629, 0, 0, 0, 0, 0, t634, -t635, 0, 0, 0, 0, 0, -qJD(5) * t206 + t601, qJD(5) * t175 - t602, qJD(5) * t176 + t601, -t637, qJD(5) * t204 + t602, qJD(5) * t13 + qJD(6) * t206 + t744 + t783; 0, 0, 0, 0, 0, 0, -t593, t633, -t631, -t628, 0, -t474 * t631 + t605, -t474 * t628 - t606, t324, t641, -t797, t298, 0, t661 - t807, t660 + t809, t91, t77, t88, t87, -t219, t704 + (t469 * t530 - t798) * qJD(3) + t58 * qJD(4) + t68 * qJD(5), t699 + (t472 * t530 + t799) * qJD(3) + t59 * qJD(4) + t66 * qJD(5), t718 + (-t469 * t531 - t800) * qJD(3) + t33 * qJD(4) + t32 * qJD(5) + t751, qJD(3) * t533 + t29 * qJD(4) - t728 + t749, t719 + (t472 * t531 - t801) * qJD(3) + t35 * qJD(4) + t39 * qJD(5) - t363 (t418 * t487 + t457 * t533) * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + t49 * qJD(6) + t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t641, -t797, t298, 0, -t604 - t807, -t603 + t809, t91, t77, t88, t87, -t219, t698 + t58 * qJD(3) + (t469 * t545 - t798) * qJD(4) + t72 * qJD(5), t695 + t59 * qJD(3) + (t472 * t545 + t799) * qJD(4) + t70 * qJD(5), t716 + t33 * qJD(3) + (t469 * t536 - t800) * qJD(4) + t38 * qJD(5) + t751, t29 * qJD(3) + qJD(4) * t532 - t727 + t749, t717 + t35 * qJD(3) + (-t472 * t536 - t801) * qJD(4) + t41 * qJD(5) - t363, t1 * qJD(3) + (pkin(9) * t532 + t725) * qJD(4) + t5 * qJD(5) + t51 * qJD(6) + t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t790, 0.2e1 * t518 - t638, -t78, -t506, t281, qJD(3) * t68 + qJD(4) * t72 - t153 - t658 + t742, qJD(2) * t175 + qJD(3) * t66 + qJD(4) * t70 + t651 + t743, qJD(2) * t176 + qJD(3) * t32 + qJD(4) * t38 - t153 + t747, qJD(5) * t544 + t24 * t623 + t594 - t720, qJD(2) * t204 + qJD(3) * t39 + qJD(4) * t41 + t415 - t651 + t746, t731 + t13 * qJD(2) + t3 * qJD(3) + t5 * qJD(4) + (-pkin(5) * t155 - qJ(6) * t154) * qJD(5) + t127 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210 * t623 - t362, -t78, -t519 - t525, qJD(3) * t49 + qJD(4) * t51 + qJD(5) * t127 + t658 + t745; 0, 0, 0, 0, -qJD(1), -t665, 0, 0, 0, 0, 0, -t632, -t629, 0, 0, 0, 0, 0, -t634, t635, 0, 0, 0, 0, 0, t788, -qJD(5) * t173 + t602, t788, t637, -qJD(5) * t203 - t602, qJD(5) * t14 + qJD(6) * t504 - t744 + t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t623 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t631, -t628, 0, 0, 0, 0, 0, -t797, t298, 0, 0, 0, 0, 0, t505, t56, t505, t156, -t56, t752 + (t284 * t457 + t707) * qJD(3) + t53 * qJD(4) + t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t797, t298, 0, 0, 0, 0, 0, t505, t56, t505, t156, -t56, t741 + t53 * qJD(3) + (pkin(9) * t284 + t706) * qJD(4) + t613; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t815, -t645 - t789, t815, 0, -t643 + t789, t190 * t623 + t427 * t486 + t734; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t815; 0, 0, 0, 0, 0, 0, t593, -t633, 0, 0, 0, -t605, t606, -t324, -t641, 0, 0, 0, -t661, -t660, t165, t193, t131, t130, t219, qJD(4) * t57 + qJD(5) * t67 - t704, qJD(4) * t60 + qJD(5) * t65 - t699, qJD(4) * t34 + qJD(5) * t31 - t718 + t750, qJD(4) * t30 + t728 + t748, qJD(4) * t36 + qJD(5) * t40 - t363 - t719, qJD(4) * t2 + qJD(5) * t4 + qJD(6) * t50 - t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t54 + t610 - t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t617, -pkin(3) * t569, t453, t437, 0, 0, 0, t458 * t649 - t558, t458 * t648 + t445, -qJD(5) * t294 + t523, t402, -qJD(5) * t293 + t667, t215 * qJD(4) + t418 * t524; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t623 * t757 (-t570 - t569) * pkin(3), t453, t437, 0, 0, 0, qJD(5) * t349 - t558 - t565, qJD(5) * t350 + t445 - t566, qJD(5) * t197 + t523 - t568, t402 - t527, qJD(5) * t196 - t567 + t667 (pkin(9) * t499 + t434 * t470) * t753 + t169 * qJD(5) + t225 * qJD(6) + t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, t228, t223, t222, -t636, qJD(4) * t349 - t508 - t597, qJD(4) * t350 - t507 + t598, qJD(4) * t197 - t528 - t597, t18, qJD(4) * t196 - t529 - t598, t169 * qJD(4) + t457 * t486 - t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t223, t295, qJD(4) * t225 - t509 + t597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, -t641, 0, 0, 0, t604, t603, t165, t193, t131, t130, t219, -qJD(3) * t57 + qJD(5) * t71 - t698, -qJD(3) * t60 + qJD(5) * t69 - t695, -qJD(3) * t34 + qJD(5) * t37 - t716 + t750, -qJD(3) * t30 + t727 + t748, -qJD(3) * t36 + qJD(5) * t42 - t363 - t717, -qJD(3) * t2 + qJD(5) * t6 + qJD(6) * t52 - t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t54 + t610 - t741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t620, pkin(3) * t570, t453, t437, 0, 0, 0, -qJD(5) * t347 + t565, -qJD(5) * t348 + t566, -qJD(5) * t194 + t450 + t568, t527, -qJD(5) * t195 + t462 + t567, -qJD(5) * t168 - qJD(6) * t224 - t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t453, t437, 0, 0, 0, -pkin(4) * t649, -pkin(4) * t648, -qJD(5) * t326 + t450, 0, -qJD(5) * t325 + t462, t524 * t434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, t228, t223, t222, -t636, -t491 - t618, -t490 + t619, -t496 - t618, t18, -t495 - t619, pkin(9) * t486 - t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t223, t295, -t488 + t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t790, t500 + t638, t92, t257 * t623 - t324 * t472, t281, -qJD(3) * t67 - qJD(4) * t71 + t659 - t742, qJD(2) * t173 - qJD(3) * t65 - qJD(4) * t69 - t743, -qJD(3) * t31 - qJD(4) * t37 + t659 - t747, t23 * t623 + t720, qJD(2) * t203 - qJD(3) * t40 - qJD(4) * t42 + t415 - t746, qJ(6) * t415 - qJD(2) * t14 - qJD(3) * t4 - qJD(4) * t6 - t731; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t793, t645, t793, 0, t643, t189 * t623 - t734; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, -t228, -t242, t640, t636, qJD(4) * t347 + t508, qJD(4) * t348 + t507, qJD(4) * t194 + t528, t20, qJD(4) * t195 + t529, qJD(4) * t168 + t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, -t228, -t242, t640, t636, t491, t490, t496, t20, t495, t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t791, -t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211 * t623 + t362, t92, t518 + t525, -qJ(6) * t650 - qJD(3) * t50 - qJD(4) * t52 - t659 - t745; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t793; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t242, -t295, qJD(4) * t224 + t509; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t242, -t295, t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t791, t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;

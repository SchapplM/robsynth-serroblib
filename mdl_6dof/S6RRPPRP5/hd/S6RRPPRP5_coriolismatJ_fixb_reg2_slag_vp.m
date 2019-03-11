% Calculate inertial parameters regressor of coriolis matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPPRP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:44:05
% EndTime: 2019-03-09 08:44:26
% DurationCPUTime: 15.08s
% Computational Cost: add. (12053->661), mult. (22718->737), div. (0->0), fcn. (23351->6), ass. (0->484)
t454 = sin(qJ(2));
t453 = cos(pkin(9));
t727 = sin(qJ(5));
t572 = t727 * t453;
t452 = sin(pkin(9));
t728 = cos(qJ(5));
t576 = t728 * t452;
t810 = t572 + t576;
t814 = t810 * t454;
t828 = -t814 / 0.2e1;
t720 = -pkin(2) - qJ(4);
t585 = pkin(8) - t720;
t391 = t585 * t452;
t540 = t585 * t453;
t309 = -t391 * t728 - t540 * t727;
t575 = t728 * t453;
t411 = t454 * t575;
t573 = t727 * t452;
t538 = t454 * t573;
t343 = t538 - t411;
t207 = t343 * t309;
t455 = cos(qJ(2));
t658 = t453 * t455;
t420 = pkin(4) * t658;
t443 = t455 * pkin(7);
t444 = t455 * pkin(3);
t648 = t443 / 0.2e1 + t444 / 0.2e1;
t308 = -t391 * t727 + t540 * t728;
t678 = t308 * t814;
t473 = -t678 / 0.2e1 + t207 / 0.2e1 + t420 / 0.2e1 + t648;
t384 = t575 - t573;
t549 = -pkin(5) * t810 + qJ(6) * t384;
t600 = t810 * qJD(6);
t550 = t720 * t455;
t531 = -pkin(1) + t550;
t657 = t454 * qJ(3);
t660 = t452 * t454;
t741 = pkin(3) + pkin(7);
t302 = t453 * (t531 - t657) + t741 * t660;
t245 = -pkin(8) * t658 + t302;
t661 = t452 * qJ(3);
t783 = t453 * t741;
t493 = t661 + t783;
t456 = (t455 * t585 + pkin(1)) * t452 + (pkin(4) + t493) * t454;
t126 = t245 * t727 - t456 * t728;
t120 = -pkin(5) * t454 + t126;
t736 = t810 / 0.2e1;
t737 = -t810 / 0.2e1;
t492 = t120 * t736 + t126 * t737;
t740 = qJ(6) / 0.2e1;
t784 = -pkin(5) / 0.2e1;
t495 = t343 * t740 - t784 * t814;
t127 = t245 * t728 + t456 * t727;
t119 = t454 * qJ(6) + t127;
t552 = t119 / 0.2e1 - t127 / 0.2e1;
t3 = t384 * t552 + t492 - t495;
t825 = t3 * qJD(1);
t827 = qJD(5) * t549 + t600 + t825;
t826 = qJD(5) * t3;
t347 = t810 * t455;
t668 = t347 * qJ(6);
t412 = t455 * t575;
t344 = t455 * t573 - t412;
t723 = t344 * pkin(5);
t824 = t723 / 0.2e1 - t668 / 0.2e1 - t473;
t535 = -t573 / 0.2e1;
t734 = t384 / 0.2e1;
t466 = t411 / 0.2e1 + (t734 + t535) * t454;
t242 = t466 * qJD(5);
t221 = t343 * t454 - t344 * t455;
t793 = t221 * qJD(1);
t823 = t793 - t242;
t524 = t454 * t734 + t538 / 0.2e1 - t411 / 0.2e1;
t241 = t524 * qJD(5);
t441 = t455 * qJD(2);
t822 = t441 * t810 + t241 + t793;
t602 = t810 * qJD(2);
t638 = qJD(1) * t344;
t267 = -t602 + t638;
t162 = t344 * t734 - t347 * t737;
t794 = t162 * qJD(5);
t821 = -t267 * t343 - t794;
t820 = t343 * t638 - t794;
t145 = -t343 * t347 - t344 * t814;
t795 = t145 * qJD(1);
t664 = t384 * t347;
t665 = t810 * t344;
t112 = -t664 + t665;
t796 = t112 * qJD(5);
t819 = t796 - t795;
t376 = t810 ^ 2;
t742 = t384 ^ 2;
t230 = t742 - t376;
t818 = t230 * qJD(5);
t817 = t796 + t795 + qJD(2) * (t343 * t384 + t810 * t814);
t340 = t344 ^ 2;
t743 = t347 ^ 2;
t177 = t743 - t340;
t508 = qJD(1) * t177 + qJD(2) * t112;
t507 = qJD(1) * t112 + qJD(2) * t230;
t816 = qJD(2) * t145 + qJD(5) * t177;
t735 = -t384 / 0.2e1;
t300 = t343 * t810;
t669 = t814 * t384;
t556 = -t669 / 0.2e1;
t64 = 0.2e1 * t556 - t300;
t522 = t343 * t344 + t347 * t814;
t762 = qJD(1) * t522;
t815 = qJD(2) * t64 + t762;
t410 = t443 + t444;
t359 = t420 + t410;
t527 = t668 - t723;
t183 = t527 + t359;
t233 = -pkin(5) * t347 - qJ(6) * t344;
t422 = pkin(4) * t452 + qJ(3);
t271 = t422 - t549;
t310 = pkin(5) * t384 + qJ(6) * t810;
t812 = -t183 * t310 / 0.2e1 - t233 * t271 / 0.2e1 + t552 * t308;
t470 = t576 / 0.2e1 + t572 / 0.2e1;
t636 = qJD(1) * t454;
t563 = t347 * t636;
t603 = t347 * qJD(5);
t811 = -qJD(2) * t524 + t563 + t603;
t544 = qJD(2) * t466 - t563;
t597 = t384 * qJD(5);
t752 = t524 * qJD(1);
t809 = t597 + t752;
t733 = t412 / 0.2e1;
t465 = t733 + (t535 + t735) * t455;
t808 = -qJD(1) * t465 - t602;
t533 = -t376 / 0.2e1 - t742 / 0.2e1;
t176 = -0.1e1 / 0.2e1 + t533;
t484 = t665 / 0.2e1 + t664 / 0.2e1;
t806 = qJD(1) * t484 + qJD(2) * t176;
t553 = t814 / 0.2e1;
t650 = t470 * t454;
t259 = t553 + t650;
t604 = t347 * qJD(4);
t805 = qJD(3) * t259 + t604;
t753 = t466 * qJD(1);
t189 = t753 + t597;
t775 = t556 + t669 / 0.2e1;
t803 = qJD(1) * t775;
t692 = t127 * t384;
t100 = t692 / 0.2e1;
t101 = -t692 / 0.2e1;
t780 = t100 + t101;
t802 = qJD(1) * t780;
t799 = qJD(3) * t775;
t778 = t376 + t742;
t798 = qJD(4) * t778;
t797 = qJD(5) * t780;
t598 = t384 * qJD(2);
t568 = t810 * t598;
t792 = -qJD(1) * t162 + t568;
t520 = -t664 - t665;
t779 = t340 + t743;
t791 = qJD(1) * t779 + qJD(2) * t520;
t637 = qJD(1) * t347;
t570 = t344 * t637;
t790 = qJD(2) * t162 - t570;
t789 = qJD(2) * t221 - t454 * t603;
t788 = -qJD(2) * t775 - t762;
t787 = qJD(1) * t520 + qJD(2) * t778;
t786 = qJD(3) * t522 + qJD(4) * t779;
t628 = qJD(4) * t344;
t782 = t454 * (qJD(3) * t343 + t628);
t446 = t452 ^ 2;
t447 = t453 ^ 2;
t415 = t446 + t447;
t450 = t454 ^ 2;
t451 = t455 ^ 2;
t417 = t451 - t450;
t258 = (t737 + t470) * t455;
t325 = t343 * t636;
t777 = qJD(2) * t258 - t242 + t325;
t729 = t455 / 0.2e1;
t257 = t455 * t470 + t729 * t810;
t776 = -qJD(2) * t257 - t241 - t325;
t670 = t814 * qJ(6);
t724 = t343 * pkin(5);
t774 = -t670 / 0.2e1 + t724 / 0.2e1;
t586 = t455 * qJD(3);
t627 = qJD(4) * t454;
t773 = t586 - t627;
t270 = t598 - t637;
t187 = -qJD(1) * t258 + t598;
t601 = t810 * qJD(5);
t261 = t828 + t650;
t611 = t261 * qJD(1);
t770 = t611 - t601;
t612 = t259 * qJD(1);
t769 = t612 + t601;
t406 = pkin(2) * t454 - qJ(3) * t455;
t375 = qJ(4) * t454 + t406;
t307 = t375 * t453 + t410 * t452;
t680 = t307 * t452;
t390 = t410 * t453;
t306 = -t375 * t452 + t390;
t681 = t306 * t453;
t766 = t680 + t681;
t523 = t308 * t384 - t309 * t810;
t758 = qJD(2) * t523;
t756 = qJD(3) * t484;
t754 = qJD(4) * t523;
t588 = t454 * qJD(2);
t154 = t484 * qJD(4);
t157 = t520 * qJD(4);
t725 = pkin(8) * t454;
t228 = pkin(4) * t455 + t390 + (-t375 - t725) * t452;
t220 = t727 * t228;
t248 = t453 * t725 + t307;
t232 = t728 * t248;
t652 = -t220 / 0.2e1 - t232 / 0.2e1;
t574 = t727 * t248;
t577 = t728 * t228;
t651 = t577 / 0.2e1 - t574 / 0.2e1;
t750 = -t309 * t344 / 0.2e1 + t308 * t347 / 0.2e1;
t749 = t126 * t735 + t127 * t736;
t748 = t119 * t736 + t120 * t735;
t463 = t384 * t729 + t455 * t535 + t733;
t569 = t814 * t636;
t747 = qJD(2) * t463 + t569;
t746 = -qJD(2) * t465 - t569;
t745 = qJD(3) * t463 - qJD(4) * t466;
t744 = qJD(3) * t465 - qJD(4) * t524;
t739 = t344 / 0.2e1;
t738 = -t347 / 0.2e1;
t732 = -t422 / 0.2e1;
t730 = -t455 / 0.2e1;
t722 = t453 * pkin(4);
t721 = t455 * pkin(5);
t719 = t157 + t797;
t696 = t126 * t309;
t74 = -t696 / 0.2e1;
t9 = t74 + t696 / 0.2e1;
t717 = qJD(1) * t9;
t716 = qJD(2) * t780;
t130 = t232 + t220;
t655 = t455 * qJ(6);
t121 = t130 + t655;
t129 = -t574 + t577;
t122 = -t129 - t721;
t358 = (-t722 - t741) * t454;
t182 = t358 - t670 + t724;
t8 = t119 * t121 + t120 * t122 + t182 * t183;
t714 = t8 * qJD(1);
t713 = t154 + t797;
t672 = t309 * t454;
t277 = -t672 / 0.2e1;
t488 = t183 * t734 + t271 * t738;
t529 = -t721 / 0.2e1 - t651;
t35 = t277 + t488 + t529;
t711 = qJD(1) * t35;
t694 = t126 * t454;
t39 = -t183 * t344 + t233 * t347 - t694;
t710 = qJD(1) * t39;
t685 = t183 * t347;
t691 = t127 * t454;
t40 = -t233 * t344 - t685 - t691;
t709 = qJD(1) * t40;
t676 = t308 * t454;
t274 = -t676 / 0.2e1;
t485 = t344 * t732 + t359 * t736;
t51 = t274 + t485 + t652;
t708 = qJD(1) * t51;
t56 = t119 * t454 + t685;
t707 = qJD(1) * t56;
t69 = -t344 * t359 - t694;
t704 = qJD(1) * t69;
t70 = -t347 * t359 - t691;
t703 = qJD(1) * t70;
t12 = -t119 * t126 + t120 * t127 + t183 * t233;
t700 = t12 * qJD(1);
t698 = t121 * t810;
t697 = t122 * t384;
t690 = t129 * t384;
t41 = t119 * t343 - t120 * t814;
t13 = t121 * t344 - t122 * t347 - t41;
t689 = t13 * qJD(1);
t688 = t130 * t810;
t14 = (t119 - t127) * t347 + (t120 - t126) * t344;
t687 = t14 * qJD(1);
t43 = -t126 * t814 + t127 * t343;
t16 = t129 * t347 + t130 * t344 - t43;
t686 = t16 * qJD(1);
t19 = -t126 * t129 + t127 * t130 + t358 * t359;
t684 = t19 * qJD(1);
t29 = t119 * t455 + t121 * t454 + t182 * t347 - t183 * t814;
t683 = t29 * qJD(1);
t30 = -t120 * t455 - t122 * t454 - t182 * t344 + t183 * t343;
t682 = t30 * qJD(1);
t675 = t308 * t455;
t671 = t309 * t455;
t37 = -t126 * t455 + t129 * t454 + t343 * t359 - t344 * t358;
t667 = t37 * qJD(1);
t38 = t127 * t455 + t130 * t454 + t347 * t358 - t359 * t814;
t666 = t38 * qJD(1);
t663 = t41 * qJD(1);
t662 = t43 * qJD(1);
t659 = t452 * t455;
t301 = -t452 * t531 + t454 * t493;
t57 = t301 * t660 - t306 * t659 + (-t302 * t454 + t307 * t455) * t453;
t654 = t57 * qJD(1);
t579 = t741 * t454;
t58 = t301 * t306 + t302 * t307 - t410 * t579;
t653 = t58 * qJD(1);
t140 = (-t301 * t452 + t302 * t453) * t454;
t647 = qJD(1) * t140;
t141 = -t301 * t659 + t302 * t658;
t646 = qJD(1) * t141;
t635 = qJD(2) * t271;
t634 = qJD(2) * t422;
t260 = t553 + t828;
t632 = qJD(3) * t260;
t631 = qJD(3) * t384;
t630 = qJD(3) * t450;
t629 = qJD(3) * t454;
t626 = qJD(5) * t308;
t625 = qJD(5) * t344;
t624 = qJD(5) * t454;
t623 = qJD(6) * t384;
t578 = t741 * t455;
t109 = t301 * t455 + (t306 + (-t578 - t410) * t453) * t454;
t622 = t109 * qJD(1);
t110 = t302 * t455 - t410 * t660 + (-t452 * t578 + t307) * t454;
t621 = t110 * qJD(1);
t223 = -t347 * t455 + t454 * t814;
t615 = t223 * qJD(1);
t262 = 0.2e1 * t828;
t239 = t262 * qJD(1);
t299 = t309 * qJD(5);
t528 = -pkin(2) * t455 - t657;
t394 = -pkin(1) + t528;
t320 = t394 * t455 + t406 * t454;
t608 = t320 * qJD(1);
t321 = -t394 * t454 + t406 * t455;
t607 = t321 * qJD(1);
t606 = t344 * qJD(6);
t345 = t415 * t455 * t454;
t605 = t345 * qJD(1);
t383 = t415 * t451;
t599 = t383 * qJD(1);
t386 = t417 * t452;
t596 = t386 * qJD(1);
t387 = t417 * t453;
t595 = t387 * qJD(1);
t551 = -t446 / 0.2e1 - t447 / 0.2e1;
t393 = -0.1e1 / 0.2e1 + t551;
t594 = t393 * qJD(2);
t593 = t415 * qJD(2);
t592 = t417 * qJD(1);
t591 = t450 * qJD(1);
t590 = t452 * qJD(2);
t589 = t453 * qJD(2);
t440 = t454 * qJD(6);
t587 = t455 * qJD(1);
t584 = -pkin(7) / 0.2e1 - pkin(3) / 0.2e1;
t583 = pkin(1) * t636;
t582 = pkin(1) * t587;
t581 = pkin(7) * t588;
t566 = t452 * t589;
t565 = t452 * t441;
t319 = t810 * t597;
t562 = t394 * t406 * qJD(1);
t561 = t394 * t636;
t560 = t452 * t591;
t559 = t453 * t591;
t419 = t454 * t441;
t418 = t454 * t587;
t276 = t672 / 0.2e1;
t548 = -0.2e1 * t452 * t658;
t547 = -t207 + t678;
t486 = -t681 / 0.2e1 - t680 / 0.2e1;
t137 = t486 + t648;
t448 = qJD(2) * qJ(3);
t546 = qJD(1) * t137 + t448;
t326 = t344 * t636;
t541 = qJD(2) * t262 + t326;
t421 = qJD(5) + t636;
t539 = t584 * t454;
t537 = t453 * t418;
t532 = t651 + t721;
t530 = t655 - t652;
t494 = t121 * t740 + t122 * t784;
t1 = (t126 / 0.2e1 - t120 / 0.2e1) * t309 + t494 + t812;
t45 = t271 * t310;
t526 = -t1 * qJD(1) + t45 * qJD(2);
t525 = qJD(2) * t9 + qJD(3) * t780;
t519 = (t300 + t669) * qJD(3);
t115 = t271 * t384 + t310 * t810;
t459 = t233 * t737 + t310 * t739 - t488;
t26 = t276 + t459 + t532;
t518 = -qJD(1) * t26 + qJD(2) * t115;
t116 = t271 * t810 - t310 * t384;
t275 = t676 / 0.2e1;
t460 = t183 * t737 + t233 * t734 + t271 * t739 + t310 * t738;
t24 = t275 + t460 + t530;
t517 = qJD(1) * t24 - qJD(2) * t116;
t476 = (-t722 / 0.2e1 + t584) * t454;
t458 = t476 + t750;
t20 = t458 + t748 + t774;
t516 = -qJD(1) * t20 + t758;
t31 = t458 + t749;
t515 = -qJD(1) * t31 + t758;
t491 = t698 / 0.2e1 - t697 / 0.2e1;
t22 = t491 + t824;
t514 = -qJD(1) * t22 + t635;
t490 = -t690 / 0.2e1 - t688 / 0.2e1;
t34 = t473 + t490;
t513 = qJD(1) * t34 + t634;
t44 = -t126 * t347 + t127 * t344;
t512 = -qJD(1) * t44 - t756;
t42 = t119 * t344 - t120 * t347;
t511 = qJD(1) * t42 + t756;
t487 = -t301 * t453 / 0.2e1 - t302 * t452 / 0.2e1;
t134 = t539 - t487;
t379 = t415 * t720;
t506 = -qJD(1) * t134 - qJD(2) * t379;
t501 = -qJD(1) * t233 - qJD(2) * t310;
t500 = -qJD(3) * t257 - qJD(4) * t262;
t499 = -qJD(3) * t258 - qJD(4) * t261;
t498 = qJD(3) * t466 - t628;
t496 = -qJD(3) * t524 + qJD(5) * t126;
t467 = -t347 * t732 + t359 * t735 + t276;
t50 = t467 + t651;
t483 = -qJD(1) * t50 + t422 * t598;
t149 = t730 - t162;
t482 = qJD(1) * t149 + t568;
t477 = t455 * t627 + t630;
t132 = qJD(2) * t261 + t326 + t625;
t474 = (-qJD(2) * t343 + t603) * t344;
t324 = t450 + t743;
t468 = qJD(1) * t324 - t347 * t598 + t624;
t464 = qJD(2) * t528 + t586;
t462 = t476 - t750;
t457 = t119 * t734 + t492 + t495;
t449 = qJ(3) * qJD(3);
t437 = pkin(7) * t441;
t429 = -t587 / 0.2e1;
t428 = t587 / 0.2e1;
t427 = t441 / 0.2e1;
t416 = t453 * t441;
t407 = t417 * qJD(2);
t405 = t421 * qJ(6);
t392 = 0.1e1 / 0.2e1 + t551;
t380 = qJD(5) * t729 + t418;
t372 = qJD(3) * t810;
t370 = t384 * qJD(4);
t298 = t347 * t623;
t238 = t262 * qJD(5);
t236 = t261 * qJD(5);
t206 = qJD(2) * t742 - t384 * t637;
t188 = -t601 + t239;
t175 = 0.1e1 / 0.2e1 + t533;
t173 = t176 * qJD(3);
t171 = t176 * qJD(4);
t170 = t175 * qJD(3);
t169 = t175 * qJD(4);
t151 = qJD(2) * t223 + t344 * t624;
t150 = t730 + t162;
t139 = (-qJD(2) * t814 - t625) * t347;
t136 = -t486 + t648;
t135 = t539 + t487;
t131 = t238 - t615;
t123 = t127 * qJD(5);
t104 = t384 * t441 + t236 + t615;
t93 = t637 * t814 + t794;
t60 = t64 * qJD(3);
t59 = t270 * t814 + t794;
t53 = -t467 + t651;
t52 = t275 - t485 + t652;
t36 = t276 - t488 + t529;
t33 = t473 - t490;
t32 = t462 - t749;
t27 = t277 - t459 + t532;
t25 = t274 - t460 + t530;
t23 = t491 - t824;
t21 = t462 - t748 + t774;
t11 = t101 + t457;
t7 = t9 * qJD(5);
t4 = t100 - t457;
t2 = t74 + t120 * t309 / 0.2e1 + t494 - t812;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, t407, 0, -t419, 0, 0, -pkin(1) * t588, -pkin(1) * t441, 0, 0, 0, 0, 0, t419, t407, -t419, 0, qJD(2) * t321 - t454 * t586, -qJD(2) * t320 + t630 (qJD(2) * t406 - t629) * t394, -t446 * t419, t548 * t588, -t386 * qJD(2), -t447 * t419, -t387 * qJD(2), t419, t109 * qJD(2) + t452 * t477, -t110 * qJD(2) + t453 * t477, -qJD(2) * t57 + qJD(3) * t345 + qJD(4) * t383, qJD(2) * t58 - qJD(3) * t140 - qJD(4) * t141, t139, -t816, t151, t474, -t789, t419, qJD(2) * t37 + qJD(5) * t70 + (qJD(3) * t814 + t604) * t454, -qJD(2) * t38 - qJD(5) * t69 - t782, qJD(2) * t16 + t786, qJD(2) * t19 + qJD(3) * t43 + qJD(4) * t44, t139, t151, t816, t419, t789, t474, t814 * t629 + qJD(2) * t30 + qJD(5) * t40 + (-t606 + t627) * t347, qJD(2) * t13 + qJD(5) * t14 + t344 * t440 + t786, qJD(2) * t29 + qJD(5) * t39 + qJD(6) * t324 + t782, qJD(2) * t8 + qJD(3) * t41 + qJD(4) * t42 + qJD(5) * t12 + qJD(6) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, t592, t441, -t418, -t588, 0, -t437 - t583, t581 - t582, 0, 0, 0, -t441, t588, t418, t592, -t418, t464, t437 + t607, -t581 - t608, pkin(7) * t464 + t562 (-t446 * t587 + t566) * t454 (qJD(1) * t548 + (-t446 + t447) * qJD(2)) * t454, t416 - t596 (-t447 * t587 - t566) * t454, -t565 - t595, t418, -t741 * t452 * t588 + t622 + ((t550 - t657) * qJD(2) + t773) * t453, -t621 + (t661 - t783) * t588 + (-qJD(2) * t550 - t773) * t452, -qJD(2) * t766 - t654, t653 + (-qJ(3) * t579 + t720 * t766) * qJD(2) + t136 * qJD(3) + t135 * qJD(4), t59, -t817, t104, t821, -t822, t380, t667 + (t343 * t422 + t358 * t810 - t675) * qJD(2) + t53 * qJD(5) + t745, -t666 + (t358 * t384 + t422 * t814 - t671) * qJD(2) + t52 * qJD(5) + t500, t686 + (t547 - t688 - t690) * qJD(2) + t60 + t719, t684 + (-t129 * t308 + t130 * t309 + t358 * t422) * qJD(2) + t33 * qJD(3) + t32 * qJD(4) + t7, t59, t104, t817, t380, t822, t821, t682 + (t182 * t810 + t271 * t343 - t675) * qJD(2) + t27 * qJD(5) + t150 * qJD(6) + t745, t689 + (t547 + t697 - t698) * qJD(2) + t60 + t157 + t4 * qJD(5) + t261 * qJD(6), t683 + (-t182 * t384 - t271 * t814 + t671) * qJD(2) + t25 * qJD(5) - t298 - t500, t714 + (t121 * t309 + t122 * t308 + t182 * t271) * qJD(2) + t23 * qJD(3) + t21 * qJD(4) + t2 * qJD(5) + t36 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t441, -t418, t591, t437 - t561, 0, 0, 0, 0, 0, 0, t416 + t560, t559 - t565, t605, qJD(2) * t136 - t647, 0, 0, 0, 0, 0, 0, t236 + t747, t776, t815, qJD(2) * t33 + t519 + t662 + t713, 0, 0, 0, 0, 0, 0, -qJD(5) * t260 + t747, t815, -t776, qJD(2) * t23 + qJD(5) * t11 + qJD(6) * t260 + t154 + t519 + t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t452 * t587 - t589) * t454 (t453 * t587 + t590) * t454, t599, qJD(2) * t135 - t646, 0, 0, 0, 0, 0, 0, -t544, -t541, t791, qJD(2) * t32 - t512, 0, 0, 0, 0, 0, 0, -t544, t791, t541, qJD(2) * t21 + t511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t790, -t508, t132, -t790, t811, t427, qJD(2) * t53 + qJD(3) * t261 - t123 + t703, qJD(2) * t52 + t496 - t704, t716, t525, t790, t132, t508, t427, -t811, -t790, qJD(2) * t27 - t123 - t632 + t709, t4 * qJD(2) + qJD(5) * t527 + t606 + t687, qJD(2) * t25 + t440 - t496 + t710, t700 + t2 * qJD(2) + t11 * qJD(3) + (-pkin(5) * t127 - qJ(6) * t126) * qJD(5) + t119 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t150 - t570, t132, t468, qJD(2) * t36 + qJD(5) * t119 + t632 + t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, -t592, 0, t418, 0, 0, t583, t582, 0, 0, 0, 0, 0, -t418, -t592, t418, 0, -t607, t608, -t562, t446 * t418, 0.2e1 * t452 * t537, t596, t447 * t418, t595, -t418, -t622, t621, t654, qJD(3) * t137 - qJD(4) * t134 - t653, t93, -t819, t131, t820, t823, -t380, -qJD(5) * t50 - t667 + t744, -qJD(5) * t51 + t499 + t666, -t686 + t719 + t799, qJD(3) * t34 - qJD(4) * t31 - t684 + t7, t93, t131, t819, -t380, -t823, t820, -qJD(5) * t26 - qJD(6) * t149 - t682 + t744, qJD(6) * t262 + t157 - t689 + t799 - t826, -qJD(5) * t24 - t298 - t499 - t683, -qJD(3) * t22 - qJD(4) * t20 - qJD(5) * t1 - qJD(6) * t35 - t714; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t449, 0, 0, 0, 0, 0, 0, qJD(3) * t452, qJD(3) * t453, t415 * qJD(4), -qJD(4) * t379 + t449, -t319, -t818, 0, t319, 0, 0, t422 * t597 + t372, -t422 * t601 + t631, t798, qJD(3) * t422 + t754, -t319, 0, t818, 0, 0, t319, qJD(5) * t115 - t384 * t600 + t372, t798, qJD(5) * t116 + qJD(6) * t742 - t631, t754 + qJD(5) * t45 + (qJD(3) - t623) * t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t448, 0, 0, 0, 0, 0, 0, t590, t589, 0, qJD(4) * t392 + t546, 0, 0, 0, 0, 0, 0, -t808, t187, t803, t169 + t513, 0, 0, 0, 0, 0, 0, -t808, t803, -t187, t169 + t514; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t593, qJD(3) * t392 + t506, 0, 0, 0, 0, 0, 0, -t752, -t611, t787, t170 + t515, 0, 0, 0, 0, 0, 0, -t752, t787, t611, t170 + t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t792, -t507, t188, t792, -t189, t429, -t299 + t483, -t422 * t602 + t626 - t708, t802, t717, -t792, t188, t507, t429, t189, t792, -t299 + t518, -t827, -t517 - t626 (-pkin(5) * t309 - qJ(6) * t308) * qJD(5) + t309 * qJD(6) + t526; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t482, t188, t206, -t271 * t598 + t299 - t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, -t591, t561, 0, 0, 0, 0, 0, 0, -t560, -t559, -t605, -qJD(2) * t137 + t647, 0, 0, 0, 0, 0, 0, t238 + t746, t777, t788, -qJD(2) * t34 - t662 + t713, 0, 0, 0, 0, 0, 0, -qJD(5) * t259 + t746, t788, -t777, qJD(2) * t22 + qJD(6) * t259 + t154 - t663 + t826; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t448, 0, 0, 0, 0, 0, 0, -t590, -t589, 0, qJD(4) * t393 - t546, 0, 0, 0, 0, 0, 0, t808, -t187, -t803, t171 - t513, 0, 0, 0, 0, 0, 0, t808, -t803, t187, t171 - t514; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t594, 0, 0, 0, 0, 0, 0, 0, 0, 0, t806, 0, 0, 0, 0, 0, 0, 0, 0, 0, t806; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, -t189, 0, t802, 0, 0, 0, 0, 0, 0, -t769, 0, t189, t827; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t769; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t452 * t418, -t537, -t599, qJD(2) * t134 + t646, 0, 0, 0, 0, 0, 0, -t811, t132, -t791, qJD(2) * t31 + t512, 0, 0, 0, 0, 0, 0, -t811, -t791, -t132, qJD(2) * t20 + qJD(5) * t233 + qJD(6) * t347 - t511; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t593, -qJD(3) * t393 - t506, 0, 0, 0, 0, 0, 0, t809, t770, -t787, -t173 - t515, 0, 0, 0, 0, 0, 0, t809, -t787, -t770, qJD(5) * t310 - t173 - t516 - t623; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t594, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t806, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t806; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, t267, 0, 0, 0, 0, 0, 0, 0, 0, t270, 0, -t267, -t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t790, t508, -t541, t790, t544, t427, qJD(2) * t50 - qJD(3) * t262 + t604 - t703, qJD(2) * t51 + t498 + t704, -t716, -t525, -t790, -t541, -t508, t427, -t544, t790, qJD(2) * t26 - t709 + t805, qJD(2) * t3 - t687, qJD(2) * t24 + t440 - t498 - t710, qJ(6) * t440 + qJD(2) * t1 - qJD(3) * t3 - qJD(4) * t233 - t700; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t792, t507, -t239, -t792, t753, t428, -t370 - t483, t708 + (qJD(4) + t634) * t810, -t802, -t717, t792, -t239, -t507, t428, -t753, -t792, -t370 - t518, t825, -qJD(4) * t810 + t517, -qJD(4) * t310 - t526; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239, t753, 0, -t802, 0, 0, 0, 0, 0, 0, t612, 0, -t753, -t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, -t267, 0, 0, 0, 0, 0, 0, 0, 0, -t270, 0, t267, t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t421, t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t149 + t570, -t541, -t468, -qJ(6) * t624 + qJD(2) * t35 - t707 - t805; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t482, -t239, -t206, t711 + (qJD(4) + t635) * t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421, -t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;

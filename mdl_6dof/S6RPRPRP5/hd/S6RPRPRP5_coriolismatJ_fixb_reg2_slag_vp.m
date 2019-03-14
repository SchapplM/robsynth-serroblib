% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPRP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:27
% EndTime: 2019-03-09 03:16:52
% DurationCPUTime: 17.46s
% Computational Cost: add. (19743->652), mult. (39893->764), div. (0->0), fcn. (46728->8), ass. (0->489)
t475 = sin(pkin(9));
t477 = cos(pkin(9));
t478 = sin(qJ(3));
t767 = cos(qJ(3));
t450 = t475 * t478 - t477 * t767;
t474 = sin(pkin(10));
t765 = sin(qJ(5));
t599 = t765 * t474;
t402 = t450 * t599;
t476 = cos(pkin(10));
t766 = cos(qJ(5));
t601 = t766 * t476;
t567 = t450 * t601;
t342 = -t567 + t402;
t704 = t342 * qJ(6);
t598 = t765 * t476;
t602 = t766 * t474;
t452 = t602 + t598;
t338 = t452 * t450;
t764 = t338 * pkin(5);
t571 = t704 / 0.2e1 + t764 / 0.2e1;
t772 = t452 / 0.2e1;
t448 = t599 - t601;
t775 = -t448 / 0.2e1;
t760 = pkin(7) + qJ(2);
t574 = t760 * t475;
t560 = t478 * t574;
t459 = t760 * t477;
t605 = t767 * t459;
t498 = t605 - t560;
t357 = t476 * t498;
t469 = -t477 * pkin(2) - pkin(1);
t604 = t767 * t475;
t684 = t478 * t477;
t454 = t604 + t684;
t553 = t450 * pkin(3) - qJ(4) * t454;
t500 = t469 + t553;
t217 = t474 * t500 + t357;
t349 = t474 * t454;
t185 = -pkin(8) * t349 + t217;
t487 = t474 * t498;
t216 = t476 * t500 - t487;
t686 = t454 * t476;
t479 = t450 * pkin(4) - pkin(8) * t686 + t216;
t104 = t185 * t766 + t479 * t765;
t93 = t450 * qJ(6) + t104;
t103 = t185 * t765 - t766 * t479;
t94 = -t450 * pkin(5) + t103;
t758 = t772 * t94 + t775 * t93;
t858 = -t571 - t758;
t384 = t459 * t478 + t767 * t574;
t311 = pkin(4) * t349 + t384;
t687 = t454 * t452;
t403 = t454 * t599;
t404 = t454 * t601;
t822 = -t403 + t404;
t552 = pkin(5) * t687 - qJ(6) * t822;
t145 = t311 + t552;
t468 = -pkin(4) * t476 - pkin(3);
t551 = pkin(5) * t448 - qJ(6) * t452;
t346 = t468 + t551;
t759 = pkin(8) + qJ(4);
t458 = t759 * t476;
t573 = t759 * t474;
t385 = t458 * t766 - t573 * t765;
t693 = t385 * t450;
t781 = t822 / 0.2e1;
t491 = t145 * t772 + t346 * t781 - t693 / 0.2e1;
t356 = t476 * t384;
t750 = qJ(4) * t450;
t762 = t454 * pkin(3);
t359 = t750 + t762;
t252 = t474 * t359 - t356;
t347 = t474 * t450;
t196 = pkin(8) * t347 + t252;
t600 = t765 * t196;
t251 = t476 * t359 + t384 * t474;
t689 = t450 * t476;
t182 = pkin(4) * t454 + pkin(8) * t689 + t251;
t603 = t766 * t182;
t677 = t603 / 0.2e1 - t600 / 0.2e1;
t857 = -t491 - t677;
t614 = t454 * qJD(3);
t708 = t687 * t454;
t712 = t338 * t450;
t171 = t708 - t712;
t839 = qJD(1) * t171;
t856 = -t448 * t614 - t839;
t834 = t687 * qJD(1);
t176 = qJD(3) * t448 + t834;
t164 = -t687 * t772 + t775 * t822;
t835 = t164 * qJD(5);
t855 = -t176 * t338 - t835;
t854 = t338 * t834 - t835;
t709 = t687 * t342;
t713 = t338 * t822;
t135 = t709 - t713;
t836 = t135 * qJD(1);
t847 = t448 * t687;
t132 = t452 * t822 - t847;
t837 = t132 * qJD(5);
t853 = t837 - t836;
t579 = t450 * t772;
t496 = t598 / 0.2e1 + t602 / 0.2e1;
t846 = t450 * t496;
t549 = t579 + t846;
t618 = t450 * qJD(1);
t845 = t822 * t618;
t852 = qJD(3) * t549 + t845;
t441 = t448 ^ 2;
t784 = t452 ^ 2;
t361 = t441 + t784;
t851 = qJD(4) * t361;
t336 = t784 - t441;
t850 = t336 * qJD(5);
t702 = t342 * t448;
t711 = t338 * t452;
t849 = t837 + t836 + qJD(3) * (t702 - t711);
t523 = qJD(1) * t132 + qJD(3) * t336;
t576 = 0.2e1 * t781;
t486 = t452 * t576 + t847;
t848 = qJD(1) * t486 + qJD(3) * t361;
t785 = t822 ^ 2;
t825 = -t687 / 0.2e1;
t710 = t687 * qJ(6);
t617 = t452 * qJD(5);
t800 = t549 * qJD(1);
t193 = -t617 - t800;
t484 = (t772 - t496) * t450;
t801 = t484 * qJD(1);
t844 = t617 + t801;
t540 = -t702 - t711;
t542 = -t709 - t713;
t812 = qJD(1) * t542;
t843 = -qJD(3) * t540 - t812;
t227 = t484 * qJD(5);
t541 = t708 + t712;
t813 = qJD(1) * t541;
t842 = t813 - t227;
t700 = t822 * t454;
t701 = t342 * t450;
t168 = -t700 + t701;
t840 = qJD(1) * t168;
t838 = qJD(5) * t687;
t226 = t549 * qJD(5);
t688 = t454 * t448;
t833 = -qJD(3) * t688 - t226 - t813;
t757 = t103 * t772 + t104 * t775;
t653 = qJD(3) * t452;
t592 = t448 * t653;
t832 = -qJD(1) * t164 + t592;
t649 = qJD(4) * t450;
t831 = -qJD(2) * t168 + t687 * t649;
t660 = qJD(1) * t822;
t595 = t687 * t660;
t830 = qJD(3) * t164 - t595;
t647 = qJD(5) * t450;
t829 = qJD(3) * t171 + t822 * t647;
t358 = pkin(5) * t452 + qJ(6) * t448;
t616 = t452 * qJD(6);
t828 = qJD(5) * t358 - t616;
t826 = t687 ^ 2;
t156 = -t785 + t826;
t526 = qJD(1) * t156 - qJD(3) * t132;
t827 = qJD(3) * t135 - qJD(5) * t156;
t782 = qJ(6) / 0.2e1;
t796 = -t618 - qJD(5);
t824 = t822 * t796;
t594 = t687 * t618;
t442 = t450 ^ 2;
t445 = t454 ^ 2;
t821 = -t445 - t442;
t610 = t445 - t442;
t558 = t601 / 0.2e1;
t208 = t454 * t558 - t404 / 0.2e1;
t820 = t208 * qJD(5) - t852;
t207 = -t403 / 0.2e1 + t404 / 0.2e1 + (t558 - t599 / 0.2e1) * t454;
t640 = t484 * qJD(3);
t819 = t207 * qJD(5) + t640 + t845;
t763 = t822 * pkin(5);
t197 = t710 + t763;
t774 = t448 / 0.2e1;
t779 = -t358 / 0.2e1;
t818 = t197 * t774 - t687 * t779;
t817 = -t251 * t474 + t252 * t476;
t674 = -t684 / 0.2e1 - t604 / 0.2e1;
t794 = t448 * t781 + t452 * t825;
t499 = -t674 - t794;
t815 = qJD(1) * t499;
t814 = qJD(1) * t540;
t811 = qJD(2) * t499;
t518 = -t674 + t794;
t810 = qJD(2) * t518;
t809 = qJD(2) * t540;
t383 = t458 * t765 + t573 * t766;
t538 = t383 * t452 - t448 * t385;
t808 = qJD(3) * t538;
t378 = t687 / 0.2e1;
t562 = t825 + t378;
t805 = qJD(3) * t562;
t803 = qJD(4) * t499;
t802 = qJD(4) * t538;
t70 = t486 * qJD(4);
t799 = t549 * qJD(4);
t798 = t562 * qJD(2);
t797 = t688 * qJD(1);
t793 = t383 * t781 + t385 * t825;
t204 = t207 * qJD(4);
t792 = qJD(2) * t549 - t204;
t203 = t208 * qJD(4);
t791 = -qJD(2) * t484 - t104 * qJD(5) + t203;
t790 = qJD(2) * t541 - t649 * t822;
t789 = qJD(4) * t518 + (-t338 * t448 + t342 * t452) * qJD(2);
t543 = t826 + t785;
t788 = qJD(1) * t543 + qJD(3) * t486;
t787 = qJD(2) * t688 - qJD(4) * t484;
t786 = qJD(2) * t542 + qJD(4) * t543;
t783 = -pkin(5) / 0.2e1;
t780 = -t346 / 0.2e1;
t778 = -t383 / 0.2e1;
t777 = t383 / 0.2e1;
t776 = -t402 / 0.2e1;
t773 = -t452 / 0.2e1;
t770 = t454 / 0.2e1;
t769 = -t468 / 0.2e1;
t733 = t103 * t448;
t493 = t733 / 0.2e1;
t6 = -t733 / 0.2e1 + t493;
t768 = t6 * qJD(5) + t70;
t761 = t454 * pkin(5);
t756 = qJD(3) * t6;
t175 = t765 * t182;
t190 = t766 * t196;
t111 = t190 + t175;
t437 = t454 * qJ(6);
t100 = t437 + t111;
t110 = -t600 + t603;
t101 = -t110 - t761;
t413 = pkin(4) * t347;
t312 = -t413 + t498;
t146 = t312 - t704 - t764;
t7 = t100 * t93 + t101 * t94 + t145 * t146;
t755 = t7 * qJD(1);
t8 = -t103 * t93 + t104 * t94 + t145 * t197;
t754 = t8 * qJD(1);
t9 = -t100 * t687 + t101 * t822 + t338 * t93 + t342 * t94;
t753 = t9 * qJD(1);
t734 = t103 * t385;
t63 = -t734 / 0.2e1;
t16 = t63 + t734 / 0.2e1;
t749 = qJD(1) * t16;
t582 = t338 * t777;
t512 = t346 * t770 - t582;
t517 = t100 * t772 + t101 * t774;
t703 = t342 * t385;
t31 = -t703 / 0.2e1 - t512 + t517;
t747 = qJD(1) * t31;
t233 = t703 / 0.2e1;
t501 = t468 * t770 + t233 - t582;
t516 = t110 * t774 + t111 * t773;
t36 = t501 + t516;
t746 = qJD(1) * t36;
t39 = -t687 * t93 + t822 * t94;
t745 = qJD(1) * t39;
t732 = t103 * t450;
t40 = t145 * t687 - t197 * t822 - t732;
t744 = qJD(1) * t40;
t726 = t145 * t822;
t729 = t104 * t450;
t41 = t197 * t687 + t726 - t729;
t743 = qJD(1) * t41;
t42 = t103 * t822 - t104 * t687;
t742 = qJD(1) * t42;
t45 = t450 * t93 - t726;
t741 = qJD(1) * t45;
t52 = -t311 * t687 + t732;
t740 = qJD(1) * t52;
t53 = t311 * t822 - t729;
t739 = qJD(1) * t53;
t695 = t384 * t454;
t722 = t217 * t476;
t723 = t216 * t474;
t74 = t695 + (-t722 + t723) * t450;
t738 = qJD(1) * t74;
t10 = (t104 - t93) * t822 + (t103 - t94) * t687;
t735 = t10 * qJD(1);
t606 = t93 / 0.2e1 - t104 / 0.2e1;
t607 = t103 / 0.2e1 - t94 / 0.2e1;
t11 = t448 * t606 + t452 * t607 + t571;
t728 = t11 * qJD(1);
t14 = t103 * t342 + t104 * t338 - t110 * t822 - t111 * t687;
t727 = t14 * qJD(1);
t17 = -t103 * t110 + t104 * t111 + t311 * t312;
t725 = t17 * qJD(1);
t21 = t100 * t450 - t145 * t342 - t146 * t822 + t454 * t93;
t724 = t21 * qJD(1);
t22 = -t101 * t450 - t145 * t338 + t146 * t687 - t454 * t94;
t721 = t22 * qJD(1);
t25 = t145 * t454 - t338 * t94 + t342 * t93;
t720 = t25 * qJD(1);
t718 = t251 * t476;
t717 = t252 * t474;
t34 = -t103 * t454 + t110 * t450 - t311 * t338 + t312 * t687;
t707 = t34 * qJD(1);
t35 = -t104 * t454 - t111 * t450 + t311 * t342 + t312 * t822;
t699 = t35 * qJD(1);
t38 = -t103 * t338 + t104 * t342 + t311 * t454;
t698 = t38 * qJD(1);
t696 = t383 * t454;
t692 = t385 * t454;
t480 = t491 - t677;
t608 = -t761 / 0.2e1;
t43 = t608 + t480;
t691 = t43 * qJD(1);
t690 = t448 * t450;
t544 = -t216 * t476 - t217 * t474;
t46 = (t717 + t718) * t454 + t544 * t450;
t685 = t46 * qJD(1);
t580 = t450 * t778;
t490 = t311 * t774 - t687 * t769 + t580;
t578 = -t175 / 0.2e1 - t190 / 0.2e1;
t48 = t490 + t578;
t683 = t48 * qJD(1);
t51 = t216 * t251 + t217 * t252 + t384 * t498;
t682 = t51 * qJD(1);
t54 = t251 * t450 - t347 * t384 + (t216 + t487) * t454;
t681 = t54 * qJD(1);
t55 = (-t217 + t357) * t454 + (-t252 - t356) * t450;
t680 = t55 * qJD(1);
t675 = t402 / 0.2e1 - t567 / 0.2e1;
t470 = t474 ^ 2;
t472 = t476 ^ 2;
t460 = t470 + t472;
t461 = t475 ^ 2 + t477 ^ 2;
t116 = t544 * t454;
t673 = qJD(1) * t116;
t167 = t700 + t701;
t666 = qJD(1) * t167;
t352 = t821 * t476;
t659 = qJD(1) * t352;
t245 = t690 / 0.2e1 + t675;
t655 = qJD(3) * t245;
t652 = qJD(3) * t468;
t651 = qJD(3) * t476;
t650 = qJD(4) * t687;
t648 = qJD(5) * t383;
t575 = -t470 / 0.2e1 - t472 / 0.2e1;
t494 = t575 * t750 - t762 / 0.2e1;
t513 = t718 / 0.2e1 + t717 / 0.2e1;
t105 = t494 - t513;
t646 = t105 * qJD(1);
t200 = -t450 * t498 + t695;
t643 = t200 * qJD(1);
t240 = (t773 - t496) * t450;
t221 = t240 * qJD(1);
t223 = t245 * qJD(1);
t246 = -t690 / 0.2e1 + t675;
t638 = t246 * qJD(1);
t247 = t776 + (t558 + t774) * t450;
t637 = t247 * qJD(1);
t636 = t247 * qJD(3);
t248 = t776 + (t558 + t775) * t450;
t635 = t248 * qJD(1);
t253 = t460 * t445;
t634 = t253 * qJD(1);
t267 = 0.2e1 * t378;
t630 = t267 * qJD(1);
t556 = t575 * t454;
t302 = t556 + t674;
t629 = t302 * qJD(1);
t305 = t610 * t474;
t628 = t305 * qJD(1);
t306 = t821 * t474;
t627 = t306 * qJD(1);
t307 = t610 * t476;
t626 = t307 * qJD(1);
t625 = t610 * qJD(1);
t624 = t687 * qJD(6);
t623 = t349 * qJD(1);
t415 = t470 * t450;
t416 = t472 * t450;
t351 = t415 + t416;
t622 = t351 * qJD(1);
t621 = t821 * qJD(1);
t360 = t385 * qJD(5);
t620 = t674 * qJD(1);
t619 = t448 * qJD(5);
t432 = t450 * qJD(3);
t430 = t450 * qJD(6);
t615 = t454 * qJD(1);
t456 = t461 * qJ(2);
t613 = t456 * qJD(1);
t612 = t460 * qJD(3);
t611 = t461 * qJD(1);
t590 = t474 * t651;
t382 = t448 * t617;
t381 = t450 * t615;
t380 = t450 * t614;
t586 = t476 * t649;
t585 = t476 * t615;
t412 = t476 * t614;
t572 = t385 * t338 + t342 * t383;
t309 = -qJD(5) * t674 + t381;
t570 = qJD(1) * t469 + qJD(2);
t566 = t450 * t585;
t565 = t474 * t381;
t561 = -t437 + t578;
t550 = t476 * t565;
t520 = t338 * t782 + t342 * t783;
t3 = -t448 * t607 + t452 * t606 + t520;
t548 = t3 * qJD(1);
t547 = qJD(1) * t6;
t515 = t145 * t779 + t197 * t780;
t519 = t100 * t782 + t101 * t783;
t1 = t383 * t606 + t385 * t607 + t515 + t519;
t115 = t346 * t358;
t546 = -t1 * qJD(1) + t115 * qJD(3);
t545 = qJD(3) * t16;
t183 = t346 * t452 + t358 * t448;
t26 = t480 - t761 + t818;
t534 = qJD(1) * t26 + qJD(3) * t183;
t184 = t346 * t448 - t358 * t452;
t482 = t145 * t774 + t197 * t773 - t687 * t780 + t779 * t822 + t580;
t29 = t482 + t561;
t533 = -qJD(1) * t29 - qJD(3) * t184;
t488 = t605 / 0.2e1 - t560 / 0.2e1;
t485 = -t413 / 0.2e1 + t488;
t481 = t485 - t793;
t23 = t481 + t858;
t532 = -qJD(1) * t23 + t808;
t32 = t481 - t757;
t531 = qJD(1) * t32 - t808;
t514 = -t723 / 0.2e1 + t722 / 0.2e1;
t113 = t488 - t514;
t457 = t460 * qJ(4);
t528 = qJD(1) * t113 - qJD(3) * t457;
t120 = pkin(5) * t576 + t710;
t527 = qJD(1) * t120 + qJD(3) * t358;
t178 = qJD(1) * t207 + t653;
t205 = t496 * t454 + t825;
t522 = -qJD(4) * t205 + qJD(5) * t103;
t489 = t311 * t773 + t693 / 0.2e1 + t822 * t769;
t47 = t489 + t677;
t511 = -qJD(1) * t47 + t452 * t652;
t150 = -t164 + t674;
t510 = qJD(1) * t150 + t592;
t505 = qJD(5) * t205 - t594;
t504 = -t594 - t838;
t502 = (-qJD(3) * t338 + qJD(5) * t822) * t687;
t270 = t442 + t785;
t492 = qJD(1) * t270 + t653 * t822 + t647;
t483 = t485 + t793;
t436 = t452 * qJD(4);
t417 = t674 * qJD(3);
t399 = t796 * qJ(6);
t301 = t556 - t674;
t298 = t822 * t616;
t249 = t579 - t846;
t225 = t245 * qJD(5);
t224 = t247 * qJD(5);
t220 = qJD(3) * t784 + t452 * t660;
t192 = t617 - t221;
t191 = -t619 - t223;
t153 = t164 + t674;
t147 = t655 + t594;
t130 = (qJD(3) * t342 - t838) * t822;
t127 = qJD(3) * t167 - t647 * t687;
t123 = t687 * t796 - t636;
t119 = t710 / 0.2e1 + t763 / 0.2e1 - t687 * t782 + t822 * t783;
t114 = t488 + t514;
t112 = -t225 - t666;
t109 = -t342 * t660 + t835;
t106 = t494 + t513;
t95 = t452 * t614 - t224 + t666;
t56 = t835 + (t653 + t660) * t342;
t50 = -t489 + t677;
t49 = -t490 + t578;
t44 = t608 + t857;
t37 = t501 - t516;
t33 = t483 + t757;
t30 = t233 + t512 + t517;
t28 = t482 - t561;
t27 = t761 + t818 - t857;
t24 = t483 + t758 - t571;
t15 = t16 * qJD(5);
t12 = -t757 - t858;
t4 = t104 * t772 + t773 * t93 + t775 * t94 + t493 + t520;
t2 = t63 + t93 * t778 + t104 * t777 + t94 * t385 / 0.2e1 - t515 + t519;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t461 * qJD(2), t456 * qJD(2), -t380, -t610 * qJD(3), 0, t380, 0, 0, t469 * t614, -t469 * t432, -qJD(2) * t821, qJD(2) * t200, -t472 * t380, 0.2e1 * t412 * t347, t307 * qJD(3), -t470 * t380, -t305 * qJD(3), t380, -qJD(2) * t306 + qJD(3) * t54 - t454 * t586, -qJD(2) * t352 + qJD(3) * t55 + t349 * t649, -qJD(3) * t46 + qJD(4) * t253, qJD(2) * t74 + qJD(3) * t51 + qJD(4) * t116, t130, -t827, t127, t502, -t829, t380, qJD(3) * t34 + qJD(5) * t53 + t790, qJD(3) * t35 + qJD(5) * t52 + t831, qJD(3) * t14 + t786, qJD(2) * t38 + qJD(3) * t17 + qJD(4) * t42, t130, t127, t827, t380, t829, t502, qJD(3) * t22 + qJD(5) * t41 - t624 * t822 + t790, qJD(3) * t9 + qJD(5) * t10 - t430 * t687 + t786, qJD(3) * t21 + qJD(5) * t40 + qJD(6) * t270 - t831, qJD(2) * t25 + qJD(3) * t7 + qJD(4) * t39 + qJD(5) * t8 + qJD(6) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t611, t613, 0, 0, 0, 0, 0, 0, 0, 0, -t621, t643, 0, 0, 0, 0, 0, 0, -t627, -t659, 0, qJD(3) * t106 + qJD(4) * t301 + t738, 0, 0, 0, 0, 0, 0, t842, -qJD(5) * t246 + t805 - t840, t812, t37 * qJD(3) + t698 + t789, 0, 0, 0, 0, 0, 0, t842, t812, -t224 + t840 + t805, t30 * qJD(3) + t12 * qJD(5) + t249 * qJD(6) + t720 + t789; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t381, -t625, -t432, t381, -t614, 0, -qJD(3) * t498 + t469 * t615, qJD(3) * t384 - t469 * t618, 0, 0 (-t472 * t615 - t590) * t450, 0.2e1 * t550 + (t415 - t416) * qJD(3), t474 * t614 + t626 (-t470 * t615 + t590) * t450, t412 - t628, t381, t681 + (t474 * t553 - t357) * qJD(3) - t347 * qJD(4), t680 + (pkin(3) * t689 - qJ(4) * t686 + t487) * qJD(3) - t586, qJD(3) * t817 - t685, t682 + t106 * qJD(2) + (-pkin(3) * t498 + qJ(4) * t817) * qJD(3) + t114 * qJD(4), t56, -t849, t95, t855, -t227 + t856, t309, t707 + (t312 * t448 - t338 * t468 - t696) * qJD(3) - t799 + t50 * qJD(5), t699 + t798 + (t312 * t452 + t342 * t468 - t692) * qJD(3) - t248 * qJD(4) + t49 * qJD(5), t727 + (-t110 * t452 - t111 * t448 + t572) * qJD(3) + t768, t725 + t37 * qJD(2) + (-t110 * t383 + t111 * t385 + t312 * t468) * qJD(3) + t33 * qJD(4) + t15, t56, t95, t849, t309, qJD(5) * t249 - t856, t855, t721 + (t146 * t448 - t338 * t346 - t696) * qJD(3) - t799 + t27 * qJD(5) + t153 * qJD(6), t753 + (-t100 * t448 + t101 * t452 + t572) * qJD(3) + t70 + t4 * qJD(5) - t247 * qJD(6), t724 + t798 + (-t146 * t452 - t342 * t346 + t692) * qJD(3) - t245 * qJD(4) + t28 * qJD(5) + t298, t755 + t30 * qJD(2) + (t100 * t385 + t101 * t383 + t146 * t346) * qJD(3) + t24 * qJD(4) + t2 * qJD(5) + t44 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t347 - t566 (t474 * t615 - t651) * t450, t634, qJD(2) * t301 + qJD(3) * t114 + t673, 0, 0, 0, 0, 0, 0, t820, -qJD(3) * t248 - t505, t788, qJD(3) * t33 + t742 + t810, 0, 0, 0, 0, 0, 0, t820, t788, t505 - t655, qJD(3) * t24 + qJD(5) * t119 - qJD(6) * t208 + t745 + t810; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t830, t526, t123, -t830, -t640 + t824, -t417, qJD(3) * t50 + t739 + t791, -qJD(2) * t246 + qJD(3) * t49 + t522 + t740, t756, t545, t830, t123, -t526, -t417, t249 * qJD(3) - t824, -t830, qJD(3) * t27 + t743 + t791, t4 * qJD(3) + qJD(5) * t552 - t624 + t735, -qJD(2) * t247 + qJD(3) * t28 + t430 - t522 + t744, t754 + t12 * qJD(2) + t2 * qJD(3) + t119 * qJD(4) + (-pkin(5) * t104 - qJ(6) * t103) * qJD(5) + t93 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t153 - t595, t123, t492, qJD(2) * t249 + qJD(3) * t44 + qJD(5) * t93 - t203 + t741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t611, -t613, 0, 0, 0, 0, 0, 0, t614, -t432, t621, -t643, 0, 0, 0, 0, 0, 0, t412 + t627, -qJD(3) * t349 + t659, t351 * qJD(3), -qJD(3) * t105 + qJD(4) * t302 - t738, 0, 0, 0, 0, 0, 0, t833, -qJD(3) * t687 - qJD(5) * t248 + t840, t843, -qJD(3) * t36 - t698 - t803, 0, 0, 0, 0, 0, 0, t833, t843, qJD(3) * t267 - t225 - t840, qJD(3) * t31 - qJD(5) * t11 - qJD(6) * t240 - t720 - t803; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t615, -t618, 0, 0, 0, 0, 0, 0, 0, 0, t585, -t623, t622, -t646, 0, 0, 0, 0, 0, 0, -t797, -t834, -t814, -t746, 0, 0, 0, 0, 0, 0, -t797, -t814, t630, t747; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t629, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t815, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t619 - t635, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, t191, -t728 - t828; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t381, t625, 0, -t381, 0, 0, -t570 * t454, t570 * t450, 0, 0, t472 * t381, -0.2e1 * t550, -t626, t470 * t381, t628, -t381, -qJD(2) * t686 - t681, qJD(2) * t349 - t680, -qJD(2) * t351 + t685, qJD(2) * t105 - qJD(4) * t113 - t682, t109, -t853, t112, t854, -t226 + t839, -t309, -qJD(5) * t47 - t707 + t787, qJD(2) * t687 - qJD(4) * t246 - qJD(5) * t48 - t699, -t727 + t768 + t809, qJD(2) * t36 - qJD(4) * t32 + t15 - t725, t109, t112, t853, -t309, -qJD(5) * t240 - t839, t854, qJD(5) * t26 - qJD(6) * t150 - t721 + t787, -qJD(5) * t3 - qJD(6) * t245 + t70 - t753 + t809, -qJD(2) * t267 - qJD(4) * t247 + qJD(5) * t29 + t298 - t724, -qJD(2) * t31 - qJD(4) * t23 - qJD(5) * t1 - qJD(6) * t43 - t755; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t615, t618, 0, 0, 0, 0, 0, 0, 0, 0, -t585, t623, -t622, t646, 0, 0, 0, 0, 0, 0, t797, t834, t814, t746, 0, 0, 0, 0, 0, 0, t797, t814, -t630, -t747; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t460 * qJD(4), t457 * qJD(4), -t382, -t850, 0, t382, 0, 0, t468 * t617, -t468 * t619, t851, t802, -t382, 0, t850, 0, 0, t382, qJD(5) * t183 - t448 * t616, t851, qJD(5) * t184 + qJD(6) * t784, qJD(5) * t115 - t346 * t616 + t802; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t612, -t528, 0, 0, 0, 0, 0, 0, -t801, -t638, t848, -t531, 0, 0, 0, 0, 0, 0, -t801, t848, -t637, t532; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t832, -t523, t191, t832, t193, t620, -t360 + t511, -t448 * t652 + t648 - t683, t547, t749, -t832, t191, t523, t620, t192, t832, -t360 + t534, qJD(5) * t551 - t448 * qJD(6) - t548, -t533 - t648 (-pkin(5) * t385 - qJ(6) * t383) * qJD(5) + t385 * qJD(6) + t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t510, t191, t220, -t346 * t653 + t360 - t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t566, -t565, -t634, -qJD(2) * t302 + qJD(3) * t113 - t673, 0, 0, 0, 0, 0, 0, t819, qJD(3) * t246 + t504, -t788, qJD(3) * t32 - t742 + t811, 0, 0, 0, 0, 0, 0, t819, -t788, -t504 + t636, qJD(3) * t23 + qJD(5) * t120 - qJD(6) * t207 - t745 + t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t629, 0, 0, 0, 0, 0, 0, 0, 0, 0, t815, 0, 0, 0, 0, 0, 0, 0, 0, 0, t815; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t612, t528, 0, 0, 0, 0, 0, 0, t844, -t619 + t638, -t848, t531, 0, 0, 0, 0, 0, 0, t844, -t848, t619 + t637, -t532 + t828; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, -t176, 0, 0, 0, 0, 0, 0, 0, 0, t178, 0, t176, t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t830, -t526, t147, t830, t852, -t417, qJD(3) * t47 - t739 + t792, qJD(2) * t248 + qJD(3) * t48 + t650 - t740, -t756, -t545, -t830, t147, t526, -t417, qJD(3) * t240 - t845, t830, -qJD(3) * t26 - t743 + t792, qJD(3) * t3 - t735, qJD(2) * t245 - qJD(3) * t29 + t430 - t650 - t744, qJ(6) * t430 + qJD(2) * t11 + qJD(3) * t1 - qJD(4) * t120 - t754; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t800, t635, 0, 0, 0, 0, 0, 0, 0, 0, t800, 0, t223, t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t832, t523, t223, -t832, t800, -t620, -t436 - t511, t683 + (qJD(4) + t652) * t448, -t547, -t749, t832, t223, -t523, -t620, t221, -t832, -t436 - t534, t548, -qJD(4) * t448 + t533, -qJD(4) * t358 - t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t176, 0, 0, 0, 0, 0, 0, 0, 0, -t178, 0, -t176, -t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t796, -t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t150 + t595, t147, -t492, -qJ(6) * t647 + qJD(2) * t240 + qJD(3) * t43 + t204 - t741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t510, t223, -t220, t691 + (qJD(3) * t346 + qJD(4)) * t452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t796, t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;
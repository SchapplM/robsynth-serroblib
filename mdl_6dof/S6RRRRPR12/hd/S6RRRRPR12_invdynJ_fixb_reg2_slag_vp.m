% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:52
% EndTime: 2019-03-09 23:46:21
% DurationCPUTime: 52.32s
% Computational Cost: add. (70902->1154), mult. (202074->1570), div. (0->0), fcn. (170741->18), ass. (0->483)
t492 = cos(qJ(2));
t724 = cos(pkin(6));
t639 = pkin(1) * t724;
t467 = t492 * t639;
t456 = qJD(1) * t467;
t489 = sin(qJ(2));
t484 = sin(pkin(6));
t723 = cos(pkin(7));
t558 = t484 * (-pkin(10) * t723 - pkin(9));
t540 = t489 * t558;
t351 = qJD(1) * t540 + t456;
t466 = t489 * t639;
t509 = t492 * t558 - t466;
t352 = t509 * qJD(1);
t483 = sin(pkin(7));
t744 = pkin(10) * t483;
t564 = pkin(2) * t489 - t492 * t744;
t661 = qJD(1) * t484;
t391 = t564 * t661;
t747 = cos(qJ(3));
t587 = t723 * t747;
t488 = sin(qJ(3));
t687 = t483 * t488;
t420 = pkin(2) * t587 - pkin(10) * t687;
t621 = t488 * t723;
t754 = t420 * qJD(3) - t747 * t351 - t352 * t621 - t391 * t687;
t524 = -t489 * t621 + t492 * t747;
t387 = t524 * t484;
t374 = qJD(1) * t387;
t629 = qJD(3) * t747;
t600 = t483 * t629;
t559 = t600 - t374;
t258 = -t352 * t483 + t723 * t391;
t522 = t488 * t492 + t489 * t587;
t386 = t522 * t484;
t373 = qJD(1) * t386;
t790 = -pkin(3) * t373 + pkin(11) * t374 - t258 + (pkin(3) * t488 - pkin(11) * t747) * t483 * qJD(3);
t633 = t489 * t661;
t603 = t483 * t633;
t789 = pkin(11) * t603 - t754;
t487 = sin(qJ(4));
t491 = cos(qJ(4));
t415 = t487 * t687 - t491 * t723;
t672 = qJD(4) * t415 + t487 * t603 - t491 * t559;
t416 = t487 * t723 + t491 * t687;
t670 = qJD(4) * t416 + t487 * t559 + t491 * t603;
t659 = qJD(3) * t488;
t630 = t483 * t659;
t755 = t373 - t630;
t638 = t483 * t747;
t422 = pkin(2) * t621 + pkin(10) * t638;
t397 = pkin(11) * t723 + t422;
t398 = (-pkin(3) * t747 - pkin(11) * t488 - pkin(2)) * t483;
t656 = qJD(4) * t491;
t657 = qJD(4) * t487;
t680 = -t397 * t657 + t398 * t656 + t487 * t790 - t789 * t491;
t282 = t491 * t397 + t487 * t398;
t679 = -qJD(4) * t282 + t789 * t487 + t491 * t790;
t614 = t724 * qJD(1);
t571 = t614 + qJD(2);
t529 = t571 * t747;
t567 = t492 * t587;
t554 = t484 * t567;
t602 = t488 * t633;
t314 = -qJD(1) * t554 - t483 * t529 + t602;
t310 = qJD(4) + t314;
t788 = -t755 * pkin(4) + qJ(5) * t672 - qJD(5) * t416 + t679;
t685 = t484 * t492;
t664 = pkin(9) * t685 + t466;
t402 = t664 * qJD(1);
t548 = t483 * t571;
t619 = t492 * t723;
t592 = t484 * t619;
t304 = t402 + (qJD(1) * t592 + t548) * pkin(10);
t514 = pkin(2) * t724 + t540;
t313 = qJD(2) * pkin(2) + qJD(1) * t514 + t456;
t556 = -pkin(2) * t492 - t489 * t744 - pkin(1);
t368 = t556 * t661;
t182 = -t488 * t304 + t313 * t587 + t368 * t638;
t634 = t747 * t489;
t523 = t488 * t619 + t634;
t516 = t523 * t484;
t531 = t488 * t548;
t316 = qJD(1) * t516 + t531;
t220 = pkin(3) * t316 + pkin(11) * t314;
t127 = -t182 * t487 + t491 * t220;
t485 = -qJ(5) - pkin(11);
t623 = qJD(4) * t485;
t787 = pkin(4) * t316 + qJD(5) * t487 + t127 + (qJ(5) * t314 - t623) * t491;
t786 = qJ(5) * t670 + qJD(5) * t415 - t680;
t128 = t491 * t182 + t487 * t220;
t707 = t314 * t487;
t785 = -qJ(5) * t707 + qJD(5) * t491 + t487 * t623 - t128;
t784 = t422 * qJD(3) - t488 * t351 + t352 * t587;
t482 = sin(pkin(13));
t722 = cos(pkin(13));
t678 = t482 * t670 + t672 * t722;
t688 = t482 * t487;
t550 = t491 * t722 - t688;
t677 = t310 * t550;
t676 = -t482 * t672 + t670 * t722;
t636 = t747 * t391;
t669 = -(-pkin(3) * t633 - t636) * t483 + t784;
t616 = t722 * t487;
t430 = t482 * t491 + t616;
t667 = t310 * t430;
t562 = qJD(3) * t587;
t519 = qJD(2) * t747 + t562;
t611 = qJDD(1) * t723;
t652 = qJD(1) * qJD(2);
t580 = t723 * t652;
t561 = t489 * t580;
t607 = t724 * qJDD(1);
t565 = t607 + qJDD(2);
t766 = t488 * t484 * t561 + qJD(3) * t602 - (qJD(3) * t529 + t488 * t565) * t483;
t497 = t766 - (qJDD(1) * t634 + t492 * (qJD(1) * t519 + t488 * t611)) * t484;
t493 = cos(qJ(1));
t618 = t493 * t724;
t746 = sin(qJ(1));
t417 = t489 * t746 - t492 * t618;
t418 = t489 * t618 + t492 * t746;
t684 = t484 * t493;
t270 = -t417 * t621 + t418 * t747 - t684 * t687;
t479 = qJ(4) + pkin(13);
t475 = sin(t479);
t476 = cos(t479);
t622 = t484 * t723;
t767 = -t417 * t483 + t493 * t622;
t205 = t270 * t476 - t475 * t767;
t269 = t417 * t587 + t418 * t488 + t638 * t684;
t486 = sin(qJ(6));
t490 = cos(qJ(6));
t783 = t205 * t486 - t269 * t490;
t782 = t205 * t490 + t269 * t486;
t731 = t482 * t788 - t786 * t722;
t725 = -t482 * t787 + t722 * t785;
t183 = t747 * t304 + (t313 * t723 + t368 * t483) * t488;
t589 = -t183 + (t657 + t707) * pkin(4);
t673 = pkin(4) * t670 + t669;
t702 = t767 * t491;
t710 = t270 * t487;
t780 = t702 + t710;
t703 = t767 * t487;
t779 = -t270 * t491 + t703;
t645 = t483 * t685;
t376 = qJD(1) * t645 - t571 * t723 - qJD(3);
t241 = t316 * t491 - t376 * t487;
t227 = -t313 * t483 + t723 * t368;
t150 = pkin(3) * t314 - pkin(11) * t316 + t227;
t156 = -t376 * pkin(11) + t183;
t94 = t491 * t150 - t156 * t487;
t79 = -qJ(5) * t241 + t94;
t72 = pkin(4) * t310 + t79;
t239 = t316 * t487 + t376 * t491;
t95 = t150 * t487 + t156 * t491;
t80 = -qJ(5) * t239 + t95;
t74 = t722 * t80;
t40 = t482 * t72 + t74;
t38 = pkin(12) * t310 + t40;
t155 = t376 * pkin(3) - t182;
t123 = t239 * pkin(4) + qJD(5) + t155;
t159 = t722 * t239 + t241 * t482;
t551 = -t482 * t239 + t241 * t722;
t63 = t159 * pkin(5) - pkin(12) * t551 + t123;
t19 = t38 * t490 + t486 * t63;
t650 = qJDD(1) * t492;
t625 = t484 * t650;
t512 = -t483 * t625 + t565 * t723 + qJDD(3);
t627 = t489 * t652;
t598 = t484 * t627;
t505 = t483 * t598 + t512;
t658 = qJD(4) * t241;
t122 = -t487 * t497 - t491 * t505 + t658;
t579 = qJD(2) * t614;
t569 = pkin(1) * t579;
t595 = pkin(1) * t607;
t641 = pkin(9) * t625 + t489 * t595 + t492 * t569;
t311 = -pkin(9) * t598 + t641;
t539 = t483 * t565;
t226 = (t539 + (t492 * t611 - t561) * t484) * pkin(10) + t311;
t455 = t492 * t595;
t527 = -t489 * t569 + t455;
t626 = t492 * t652;
t651 = qJDD(1) * t489;
t549 = -t626 - t651;
t528 = t549 * pkin(9);
t230 = t565 * pkin(2) + ((-t489 * t611 - t492 * t580) * pkin(10) + t528) * t484 + t527;
t536 = t564 * qJD(2);
t277 = (qJD(1) * t536 + qJDD(1) * t556) * t484;
t586 = qJD(3) * t621;
t557 = t488 * t226 - t230 * t587 - t277 * t638 + t304 * t629 + t313 * t586 + t368 * t630;
t78 = -pkin(3) * t505 + t557;
t54 = t122 * pkin(4) + qJDD(5) + t78;
t121 = t316 * t657 + t376 * t656 - t487 * t505 + t491 * t497;
t69 = -t121 * t482 + t722 * t122;
t70 = -t121 * t722 - t482 * t122;
t14 = t69 * pkin(5) - t70 * pkin(12) + t54;
t500 = qJD(2) * t522 + qJD(3) * t523;
t199 = qJD(3) * t531 - qJDD(1) * t554 + t484 * (qJD(1) * t500 + t488 * t651) - t747 * t539;
t196 = qJDD(4) + t199;
t89 = t747 * t226 + t230 * t621 + t277 * t687 - t304 * t659 + t313 * t562 + t368 * t600;
t77 = pkin(11) * t505 + t89;
t173 = -t483 * t230 + t723 * t277;
t85 = t199 * pkin(3) + pkin(11) * t497 + t173;
t27 = -qJD(4) * t95 - t487 * t77 + t491 * t85;
t17 = pkin(4) * t196 + qJ(5) * t121 - qJD(5) * t241 + t27;
t555 = -t150 * t656 + t156 * t657 - t487 * t85 - t491 * t77;
t21 = -qJ(5) * t122 - qJD(5) * t239 - t555;
t8 = t482 * t17 + t722 * t21;
t6 = pkin(12) * t196 + t8;
t2 = -qJD(6) * t19 + t490 * t14 - t486 * t6;
t764 = qJD(6) + t159;
t778 = t19 * t764 + t2;
t575 = t38 * t486 - t490 * t63;
t1 = -t575 * qJD(6) + t486 * t14 + t490 * t6;
t777 = t575 * t764 + t1;
t776 = pkin(12) * t755 - t731;
t775 = t676 * pkin(5) + t678 * pkin(12) + t673;
t774 = -pkin(12) * t316 + t725;
t773 = t667 * pkin(5) - pkin(12) * t677 + t589;
t772 = t159 * t551;
t588 = t724 * t746;
t525 = t493 * t489 + t492 * t588;
t357 = t525 * t483 + t746 * t622;
t768 = -t94 * t310 - t555;
t608 = t764 * t490;
t68 = qJDD(6) + t69;
t765 = -t486 * t68 - t608 * t764;
t762 = -t270 * t475 - t476 * t767;
t478 = t484 ^ 2;
t761 = 0.2e1 * t478;
t738 = t786 * t482 + t722 * t788;
t760 = -t95 * t310 - t27;
t726 = t482 * t785 + t722 * t787;
t134 = t310 * t486 + t490 * t551;
t609 = t764 * t486;
t759 = t134 * t609;
t419 = -t489 * t588 + t492 * t493;
t637 = t484 * t746;
t752 = -t483 * t637 + t525 * t723;
t273 = t419 * t488 + t747 * t752;
t617 = t724 * t483;
t568 = t747 * t617;
t683 = t488 * t489;
t345 = t484 * t683 - t554 - t568;
t544 = g(1) * t273 + g(2) * t269 + g(3) * t345;
t756 = t544 * t475;
t350 = t467 + t514;
t686 = t484 * t489;
t646 = t483 * t686;
t665 = pkin(2) * t685 + pkin(10) * t646;
t380 = -pkin(1) * t484 - t665;
t248 = -t350 * t483 + t723 * t380;
t593 = t488 * t617;
t346 = t593 + t516;
t178 = pkin(3) * t345 - pkin(11) * t346 + t248;
t339 = (t592 + t617) * pkin(10) + t664;
t201 = t747 * t339 + t350 * t621 + t380 * t687;
t412 = -t723 * t724 + t645;
t186 = -pkin(11) * t412 + t201;
t114 = t487 * t178 + t491 * t186;
t668 = t483 * t636 + t784;
t751 = (qJDD(2) + 0.2e1 * t607) * t484;
t581 = g(1) * t419 + g(2) * t418;
t406 = t664 * qJD(2);
t494 = qJD(1) ^ 2;
t745 = pkin(4) * t487;
t256 = qJD(3) * t568 + ((t567 - t683) * qJD(3) + t524 * qJD(2)) * t484;
t699 = t412 * t491;
t704 = t346 * t487;
t266 = t699 + t704;
t660 = qJD(2) * t489;
t631 = t484 * t660;
t601 = t483 * t631;
t172 = -qJD(4) * t266 + t256 * t491 + t487 * t601;
t255 = qJD(3) * t593 + t484 * t500;
t267 = t346 * t491 - t412 * t487;
t457 = qJD(2) * t467;
t353 = qJD(2) * t540 + t457;
t354 = t509 * qJD(2);
t392 = t484 * t536;
t129 = -t339 * t659 + t350 * t562 + t747 * t353 + t354 * t621 + t380 * t600 + t392 * t687;
t125 = pkin(11) * t601 + t129;
t259 = -t354 * t483 + t723 * t392;
t138 = pkin(3) * t255 - pkin(11) * t256 + t259;
t53 = -qJD(4) * t114 - t125 * t487 + t491 * t138;
t33 = pkin(4) * t255 - qJ(5) * t172 - qJD(5) * t267 + t53;
t171 = qJD(4) * t267 + t256 * t487 - t491 * t601;
t52 = t491 * t125 + t487 * t138 + t178 * t656 - t186 * t657;
t36 = -qJ(5) * t171 - qJD(5) * t266 + t52;
t12 = t482 * t33 + t722 * t36;
t281 = -t487 * t397 + t491 * t398;
t235 = -pkin(4) * t638 - t416 * qJ(5) + t281;
t245 = -qJ(5) * t415 + t282;
t164 = t482 * t235 + t722 * t245;
t152 = -pkin(12) * t638 + t164;
t306 = t415 * t722 + t416 * t482;
t307 = -t482 * t415 + t416 * t722;
t396 = -pkin(3) * t723 - t420;
t321 = t415 * pkin(4) + t396;
t184 = t306 * pkin(5) - t307 * pkin(12) + t321;
t105 = -t152 * t486 + t184 * t490;
t741 = qJD(6) * t105 + t775 * t486 - t490 * t776;
t106 = t152 * t490 + t184 * t486;
t740 = -qJD(6) * t106 + t486 * t776 + t775 * t490;
t113 = t491 * t178 - t186 * t487;
t88 = pkin(4) * t345 - qJ(5) * t267 + t113;
t93 = -qJ(5) * t266 + t114;
t51 = t482 * t88 + t722 * t93;
t739 = pkin(5) * t755 - t738;
t654 = qJD(6) * t490;
t655 = qJD(6) * t486;
t46 = -t486 * t196 - t310 * t654 - t490 * t70 + t551 * t655;
t737 = t46 * t486;
t615 = -t490 * t196 + t486 * t70;
t47 = qJD(6) * t134 + t615;
t736 = t47 * t490;
t735 = t482 * t80;
t132 = -t490 * t310 + t486 * t551;
t730 = -t132 * t654 - t486 * t47;
t474 = pkin(4) * t491 + pkin(3);
t327 = -pkin(5) * t550 - pkin(12) * t430 - t474;
t450 = t485 * t491;
t367 = -t450 * t722 + t485 * t688;
t236 = t327 * t490 - t367 * t486;
t729 = qJD(6) * t236 + t486 * t773 + t490 * t774;
t237 = t327 * t486 + t367 * t490;
t728 = -qJD(6) * t237 - t486 * t774 + t490 * t773;
t727 = t316 * pkin(5) + t726;
t721 = t132 * t551;
t720 = t132 * t159;
t719 = t134 * t132;
t718 = t134 * t551;
t717 = t551 ^ 2;
t716 = t551 * t310;
t715 = t159 ^ 2;
t714 = t159 * t310;
t713 = t239 * t310;
t712 = t241 * t239;
t711 = t241 * t310;
t274 = t419 * t747 - t488 * t752;
t709 = t274 * t487;
t708 = t310 * t316;
t706 = t316 * t314;
t705 = t345 * t486;
t701 = t357 * t487;
t700 = t357 * t491;
t697 = t418 * t483;
t696 = t419 * t483;
t695 = t430 * t486;
t694 = t430 * t490;
t693 = t475 * t483;
t692 = t476 * t483;
t691 = t476 * t486;
t690 = t476 * t490;
t689 = t478 * t494;
t275 = t486 * t307 + t490 * t638;
t682 = qJD(6) * t275 + t486 * t755 + t490 * t678;
t605 = t486 * t638;
t681 = -qJD(6) * t605 + t307 * t654 - t486 * t678 + t490 * t755;
t675 = -t269 * t474 - t270 * t485;
t674 = -t273 * t474 - t274 * t485;
t671 = -t345 * t474 - t346 * t485;
t663 = t493 * pkin(1) + pkin(9) * t637;
t480 = t489 ^ 2;
t481 = t492 ^ 2;
t662 = t480 - t481;
t647 = t492 * t689;
t300 = -t417 * t488 + t418 * t587;
t301 = -t417 * t747 - t418 * t621;
t407 = t417 * pkin(2);
t643 = -t300 * t485 + t301 * t474 - t407;
t302 = t419 * t587 - t488 * t525;
t303 = -t419 * t621 - t525 * t747;
t409 = t525 * pkin(2);
t642 = -t302 * t485 + t303 * t474 - t409;
t635 = t747 * t392;
t632 = t483 * t660;
t628 = pkin(1) * t761;
t7 = t722 * t17 - t482 * t21;
t613 = -t490 * t316 - t486 * t677;
t612 = t316 * t486 - t490 * t677;
t610 = t491 * t310;
t606 = t489 * t647;
t599 = t489 * t626;
t597 = pkin(10) * t697 - t407;
t596 = pkin(10) * t696 - t409;
t590 = -pkin(1) * t746 + pkin(9) * t684;
t585 = t484 * t494 * t724;
t584 = -pkin(5) * t476 - pkin(12) * t475;
t208 = t274 * t475 - t357 * t476;
t583 = -g(1) * t762 - g(2) * t208;
t582 = -g(1) * t269 + g(2) * t273;
t49 = pkin(12) * t345 + t51;
t200 = -t488 * t339 + t350 * t587 + t380 * t638;
t185 = t412 * pkin(3) - t200;
t139 = t266 * pkin(4) + t185;
t187 = t266 * t722 + t267 * t482;
t188 = -t482 * t266 + t267 * t722;
t73 = t187 * pkin(5) - t188 * pkin(12) + t139;
t25 = t486 * t73 + t49 * t490;
t24 = -t486 * t49 + t490 * t73;
t574 = -t386 * t485 + t387 * t474 + t646 * t745 + t665;
t145 = t188 * t490 + t705;
t570 = 0.2e1 * t614 + qJD(2);
t563 = t490 * t68 + (-t159 * t486 - t655) * t764;
t560 = g(1) * t493 + g(2) * t746;
t11 = t33 * t722 - t482 * t36;
t39 = t72 * t722 - t735;
t50 = -t482 * t93 + t722 * t88;
t553 = -pkin(11) * t196 + t155 * t310;
t37 = -t310 * pkin(5) - t39;
t471 = pkin(4) * t482 + pkin(12);
t552 = t37 * t764 - t471 * t68;
t163 = t235 * t722 - t482 * t245;
t253 = -t346 * t475 - t412 * t476;
t547 = g(1) * t208 - g(2) * t762 - g(3) * t253;
t209 = t274 * t476 + t357 * t475;
t254 = t346 * t476 - t412 * t475;
t546 = -g(1) * t209 - g(2) * t205 - g(3) * t254;
t231 = t301 * t475 - t418 * t692;
t233 = t303 * t475 - t419 * t692;
t294 = t387 * t475 - t476 * t646;
t545 = -g(1) * t233 - g(2) * t231 - g(3) * t294;
t543 = g(1) * t274 + g(2) * t270 + g(3) * t346;
t542 = -g(1) * t302 - g(2) * t300 - g(3) * t386;
t541 = g(1) * t303 + g(2) * t301 + g(3) * t387;
t535 = t430 * t654 - t613;
t534 = -t430 * t655 - t612;
t532 = t544 - t78;
t530 = -g(3) * t686 - t581;
t520 = -t418 * pkin(2) + pkin(10) * t767 + t590;
t518 = -t339 * t629 - t350 * t586 - t488 * t353 + t354 * t587 - t380 * t630;
t515 = pkin(11) * qJD(4) * t310 - t532;
t5 = -pkin(5) * t196 - t7;
t513 = qJD(6) * t471 * t764 + t5 - t547;
t508 = pkin(4) * t703 + t269 * t485 - t270 * t474 + t520;
t507 = t581 * (pkin(10) + t745) * t483;
t506 = t419 * pkin(2) + pkin(10) * t357 + t663;
t503 = pkin(4) * t701 - t273 * t485 + t274 * t474 + t506;
t502 = t483 * t505;
t126 = (-pkin(3) * t631 - t635) * t483 - t518;
t84 = t171 * pkin(4) + t126;
t496 = t483 * t497;
t472 = -pkin(4) * t722 - pkin(5);
t421 = -pkin(9) * t686 + t467;
t405 = -pkin(9) * t631 + t457;
t400 = -pkin(9) * t633 + t456;
t393 = pkin(4) * t699;
t366 = -t450 * t482 - t485 * t616;
t338 = pkin(4) * t700;
t336 = pkin(4) * t702;
t334 = t345 * t490;
t312 = t484 * t528 + t527;
t295 = t387 * t476 + t475 * t646;
t276 = t490 * t307 - t605;
t234 = t303 * t476 + t419 * t693;
t232 = t301 * t476 + t418 * t693;
t215 = t274 * t491 + t701;
t214 = t700 - t709;
t151 = pkin(5) * t638 - t163;
t147 = t209 * t490 + t273 * t486;
t146 = -t209 * t486 + t273 * t490;
t144 = t188 * t486 - t334;
t130 = t483 * t635 + t518;
t117 = -t196 * t638 - t310 * t755;
t116 = t196 * t345 + t255 * t310;
t109 = -t482 * t171 + t172 * t722;
t108 = t171 * t722 + t172 * t482;
t81 = pkin(4) * t241 + pkin(5) * t551 + pkin(12) * t159;
t65 = qJD(6) * t145 + t109 * t486 - t255 * t490;
t64 = -t490 * t109 + t188 * t655 - t255 * t486 - t345 * t654;
t48 = -t345 * pkin(5) - t50;
t44 = t722 * t79 - t735;
t43 = t482 * t79 + t74;
t34 = t108 * pkin(5) - t109 * pkin(12) + t84;
t23 = t44 * t490 + t486 * t81;
t22 = -t44 * t486 + t490 * t81;
t10 = pkin(12) * t255 + t12;
t9 = -t255 * pkin(5) - t11;
t4 = -qJD(6) * t25 - t10 * t486 + t34 * t490;
t3 = qJD(6) * t24 + t10 * t490 + t34 * t486;
t13 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t746 - g(2) * t493, t560, 0, 0 (qJDD(1) * t480 + 0.2e1 * t599) * t478 (t489 * t650 - t652 * t662) * t761, qJD(2) * t570 * t685 + t489 * t751 (qJDD(1) * t481 - 0.2e1 * t599) * t478, t492 * t751 - t570 * t631, t565 * t724, -t406 * t571 + t421 * t565 + t312 * t724 + g(1) * t418 - g(2) * t419 + (-t627 + t650) * t628, -g(1) * t417 + g(2) * t525 - t311 * t724 - t405 * t571 + t549 * t628 - t565 * t664 ((-t400 * qJD(2) + qJDD(1) * t664 + t311 + (-qJD(2) * t421 + t405) * qJD(1)) * t492 + (-t402 * qJD(2) - qJDD(1) * t421 - t312) * t489 - t560) * t484, t478 * qJDD(1) * pkin(1) ^ 2 - g(1) * t590 - g(2) * t663 + t311 * t664 + t312 * t421 - t400 * t406 + t402 * t405, t316 * t256 - t346 * t497, -t346 * t199 - t316 * t255 - t256 * t314 + t345 * t497, -t256 * t376 + t346 * t512 + t766 * t412 + (t316 * t632 - t523 * t412 * qJDD(1) + (-t412 * t492 * t519 + t346 * t632) * qJD(1)) * t484, t199 * t345 + t255 * t314, t199 * t412 + t255 * t376 - t314 * t601 - t345 * t505, -t376 * t601 - t412 * t505, g(1) * t270 - g(2) * t274 - t130 * t376 + t173 * t345 + t182 * t601 + t248 * t199 + t200 * t505 + t227 * t255 + t259 * t314 + t412 * t557, t129 * t376 + t173 * t346 - t183 * t601 - t201 * t505 + t227 * t256 - t248 * t497 + t259 * t316 + t89 * t412 + t582, -g(1) * t767 - g(2) * t357 - t129 * t314 - t130 * t316 - t182 * t256 - t183 * t255 - t201 * t199 + t200 * t497 - t89 * t345 + t346 * t557, -g(1) * t520 - g(2) * t506 + t183 * t129 + t182 * t130 + t173 * t248 - t200 * t557 + t89 * t201 + t227 * t259, -t121 * t267 + t172 * t241, t121 * t266 - t122 * t267 - t171 * t241 - t172 * t239, -t121 * t345 + t172 * t310 + t196 * t267 + t241 * t255, t122 * t266 + t171 * t239, -t122 * t345 - t171 * t310 - t196 * t266 - t239 * t255, t116, -g(1) * t779 - g(2) * t215 + t113 * t196 + t185 * t122 + t126 * t239 + t155 * t171 + t94 * t255 + t78 * t266 + t27 * t345 + t53 * t310, -g(1) * t780 - g(2) * t214 - t114 * t196 - t185 * t121 + t126 * t241 + t155 * t172 - t95 * t255 + t78 * t267 - t52 * t310 + t555 * t345, t113 * t121 - t114 * t122 - t171 * t95 - t172 * t94 - t239 * t52 - t241 * t53 + t266 * t555 - t267 * t27 - t582, -t555 * t114 + t95 * t52 + t27 * t113 + t94 * t53 + t78 * t185 + t155 * t126 - g(1) * (-pkin(3) * t270 - pkin(11) * t269 + t520) - g(2) * (t274 * pkin(3) + t273 * pkin(11) + t506) t109 * t551 + t188 * t70, -t108 * t551 - t109 * t159 - t187 * t70 - t188 * t69, t109 * t310 + t188 * t196 + t255 * t551 + t345 * t70, t108 * t159 + t187 * t69, -t108 * t310 - t159 * t255 - t187 * t196 - t345 * t69, t116, g(1) * t205 - g(2) * t209 + t108 * t123 + t11 * t310 + t139 * t69 + t159 * t84 + t187 * t54 + t196 * t50 + t255 * t39 + t345 * t7, t109 * t123 - t12 * t310 + t139 * t70 + t188 * t54 - t196 * t51 - t255 * t40 - t345 * t8 + t551 * t84 - t583, -t108 * t40 - t109 * t39 - t11 * t551 - t12 * t159 - t187 * t8 - t188 * t7 - t50 * t70 - t51 * t69 - t582, -g(1) * t508 - g(2) * t503 + t39 * t11 + t40 * t12 + t123 * t84 + t54 * t139 + t7 * t50 + t8 * t51, -t134 * t64 - t145 * t46, t132 * t64 - t134 * t65 + t144 * t46 - t145 * t47, t108 * t134 + t145 * t68 - t187 * t46 - t64 * t764, t132 * t65 + t144 * t47, -t108 * t132 - t144 * t68 - t187 * t47 - t65 * t764, t108 * t764 + t187 * t68, g(1) * t782 - g(2) * t147 - t575 * t108 + t9 * t132 + t5 * t144 + t2 * t187 + t24 * t68 + t37 * t65 + t4 * t764 + t48 * t47, -g(1) * t783 - g(2) * t146 - t1 * t187 - t19 * t108 + t9 * t134 + t5 * t145 - t25 * t68 - t3 * t764 - t37 * t64 - t48 * t46, -t1 * t144 - t132 * t3 - t134 * t4 - t145 * t2 - t19 * t65 + t24 * t46 - t25 * t47 - t575 * t64 + t583, t1 * t25 + t19 * t3 + t2 * t24 - t575 * t4 + t5 * t48 + t37 * t9 - g(1) * (-pkin(5) * t205 + pkin(12) * t762 + t508) - g(2) * (t209 * pkin(5) + t208 * pkin(12) + t503); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t606, t662 * t689, t484 * t651 - t492 * t585, t606, t489 * t585 + t625, t565, t455 + t402 * t571 + g(1) * t525 + g(2) * t417 + (-t579 + t689) * t489 * pkin(1) + (-g(3) * t492 + t528) * t484, pkin(1) * t647 + t400 * t571 + (pkin(9) * t652 + g(3)) * t686 + t581 - t641, 0, 0, t316 * t559 - t488 * t496, -t199 * t687 - t314 * t559 + t316 * t755 - t496 * t747, -t316 * t603 - t376 * t559 + t488 * t502 - t497 * t723, -t199 * t638 - t314 * t755, -t199 * t723 + t314 * t603 - t376 * t755 + t502 * t747, t376 * t603 + t505 * t723, t420 * t512 - t557 * t723 - t258 * t314 - t227 * t373 + t668 * t376 + (t227 * t659 - t173 * t747 - pkin(2) * t199 + (qJD(2) * t420 - t182) * t633) * t483 - t541, pkin(2) * t496 + t173 * t687 + t183 * t603 + t227 * t559 - t258 * t316 + t376 * t754 - t422 * t505 - t723 * t89 - t542, -g(1) * t696 - g(2) * t697 - g(3) * t646 - t182 * t559 + t183 * t755 - t422 * t199 - t314 * t754 + t316 * t668 + t420 * t497 + t557 * t687 + t638 * t89, -t173 * t483 * pkin(2) - g(1) * t596 - g(2) * t597 - g(3) * t665 - t668 * t182 + t183 * t754 - t227 * t258 - t420 * t557 + t89 * t422, -t121 * t416 - t241 * t672, t121 * t415 - t122 * t416 + t239 * t672 - t241 * t670, t121 * t638 + t416 * t196 - t241 * t755 - t310 * t672, t122 * t415 + t239 * t670, t122 * t638 - t415 * t196 + t239 * t755 - t310 * t670, t117, t396 * t122 + t281 * t196 - t94 * t373 + t78 * t415 - t541 * t491 + (-t27 * t747 + t487 * t530 + t659 * t94) * t483 + t679 * t310 + t669 * t239 + t670 * t155, -t396 * t121 - t282 * t196 + t95 * t373 + t78 * t416 + t541 * t487 + (t491 * t530 - t555 * t747 - t659 * t95) * t483 - t680 * t310 + t669 * t241 - t672 * t155, t121 * t281 - t122 * t282 - t239 * t680 - t241 * t679 - t27 * t416 + t415 * t555 - t670 * t95 + t672 * t94 + t542, -t555 * t282 + t27 * t281 + t78 * t396 - g(1) * (pkin(3) * t303 + pkin(11) * t302 + t596) - g(2) * (pkin(3) * t301 + pkin(11) * t300 + t597) - g(3) * (pkin(3) * t387 + pkin(11) * t386 + t665) + t680 * t95 + t679 * t94 + t669 * t155, t307 * t70 - t551 * t678, t159 * t678 - t306 * t70 - t307 * t69 - t551 * t676, t307 * t196 - t310 * t678 - t551 * t755 - t638 * t70, t159 * t676 + t306 * t69, t159 * t755 - t306 * t196 - t310 * t676 + t638 * t69, t117, -g(1) * t234 - g(2) * t232 - g(3) * t295 + t123 * t676 + t159 * t673 + t163 * t196 + t54 * t306 + t310 * t738 + t321 * t69 - t39 * t755 - t638 * t7, -t123 * t678 - t164 * t196 + t54 * t307 - t310 * t731 + t321 * t70 + t40 * t755 + t551 * t673 + t638 * t8 - t545, -t159 * t731 - t163 * t70 - t164 * t69 - t306 * t8 - t307 * t7 + t39 * t678 - t40 * t676 - t551 * t738 + t542, -g(1) * t642 - g(2) * t643 - g(3) * t574 + t123 * t673 + t7 * t163 + t8 * t164 + t54 * t321 + t39 * t738 + t40 * t731 - t507, -t134 * t682 - t276 * t46, t132 * t682 - t134 * t681 + t275 * t46 - t276 * t47, t134 * t676 + t276 * t68 - t306 * t46 - t682 * t764, t132 * t681 + t275 * t47, -t132 * t676 - t275 * t68 - t306 * t47 - t681 * t764, t306 * t68 + t676 * t764, t105 * t68 + t2 * t306 + t151 * t47 + t5 * t275 - g(1) * (t234 * t490 + t302 * t486) - g(2) * (t232 * t490 + t300 * t486) - g(3) * (t295 * t490 + t386 * t486) + t681 * t37 - t676 * t575 + t740 * t764 + t739 * t132, -t106 * t68 - t1 * t306 - t151 * t46 + t5 * t276 - g(1) * (-t234 * t486 + t302 * t490) - g(2) * (-t232 * t486 + t300 * t490) - g(3) * (-t295 * t486 + t386 * t490) - t682 * t37 - t676 * t19 - t741 * t764 + t739 * t134, -t1 * t275 + t105 * t46 - t106 * t47 - t132 * t741 - t134 * t740 - t19 * t681 - t2 * t276 - t575 * t682 + t545, t1 * t106 + t2 * t105 + t5 * t151 - g(1) * (pkin(5) * t234 + pkin(12) * t233 + t642) - g(2) * (pkin(5) * t232 + pkin(12) * t231 + t643) - g(3) * (pkin(5) * t295 + pkin(12) * t294 + t574) - t507 + t739 * t37 + t741 * t19 - t740 * t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t706, -t314 ^ 2 + t316 ^ 2, -t314 * t376 - t497, -t706, -t316 * t376 - t199, t505, -t183 * t376 - t227 * t316 + t544 - t557, -t182 * t376 + t227 * t314 + t543 - t89, 0, 0, -t121 * t487 + t241 * t610 (-t121 - t713) * t491 + (-t122 - t711) * t487, t196 * t487 - t241 * t316 + t310 * t610, -t122 * t491 + t487 * t713, -t310 ^ 2 * t487 + t196 * t491 + t239 * t316, -t708, -pkin(3) * t122 - t127 * t310 - t183 * t239 - t316 * t94 + t487 * t553 - t491 * t515, pkin(3) * t121 + t128 * t310 - t183 * t241 + t316 * t95 + t487 * t515 + t491 * t553, t127 * t241 + t128 * t239 + ((-t122 + t658) * pkin(11) + t768) * t491 + ((qJD(4) * t239 - t121) * pkin(11) + t760) * t487 - t543, -t94 * t127 - t95 * t128 - t155 * t183 + t532 * pkin(3) + (-t555 * t491 - t27 * t487 + (-t487 * t95 - t491 * t94) * qJD(4) - t543) * pkin(11), t430 * t70 + t551 * t677, -t159 * t677 - t430 * t69 + t550 * t70 - t551 * t667, t196 * t430 + t310 * t677 - t316 * t551, t159 * t667 - t550 * t69, t159 * t316 + t196 * t550 - t310 * t667, -t708, t123 * t667 + t159 * t589 - t196 * t366 - t310 * t726 - t316 * t39 - t474 * t69 + t476 * t544 - t54 * t550, t123 * t677 - t196 * t367 - t310 * t725 + t316 * t40 + t430 * t54 - t474 * t70 + t551 * t589 - t756, -t159 * t725 + t366 * t70 - t367 * t69 - t39 * t677 - t40 * t667 - t430 * t7 + t550 * t8 + t551 * t726 - t543, -g(1) * t674 - g(2) * t675 - g(3) * t671 + t123 * t589 - t7 * t366 + t8 * t367 - t39 * t726 + t40 * t725 - t54 * t474, t134 * t534 - t46 * t694, t613 * t134 + t612 * t132 + (t737 - t736 + (t132 * t486 - t134 * t490) * qJD(6)) * t430, t134 * t667 + t46 * t550 + t534 * t764 + t68 * t694, t132 * t535 + t47 * t695, -t132 * t667 + t47 * t550 - t535 * t764 - t68 * t695, -t550 * t68 + t667 * t764, t236 * t68 - t2 * t550 + t366 * t47 + t5 * t695 - g(1) * (-t273 * t690 + t274 * t486) - g(2) * (-t269 * t690 + t270 * t486) - g(3) * (-t345 * t690 + t346 * t486) - t667 * t575 + t728 * t764 + t727 * t132 + t535 * t37, -t237 * t68 + t1 * t550 - t366 * t46 + t5 * t694 - g(1) * (t273 * t691 + t274 * t490) - g(2) * (t269 * t691 + t270 * t490) - g(3) * (t345 * t691 + t346 * t490) - t667 * t19 - t729 * t764 + t727 * t134 + t534 * t37, t236 * t46 - t237 * t47 + t613 * t19 - t612 * t575 - t728 * t134 - t729 * t132 + t756 + (-t1 * t486 - t2 * t490 + (-t19 * t490 - t486 * t575) * qJD(6)) * t430, t1 * t237 + t2 * t236 + t5 * t366 - g(1) * (t273 * t584 + t674) - g(2) * (t269 * t584 + t675) - g(3) * (t345 * t584 + t671) + t727 * t37 + t729 * t19 - t728 * t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t712, -t239 ^ 2 + t241 ^ 2, -t121 + t713, -t712, t711 - t122, t196, -g(1) * t214 + g(2) * t780 + g(3) * t266 - t155 * t241 - t760, g(1) * t215 - g(2) * t779 + g(3) * t267 + t155 * t239 - t768, 0, 0, t772, -t715 + t717, t70 + t714, -t772, -t69 + t716, t196, -t123 * t551 + t43 * t310 + (-t159 * t241 + t196 * t722) * pkin(4) + t547 + t7, t123 * t159 + t310 * t44 + (-t196 * t482 - t241 * t551) * pkin(4) - t546 - t8 (-t482 * t69 - t70 * t722) * pkin(4) + (t40 - t43) * t551 + (t44 - t39) * t159, -g(1) * t338 + g(2) * t336 + g(3) * t393 + t39 * t43 - t40 * t44 + (-t123 * t241 + t8 * t482 + t487 * t543 + t7 * t722) * pkin(4), t134 * t608 - t737 (-t46 - t720) * t490 - t759 + t730, -t718 - t765, t132 * t609 - t736, t563 + t721, -t764 * t551, -t132 * t43 - t22 * t764 + t47 * t472 + t486 * t552 - t490 * t513 + t551 * t575, -t134 * t43 + t19 * t551 + t23 * t764 - t46 * t472 + t486 * t513 + t490 * t552, t132 * t23 + t134 * t22 + (t159 * t575 - t47 * t471 + t1 + (t134 * t471 + t575) * qJD(6)) * t490 + (-t159 * t19 - t46 * t471 - t2 + (t132 * t471 - t19) * qJD(6)) * t486 + t546, t5 * t472 - t19 * t23 + t575 * t22 - t37 * t43 - g(1) * (-pkin(4) * t709 - pkin(5) * t208 + pkin(12) * t209 + t338) - g(2) * (-pkin(4) * t710 + pkin(5) * t762 + pkin(12) * t205 - t336) - g(3) * (-pkin(4) * t704 + pkin(5) * t253 + pkin(12) * t254 - t393) + (t1 * t490 - t2 * t486 + (-t19 * t486 + t490 * t575) * qJD(6)) * t471; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 + t716, t70 - t714, -t715 - t717, t159 * t40 + t39 * t551 + t54 - t544, 0, 0, 0, 0, 0, 0, t563 - t721, -t718 + t765 (t46 - t720) * t490 + t759 + t730, -t551 * t37 + t486 * t777 + t490 * t778 - t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t719, -t132 ^ 2 + t134 ^ 2, t132 * t764 - t46, -t719, -t615 + (-qJD(6) + t764) * t134, t68, -t37 * t134 - g(1) * t146 + g(2) * t783 - g(3) * (-t254 * t486 + t334) + t778, t37 * t132 + g(1) * t147 + g(2) * t782 - g(3) * (-t254 * t490 - t705) - t777, 0, 0;];
tau_reg  = t13;

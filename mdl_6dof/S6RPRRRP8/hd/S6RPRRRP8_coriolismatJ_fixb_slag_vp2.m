% Calculate matrix of centrifugal and coriolis load on the joints for
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
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:23:01
% EndTime: 2019-03-09 06:23:23
% DurationCPUTime: 13.20s
% Computational Cost: add. (18030->609), mult. (33293->771), div. (0->0), fcn. (33746->6), ass. (0->354)
t768 = Ifges(6,1) + Ifges(7,1);
t390 = sin(qJ(5));
t393 = cos(qJ(5));
t528 = t393 * pkin(5) + t390 * qJ(6);
t357 = -pkin(4) - t528;
t671 = mrSges(7,3) * t393;
t673 = mrSges(7,1) * t390;
t539 = -t671 + t673;
t293 = t357 * t539;
t454 = -t293 / 0.2e1;
t650 = t393 * mrSges(7,1);
t653 = t390 * mrSges(7,3);
t359 = -t650 - t653;
t643 = qJ(6) * t393;
t527 = pkin(5) * t390 - t643;
t294 = t527 * t359;
t666 = Ifges(7,5) * t393;
t530 = Ifges(7,3) * t390 + t666;
t668 = Ifges(6,4) * t393;
t534 = -Ifges(6,2) * t390 + t668;
t687 = -t393 / 0.2e1;
t689 = -t390 / 0.2e1;
t667 = Ifges(7,5) * t390;
t529 = -Ifges(7,3) * t393 + t667;
t669 = Ifges(6,4) * t390;
t533 = Ifges(6,2) * t393 + t669;
t688 = t390 / 0.2e1;
t765 = t390 * t768 - t666 + t668;
t713 = t529 * t689 + t533 * t688 + t687 * t765;
t759 = t393 / 0.2e1;
t536 = Ifges(7,1) * t393 + t667;
t538 = Ifges(6,1) * t393 - t669;
t761 = t536 + t538;
t767 = -t530 * t759 - t534 * t687 - t689 * t761 + t294 - t713;
t386 = t390 ^ 2;
t387 = t393 ^ 2;
t736 = t386 + t387;
t728 = Ifges(7,4) + Ifges(6,5);
t685 = cos(qJ(4));
t595 = t685 * pkin(3);
t331 = -t595 + t357;
t375 = -t595 - pkin(4);
t672 = mrSges(6,2) * t393;
t674 = mrSges(6,1) * t390;
t540 = t672 + t674;
t497 = t375 * t540;
t498 = t357 * t527;
t694 = t331 / 0.2e1;
t700 = m(7) / 0.2e1;
t750 = -t540 / 0.2e1;
t766 = t454 - (t331 * t527 + t498) * t700 - pkin(4) * t750 - t539 * t694 - t497 / 0.2e1 - t767;
t391 = sin(qJ(4));
t392 = sin(qJ(3));
t686 = cos(qJ(3));
t351 = t391 * t392 - t685 * t686;
t352 = t391 * t686 + t685 * t392;
t172 = -Ifges(7,6) * t351 - t352 * t530;
t173 = -Ifges(6,6) * t351 - t352 * t534;
t501 = t351 * (Ifges(6,5) * t390 + Ifges(6,6) * t393);
t693 = -t351 / 0.2e1;
t174 = -Ifges(7,4) * t351 - t352 * t536;
t175 = -Ifges(6,5) * t351 - t352 * t538;
t738 = t175 + t174;
t414 = t172 * t687 + t173 * t759 - t501 / 0.2e1 + (Ifges(7,4) * t390 - Ifges(7,6) * t393) * t693 + Ifges(5,6) * t351 + t738 * t688 + (-Ifges(5,5) + t713) * t352;
t651 = t393 * mrSges(6,1);
t654 = t390 * mrSges(6,2);
t360 = -t651 + t654;
t394 = -pkin(1) - pkin(7);
t355 = (-pkin(8) + t394) * t392;
t570 = t686 * t394;
t356 = -t686 * pkin(8) + t570;
t726 = t685 * t355 + t391 * t356;
t739 = t726 * t360;
t744 = t726 * mrSges(5,1);
t253 = t355 * t391 - t685 * t356;
t745 = t253 * mrSges(5,2);
t620 = t352 * t390;
t478 = -pkin(5) * t620 + t726;
t619 = t352 * t393;
t758 = qJ(6) * t619 + t478;
t762 = t758 * t359;
t764 = t414 + t739 + t745 - t744 + t762;
t763 = t745 / 0.2e1 + t762 / 0.2e1;
t760 = t736 * t685;
t746 = mrSges(7,2) + mrSges(6,3);
t723 = t728 * t393 + (-Ifges(6,6) + Ifges(7,6)) * t390;
t740 = t723 * t352;
t683 = pkin(3) * t391;
t374 = pkin(9) + t683;
t608 = t736 * t374 * t351;
t623 = t351 * t390;
t655 = t352 * mrSges(7,3);
t509 = mrSges(7,2) * t623 + t655;
t594 = mrSges(6,3) * t623;
t510 = -mrSges(6,2) * t352 + t594;
t757 = t509 + t510;
t372 = t392 * pkin(3) + qJ(2);
t606 = -t352 * mrSges(5,1) + t351 * mrSges(5,2);
t756 = m(5) * t372 - t606;
t755 = t739 / 0.2e1 - t744 / 0.2e1;
t624 = t351 * t387;
t625 = t351 * t386;
t754 = t746 * (-t624 / 0.2e1 - t625 / 0.2e1);
t752 = 0.2e1 * t700;
t502 = t351 * t539;
t455 = t502 / 0.2e1;
t622 = t351 * t393;
t560 = -t622 / 0.2e1;
t559 = t622 / 0.2e1;
t748 = -t623 / 0.2e1;
t743 = t253 * t390;
t742 = t253 * t391;
t629 = t253 * t726;
t741 = t393 * t253;
t324 = t351 * t683;
t326 = t352 * t595;
t737 = t326 + t324;
t725 = t360 + t359;
t735 = t760 * pkin(3);
t734 = -t624 - t625;
t679 = t352 * pkin(9);
t252 = -pkin(4) * t351 + t679;
t596 = t686 * pkin(3);
t237 = t596 + t252;
t114 = t390 * t237 - t741;
t637 = t114 * t393;
t113 = t237 * t393 + t743;
t638 = t113 * t390;
t518 = t637 - t638;
t681 = pkin(5) * t351;
t91 = -t113 + t681;
t645 = t91 * t390;
t332 = t351 * qJ(6);
t89 = -t332 + t114;
t647 = t89 * t393;
t526 = t645 + t647;
t241 = mrSges(7,2) * t620 - mrSges(7,3) * t351;
t614 = t393 * t241;
t590 = mrSges(7,2) * t619;
t656 = t351 * mrSges(7,1);
t240 = -t590 + t656;
t616 = t390 * t240;
t719 = t614 / 0.2e1 + t616 / 0.2e1;
t487 = t540 * t352;
t733 = m(6) * t726 - t487;
t626 = t726 * t390;
t682 = pkin(4) * t352;
t236 = pkin(9) * t351 + t372 + t682;
t630 = t236 * t393;
t111 = -t626 + t630;
t613 = t393 * t726;
t617 = t390 * t236;
t112 = t613 + t617;
t732 = (t111 * t390 - t112 * t393 + t726) * t351;
t730 = m(7) * t331;
t729 = mrSges(6,1) + mrSges(7,1);
t676 = Ifges(6,4) - Ifges(7,5);
t727 = t676 * t393;
t621 = t352 * qJ(6);
t86 = t112 + t621;
t382 = t392 * mrSges(4,1);
t724 = -t686 * mrSges(4,2) - t382;
t126 = t390 * t252 - t741;
t94 = -t332 + t126;
t125 = t252 * t393 + t743;
t95 = -t125 + t681;
t525 = t390 * t95 + t393 * t94;
t583 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t718 = -Ifges(5,1) + Ifges(7,2) + Ifges(6,3);
t589 = mrSges(7,2) * t622;
t512 = -t352 * mrSges(7,1) - t589;
t717 = t390 * t512 + t393 * t757;
t706 = t351 ^ 2;
t716 = t352 ^ 2 / 0.2e1 + t706 / 0.2e1;
t381 = m(7) * qJ(6) + mrSges(7,3);
t511 = mrSges(6,2) * t351 + mrSges(6,3) * t620;
t469 = t393 * t511;
t514 = -mrSges(6,1) * t351 + mrSges(6,3) * t619;
t476 = t390 * t514;
t715 = -t469 / 0.2e1 + t476 / 0.2e1 - t719;
t470 = t393 * t512;
t593 = mrSges(6,3) * t622;
t513 = mrSges(6,1) * t352 + t593;
t471 = t393 * t513;
t472 = t390 * t509;
t473 = t390 * t510;
t714 = -t470 / 0.2e1 + t471 / 0.2e1 + t472 / 0.2e1 + t473 / 0.2e1;
t695 = pkin(5) * t95;
t712 = (qJ(6) * t94 - t695) * t700 - t126 * mrSges(6,2) / 0.2e1 + t125 * mrSges(6,1) / 0.2e1 - t95 * mrSges(7,1) / 0.2e1 + t94 * mrSges(7,3) / 0.2e1;
t662 = t114 * mrSges(6,2);
t663 = t113 * mrSges(6,1);
t677 = t91 * mrSges(7,1);
t678 = t89 * mrSges(7,3);
t711 = (-pkin(5) * t91 + qJ(6) * t89) * t700 - t662 / 0.2e1 + t663 / 0.2e1 - t677 / 0.2e1 + t678 / 0.2e1;
t644 = t112 - t86;
t680 = pkin(5) * t352;
t701 = -m(7) / 0.2e1;
t710 = ((t111 - t630 - t680) * t393 + (t613 + t644) * t390) * t701 + t714;
t702 = m(6) / 0.2e1;
t708 = (-pkin(4) * t726 + pkin(9) * t518) * t702 + (pkin(9) * t526 + t357 * t758) * t700 + (t637 / 0.2e1 - t638 / 0.2e1) * mrSges(6,3) + (t647 / 0.2e1 + t645 / 0.2e1) * mrSges(7,2) + t763;
t134 = -t351 * t527 + t253;
t427 = Ifges(6,5) * t352 - t351 * t538;
t417 = t393 * t427;
t428 = Ifges(7,4) * t352 - t351 * t536;
t418 = t393 * t428;
t425 = Ifges(7,6) * t352 - t351 * t530;
t419 = t390 * t425;
t426 = Ifges(6,6) * t352 - t351 * t534;
t420 = t390 * t426;
t456 = t111 + t680;
t429 = t456 * t393;
t484 = t528 * t351;
t707 = t527 * t455 - t417 / 0.4e1 - t418 / 0.4e1 - t419 / 0.4e1 + t420 / 0.4e1 + mrSges(7,2) * t429 / 0.2e1 + t359 * t484 / 0.2e1 + t529 * t559 + t533 * t560 + t253 * t750 - t134 * t539 / 0.2e1 - t740 / 0.4e1 + t761 * t622 / 0.4e1 + t765 * t748 + (t530 / 0.4e1 - t534 / 0.4e1) * t623;
t704 = 2 * qJD(4);
t703 = -m(6) / 0.2e1;
t699 = pkin(4) / 0.2e1;
t698 = m(7) * pkin(3);
t692 = t351 / 0.2e1;
t690 = t352 / 0.2e1;
t684 = m(7) * t393;
t475 = t390 * t513;
t485 = t539 * t352;
t400 = t475 * t693 + (t485 + t487) * t351 + t717 * t692 + t715 * t352;
t430 = t456 * t390;
t11 = ((t134 + t526) * t352 + (-t86 * t393 + t430 + t758) * t351) * t701 + ((t253 + t518) * t352 + t732) * t703 + t400;
t635 = t126 * t393;
t636 = t125 * t390;
t517 = t635 - t636;
t13 = ((t134 + t525) * t352 + (-(t626 - t680) * t390 + (t617 - t86 + t621) * t393 + t478) * t351) * t701 + ((t253 + t517) * t352 + t732) * t703 + t400;
t675 = -t11 * qJD(3) - t13 * qJD(4);
t670 = Ifges(5,4) * t351;
t665 = Ifges(7,2) * t351;
t664 = Ifges(6,3) * t351;
t658 = t726 * mrSges(5,3);
t383 = t393 * mrSges(7,2);
t648 = t86 * t390;
t321 = Ifges(7,6) * t622;
t483 = t528 * t706;
t486 = t351 * t359;
t488 = t360 * t351;
t561 = t623 / 0.2e1;
t9 = -t351 * mrSges(7,2) * t430 + t253 * t488 + t425 * t560 + t426 * t559 + t539 * t483 + t86 * t589 + (Ifges(7,4) * t623 - t321 + t501) * t690 + (t428 + t427) * t561 + (-m(7) * t484 + t486) * t134 + (-m(7) * t456 + t512 - t513 + t593) * t112 + (m(7) * t86 - t594 + t757) * t111 + t713 * t706;
t646 = t9 * qJD(1);
t639 = t112 * t390;
t19 = t470 - t472 - t473 - t471 - mrSges(3,3) - m(7) * (t429 + t648) - m(6) * (t111 * t393 + t639) + (-m(4) - m(3)) * qJ(2) + t724 - t756;
t642 = qJD(1) * t19;
t36 = m(7) * (t134 * t622 + t352 * t86) + t352 * t509 - t393 * t706 * t539;
t641 = qJD(1) * t36;
t640 = t11 * qJD(1);
t634 = t13 * qJD(1);
t633 = t134 * t390;
t399 = -t483 * t701 - t725 * t706 / 0.2e1 + ((t644 * t390 + (t111 - t456) * t393) * t701 + t714 + t754) * t352;
t406 = t528 * t700 - t654 / 0.2e1 + t653 / 0.2e1 + t651 / 0.2e1 + t650 / 0.2e1;
t14 = t399 + t406;
t631 = t14 * qJD(1);
t612 = t393 * t374;
t40 = 0.4e1 * (m(7) / 0.4e1 + m(6) / 0.4e1) * (0.1e1 - t736) * t352 * t351;
t611 = t40 * qJD(2);
t607 = t734 * pkin(9);
t605 = t735 * pkin(9);
t138 = (0.1e1 / 0.2e1 + t716) * t684;
t604 = qJD(1) * t138;
t603 = qJD(5) * t393;
t601 = qJD(3) + qJD(4);
t600 = t91 * t700;
t599 = t95 * t700;
t229 = m(7) * t623;
t506 = m(7) * t527;
t508 = -t672 / 0.2e1 - t674 / 0.2e1;
t50 = (t506 - t671 / 0.2e1 + t673 / 0.2e1 - t508) * t351 + (t539 + t540) * t692;
t598 = t50 * qJD(5) - t229 * qJD(6) + t611;
t592 = mrSges(6,3) * t636;
t591 = mrSges(6,3) * t635;
t584 = -Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t581 = mrSges(7,2) * t689;
t576 = -t383 / 0.2e1;
t575 = t383 / 0.2e1;
t572 = t390 * t685;
t571 = t393 * t685;
t567 = -t639 / 0.2e1;
t566 = t639 / 0.2e1;
t558 = -t619 / 0.2e1;
t552 = mrSges(7,2) * qJ(6) + Ifges(6,6);
t547 = -pkin(5) * mrSges(7,2) + t728;
t546 = -t595 / 0.2e1;
t545 = t595 / 0.2e1;
t523 = t736 * t326 + t324 - t608;
t338 = Ifges(5,4) * t352;
t477 = -t111 * t514 - t112 * t511 + t134 * t485 + t456 * t240 - t86 * t241 + t253 * t487 - t351 * t658 - t372 * (-mrSges(5,1) * t351 - mrSges(5,2) * t352);
t2 = (-t758 * t539 + (-mrSges(5,3) - t540) * t726 + (t113 * t393 + t114 * t390) * mrSges(6,3) + (t390 * t89 - t393 * t91) * mrSges(7,2)) * t351 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t392) * t392 + (-t662 + t663 - t677 + t678) * t352 + m(6) * (t111 * t113 + t112 * t114 + t629) + t172 * t748 + (-Ifges(5,2) * t352 - t670) * t692 - t477 + m(7) * (t134 * t758 - t456 * t91 + t86 * t89) + t173 * t561 + t738 * t560 + (-t723 * t351 + t718 * t352 + t670) * t693 + (t420 - t664 - t665 - t740) * t690 + ((-Ifges(4,1) + Ifges(4,2)) * t392 + qJ(2) * mrSges(4,1) - Ifges(4,4) * t686) * t686 - (-0.2e1 * t338 + t417 + t418 + t419 + (-Ifges(5,1) + Ifges(5,2)) * t351) * t352 / 0.2e1 + t756 * t596;
t522 = t2 * qJD(1) - t11 * qJD(2);
t503 = t351 * t540;
t3 = (-m(7) * (t134 * t643 - t695) - t338 + t740 + (Ifges(5,2) + (-Ifges(6,1) / 0.2e1 - Ifges(7,1) / 0.2e1) * t387 + ((-Ifges(6,2) / 0.2e1 - Ifges(7,3) / 0.2e1) * t390 + t727) * t390 + t718) * t351) * t352 + (t658 + t670 + (t175 / 0.2e1 + t174 / 0.2e1 + t584 * t351) * t393 + (-t173 / 0.2e1 + t172 / 0.2e1 + t583 * t351) * t390) * t351 - m(6) * (t111 * t125 + t112 * t126 + t629) + t726 * t503 + t477 + t758 * t502 - m(7) * (-t111 * t95 + t134 * t478 + t86 * t94) - t125 * t513 - t95 * t512 - t126 * t510 - t94 * t509;
t521 = -t3 * qJD(1) - t13 * qJD(2);
t42 = t497 + (t506 + t539) * t331 + t767;
t507 = -pkin(5) * t240 / 0.2e1 + qJ(6) * t241 / 0.2e1;
t395 = t707 + t507 + t583 * t620 + t728 * t558 + (t648 / 0.2e1 + t567) * mrSges(7,2) - t665 / 0.2e1 - t664 / 0.2e1 + t111 * t576;
t448 = -t486 / 0.2e1;
t482 = t527 * t134;
t397 = (-t331 * t484 + t482) * t701 + t331 * t448 - t375 * t488 / 0.2e1 + t710 * t374 - t746 * t608 / 0.2e1;
t5 = t395 + t397 + t711;
t520 = -t5 * qJD(1) + t42 * qJD(3);
t239 = (t359 + t730) * t390;
t408 = (mrSges(7,1) / 0.2e1 + t539 * t689 + t359 * t687) * t351 - t590;
t416 = m(7) * (-t633 + (t331 * t351 + t374 * t352) * t393);
t30 = t600 - t416 / 0.2e1 + t408;
t516 = qJD(1) * t30 + qJD(3) * t239;
t52 = t655 + 0.2e1 * (t621 / 0.2e1 + t617 / 0.4e1 + t613 / 0.4e1 - t112 / 0.4e1) * m(7);
t515 = qJD(1) * t52 + qJD(5) * t381;
t466 = t331 * t485;
t449 = t487 / 0.2e1;
t260 = t331 * t352;
t303 = t375 * t352;
t405 = (t260 + t523) * t700 + (t303 + t523) * t702;
t411 = (t607 - t682) * t702 + (t352 * t357 + t607) * t700;
t27 = t405 - t411;
t396 = (t375 * t726 + (-t111 * t572 + t112 * t571 + t742) * pkin(3)) * t703 + (t331 * t758 + (t391 * t134 - t685 * t430 + t86 * t571) * pkin(3)) * t701 + t592 / 0.2e1 - t591 / 0.2e1 + t466 / 0.2e1 + t375 * t449 + t94 * t576 + t95 * t581 + t475 * t545 + (t503 / 0.2e1 + t455) * t683 + t717 * t546 + (t517 * t703 + t525 * t701 + t715) * t374 - t755 - t763;
t442 = pkin(9) * t469;
t443 = pkin(9) * t476;
t4 = t396 - t443 / 0.2e1 + t442 / 0.2e1 + t352 * t454 + pkin(4) * t449 + t719 * pkin(9) + t708 + t755;
t412 = -mrSges(5,2) * t595 + (-mrSges(5,1) + t725) * t683 + t746 * t735;
t431 = t760 * t374;
t59 = -(t331 * t391 + t431) * t698 - m(6) * (t375 * t391 + t431) * pkin(3) - t412;
t447 = -t4 * qJD(1) + t27 * qJD(2) - t59 * qJD(3);
t403 = (-pkin(5) * t572 + qJ(6) * t571) * t698 / 0.2e1 + t545 * t671 + (t729 * t390 + t672) * t546;
t34 = t403 + t766;
t51 = -m(7) * t498 - t293 - t294 + (pkin(4) * mrSges(6,2) - t727) * t393 + (pkin(4) * mrSges(6,1) + t676 * t390 + (Ifges(6,2) + Ifges(7,3) - t768) * t393) * t390;
t398 = (-t357 * t484 + t482) * t701 + t488 * t699 + t357 * t448 + (t710 + t754) * pkin(9);
t8 = t395 + t398 + t712;
t446 = -t8 * qJD(1) - t34 * qJD(3) - t51 * qJD(4);
t445 = t725 * t352 + t734 * t746 + t606;
t136 = (t359 + (t545 + t694 + t357 / 0.2e1) * m(7)) * t390;
t259 = (m(7) * t357 + t359) * t390;
t422 = m(7) * (-t633 + (t357 * t351 + t679) * t393);
t32 = t599 - t422 / 0.2e1 + t408;
t444 = qJD(1) * t32 + qJD(3) * t136 + qJD(4) * t259;
t423 = t390 * t455 + t359 * t559 + t352 * t575 + t656 / 0.2e1 + mrSges(7,2) * t558;
t409 = t547 * t603 + (Ifges(7,6) - t552) * qJD(5) * t390;
t407 = (-m(7) * t528 + t725) * qJD(5);
t401 = -t707 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t351 + t507 + (t390 * t583 + t393 * t584) * t352 + mrSges(7,2) * t566 + t111 * t575 + t86 * t581 + (t566 + t567) * mrSges(6,3);
t358 = pkin(9) * t684 + t383;
t330 = m(7) * t612 + t383;
t139 = (-0.1e1 / 0.2e1 + t716) * t684;
t137 = m(7) * (-t331 - t357) * t688 + (m(7) * t545 - t359) * t390;
t41 = t86 * t752 + t509;
t35 = t403 - t766;
t33 = t422 / 0.2e1 + t599 + t423;
t31 = t416 / 0.2e1 + t600 + t423;
t17 = t405 + t411 + t445;
t15 = -t399 + t406;
t7 = t401 - t398 + t712;
t6 = t401 - t397 + t711;
t1 = (t360 / 0.2e1 - mrSges(5,1) / 0.2e1) * t726 + (-t351 * t508 + t719) * pkin(9) + (t540 * t699 + t454) * t352 - t396 + t414 + t708;
t10 = [-qJD(2) * t19 + qJD(3) * t2 - qJD(4) * t3 + qJD(5) * t9 + qJD(6) * t36, qJD(5) * t15 + qJD(6) * t139 - t642 + t675, t1 * qJD(4) + t6 * qJD(5) + t31 * qJD(6) + t522 + (t758 * t730 + m(5) * (-t685 * t726 - t742) * pkin(3) - t394 * t382 - mrSges(4,2) * t570 + t241 * t612 - Ifges(4,6) * t686 - t466 - Ifges(4,5) * t392 + (m(6) * t518 + m(7) * t526 + t469 - t476 + t616) * t374 + t733 * t375 + t518 * mrSges(6,3) + t737 * mrSges(5,3) + t526 * mrSges(7,2) + t764) * qJD(3), t1 * qJD(3) + (t442 - t443 + t591 - t592 + (m(7) * t758 - t485) * t357 - t733 * pkin(4) + (m(6) * t517 + m(7) * t525 + t614 + t616) * pkin(9) + t525 * mrSges(7,2) + t764) * qJD(4) + t7 * qJD(5) + t33 * qJD(6) + t521, t15 * qJD(2) + t6 * qJD(3) + t7 * qJD(4) + t41 * qJD(6) + t646 + (-t321 + (-m(7) * pkin(5) - t729) * t112 + (-mrSges(6,2) + t381) * t111 + (t390 * t547 + t393 * t552) * t351) * qJD(5), qJD(2) * t139 + qJD(3) * t31 + qJD(4) * t33 + qJD(5) * t41 + t641; -qJD(5) * t14 + qJD(6) * t138 + t642 + t675, t601 * t40, t17 * qJD(4) + t598 - t640 + (t445 + 0.2e1 * (t303 - t608) * t702 + (t260 - t608) * t752 - m(5) * t737 + t724) * qJD(3), t17 * qJD(3) + qJD(4) * t445 + t411 * t704 + t598 - t634, -t631 + t601 * t50 + (qJD(6) * t684 + t407) * t352, m(7) * t352 * t603 - t229 * t601 + t604; -qJD(4) * t4 - qJD(5) * t5 - qJD(6) * t30 - t522, qJD(4) * t27 - t611 + t640, -qJD(4) * t59 + qJD(5) * t42 - qJD(6) * t239, t412 * qJD(4) + t35 * qJD(5) + t137 * qJD(6) + ((t357 * t683 + t605) * t700 + (-pkin(4) * t683 + t605) * t702) * t704 + t447, t35 * qJD(4) + t330 * qJD(6) + t374 * t407 + t409 + t520, qJD(4) * t137 + qJD(5) * t330 - t516; qJD(3) * t4 - qJD(5) * t8 - qJD(6) * t32 - t521, -qJD(3) * t27 - t611 + t634, -qJD(5) * t34 - qJD(6) * t136 - t447, -qJD(5) * t51 - qJD(6) * t259, pkin(9) * t407 + t358 * qJD(6) + t409 + t446, qJD(5) * t358 - t444; qJD(2) * t14 + qJD(3) * t5 + qJD(4) * t8 + qJD(6) * t52 - t646, t631, qJD(4) * t34 - t520, -t446, t381 * qJD(6), t515; -qJD(2) * t138 + qJD(3) * t30 + qJD(4) * t32 - qJD(5) * t52 - t641, -t604, qJD(4) * t136 + t516, t444, -t515, 0;];
Cq  = t10;

% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:02:02
% EndTime: 2019-03-09 22:02:34
% DurationCPUTime: 17.92s
% Computational Cost: add. (36090->695), mult. (69640->867), div. (0->0), fcn. (79522->8), ass. (0->402)
t429 = sin(qJ(2));
t725 = -pkin(8) - pkin(7);
t402 = t725 * t429;
t432 = cos(qJ(2));
t403 = t725 * t432;
t428 = sin(qJ(3));
t431 = cos(qJ(3));
t339 = t402 * t428 - t403 * t431;
t391 = -t428 * t429 + t431 * t432;
t280 = pkin(9) * t391 + t339;
t427 = sin(qJ(4));
t703 = cos(qJ(4));
t392 = -t428 * t432 - t431 * t429;
t769 = t402 * t431 + t428 * t403;
t795 = t392 * pkin(9) + t769;
t811 = t703 * t280 + t427 * t795;
t847 = t811 * mrSges(6,2);
t848 = t811 * mrSges(5,1);
t810 = -t427 * t280 + t703 * t795;
t849 = t810 * mrSges(6,3);
t850 = t810 * mrSges(5,2);
t426 = sin(qJ(6));
t420 = t426 * mrSges(7,1);
t430 = cos(qJ(6));
t421 = t430 * mrSges(7,2);
t768 = t421 + t420;
t488 = t427 * t391 - t392 * t703;
t837 = -t488 * pkin(5) + t810;
t869 = t837 * t768;
t874 = t847 - t848 + t849 - t850 + t869;
t322 = t703 * t391 + t392 * t427;
t418 = -pkin(2) * t432 - pkin(1);
t352 = -pkin(3) * t391 + t418;
t466 = -qJ(5) * t488 + t352;
t165 = -pkin(4) * t322 + t466;
t642 = t430 * mrSges(7,1);
t647 = t426 * mrSges(7,2);
t513 = t642 - t647;
t203 = t513 * t322;
t425 = t430 ^ 2;
t680 = Ifges(7,6) * t430;
t682 = Ifges(7,5) * t426;
t509 = t680 + t682;
t726 = pkin(4) + pkin(10);
t91 = -t322 * t726 + t466;
t56 = -t426 * t91 - t430 * t837;
t57 = -t426 * t837 + t430 * t91;
t586 = -Ifges(5,4) / 0.2e1 - Ifges(6,6) / 0.2e1;
t684 = Ifges(7,4) * t430;
t687 = Ifges(7,1) * t426;
t716 = -t488 / 0.2e1;
t717 = -t322 / 0.2e1;
t730 = Ifges(7,3) / 0.2e1;
t731 = Ifges(7,2) / 0.2e1;
t773 = t513 * t488;
t641 = t430 * mrSges(7,3);
t791 = -mrSges(7,2) * t322 + t488 * t641;
t646 = t426 * mrSges(7,3);
t792 = mrSges(7,1) * t322 - t488 * t646;
t823 = Ifges(5,4) + Ifges(6,6);
t836 = pkin(5) * t322 + t811;
t873 = t837 * t203 + t56 * t792 + t57 * t791 - t773 * t836 + (-t322 * (t425 * t731 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + (t684 + t687 / 0.2e1) * t426 - t730 - Ifges(5,1) / 0.2e1 - Ifges(6,2) / 0.2e1) + (t509 + t586) * t488 + t352 * mrSges(5,1) - t165 * mrSges(6,2) + (Ifges(5,2) + Ifges(6,3)) * t717 + t823 * t716) * t488;
t314 = Ifges(5,6) * t488;
t315 = Ifges(6,5) * t488;
t316 = Ifges(5,5) * t322;
t317 = Ifges(6,4) * t322;
t396 = -Ifges(7,2) * t426 + t684;
t603 = t430 * t396;
t685 = Ifges(7,4) * t426;
t397 = Ifges(7,1) * t430 - t685;
t610 = t426 * t397;
t489 = t610 / 0.2e1 + t603 / 0.2e1;
t511 = t684 + t687;
t790 = Ifges(7,5) * t322 + t488 * t511;
t608 = t430 * t790;
t510 = Ifges(7,2) * t430 + t685;
t789 = Ifges(7,6) * t322 + t488 * t510;
t615 = t426 * t789;
t395 = Ifges(7,5) * t430 - Ifges(7,6) * t426;
t623 = t322 * t395;
t457 = -t615 / 0.2e1 + t608 / 0.2e1 + t623 / 0.2e1 + t315 + t316 + Ifges(4,6) * t392 + Ifges(4,5) * t391 - t314 - t317 + t489 * t488;
t802 = t769 * mrSges(4,2);
t803 = t339 * mrSges(4,1);
t870 = t457 - t802 - t803 + t874;
t854 = -t847 / 0.2e1 + t848 / 0.2e1 - t849 / 0.2e1 + t850 / 0.2e1;
t746 = -t869 / 0.2e1 + t854;
t864 = qJ(5) * t837;
t561 = t703 * t428;
t699 = pkin(2) * t431;
t565 = pkin(3) + t699;
t373 = pkin(2) * t561 + t427 * t565;
t365 = qJ(5) + t373;
t863 = t365 * t837;
t695 = t427 * pkin(3);
t414 = qJ(5) + t695;
t862 = t414 * t837;
t860 = t836 * t837;
t771 = m(6) * t165 + mrSges(6,2) * t322 - mrSges(6,3) * t488;
t706 = t430 / 0.2e1;
t715 = t488 / 0.2e1;
t718 = t322 / 0.2e1;
t751 = t823 * t718;
t852 = t790 / 0.2e1;
t868 = t426 * t852 + t706 * t789 - t352 * mrSges(5,2) + t165 * mrSges(6,3) - t509 * t717 - t751 - (Ifges(5,1) + Ifges(6,2) + Ifges(7,3)) * t715;
t865 = m(5) * t352;
t845 = t426 * t836;
t843 = t430 * t836;
t734 = m(5) * pkin(3);
t832 = t427 * t810 - t703 * t811;
t861 = t832 * t734;
t596 = t703 * pkin(3);
t417 = -t596 - pkin(4);
t833 = t810 * t414 + t417 * t811;
t764 = -t610 / 0.4e1 - t603 / 0.4e1;
t858 = -t488 * t764 + t608 / 0.4e1;
t739 = m(6) / 0.2e1;
t831 = (-pkin(4) * t811 + qJ(5) * t810) * t739;
t853 = t418 * (-mrSges(4,1) * t392 + mrSges(4,2) * t391) + (-t322 * t586 - t868) * t322 + t873;
t804 = t322 * mrSges(6,1);
t613 = t426 * t791;
t196 = t613 / 0.2e1;
t607 = t430 * t792;
t840 = t196 + t607 / 0.2e1;
t800 = qJ(5) * t773;
t839 = -t314 / 0.2e1 + t315 / 0.2e1 + t316 / 0.2e1 - t317 / 0.2e1 - t615 / 0.4e1 + t623 / 0.4e1 - t800 / 0.2e1 + t858;
t635 = qJ(5) * t322;
t131 = t488 * t726 - t635;
t68 = -t131 * t426 + t843;
t638 = t430 * t68;
t69 = t131 * t430 + t845;
t643 = t426 * t69;
t505 = t638 + t643;
t736 = m(7) / 0.2e1;
t838 = t505 * t736 + t840;
t700 = pkin(2) * t428;
t407 = t427 * t700;
t372 = -t703 * t565 + t407;
t367 = -pkin(4) + t372;
t830 = t365 * t810 + t367 * t811;
t829 = t372 * t811 + t373 * t810;
t492 = t647 / 0.2e1 - t642 / 0.2e1;
t462 = (-mrSges(6,1) / 0.2e1 + t492) * t322;
t576 = -t804 / 0.2e1;
t548 = -t607 / 0.2e1;
t765 = t548 - t613 / 0.2e1;
t828 = -t462 + t576 + t765;
t663 = t488 * mrSges(6,1);
t585 = mrSges(6,3) + t768;
t820 = t367 * t804;
t424 = t426 ^ 2;
t602 = t424 + t425;
t782 = mrSges(7,3) * t602;
t816 = (mrSges(5,1) + t782) * t695;
t754 = t489 + t511 * t706 - t426 * t510 / 0.2e1;
t688 = mrSges(7,3) * t425;
t689 = mrSges(7,3) * t424;
t761 = -t689 / 0.2e1 - t688 / 0.2e1;
t815 = pkin(4) * t488 - t635;
t624 = t322 * t430;
t210 = -mrSges(7,2) * t488 - mrSges(7,3) * t624;
t604 = t430 * t210;
t625 = t322 * t426;
t207 = mrSges(7,1) * t488 + mrSges(7,3) * t625;
t614 = t426 * t207;
t491 = t614 / 0.2e1 - t604 / 0.2e1;
t766 = t761 * t322;
t749 = -t766 - t491;
t812 = t802 / 0.2e1 + t803 / 0.2e1;
t805 = t602 * t736;
t798 = -t820 / 0.2e1;
t797 = t414 * t663;
t381 = t699 * t703 - t407;
t340 = t381 * t768;
t380 = (t427 * t431 + t561) * pkin(2);
t375 = t380 * mrSges(6,2);
t376 = t380 * mrSges(5,1);
t793 = -(mrSges(4,1) * t428 + mrSges(4,2) * t431) * pkin(2) + t340 + t375 - t376;
t738 = m(6) / 0.4e1;
t787 = t805 + 0.2e1 * t738;
t598 = m(7) / 0.4e1 + t738;
t786 = 0.4e1 * t598;
t727 = -m(7) - m(6);
t785 = t395 / 0.4e1;
t709 = -t414 / 0.2e1;
t694 = mrSges(5,2) - mrSges(6,3);
t783 = mrSges(6,2) - mrSges(5,1);
t536 = t424 / 0.2e1 + t425 / 0.2e1;
t517 = mrSges(7,3) * t536;
t778 = -t797 / 0.2e1;
t698 = pkin(3) * t392;
t173 = -t698 + t815;
t423 = t429 * pkin(2);
t166 = t173 + t423;
t319 = t488 * pkin(10);
t94 = t166 + t319;
t58 = -t426 * t94 + t843;
t640 = t430 * t58;
t59 = t430 * t94 + t845;
t645 = t426 * t59;
t758 = -t645 - t640;
t774 = t758 * t736;
t767 = t805 + t739;
t107 = t173 + t319;
t62 = -t107 * t426 + t843;
t63 = t107 * t430 + t845;
t732 = -mrSges(7,2) / 0.2e1;
t733 = mrSges(7,1) / 0.2e1;
t760 = t62 * t733 + t63 * t732;
t729 = -t59 / 0.2e1;
t759 = mrSges(7,2) * t729 + t58 * t733;
t757 = -t688 - t689;
t755 = t736 + t739;
t606 = t430 * t207;
t612 = t426 * t210;
t753 = t612 / 0.2e1 + t606 / 0.2e1 + mrSges(6,1) * t715;
t752 = (-mrSges(5,2) + t585) * t703;
t374 = -t727 * qJ(5) + t585;
t569 = t641 / 0.2e1;
t570 = -t641 / 0.2e1;
t681 = Ifges(7,6) * t488;
t683 = Ifges(7,5) * t488;
t708 = -t426 / 0.4e1;
t723 = t836 / 0.2e1;
t750 = t513 * t723 + (-t322 * t511 + t683) * t708 - t430 * (-t322 * t510 + t681) / 0.4e1 + (-t397 / 0.2e1 + t510 / 0.4e1) * t624 + (t569 + t570) * t57 + (0.2e1 * t396 + t511) * t625 / 0.4e1 + (t682 / 0.2e1 + t680 / 0.2e1 - t509 / 0.4e1) * t488;
t748 = t869 / 0.2e1 + (-t645 / 0.2e1 - t640 / 0.2e1) * mrSges(7,3);
t743 = 0.2e1 * m(7);
t742 = 0.2e1 * qJ(5);
t741 = 2 * qJD(4);
t740 = -m(6) / 0.2e1;
t737 = -m(7) / 0.2e1;
t728 = t68 / 0.2e1;
t724 = -qJ(5) / 0.2e1;
t722 = t836 / 0.4e1;
t721 = t811 / 0.2e1;
t719 = t207 / 0.2e1;
t362 = -pkin(10) + t367;
t714 = -t362 / 0.2e1;
t713 = -t365 / 0.2e1;
t712 = t372 / 0.2e1;
t711 = -t373 / 0.2e1;
t710 = -t380 / 0.2e1;
t705 = t726 / 0.2e1;
t704 = -t726 / 0.2e1;
t701 = m(6) * t811;
t690 = mrSges(5,3) * t322;
t686 = Ifges(4,4) * t392;
t217 = -mrSges(5,1) * t322 + mrSges(5,2) * t488;
t514 = Ifges(4,4) * t391 + (-Ifges(4,1) + Ifges(4,2)) * t392;
t1 = m(7) * (t56 * t58 + t57 * t59 + t860) + m(4) * t418 * t423 + (-mrSges(4,1) * t423 + t514) * t391 + (-mrSges(4,2) * t423 - t686) * t392 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t429) * t429 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t432 + (Ifges(3,1) - Ifges(3,2)) * t429) * t432 + t58 * t207 + t59 * t210 + t853 + (t865 + t217) * (t423 - t698) + t771 * t166;
t679 = t1 * qJD(1);
t678 = t836 * mrSges(7,1);
t677 = t836 * mrSges(7,2);
t2 = m(7) * (t56 * t62 + t57 * t63 + t860) - t698 * t865 + (-pkin(3) * t217 - t686) * t392 + t514 * t391 + t62 * t207 + t63 * t210 + t771 * t173 + t853;
t664 = t2 * qJD(1);
t662 = t488 * mrSges(5,3);
t655 = t365 * mrSges(6,1);
t653 = t372 * mrSges(5,2);
t652 = t372 * mrSges(6,3);
t651 = t381 * mrSges(5,2);
t650 = t381 * mrSges(6,3);
t648 = t417 * mrSges(6,1);
t644 = t426 * t63;
t639 = t430 * t62;
t9 = t68 * t207 + t69 * t210 + m(7) * (t56 * t68 + t57 * t69 + t860) - (-t751 + t868) * t322 + t771 * t815 + t873;
t637 = t9 * qJD(1);
t634 = qJ(5) * t372;
t633 = qJ(5) * t381;
t24 = (-t604 + m(7) * (t426 * t56 - t430 * t57) + t614 - t771) * t488;
t632 = qJD(1) * t24;
t299 = Ifges(7,5) * t624;
t10 = -t299 * t715 + t56 * t210 - t57 * t207 - ((t677 - t56 * mrSges(7,3) + t683 / 0.2e1 - Ifges(7,4) * t624) * t430 + (t678 - t57 * mrSges(7,3) - t681 - (-t685 + (Ifges(7,1) - Ifges(7,2)) * t430) * t322) * t426) * t322;
t631 = t10 * qJD(1);
t477 = (t421 / 0.2e1 + t420 / 0.2e1) * t488;
t496 = t322 * t517;
t27 = t477 - t496 + t491;
t626 = t27 * qJD(1);
t622 = t365 * t773;
t621 = t365 * t372;
t620 = t365 * t381;
t617 = t414 * t773;
t616 = t414 * t372;
t601 = 0.2e1 * t695;
t599 = mrSges(6,2) * t695;
t597 = mrSges(7,3) * t644;
t595 = t372 * t690;
t594 = t373 * t662;
t593 = m(7) * t723;
t587 = mrSges(6,3) / 0.2e1 - mrSges(5,2) / 0.2e1;
t580 = t322 * t730;
t579 = -t678 / 0.2e1;
t578 = t677 / 0.2e1;
t577 = t663 / 0.2e1;
t302 = t804 / 0.2e1;
t572 = -t646 / 0.2e1;
t334 = t372 * t768;
t363 = t373 * mrSges(6,2);
t364 = t373 * mrSges(5,1);
t567 = t334 + t364 - t363;
t564 = t703 * t365;
t562 = t703 * t414;
t549 = t426 * t705;
t541 = t430 * t704;
t537 = 0.2e1 * t715;
t533 = t602 * t373;
t532 = t602 * t380;
t408 = -pkin(10) + t417;
t531 = t602 * t408;
t530 = t602 * t427;
t529 = t602 * t726;
t528 = t662 * t695;
t526 = -t596 / 0.2e1;
t525 = t596 / 0.2e1;
t516 = t596 * t690;
t508 = t426 * t57 + t430 * t56;
t506 = t639 + t644;
t487 = qJ(5) * t513;
t504 = t487 / 0.2e1 - t754;
t45 = -t694 * t372 + mrSges(7,3) * t533 - m(7) * (t362 * t533 - t621) - m(6) * (t367 * t373 - t621) + t567;
t438 = ((t367 - t372) * t811 - (-t365 + t373) * t810) * t740 + (t362 * t505 - t372 * t836 + t373 * t508 + t863) * t737 - t773 * t713 + t203 * t712;
t453 = pkin(4) * t576 - t196 * t726 + t541 * t792 + t663 * t724 + t839;
t439 = t453 + t831 + (t726 * t758 + t864) * t736;
t5 = t438 + (Ifges(5,6) / 0.2e1 - Ifges(6,5) / 0.2e1 + (t365 / 0.2e1 + t711) * mrSges(6,1) + t764) * t488 - (t785 - Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 + (t367 / 0.2e1 - t372 / 0.2e1) * mrSges(6,1)) * t322 + (t789 / 0.4e1 + t210 * t711 + t791 * t714 + (t69 / 0.2e1 + t729) * mrSges(7,3)) * t426 + (-t790 / 0.4e1 + t207 * t711 + t792 * t714 + (t728 - t58 / 0.2e1) * mrSges(7,3)) * t430 + t439 + t783 * (t721 - t811 / 0.2e1);
t503 = -t5 * qJD(1) - t45 * qJD(2);
t46 = t694 * t381 + mrSges(7,3) * t532 - m(7) * (t362 * t532 + t620) - m(6) * (t367 * t380 + t620) - m(5) * (t372 * t380 + t373 * t381) - t793;
t499 = -t380 * t810 + t381 * t811;
t434 = -m(5) * (t499 + t829) / 0.2e1 + (t499 + t830) * t740 + (t380 * t508 + t381 * t836 + t863) * t737 + t622 / 0.2e1 - t381 * t203 / 0.2e1 + t365 * t577 + t798 - t595 / 0.2e1 + t594 / 0.2e1 + (t612 + t606) * t710 + (t506 * t737 + t765) * t362 + (mrSges(5,3) + mrSges(6,1)) * (t381 * t717 + t488 * t710) + t812;
t436 = t833 * t739 + t736 * t862 - t617 / 0.2e1 + t861 / 0.2e1 + t778 + t417 * t302 - t528 / 0.2e1 - t516 / 0.2e1 + (t840 - t774) * t408 + t748 - t812 - t854;
t445 = t62 * t569 + t597 / 0.2e1 + t746;
t6 = t434 + t436 + t445;
t502 = -t6 * qJD(1) - t46 * qJD(2);
t461 = (-Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.2e1) * t426 + t396 / 0.4e1 + 0.3e1 / 0.2e1 * t684;
t464 = Ifges(7,6) * t537 + t579;
t465 = Ifges(7,5) * t537 + t578;
t482 = (t731 - Ifges(7,1) / 0.4e1) * t430 - t397 / 0.4e1;
t11 = t580 - t362 * t496 + (t210 * t714 - (mrSges(7,2) * t713 + t482) * t322 + t464) * t430 + (t362 * t719 - (mrSges(7,1) * t713 + t461) * t322 + t465) * t426 + t759;
t486 = t365 * t513;
t183 = t486 - t754;
t498 = -t11 * qJD(1) + t183 * qJD(2);
t20 = (t722 - t645 / 0.4e1 - t640 / 0.4e1) * t743 + t828;
t229 = t365 * t786 + t585;
t497 = -qJD(1) * t20 - qJD(2) * t229;
t494 = mrSges(7,1) * t728 + t69 * t732;
t485 = t414 * t513;
t481 = t62 * t570 + t63 * t572 - t746;
t479 = t492 * t373;
t478 = t492 * t380;
t476 = -t487 / 0.2e1;
t475 = -t486 / 0.2e1;
t472 = t768 * t717;
t106 = t599 + (m(7) * (t408 * t530 + t562) + m(6) * (t417 * t427 + t562) + t752) * pkin(3) - t816;
t437 = -t334 / 0.2e1 + t363 / 0.2e1 - t364 / 0.2e1 + (-t616 + (t367 * t427 + t564) * pkin(3)) * t739 + (-t616 + (t362 * t530 + t564) * pkin(3)) * t736 + t653 / 0.2e1 - t652 / 0.2e1 + t599 / 0.2e1 + mrSges(5,2) * t526 + t585 * t525 + (t417 * t739 + t531 * t736 + t761) * t373 - t816 / 0.2e1;
t447 = t340 / 0.2e1 + t375 / 0.2e1 - t376 / 0.2e1 + (-pkin(4) * t380 + t633) * t739 + (-t380 * t529 + t633) * t736;
t29 = -t437 - t651 / 0.2e1 + t650 / 0.2e1 + t447 + t761 * t380;
t452 = Ifges(6,4) * t717 + Ifges(5,5) * t718 + Ifges(6,5) * t715 + Ifges(5,6) * t716 + t322 * t785 + t68 * t570 + t69 * t572 + t789 * t708 - t746 + t858;
t435 = t203 * t525 + t452 - t526 * t804 + (t862 + (t427 * t508 + t703 * t836) * pkin(3)) * t736 + (-pkin(3) * t832 + t833) * t739 + t773 * t709 + t778 + t648 * t718 + t838 * t408 + t753 * t695;
t448 = t831 + (-t506 * t726 + t864) * t736;
t8 = pkin(4) * t302 + qJ(5) * t577 - t548 * t726 + t549 * t791 + t435 + t445 - t448 - t839;
t470 = t8 * qJD(1) - t29 * qJD(2) + t106 * qJD(3);
t469 = t373 * t787;
t13 = t580 - t408 * t496 + (-t408 * t210 / 0.2e1 - (mrSges(7,2) * t709 + t482) * t322 + t464) * t430 + (t408 * t719 - (mrSges(7,1) * t709 + t461) * t322 + t465) * t426 + t760;
t213 = t485 - t754;
t442 = -t485 / 0.2e1 + t754;
t82 = t475 - t478 + t442;
t468 = -t13 * qJD(1) - t82 * qJD(2) + t213 * qJD(3);
t22 = (t722 - t644 / 0.4e1 - t639 / 0.4e1) * t743 + t828;
t324 = t414 * t786 + t585;
t441 = t585 + t755 * (t742 + t380 + t601);
t92 = t380 * t767 - t441;
t467 = -qJD(1) * t22 + qJD(2) * t92 - qJD(3) * t324;
t463 = t492 * t695;
t132 = t476 - t463 + t442;
t15 = -(-Ifges(7,3) / 0.2e1 - t726 * t517) * t322 + (t681 + t579 + t210 * t705 - (mrSges(7,2) * t724 + t482) * t322) * t430 + (t683 + t578 + t207 * t704 - (mrSges(7,1) * t724 + t461) * t322) * t426 + t494;
t224 = t487 - t754;
t84 = t476 + t475 - t479 + t754;
t459 = t15 * qJD(1) + t84 * qJD(2) + t132 * qJD(3) - t224 * qJD(4);
t443 = -t585 + (t737 + t740) * (t742 + t373);
t104 = t469 + t443;
t240 = (-0.1e1 / 0.2e1 + t536) * m(7) * t695 - t374;
t26 = -t492 * t322 + (t722 - t643 / 0.4e1 - t638 / 0.4e1) * t743 - t840;
t458 = qJD(1) * t26 - qJD(2) * t104 - qJD(3) * t240 + qJD(4) * t374;
t450 = t302 - t462 + t701 / 0.2e1 + t593 + t840;
t446 = t580 + t750;
t366 = t485 / 0.2e1;
t326 = t486 / 0.2e1;
t241 = t767 * t695 + t585 + t598 * (0.4e1 * qJ(5) + t601);
t133 = t366 - t463 + t504;
t105 = t469 - t443;
t93 = t380 * t787 + t441;
t85 = t326 - t479 + t504;
t83 = t326 + t366 - t478 - t754;
t30 = -t380 * t517 + t381 * t587 + t437 + t447;
t28 = t477 + t749;
t25 = t593 + t701 - (-mrSges(6,1) + t492) * t322 + t838;
t23 = t506 * t736 + t739 * t811 + t450;
t21 = m(6) * t721 + t450 - t774;
t16 = Ifges(7,3) * t718 + qJ(5) * t472 + t207 * t549 + t210 * t541 + t726 * t766 + t494 + t750;
t14 = t408 * t749 + t414 * t472 + t446 + t760;
t12 = t362 * t749 + t365 * t472 + t446 + t759;
t7 = t453 + t448 + t435 + t481;
t4 = -t438 + t452 + (mrSges(6,2) / 0.2e1 - mrSges(5,1) / 0.2e1) * t811 + t587 * t810 + t655 * t716 - t798 - t804 * t712 + t439 + t840 * t362 + t753 * t373 + t748;
t3 = -t434 + t436 + t457 + t481;
t17 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t9 + qJD(5) * t24 + qJD(6) * t10, t3 * qJD(3) + t4 * qJD(4) + t21 * qJD(5) + t12 * qJD(6) + t679 + (m(4) * (-t339 * t431 + t428 * t769) * pkin(2) - t622 - t594 - t488 * t655 + (-mrSges(3,1) * t432 + mrSges(3,2) * t429) * pkin(7) + (-t391 * t699 + t392 * t700) * mrSges(4,3) + t362 * t613 + t362 * t607 + t820 + t595 + Ifges(3,5) * t432 - Ifges(3,6) * t429 + 0.2e1 * (-t362 * t758 + t863) * t736 + m(5) * t829 + 0.2e1 * t830 * t739 + t758 * mrSges(7,3) + t870) * qJD(2), t664 + t3 * qJD(2) + (t861 + m(6) * t833 - t528 - t516 + m(7) * t862 - t617 - t797 + t322 * t648 - t597 - mrSges(7,3) * t639 + (m(7) * t506 + t607 + t613) * t408 + t870) * qJD(3) + t7 * qJD(4) + t23 * qJD(5) + t14 * qJD(6), t637 + t4 * qJD(2) + t7 * qJD(3) + t25 * qJD(5) + t16 * qJD(6) + ((-t505 * t726 + t864) * t736 + t831) * t741 + (-t800 + (-t68 * mrSges(7,3) - t726 * t792 + t852) * t430 + (-t789 / 0.2e1 - t69 * mrSges(7,3) - t726 * t791) * t426 - (-t395 / 0.2e1 - Ifges(5,5) + Ifges(6,4) + pkin(4) * mrSges(6,1)) * t322 + (-qJ(5) * mrSges(6,1) + Ifges(6,5) - Ifges(5,6) + t489) * t488 + t874) * qJD(4), qJD(2) * t21 + qJD(3) * t23 + qJD(4) * t25 + qJD(6) * t28 + t632, t631 + t12 * qJD(2) + t14 * qJD(3) + t16 * qJD(4) + t28 * qJD(5) + (-mrSges(7,1) * t57 - mrSges(7,2) * t56 + Ifges(7,6) * t625 - t299) * qJD(6); -qJD(3) * t6 - qJD(4) * t5 + qJD(5) * t20 - qJD(6) * t11 - t679, -qJD(3) * t46 - qJD(4) * t45 + qJD(5) * t229 + qJD(6) * t183 (t650 - t651 + (m(6) * t417 + m(7) * t531 - t703 * t734 + t757) * t380 + (-t414 * t727 + t427 * t734) * t381 + t793) * qJD(3) + t30 * qJD(4) + t93 * qJD(5) + t83 * qJD(6) + t502, t30 * qJD(3) + (t373 * t757 - t567 - t652 + t653) * qJD(4) + t105 * qJD(5) + t85 * qJD(6) + ((-t373 * t529 - t634) * t736 + (-pkin(4) * t373 - t634) * t739) * t741 + t503, qJD(3) * t93 + qJD(4) * t105 - t497, t83 * qJD(3) + t85 * qJD(4) + (-t362 * t768 - t509) * qJD(6) + t498; qJD(2) * t6 + qJD(4) * t8 + qJD(5) * t22 - qJD(6) * t13 - t664, -qJD(4) * t29 - qJD(5) * t92 - qJD(6) * t82 - t502, qJD(4) * t106 + qJD(5) * t324 + qJD(6) * t213, t241 * qJD(5) + t133 * qJD(6) + t755 * qJ(5) * t596 * t741 + ((-m(6) * pkin(4) - m(7) * t529 - t782 + t783) * t427 + t752) * qJD(4) * pkin(3) + t470, qJD(4) * t241 - t467, t133 * qJD(4) + (-t408 * t768 - t509) * qJD(6) + t468; qJD(2) * t5 - qJD(3) * t8 + qJD(5) * t26 - qJD(6) * t15 - t637, qJD(3) * t29 - qJD(5) * t104 - qJD(6) * t84 - t503, -qJD(5) * t240 - qJD(6) * t132 - t470, qJD(5) * t374 + qJD(6) * t224, t458 ((mrSges(7,2) * t726 - Ifges(7,6)) * t430 + (mrSges(7,1) * t726 - Ifges(7,5)) * t426) * qJD(6) - t459; -qJD(2) * t20 - qJD(3) * t22 - qJD(4) * t26 - qJD(6) * t27 - t632, qJD(3) * t92 + qJD(4) * t104 + t497, qJD(4) * t240 + t467, -t458, 0, -qJD(6) * t768 - t626; qJD(2) * t11 + qJD(3) * t13 + qJD(4) * t15 + qJD(5) * t27 - t631, qJD(3) * t82 + qJD(4) * t84 - t498, qJD(4) * t132 - t468, t459, t626, 0;];
Cq  = t17;

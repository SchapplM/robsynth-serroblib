% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:32
% EndTime: 2019-12-31 22:15:07
% DurationCPUTime: 16.91s
% Computational Cost: add. (18922->900), mult. (47285->1220), div. (0->0), fcn. (47548->8), ass. (0->432)
t884 = Ifges(6,4) + Ifges(5,5);
t591 = cos(qJ(3));
t587 = sin(qJ(4));
t577 = Ifges(6,6) * t587;
t590 = cos(qJ(4));
t581 = Ifges(6,4) * t590;
t874 = t581 + t577;
t579 = Ifges(5,5) * t590;
t776 = Ifges(5,6) * t587;
t875 = t579 - t776;
t885 = t874 / 0.4e1 + t875 / 0.4e1;
t893 = t885 * t591;
t665 = t590 * mrSges(5,1) - t587 * mrSges(5,2);
t892 = -mrSges(4,1) - t665;
t586 = sin(pkin(5));
t592 = cos(qJ(2));
t737 = t586 * t592;
t709 = t587 * t737;
t589 = sin(qJ(2));
t738 = t586 * t589;
t438 = -t590 * t738 + t591 * t709;
t727 = t591 * t592;
t439 = (t587 * t589 + t590 * t727) * t586;
t871 = (Ifges(5,1) + Ifges(6,1)) * t439 + (-Ifges(5,4) + Ifges(6,5)) * t438;
t859 = m(6) / 0.2e1;
t891 = 0.2e1 * t859;
t861 = m(5) / 0.2e1;
t890 = 0.2e1 * t861;
t588 = sin(qJ(3));
t731 = t588 * t592;
t710 = t586 * t731;
t889 = t884 * t710 + t871;
t750 = cos(pkin(5));
t705 = pkin(1) * t750;
t496 = -pkin(7) * t738 + t592 * t705;
t497 = (pkin(2) * t589 - pkin(8) * t592) * t586;
t337 = t591 * t496 + t588 * t497;
t301 = pkin(9) * t738 + t337;
t498 = pkin(7) * t737 + t589 * t705;
t794 = pkin(9) * t591;
t557 = pkin(3) * t588 - t794;
t350 = t557 * t737 + t498;
t155 = -t587 * t301 + t590 * t350;
t123 = -pkin(4) * t710 - t155;
t746 = t123 * t587;
t156 = t590 * t301 + t587 * t350;
t122 = qJ(5) * t710 + t156;
t747 = t122 * t590;
t888 = t746 + t747;
t732 = t588 * t590;
t570 = Ifges(6,5) * t732;
t734 = t587 * t588;
t456 = -Ifges(6,6) * t591 + Ifges(6,3) * t734 + t570;
t775 = Ifges(5,6) * t591;
t582 = Ifges(5,4) * t590;
t873 = -Ifges(5,2) * t587 + t582;
t462 = t588 * t873 - t775;
t689 = t456 / 0.2e1 - t462 / 0.2e1;
t728 = t590 * t591;
t530 = -pkin(3) * t591 - pkin(9) * t588 - pkin(2);
t735 = t587 * t530;
t441 = pkin(8) * t728 + t735;
t765 = t441 * mrSges(5,3);
t887 = -t765 + t689;
t494 = t588 * t750 + t591 * t738;
t380 = t494 * t587 + t590 * t737;
t381 = t494 * t590 - t709;
t223 = pkin(4) * t381 + qJ(5) * t380;
t224 = mrSges(6,1) * t381 + mrSges(6,3) * t380;
t228 = -Ifges(5,5) * t380 - Ifges(5,6) * t381;
t229 = -Ifges(6,4) * t380 + Ifges(6,6) * t381;
t795 = pkin(8) * t590;
t677 = -qJ(5) + t795;
t408 = t591 * t677 + t735;
t706 = pkin(8) * t587 + pkin(4);
t740 = t530 * t590;
t409 = t591 * t706 - t740;
t733 = t587 * t591;
t440 = -pkin(8) * t733 + t740;
t533 = pkin(4) * t587 - qJ(5) * t590;
t645 = pkin(8) + t533;
t471 = t645 * t588;
t493 = t588 * t738 - t591 * t750;
t650 = pkin(4) * t590 + qJ(5) * t587;
t501 = t650 * t588;
t663 = t590 * mrSges(6,1) + t587 * mrSges(6,3);
t502 = t663 * t588;
t503 = t665 * t588;
t569 = Ifges(6,6) * t732;
t578 = Ifges(6,5) * t587;
t870 = t590 * Ifges(6,3) - t578;
t508 = t870 * t588;
t758 = t587 * Ifges(5,4);
t543 = t590 * Ifges(5,2) + t758;
t509 = t543 * t588;
t876 = Ifges(6,1) * t590 + t578;
t464 = -Ifges(6,4) * t591 + t588 * t876;
t550 = Ifges(5,1) * t590 - t758;
t466 = -Ifges(5,5) * t591 + t550 * t588;
t684 = -t466 / 0.4e1 - t464 / 0.4e1;
t635 = -t508 / 0.4e1 - t509 / 0.4e1 - t684;
t510 = -Ifges(6,1) * t734 + t570;
t869 = -t587 * Ifges(5,1) - t582;
t511 = t869 * t588;
t687 = t462 / 0.4e1 - t456 / 0.4e1;
t636 = t510 / 0.4e1 + t511 / 0.4e1 - t687;
t743 = t493 * qJ(5);
t447 = pkin(8) * t750 + t498;
t448 = (-pkin(2) * t592 - pkin(8) * t589 - pkin(1)) * t586;
t292 = t591 * t447 + t588 * t448;
t262 = -pkin(9) * t737 + t292;
t730 = t590 * t262;
t615 = -pkin(2) * t750 - t496;
t254 = t493 * pkin(3) - t494 * pkin(9) + t615;
t736 = t587 * t254;
t92 = t730 + t736;
t68 = t92 + t743;
t517 = -mrSges(5,1) * t591 - mrSges(5,3) * t732;
t518 = mrSges(6,1) * t591 + mrSges(6,2) * t732;
t818 = t518 / 0.2e1;
t682 = -t517 / 0.2e1 + t818;
t576 = t591 * mrSges(6,3);
t522 = -mrSges(6,2) * t734 - t576;
t816 = t522 / 0.2e1;
t723 = mrSges(5,3) * t734;
t515 = mrSges(5,2) * t591 - t723;
t820 = t515 / 0.2e1;
t683 = t820 + t816;
t91 = t254 * t590 - t262 * t587;
t70 = -pkin(4) * t493 - t91;
t766 = t440 * mrSges(5,3);
t769 = t409 * mrSges(6,2);
t770 = t408 * mrSges(6,2);
t534 = mrSges(6,1) * t587 - mrSges(6,3) * t590;
t504 = t534 * t588;
t821 = t504 / 0.2e1;
t824 = t471 / 0.2e1;
t291 = -t588 * t447 + t448 * t591;
t261 = pkin(3) * t737 - t291;
t840 = t261 / 0.2e1;
t226 = mrSges(6,1) * t380 - mrSges(6,3) * t381;
t841 = t226 / 0.2e1;
t651 = pkin(4) * t380 - qJ(5) * t381;
t93 = t261 + t651;
t848 = t93 / 0.2e1;
t886 = (-t229 / 0.4e1 - t228 / 0.4e1) * t591 + (t766 / 0.2e1 - t769 / 0.2e1 - t635) * t380 + (-t765 / 0.2e1 - t770 / 0.2e1 + t636) * t381 + (t223 * t471 + t408 * t91 + t409 * t92 + t440 * t68 + t441 * t70 + t501 * t93) * t859 + t223 * t821 + t503 * t840 + t224 * t824 + t493 * t569 / 0.4e1 + t501 * t841 + t502 * t848 + t683 * t91 + t682 * t92;
t883 = Ifges(6,2) + Ifges(5,3);
t882 = t496 * mrSges(3,2);
t536 = mrSges(4,1) * t588 + mrSges(4,2) * t591;
t879 = t615 * t536;
t764 = t493 * mrSges(6,3);
t283 = -t380 * mrSges(6,2) + t764;
t284 = -mrSges(5,2) * t493 - t380 * mrSges(5,3);
t878 = t284 + t283;
t528 = -pkin(3) - t650;
t877 = m(6) * t528 - t663;
t583 = Ifges(4,4) * t591;
t546 = -t588 * Ifges(4,2) + t583;
t774 = Ifges(6,6) * t438;
t777 = Ifges(5,6) * t438;
t781 = Ifges(5,5) * t439;
t782 = Ifges(6,4) * t439;
t872 = t774 + t782 - t777 + t781;
t754 = t590 * Ifges(6,5);
t547 = t587 * Ifges(6,1) - t754;
t868 = -t547 / 0.4e1 + t869 / 0.4e1;
t686 = t464 / 0.2e1 + t466 / 0.2e1;
t712 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t371 = Ifges(6,5) * t381;
t165 = t493 * Ifges(6,6) + Ifges(6,3) * t380 + t371;
t763 = t493 * Ifges(5,6);
t783 = Ifges(5,4) * t381;
t166 = -Ifges(5,2) * t380 + t763 + t783;
t697 = t165 / 0.2e1 - t166 / 0.2e1;
t867 = t870 / 0.4e1 + t543 / 0.4e1;
t575 = m(6) * qJ(5) + mrSges(6,3);
t866 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t865 = -mrSges(5,2) + t575;
t679 = t547 / 0.2e1 - t869 / 0.2e1;
t681 = -t870 / 0.2e1 - t543 / 0.2e1;
t864 = t681 * t587 + t679 * t590;
t863 = 0.2e1 * m(6);
t862 = t586 ^ 2;
t860 = -m(6) / 0.2e1;
t858 = -pkin(3) / 0.2e1;
t856 = -mrSges(5,1) / 0.2e1;
t855 = -mrSges(6,1) / 0.2e1;
t854 = mrSges(5,2) / 0.2e1;
t852 = -mrSges(5,3) / 0.2e1;
t851 = -mrSges(6,3) / 0.2e1;
t850 = t91 / 0.2e1;
t849 = t92 / 0.2e1;
t847 = pkin(2) * mrSges(4,1);
t846 = pkin(2) * mrSges(4,2);
t845 = pkin(4) * mrSges(6,2);
t844 = pkin(8) * mrSges(4,2);
t843 = -qJ(5) / 0.2e1;
t344 = pkin(3) * t494 + pkin(9) * t493;
t140 = -t291 * t587 + t344 * t590;
t100 = -pkin(4) * t494 - t140;
t842 = t100 / 0.2e1;
t839 = t283 / 0.2e1;
t285 = mrSges(5,1) * t493 - t381 * mrSges(5,3);
t838 = -t285 / 0.2e1;
t837 = t292 / 0.2e1;
t535 = mrSges(5,1) * t587 + mrSges(5,2) * t590;
t306 = t535 * t493;
t836 = -t306 / 0.2e1;
t835 = -t381 / 0.2e1;
t834 = t409 / 0.2e1;
t739 = t557 * t590;
t416 = -t588 * t706 - t739;
t833 = t416 / 0.2e1;
t831 = t440 / 0.2e1;
t830 = t441 / 0.2e1;
t514 = t587 * t557;
t444 = -pkin(8) * t732 + t514;
t829 = t444 / 0.2e1;
t718 = mrSges(6,2) * t728;
t757 = t588 * mrSges(6,1);
t520 = t718 - t757;
t817 = t520 / 0.2e1;
t815 = t528 / 0.2e1;
t814 = t663 / 0.2e1;
t813 = -t663 / 0.2e1;
t812 = t533 / 0.2e1;
t811 = t534 / 0.2e1;
t810 = t535 / 0.2e1;
t551 = t588 * Ifges(4,1) + t583;
t807 = t551 / 0.2e1;
t806 = -t586 / 0.2e1;
t805 = -t587 / 0.2e1;
t804 = -t587 / 0.4e1;
t802 = -t590 / 0.4e1;
t800 = t590 / 0.4e1;
t799 = -t591 / 0.2e1;
t796 = pkin(8) * t588;
t793 = t68 * mrSges(6,2);
t792 = t70 * mrSges(6,2);
t791 = t91 * mrSges(5,3);
t790 = t92 * mrSges(5,3);
t789 = mrSges(5,1) * t380;
t788 = mrSges(5,2) * t381;
t786 = Ifges(3,4) * t589;
t785 = Ifges(3,4) * t592;
t784 = Ifges(4,4) * t494;
t580 = Ifges(4,5) * t591;
t780 = Ifges(6,5) * t380;
t778 = Ifges(6,2) * t588;
t773 = Ifges(4,3) * t589;
t772 = Ifges(5,3) * t588;
t771 = qJ(5) * mrSges(6,2);
t768 = t438 * mrSges(6,2);
t767 = t439 * mrSges(6,2);
t762 = t494 * mrSges(6,1);
t761 = t494 * mrSges(4,3);
t336 = -t588 * t496 + t497 * t591;
t300 = -pkin(3) * t738 - t336;
t130 = pkin(4) * t438 - qJ(5) * t439 + t300;
t167 = Ifges(6,1) * t381 + t493 * Ifges(6,4) + t780;
t374 = Ifges(5,4) * t380;
t168 = Ifges(5,1) * t381 + t493 * Ifges(5,5) - t374;
t286 = -mrSges(6,1) * t493 + t381 * mrSges(6,2);
t756 = t588 * Ifges(4,4);
t552 = Ifges(4,1) * t591 - t756;
t652 = Ifges(6,5) * t439 + Ifges(6,3) * t438;
t616 = Ifges(6,6) * t710 + t652;
t659 = Ifges(5,4) * t439 - Ifges(5,2) * t438;
t617 = Ifges(5,6) * t710 + t659;
t629 = mrSges(6,3) * t710 - t768;
t630 = -mrSges(5,2) * t710 - t438 * mrSges(5,3);
t631 = -mrSges(6,1) * t710 + t767;
t632 = mrSges(5,1) * t710 - t439 * mrSges(5,3);
t639 = -Ifges(4,6) * t737 + t784;
t640 = mrSges(4,2) * t737 - t493 * mrSges(4,3);
t641 = -mrSges(4,1) * t737 - t761;
t655 = Ifges(3,5) * t592 - Ifges(3,6) * t589;
t664 = mrSges(6,1) * t438 - mrSges(6,3) * t439;
t666 = mrSges(5,1) * t438 + mrSges(5,2) * t439;
t667 = t788 + t789;
t482 = Ifges(4,4) * t493;
t722 = Ifges(4,5) * t737;
t669 = -t482 - t722;
t701 = t737 / 0.2e1;
t674 = t588 * t701;
t702 = -t737 / 0.2e1;
t675 = t588 * t702;
t703 = t738 / 0.2e1;
t5 = (-t616 / 0.2e1 + t617 / 0.2e1 + (-Ifges(5,6) + Ifges(6,6)) * t675) * t380 + t494 * (Ifges(4,5) * t589 + t552 * t592) * t806 + (pkin(1) * (mrSges(3,1) * t589 + mrSges(3,2) * t592) + t592 * (t773 + t592 * (-Ifges(4,6) * t588 + t580)) / 0.2e1) * t862 + t889 * t835 - (Ifges(4,5) * t494 - Ifges(4,6) * t493 - Ifges(4,3) * t737) * t738 / 0.2e1 - (t883 * t710 + t872) * t493 / 0.2e1 + (t884 * t381 + t883 * t493) * t675 + (t498 * mrSges(3,1) + Ifges(3,5) * t702 + Ifges(3,6) * t703 + t655 * t806 + t882) * t750 - t737 * t879 + ((Ifges(3,1) * t589 + t785) * t702 + (Ifges(3,2) * t592 + t786) * t703 - t291 * (mrSges(4,1) * t589 - mrSges(4,3) * t727) - t292 * (-mrSges(4,2) * t589 - mrSges(4,3) * t731) + t493 * (Ifges(4,6) * t589 + t546 * t592) / 0.2e1) * t586 - t697 * t438 + t591 * (Ifges(4,1) * t494 + t669) * t702 - t93 * t664 - t261 * t666 - t300 * t667 - t337 * t640 - t336 * t641 - t498 * (mrSges(4,1) * t493 + mrSges(4,2) * t494) - t68 * t629 - t92 * t630 - t70 * t631 - t91 * t632 + (-Ifges(4,2) * t493 + t639) * t674 - t122 * t283 - t156 * t284 - t155 * t285 - t123 * t286 - t130 * t226 - m(4) * (t291 * t336 + t292 * t337 + t498 * t615) - m(5) * (t155 * t91 + t156 * t92 + t261 * t300) - m(6) * (t122 * t68 + t123 * t70 + t130 * t93) - (t592 * (-Ifges(3,2) * t589 + t785) + t589 * (Ifges(3,1) * t592 - t786)) * t862 / 0.2e1 - (t168 + t167) * t439 / 0.2e1;
t760 = t5 * qJD(1);
t141 = t590 * t291 + t587 * t344;
t146 = -t493 * t533 + t292;
t305 = t534 * t493;
t742 = t493 * t587;
t331 = -mrSges(5,2) * t494 + mrSges(5,3) * t742;
t741 = t493 * t590;
t332 = mrSges(5,1) * t494 + mrSges(5,3) * t741;
t721 = mrSges(6,2) * t741;
t333 = -t721 - t762;
t334 = mrSges(6,2) * t742 + mrSges(6,3) * t494;
t670 = t712 * t493;
t715 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t671 = t715 * t493;
t245 = Ifges(6,4) * t494 - t493 * t876;
t246 = Ifges(5,5) * t494 - t493 * t550;
t692 = -t246 / 0.2e1 - t245 / 0.2e1;
t538 = Ifges(6,3) * t587 + t754;
t243 = Ifges(6,6) * t494 - t493 * t538;
t244 = Ifges(5,6) * t494 - t493 * t873;
t693 = t244 / 0.2e1 - t243 / 0.2e1;
t695 = t167 / 0.2e1 + t168 / 0.2e1;
t98 = qJ(5) * t494 + t141;
t6 = -t68 * t334 + t261 * t306 - t92 * t331 - t91 * t332 - t141 * t284 - t140 * t285 - t100 * t286 + t93 * t305 - t98 * t283 - t70 * t333 - t146 * t226 - t615 * (mrSges(4,1) * t494 - mrSges(4,2) * t493) - t291 * t640 + t692 * t381 + t693 * t380 - m(5) * (t140 * t91 + t141 * t92) - m(6) * (t100 * t70 + t146 * t93 + t68 * t98) + (-t712 * t380 - t715 * t381 + t639) * t494 + (-t291 * mrSges(4,3) + (t671 + t695) * t590 + (t670 + t697) * t587 + (-Ifges(4,2) + Ifges(4,1) - t883) * t494 + t669) * t493 + (-m(5) * t261 + t641 - t667 + t761) * t292;
t752 = t6 * qJD(1);
t225 = mrSges(5,1) * t381 - mrSges(5,2) * t380;
t227 = Ifges(6,3) * t381 - t780;
t230 = -Ifges(5,2) * t381 - t374;
t231 = -Ifges(6,1) * t380 + t371;
t232 = -Ifges(5,1) * t380 - t783;
t9 = t93 * t224 + t261 * t225 + (t229 / 0.2e1 + t228 / 0.2e1) * t493 + (-t793 + t232 / 0.2e1 - t790 + t231 / 0.2e1 + t697) * t381 + (t227 / 0.2e1 - t230 / 0.2e1 + t791 - t792 - t695) * t380 + (m(6) * t93 + t226) * t223 + (m(6) * t70 - t285 + t286) * t92 + (m(6) * t68 + t878) * t91;
t751 = t9 * qJD(1);
t27 = m(6) * (-t381 * t93 + t493 * t68) + t493 * t283 - t381 * t226;
t748 = qJD(1) * t27;
t745 = t155 * t587;
t744 = t156 * t590;
t729 = t590 * t663;
t726 = qJD(3) * t587;
t725 = qJD(3) * t590;
t717 = m(6) * t833;
t716 = t852 - mrSges(6,2) / 0.2e1;
t714 = -Ifges(5,5) / 0.4e1 - Ifges(6,4) / 0.4e1;
t713 = -Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1;
t711 = Ifges(6,6) / 0.4e1 - Ifges(5,6) / 0.4e1;
t708 = -t775 / 0.2e1;
t707 = t577 / 0.2e1;
t700 = -t590 * t226 / 0.2e1;
t698 = -t592 * t580 / 0.4e1;
t696 = -t166 / 0.4e1 + t165 / 0.4e1;
t694 = t168 / 0.4e1 + t167 / 0.4e1;
t691 = t284 / 0.2e1 + t839;
t690 = t286 / 0.2e1 + t838;
t457 = t588 * Ifges(6,6) + t538 * t591;
t463 = Ifges(5,6) * t588 + t591 * t873;
t688 = t457 / 0.2e1 - t463 / 0.2e1;
t465 = Ifges(6,4) * t588 + t591 * t876;
t467 = Ifges(5,5) * t588 + t550 * t591;
t685 = t465 / 0.2e1 + t467 / 0.2e1;
t539 = t587 * Ifges(5,5) + t590 * Ifges(5,6);
t541 = t587 * Ifges(6,4) - t590 * Ifges(6,6);
t676 = t539 / 0.2e1 + t541 / 0.2e1 - Ifges(4,6);
t672 = t715 * t591;
t458 = -t591 * Ifges(5,3) + t588 * t875;
t460 = -t591 * Ifges(6,2) + t588 * t874;
t545 = t591 * Ifges(4,2) + t756;
t668 = t458 / 0.2e1 + t460 / 0.2e1 - t545 / 0.2e1;
t410 = -t588 * t677 + t514;
t443 = pkin(8) * t734 + t739;
t472 = t645 * t591;
t505 = t535 * t588;
t506 = t534 * t591;
t507 = t535 * t591;
t516 = -mrSges(5,2) * t588 - mrSges(5,3) * t733;
t519 = mrSges(5,1) * t588 - mrSges(5,3) * t728;
t521 = -mrSges(6,2) * t733 + mrSges(6,3) * t588;
t593 = (t100 * t409 + t146 * t471 + t408 * t98 + t410 * t68 + t416 * t70 + t472 * t93) * t859 + t507 * t840 + t472 * t841 + t506 * t848 + t516 * t849 + t519 * t850 + t100 * t818 + t141 * t820 + t146 * t821 - t305 * t824 + t284 * t829 + t331 * t830 + t332 * t831 + t286 * t833 + t333 * t834 + t505 * t837 + t410 * t839 + t98 * t816 + t70 * t817 + t879 / 0.2e1 + (t467 / 0.4e1 + t465 / 0.4e1) * t381 + (t457 / 0.4e1 - t463 / 0.4e1) * t380 + t140 * t517 / 0.2e1 + t68 * t521 / 0.2e1 + t443 * t285 / 0.2e1 + t408 * t334 / 0.2e1;
t595 = -m(5) * (-pkin(3) * t300 + (t744 - t745) * pkin(9)) / 0.2e1 + (t888 * pkin(9) + t130 * t528) * t860 + pkin(3) * t666 / 0.2e1 + t130 * t814 + t300 * t665 / 0.2e1 - t336 * mrSges(4,1) / 0.2e1 + t337 * mrSges(4,2) / 0.2e1 - t528 * t664 / 0.2e1 + t868 * t439 + t867 * t438;
t644 = Ifges(4,2) / 0.4e1 - Ifges(4,1) / 0.4e1 + Ifges(6,2) / 0.4e1 + Ifges(5,3) / 0.4e1;
t598 = t696 * t587 + t694 * t590 - t644 * t494 + (mrSges(4,1) * t701 + t789 / 0.2e1 + t788 / 0.2e1 + m(5) * t840) * pkin(8) - t482 / 0.2e1;
t606 = -t714 * t381 + t711 * t380 + (-t244 / 0.4e1 + t243 / 0.4e1) * t587 + (t246 / 0.4e1 + t245 / 0.4e1) * t590;
t614 = t140 * t440 + t141 * t441 + t443 * t91 + t444 * t92;
t620 = -t545 / 0.4e1 + t552 / 0.4e1 + t460 / 0.4e1 + t458 / 0.4e1 - t847 / 0.2e1;
t459 = t591 * t875 + t772;
t461 = t591 * t874 + t778;
t621 = t459 / 0.4e1 + t461 / 0.4e1 - t546 / 0.4e1 - t551 / 0.4e1 + t846 / 0.2e1;
t622 = t590 * t629;
t623 = t590 * t630;
t624 = t587 * t631;
t625 = t587 * t632;
t1 = t595 + t593 + (t587 * t687 + t588 * t644 + t590 * t684 + t621) * t493 + (pkin(8) * t836 + (-t844 / 0.2e1 - t541 / 0.4e1 - t539 / 0.4e1 + t711 * t590 + t714 * t587) * t737 + t606) * t588 + (-t773 / 0.2e1 + t698) * t586 + (t292 * t796 + t614) * t861 + Ifges(4,6) * t710 + t652 * t800 + t659 * t802 + (-t756 / 0.2e1 + t620) * t494 + (-t744 / 0.2e1 + t745 / 0.2e1) * mrSges(5,3) + (-t747 / 0.2e1 - t746 / 0.2e1) * mrSges(6,2) + (-t623 / 0.2e1 - t622 / 0.2e1 - t624 / 0.2e1 + t625 / 0.2e1) * pkin(9) + (-0.3e1 / 0.4e1 * t722 + (t587 * t711 - t590 * t714) * t493 + t598) * t591 + t871 * t804;
t12 = -pkin(2) * t536 + t408 * t521 + t409 * t520 + t410 * t522 + t416 * t518 + t440 * t519 + t441 * t516 + t443 * t517 + t444 * t515 + t471 * t506 + t472 * t504 + m(5) * (t440 * t443 + t441 * t444) + m(6) * (t408 * t410 + t409 * t416 + t471 * t472) + (pkin(8) * t505 + t546 / 0.2e1 + t807 - t459 / 0.2e1 - t461 / 0.2e1 + t686 * t590 + t689 * t587) * t591 + (t552 / 0.2e1 + t685 * t590 + t688 * t587 + (m(5) * t591 * pkin(8) + t507) * pkin(8) + t668) * t588;
t649 = t1 * qJD(1) + t12 * qJD(2);
t15 = -t503 * t796 - t501 * t504 + t591 * t569 / 0.2e1 + ((-t510 / 0.2e1 - t511 / 0.2e1 + t770 + t708 - t887) * t590 + (-t508 / 0.2e1 - t509 / 0.2e1 + t769 - t672 + t686) * t587) * t588 + (-m(6) * t501 - t502) * t471 + (-m(6) * t409 + t517 - t518) * t441 + (-m(6) * t408 - t515 - t522 - t723) * t440;
t637 = -t227 / 0.4e1 + t230 / 0.4e1 + t694;
t638 = t231 / 0.4e1 + t232 / 0.4e1 + t696;
t597 = pkin(8) * t225 / 0.2e1 + (-t763 / 0.4e1 - t790 / 0.2e1 - t793 / 0.2e1 + t638) * t590 + (t791 / 0.2e1 - t792 / 0.2e1 + t714 * t493 - t637) * t587;
t607 = (-pkin(4) * t123 + qJ(5) * t122) * t860 + t122 * t851 + t123 * mrSges(6,1) / 0.2e1 + t155 * t856 + t156 * t854;
t4 = (t845 / 0.2e1 - t715) * t439 + (t771 / 0.2e1 - t712) * t438 + ((mrSges(6,3) * t843 + pkin(4) * t855 + t713) * t737 + t597) * t588 + t690 * t441 + t691 * t440 + t607 + t886;
t648 = t4 * qJD(1) - t15 * qJD(2);
t647 = -t336 * t588 + t337 * t591;
t125 = t504 * t732 - m(6) * (-t408 * t591 - t471 * t732) + t591 * t522;
t603 = (-t381 * t471 + t408 * t493 - t591 * t68 - t732 * t93) * t859 + t504 * t835 + t493 * t816 + t283 * t799;
t642 = t123 * t860 - t767 / 0.2e1;
t19 = (mrSges(6,1) * t701 + t700) * t588 + t603 + t642;
t646 = qJD(1) * t19 - qJD(2) * t125;
t643 = m(6) * t842 - t762 / 0.2e1;
t634 = t876 / 0.4e1 + t550 / 0.4e1 - t867;
t633 = t538 / 0.4e1 - t873 / 0.4e1 + t868;
t599 = pkin(8) * t810 + t633 * t587 + t634 * t590 + (mrSges(6,2) + mrSges(5,3)) * pkin(9) * (-t590 ^ 2 / 0.2e1 - t587 ^ 2 / 0.2e1);
t600 = -t893 + (t471 * t533 + t501 * t528) * t859 + t503 * t858 + t471 * t811 + t501 * t813 + t502 * t815 + t504 * t812;
t601 = (-pkin(4) * t416 + qJ(5) * t410) * t860 + pkin(4) * t817 + t521 * t843 + t410 * t851 + mrSges(6,1) * t833 + t443 * t856 + mrSges(5,2) * t829;
t604 = (t830 - t408 / 0.2e1) * mrSges(6,2) + ((-t408 + t441) * t859 - t683) * pkin(9) + t636;
t605 = (t831 + t834) * mrSges(6,2) + ((t409 + t440) * t859 + t682) * pkin(9) + t635;
t10 = (-t672 + t605) * t590 + (-t591 * t712 + t604) * t587 + (t599 + t713) * t588 + t600 + t601;
t39 = -pkin(3) * t535 - t533 * t663 + (m(6) * t533 + t534) * t528 + (-t538 / 0.2e1 + t873 / 0.2e1 + t679) * t590 + (t876 / 0.2e1 + t550 / 0.2e1 + t681) * t587;
t596 = t885 * t493 + t634 * t381 + t633 * t380 + (t223 * t528 + t533 * t93) * t859 + t225 * t858 + t223 * t813 + t261 * t810 + t224 * t815 + t226 * t812 + t93 * t811;
t602 = (-pkin(4) * t100 + qJ(5) * t98) * t860 + pkin(4) * t333 / 0.2e1 + t334 * t843 + mrSges(6,1) * t842 + t140 * t856 + t141 * t854 + t98 * t851;
t610 = (-t68 + t92) * t859 + t716 * t380 - t691;
t611 = (t70 + t91) * t859 + t716 * t381 + t690;
t612 = (t849 - t68 / 0.2e1) * mrSges(6,2) + t638;
t613 = (t850 + t70 / 0.2e1) * mrSges(6,2) + t637;
t8 = t596 + t713 * t494 + (pkin(9) * t611 + t613 + t671) * t590 + (pkin(9) * t610 + t612 + t670) * t587 + t602;
t628 = t8 * qJD(1) + t10 * qJD(2) + t39 * qJD(3);
t608 = (pkin(9) * t741 - t381 * t528 - t587 * t93) * t859 + t381 * t814 + t226 * t805;
t26 = t608 - t643 + t721;
t384 = t877 * t587;
t609 = (-t471 * t587 + (-t528 * t588 - t794) * t590) * t859 + t504 * t805;
t80 = t718 + (t855 - t729 / 0.2e1) * t588 + t717 - t609;
t627 = -qJD(1) * t26 + qJD(2) * t80 + qJD(3) * t384;
t241 = -t576 + (t735 / 0.4e1 - t441 / 0.4e1 + (t795 / 0.4e1 + t843) * t591) * t863;
t30 = t764 + (t743 / 0.2e1 + t736 / 0.4e1 + t730 / 0.4e1 - t92 / 0.4e1) * t863;
t626 = qJD(1) * t30 + qJD(2) * t241 + qJD(4) * t575;
t529 = (m(6) * pkin(9) + mrSges(6,2)) * t590;
t198 = (t735 + (-0.2e1 * qJ(5) + t795) * t591) * t859 + m(6) * t830 + t522;
t90 = t588 * t729 / 0.2e1 - t757 / 0.2e1 + t717 + t609;
t29 = (t92 + 0.2e1 * t743) * t859 + m(6) * t849 + t283;
t25 = t608 + t643;
t18 = mrSges(6,1) * t675 + t588 * t700 + t603 - t642;
t11 = t605 * t590 - t601 + t591 * t707 + t778 / 0.2e1 + t772 / 0.2e1 + t599 * t588 + t600 + t884 * t728 / 0.2e1 + (t604 + t708) * t587;
t7 = t612 * t587 + t613 * t590 - t602 + t596 + (t587 * t610 + t590 * t611) * pkin(9) + t883 * t494 / 0.2e1 - t712 * t742 - t884 * t741 / 0.2e1;
t3 = -t607 - pkin(4) * t631 / 0.2e1 + qJ(5) * t629 / 0.2e1 + t286 * t830 + t441 * t838 + t782 / 0.2e1 + t781 / 0.2e1 + t774 / 0.2e1 - t777 / 0.2e1 + t597 * t588 + t878 * t831 + t883 * t674 + t886;
t2 = t614 * t861 + t745 * t852 + t616 * t802 + t701 * t580 + t617 * t800 + Ifges(4,3) * t703 + t586 * t698 + (Ifges(4,6) * t701 - t784 / 0.2e1 + (m(5) * t837 + mrSges(4,2) * t702 + t836) * pkin(8) + t644 * t493 + t606) * t588 - t595 + t593 + mrSges(5,3) * t744 / 0.2e1 + (-t722 / 0.4e1 + t598) * t591 + (t456 * t804 + t621 + (t466 + t464) * t802 + t893) * t493 + Ifges(4,6) * t675 - pkin(9) * t625 / 0.2e1 + t620 * t494 + t888 * mrSges(6,2) / 0.2e1 + (t462 * t493 + t889) * t587 / 0.4e1 + (t541 + t539) * t710 / 0.4e1 + (t624 + t623 + t622) * pkin(9) / 0.2e1;
t13 = [-qJD(2) * t5 - qJD(3) * t6 + qJD(4) * t9 + qJD(5) * t27, t2 * qJD(3) + t3 * qJD(4) + t18 * qJD(5) - t760 + (-t882 + t122 * t522 + t123 * t518 + t130 * t504 + t155 * t517 + t156 * t515 + t300 * t505 - t408 * t768 + t409 * t767 + t471 * t664 + t872 * t799 + (-mrSges(4,1) * t591 - mrSges(3,1)) * t498 + (-t766 + t686) * t439 + t887 * t438 + (t587 * t652 / 0.2e1 + t659 * t805 + pkin(8) * t666 + t498 * mrSges(4,2) + t871 * t590 / 0.2e1) * t588 + m(4) * (-pkin(2) * t498 + pkin(8) * t647) + (t155 * t440 + t156 * t441 + t300 * t796) * t890 + (t122 * t408 + t123 * t409 + t130 * t471) * t891 + (((Ifges(4,6) - t844) * t589 + (-t846 + t807 + t583 / 0.2e1) * t592) * t591 + ((-mrSges(4,1) * pkin(8) + Ifges(4,5)) * t589 + (t440 * mrSges(5,1) - t441 * mrSges(5,2) - t409 * mrSges(6,1) + t408 * mrSges(6,3) - t847 + (t581 / 0.2e1 + t579 / 0.2e1 + t707 - t776 / 0.2e1 - Ifges(4,4) / 0.2e1) * t588 + (-Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + t713) * t591 + t668) * t592) * t588 + t655) * t586 + t647 * mrSges(4,3)) * qJD(2), -t752 + t2 * qJD(2) + t7 * qJD(4) + t25 * qJD(5) + (t98 * mrSges(6,2) + t141 * mrSges(5,3) + (t331 + t334) * pkin(9) + t693) * t725 + (t100 * mrSges(6,2) - t140 * mrSges(5,3) + (-t332 + t333) * pkin(9) - t692) * t726 + (-t291 * mrSges(4,2) + pkin(3) * t306 - t146 * t663 - t528 * t305 + (t146 * t528 + (t100 * t587 + t590 * t98) * pkin(9)) * t891 + (-t140 * t587 + t141 * t590) * pkin(9) * t890 + t676 * t494 + (-Ifges(4,5) - t864) * t493 + (-pkin(3) * t890 + t892) * t292) * qJD(3), t751 + t3 * qJD(2) + t7 * qJD(3) + (t651 * mrSges(6,2) + t865 * t91 + t866 * t92 + t228 + t229) * qJD(4) + t29 * qJD(5), qJD(2) * t18 + qJD(3) * t25 + qJD(4) * t29 + t748; qJD(3) * t1 + qJD(4) * t4 + qJD(5) * t19 + t760, qJD(3) * t12 - qJD(4) * t15 - qJD(5) * t125, t11 * qJD(4) + t90 * qJD(5) + (t410 * mrSges(6,2) + t444 * mrSges(5,3) - t688) * t725 + (t416 * mrSges(6,2) - t443 * mrSges(5,3) + t685) * t726 + t649 + (-pkin(3) * t507 + t528 * t506 + t580 + ((t516 + t521) * t590 + (-t519 + t520) * t587 + m(5) * (-t443 * t587 + t444 * t590) + m(6) * (t410 * t590 + t416 * t587)) * pkin(9) + ((-m(5) * pkin(3) + t892) * pkin(8) + t864) * t591 + (t676 + t844) * t588 + t877 * t472) * qJD(3), t11 * qJD(3) + t198 * qJD(5) + t648 + (t569 + ((-Ifges(5,6) - t771) * t590 + (t845 - t884) * t587) * t588 + t866 * t441 + t865 * t440) * qJD(4), qJD(3) * t90 + qJD(4) * t198 + t646; -qJD(2) * t1 + qJD(4) * t8 + qJD(5) * t26 + t752, qJD(4) * t10 - qJD(5) * t80 - t649, qJD(4) * t39 - qJD(5) * t384, t529 * qJD(5) + t628 + (t875 + t874 + (-m(6) * t650 - t663 - t665) * pkin(9) - t650 * mrSges(6,2)) * qJD(4), qJD(4) * t529 - t627; -qJD(2) * t4 - qJD(3) * t8 + qJD(5) * t30 - t751, -qJD(3) * t10 + qJD(5) * t241 - t648, -t628, t575 * qJD(5), t626; -qJD(2) * t19 - qJD(3) * t26 - qJD(4) * t30 - t748, qJD(3) * t80 - qJD(4) * t241 - t646, t627, -t626, 0;];
Cq = t13;

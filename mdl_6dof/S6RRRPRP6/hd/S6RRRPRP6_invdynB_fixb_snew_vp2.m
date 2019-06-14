% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 08:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:01:30
% EndTime: 2019-05-07 08:01:58
% DurationCPUTime: 25.53s
% Computational Cost: add. (410952->376), mult. (906533->475), div. (0->0), fcn. (724472->12), ass. (0->154)
t800 = Ifges(6,1) + Ifges(7,1);
t793 = Ifges(6,4) + Ifges(7,4);
t792 = Ifges(6,5) + Ifges(7,5);
t799 = Ifges(6,2) + Ifges(7,2);
t798 = Ifges(6,6) + Ifges(7,6);
t797 = Ifges(6,3) + Ifges(7,3);
t750 = sin(pkin(6));
t755 = sin(qJ(2));
t759 = cos(qJ(2));
t777 = qJD(1) * qJD(2);
t739 = (-qJDD(1) * t759 + t755 * t777) * t750;
t780 = qJD(1) * t750;
t737 = (-pkin(2) * t759 - pkin(9) * t755) * t780;
t752 = cos(pkin(6));
t746 = t752 * qJD(1) + qJD(2);
t744 = t746 ^ 2;
t745 = t752 * qJDD(1) + qJDD(2);
t779 = qJD(1) * t759;
t756 = sin(qJ(1));
t760 = cos(qJ(1));
t742 = t756 * g(1) - t760 * g(2);
t761 = qJD(1) ^ 2;
t796 = pkin(8) * t750;
t734 = qJDD(1) * pkin(1) + t761 * t796 + t742;
t743 = -t760 * g(1) - t756 * g(2);
t735 = -t761 * pkin(1) + qJDD(1) * t796 + t743;
t788 = t752 * t755;
t781 = t734 * t788 + t759 * t735;
t694 = -t744 * pkin(2) + t745 * pkin(9) + (-g(3) * t755 + t737 * t779) * t750 + t781;
t738 = (qJDD(1) * t755 + t759 * t777) * t750;
t795 = t752 * g(3);
t695 = t739 * pkin(2) - t738 * pkin(9) - t795 + (-t734 + (pkin(2) * t755 - pkin(9) * t759) * t746 * qJD(1)) * t750;
t754 = sin(qJ(3));
t758 = cos(qJ(3));
t660 = -t754 * t694 + t758 * t695;
t774 = t755 * t780;
t727 = t758 * t746 - t754 * t774;
t708 = t727 * qJD(3) + t758 * t738 + t754 * t745;
t728 = t754 * t746 + t758 * t774;
t731 = qJDD(3) + t739;
t773 = t750 * t779;
t741 = qJD(3) - t773;
t649 = (t727 * t741 - t708) * qJ(4) + (t727 * t728 + t731) * pkin(3) + t660;
t661 = t758 * t694 + t754 * t695;
t707 = -t728 * qJD(3) - t754 * t738 + t758 * t745;
t718 = t741 * pkin(3) - t728 * qJ(4);
t726 = t727 ^ 2;
t652 = -t726 * pkin(3) + t707 * qJ(4) - t741 * t718 + t661;
t749 = sin(pkin(11));
t751 = cos(pkin(11));
t715 = t749 * t727 + t751 * t728;
t643 = -0.2e1 * qJD(4) * t715 + t751 * t649 - t749 * t652;
t794 = -mrSges(6,2) - mrSges(7,2);
t790 = t750 * t755;
t789 = t750 * t759;
t787 = t752 * t759;
t710 = -g(3) * t790 + t781;
t732 = t746 * mrSges(3,1) - mrSges(3,3) * t774;
t736 = (-mrSges(3,1) * t759 + mrSges(3,2) * t755) * t780;
t714 = t751 * t727 - t749 * t728;
t644 = 0.2e1 * qJD(4) * t714 + t749 * t649 + t751 * t652;
t682 = t751 * t707 - t749 * t708;
t688 = -t714 * mrSges(5,1) + t715 * mrSges(5,2);
t700 = t741 * mrSges(5,1) - t715 * mrSges(5,3);
t689 = -t714 * pkin(4) - t715 * pkin(10);
t740 = t741 ^ 2;
t642 = -t740 * pkin(4) + t731 * pkin(10) + t714 * t689 + t644;
t709 = -g(3) * t789 + t734 * t787 - t755 * t735;
t693 = -t745 * pkin(2) - t744 * pkin(9) + t737 * t774 - t709;
t653 = -t707 * pkin(3) - t726 * qJ(4) + t728 * t718 + qJDD(4) + t693;
t683 = t749 * t707 + t751 * t708;
t647 = (-t714 * t741 - t683) * pkin(10) + (t715 * t741 - t682) * pkin(4) + t653;
t753 = sin(qJ(5));
t757 = cos(qJ(5));
t637 = -t753 * t642 + t757 * t647;
t697 = -t753 * t715 + t757 * t741;
t658 = t697 * qJD(5) + t757 * t683 + t753 * t731;
t698 = t757 * t715 + t753 * t741;
t672 = -t697 * mrSges(7,1) + t698 * mrSges(7,2);
t673 = -t697 * mrSges(6,1) + t698 * mrSges(6,2);
t713 = qJD(5) - t714;
t675 = -t713 * mrSges(6,2) + t697 * mrSges(6,3);
t681 = qJDD(5) - t682;
t634 = -0.2e1 * qJD(6) * t698 + (t697 * t713 - t658) * qJ(6) + (t697 * t698 + t681) * pkin(5) + t637;
t674 = -t713 * mrSges(7,2) + t697 * mrSges(7,3);
t776 = m(7) * t634 + t681 * mrSges(7,1) + t713 * t674;
t626 = m(6) * t637 + t681 * mrSges(6,1) + t713 * t675 + (-t672 - t673) * t698 + (-mrSges(6,3) - mrSges(7,3)) * t658 + t776;
t638 = t757 * t642 + t753 * t647;
t657 = -t698 * qJD(5) - t753 * t683 + t757 * t731;
t676 = t713 * pkin(5) - t698 * qJ(6);
t696 = t697 ^ 2;
t636 = -t696 * pkin(5) + t657 * qJ(6) + 0.2e1 * qJD(6) * t697 - t713 * t676 + t638;
t775 = m(7) * t636 + t657 * mrSges(7,3) + t697 * t672;
t677 = t713 * mrSges(7,1) - t698 * mrSges(7,3);
t782 = -t713 * mrSges(6,1) + t698 * mrSges(6,3) - t677;
t631 = m(6) * t638 + t657 * mrSges(6,3) + t697 * t673 + t794 * t681 + t782 * t713 + t775;
t769 = -t753 * t626 + t757 * t631;
t622 = m(5) * t644 - t731 * mrSges(5,2) + t682 * mrSges(5,3) + t714 * t688 - t741 * t700 + t769;
t699 = -t741 * mrSges(5,2) + t714 * mrSges(5,3);
t641 = -t731 * pkin(4) - t740 * pkin(10) + t715 * t689 - t643;
t639 = -t657 * pkin(5) - t696 * qJ(6) + t698 * t676 + qJDD(6) + t641;
t767 = m(7) * t639 - t657 * mrSges(7,1) - t697 * t674;
t763 = -m(6) * t641 + t657 * mrSges(6,1) + t794 * t658 + t697 * t675 + t782 * t698 - t767;
t628 = m(5) * t643 + t731 * mrSges(5,1) - t683 * mrSges(5,3) - t715 * t688 + t741 * t699 + t763;
t615 = t749 * t622 + t751 * t628;
t716 = -t727 * mrSges(4,1) + t728 * mrSges(4,2);
t717 = -t741 * mrSges(4,2) + t727 * mrSges(4,3);
t613 = m(4) * t660 + t731 * mrSges(4,1) - t708 * mrSges(4,3) - t728 * t716 + t741 * t717 + t615;
t719 = t741 * mrSges(4,1) - t728 * mrSges(4,3);
t770 = t751 * t622 - t749 * t628;
t614 = m(4) * t661 - t731 * mrSges(4,2) + t707 * mrSges(4,3) + t727 * t716 - t741 * t719 + t770;
t771 = -t754 * t613 + t758 * t614;
t605 = m(3) * t710 - t745 * mrSges(3,2) - t739 * mrSges(3,3) - t746 * t732 + t736 * t773 + t771;
t608 = t758 * t613 + t754 * t614;
t723 = -t750 * t734 - t795;
t733 = -t746 * mrSges(3,2) + mrSges(3,3) * t773;
t607 = m(3) * t723 + t739 * mrSges(3,1) + t738 * mrSges(3,2) + (t732 * t755 - t733 * t759) * t780 + t608;
t624 = t757 * t626 + t753 * t631;
t764 = m(5) * t653 - t682 * mrSges(5,1) + t683 * mrSges(5,2) - t714 * t699 + t715 * t700 + t624;
t762 = -m(4) * t693 + t707 * mrSges(4,1) - t708 * mrSges(4,2) + t727 * t717 - t728 * t719 - t764;
t619 = m(3) * t709 + t745 * mrSges(3,1) - t738 * mrSges(3,3) + t746 * t733 - t736 * t774 + t762;
t595 = t605 * t788 - t750 * t607 + t619 * t787;
t593 = m(2) * t742 + qJDD(1) * mrSges(2,1) - t761 * mrSges(2,2) + t595;
t600 = t759 * t605 - t755 * t619;
t599 = m(2) * t743 - t761 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t600;
t786 = t760 * t593 + t756 * t599;
t785 = t697 * t798 + t698 * t792 + t713 * t797;
t784 = -t697 * t799 - t698 * t793 - t713 * t798;
t783 = t793 * t697 + t698 * t800 + t792 * t713;
t594 = t605 * t790 + t752 * t607 + t619 * t789;
t772 = -t756 * t593 + t760 * t599;
t616 = -mrSges(6,1) * t641 + mrSges(6,3) * t638 - mrSges(7,1) * t639 + mrSges(7,3) * t636 - pkin(5) * t767 + qJ(6) * t775 + (-qJ(6) * t677 + t783) * t713 + (-pkin(5) * t677 - t785) * t698 + (-qJ(6) * mrSges(7,2) + t798) * t681 + (-pkin(5) * mrSges(7,2) + t793) * t658 + t799 * t657;
t632 = -t658 * mrSges(7,3) - t698 * t672 + t776;
t623 = mrSges(6,2) * t641 + mrSges(7,2) * t639 - mrSges(6,3) * t637 - mrSges(7,3) * t634 - qJ(6) * t632 + t793 * t657 + t658 * t800 + t792 * t681 + t785 * t697 + t784 * t713;
t684 = Ifges(5,5) * t715 + Ifges(5,6) * t714 + Ifges(5,3) * t741;
t685 = Ifges(5,4) * t715 + Ifges(5,2) * t714 + Ifges(5,6) * t741;
t601 = mrSges(5,2) * t653 - mrSges(5,3) * t643 + Ifges(5,1) * t683 + Ifges(5,4) * t682 + Ifges(5,5) * t731 - pkin(10) * t624 - t753 * t616 + t757 * t623 + t714 * t684 - t741 * t685;
t686 = Ifges(5,1) * t715 + Ifges(5,4) * t714 + Ifges(5,5) * t741;
t609 = -mrSges(5,1) * t653 - mrSges(6,1) * t637 - mrSges(7,1) * t634 + mrSges(6,2) * t638 + mrSges(7,2) * t636 + mrSges(5,3) * t644 + Ifges(5,4) * t683 + Ifges(5,2) * t682 + Ifges(5,6) * t731 - pkin(4) * t624 - pkin(5) * t632 - t715 * t684 + t741 * t686 + t784 * t698 + t783 * t697 - t797 * t681 - t792 * t658 - t798 * t657;
t701 = Ifges(4,5) * t728 + Ifges(4,6) * t727 + Ifges(4,3) * t741;
t703 = Ifges(4,1) * t728 + Ifges(4,4) * t727 + Ifges(4,5) * t741;
t591 = -mrSges(4,1) * t693 + mrSges(4,3) * t661 + Ifges(4,4) * t708 + Ifges(4,2) * t707 + Ifges(4,6) * t731 - pkin(3) * t764 + qJ(4) * t770 + t749 * t601 + t751 * t609 - t728 * t701 + t741 * t703;
t702 = Ifges(4,4) * t728 + Ifges(4,2) * t727 + Ifges(4,6) * t741;
t596 = mrSges(4,2) * t693 - mrSges(4,3) * t660 + Ifges(4,1) * t708 + Ifges(4,4) * t707 + Ifges(4,5) * t731 - qJ(4) * t615 + t751 * t601 - t749 * t609 + t727 * t701 - t741 * t702;
t720 = Ifges(3,3) * t746 + (Ifges(3,5) * t755 + Ifges(3,6) * t759) * t780;
t721 = Ifges(3,6) * t746 + (Ifges(3,4) * t755 + Ifges(3,2) * t759) * t780;
t589 = mrSges(3,2) * t723 - mrSges(3,3) * t709 + Ifges(3,1) * t738 - Ifges(3,4) * t739 + Ifges(3,5) * t745 - pkin(9) * t608 - t754 * t591 + t758 * t596 + t720 * t773 - t746 * t721;
t722 = Ifges(3,5) * t746 + (Ifges(3,1) * t755 + Ifges(3,4) * t759) * t780;
t590 = (-Ifges(4,3) - Ifges(5,3)) * t731 - pkin(2) * t608 - t757 * t616 - t753 * t623 + Ifges(3,6) * t745 + t746 * t722 + Ifges(3,4) * t738 - Ifges(3,2) * t739 - t728 * t702 - mrSges(3,1) * t723 + t727 * t703 + t714 * t686 - t715 * t685 - Ifges(4,6) * t707 - Ifges(4,5) * t708 + mrSges(3,3) * t710 - Ifges(5,6) * t682 - Ifges(5,5) * t683 - mrSges(4,1) * t660 + mrSges(4,2) * t661 + mrSges(5,2) * t644 - mrSges(5,1) * t643 - t720 * t774 - pkin(10) * t769 - pkin(3) * t615 - pkin(4) * t763;
t765 = pkin(8) * t600 + t589 * t755 + t590 * t759;
t588 = Ifges(3,5) * t738 - Ifges(3,6) * t739 + Ifges(3,3) * t745 + mrSges(3,1) * t709 - mrSges(3,2) * t710 + t754 * t596 + t758 * t591 + pkin(2) * t762 + pkin(9) * t771 + (t721 * t755 - t722 * t759) * t780;
t587 = -mrSges(2,2) * g(3) - mrSges(2,3) * t742 + Ifges(2,5) * qJDD(1) - t761 * Ifges(2,6) + t759 * t589 - t755 * t590 + (-t594 * t750 - t595 * t752) * pkin(8);
t586 = mrSges(2,1) * g(3) + mrSges(2,3) * t743 + t761 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t594 - t750 * t588 + t765 * t752;
t1 = [-m(1) * g(1) + t772; -m(1) * g(2) + t786; (-m(1) - m(2)) * g(3) + t594; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t786 - t756 * t586 + t760 * t587; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t772 + t760 * t586 + t756 * t587; -mrSges(1,1) * g(2) + mrSges(2,1) * t742 + mrSges(1,2) * g(1) - mrSges(2,2) * t743 + Ifges(2,3) * qJDD(1) + pkin(1) * t595 + t752 * t588 + t765 * t750;];
tauB  = t1;

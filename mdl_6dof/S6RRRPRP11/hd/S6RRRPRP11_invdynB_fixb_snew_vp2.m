% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:08:36
% EndTime: 2019-05-07 09:08:53
% DurationCPUTime: 9.64s
% Computational Cost: add. (145862->364), mult. (312851->432), div. (0->0), fcn. (237886->10), ass. (0->150)
t839 = Ifges(4,1) + Ifges(5,2);
t838 = Ifges(6,1) + Ifges(7,1);
t827 = Ifges(6,4) + Ifges(7,4);
t826 = Ifges(4,5) - Ifges(5,4);
t825 = Ifges(6,5) + Ifges(7,5);
t837 = Ifges(4,2) + Ifges(5,3);
t836 = Ifges(6,2) + Ifges(7,2);
t824 = Ifges(4,6) - Ifges(5,5);
t823 = -Ifges(5,6) - Ifges(4,4);
t822 = Ifges(6,6) + Ifges(7,6);
t835 = Ifges(4,3) + Ifges(5,1);
t834 = Ifges(6,3) + Ifges(7,3);
t777 = cos(pkin(6));
t773 = qJD(1) * t777 + qJD(2);
t779 = sin(qJ(3));
t780 = sin(qJ(2));
t776 = sin(pkin(6));
t805 = qJD(1) * t776;
t798 = t780 * t805;
t831 = cos(qJ(3));
t749 = t779 * t773 + t798 * t831;
t783 = cos(qJ(2));
t802 = qJD(1) * qJD(2);
t761 = (qJDD(1) * t780 + t783 * t802) * t776;
t772 = qJDD(1) * t777 + qJDD(2);
t721 = qJD(3) * t749 + t761 * t779 - t772 * t831;
t748 = -t773 * t831 + t779 * t798;
t804 = qJD(1) * t783;
t797 = t776 * t804;
t767 = -qJD(3) + t797;
t778 = sin(qJ(5));
t782 = cos(qJ(5));
t732 = t748 * t778 - t767 * t782;
t762 = (-qJDD(1) * t783 + t780 * t802) * t776;
t754 = qJDD(3) + t762;
t680 = -qJD(5) * t732 + t721 * t782 - t754 * t778;
t833 = (mrSges(6,1) + mrSges(7,1)) * t680;
t760 = (-pkin(2) * t783 - pkin(9) * t780) * t805;
t771 = t773 ^ 2;
t781 = sin(qJ(1));
t784 = cos(qJ(1));
t769 = t781 * g(1) - g(2) * t784;
t785 = qJD(1) ^ 2;
t830 = pkin(8) * t776;
t757 = qJDD(1) * pkin(1) + t785 * t830 + t769;
t770 = -g(1) * t784 - g(2) * t781;
t758 = -pkin(1) * t785 + qJDD(1) * t830 + t770;
t817 = t777 * t780;
t806 = t757 * t817 + t783 * t758;
t695 = -t771 * pkin(2) + t772 * pkin(9) + (-g(3) * t780 + t760 * t804) * t776 + t806;
t829 = t777 * g(3);
t696 = t762 * pkin(2) - t761 * pkin(9) - t829 + (-t757 + (pkin(2) * t780 - pkin(9) * t783) * t773 * qJD(1)) * t776;
t676 = t831 * t695 + t779 * t696;
t725 = pkin(3) * t748 - qJ(4) * t749;
t766 = t767 ^ 2;
t672 = pkin(3) * t766 - t754 * qJ(4) + 0.2e1 * qJD(4) * t767 + t748 * t725 - t676;
t731 = t748 * t782 + t767 * t778;
t746 = qJD(5) + t749;
t703 = -mrSges(7,2) * t746 + mrSges(7,3) * t731;
t704 = -mrSges(6,2) * t746 + mrSges(6,3) * t731;
t832 = -t833 - (t703 + t704) * t731;
t821 = t731 * t703;
t820 = t748 * t767;
t819 = t776 * t780;
t818 = t776 * t783;
t816 = t777 * t783;
t724 = -g(3) * t819 + t806;
t755 = mrSges(3,1) * t773 - mrSges(3,3) * t798;
t759 = (-mrSges(3,1) * t783 + mrSges(3,2) * t780) * t805;
t675 = -t779 * t695 + t696 * t831;
t722 = -t748 * qJD(3) + t761 * t831 + t779 * t772;
t726 = mrSges(4,1) * t748 + mrSges(4,2) * t749;
t733 = mrSges(4,2) * t767 - mrSges(4,3) * t748;
t735 = mrSges(5,1) * t748 + mrSges(5,3) * t767;
t673 = -t754 * pkin(3) - t766 * qJ(4) + t749 * t725 + qJDD(4) - t675;
t667 = (t748 * t749 - t754) * pkin(10) + (t722 - t820) * pkin(4) + t673;
t737 = pkin(4) * t749 + pkin(10) * t767;
t747 = t748 ^ 2;
t723 = -g(3) * t818 + t757 * t816 - t780 * t758;
t694 = -t772 * pkin(2) - t771 * pkin(9) + t760 * t798 - t723;
t786 = (-t722 - t820) * qJ(4) + t694 + (-t767 * pkin(3) - 0.2e1 * qJD(4)) * t749;
t671 = -t747 * pkin(4) - t749 * t737 + (pkin(3) + pkin(10)) * t721 + t786;
t661 = t782 * t667 - t778 * t671;
t681 = qJD(5) * t731 + t721 * t778 + t754 * t782;
t698 = -mrSges(7,1) * t731 + mrSges(7,2) * t732;
t699 = -mrSges(6,1) * t731 + mrSges(6,2) * t732;
t718 = qJDD(5) + t722;
t658 = -0.2e1 * qJD(6) * t732 + (t731 * t746 - t681) * qJ(6) + (t731 * t732 + t718) * pkin(5) + t661;
t801 = m(7) * t658 + t718 * mrSges(7,1) + t746 * t703;
t653 = m(6) * t661 + t718 * mrSges(6,1) + t746 * t704 + (-t698 - t699) * t732 + (-mrSges(6,3) - mrSges(7,3)) * t681 + t801;
t662 = t778 * t667 + t782 * t671;
t706 = mrSges(7,1) * t746 - mrSges(7,3) * t732;
t707 = mrSges(6,1) * t746 - mrSges(6,3) * t732;
t705 = pkin(5) * t746 - qJ(6) * t732;
t730 = t731 ^ 2;
t660 = -pkin(5) * t730 + qJ(6) * t680 + 0.2e1 * qJD(6) * t731 - t705 * t746 + t662;
t800 = m(7) * t660 + t680 * mrSges(7,3) + t731 * t698;
t655 = m(6) * t662 + t680 * mrSges(6,3) + t731 * t699 + (-t706 - t707) * t746 + (-mrSges(6,2) - mrSges(7,2)) * t718 + t800;
t648 = t782 * t653 + t778 * t655;
t727 = -mrSges(5,2) * t748 - mrSges(5,3) * t749;
t789 = -m(5) * t673 - t722 * mrSges(5,1) - t749 * t727 - t648;
t645 = m(4) * t675 - t722 * mrSges(4,3) - t749 * t726 + (-t733 + t735) * t767 + (mrSges(4,1) - mrSges(5,2)) * t754 + t789;
t734 = -mrSges(4,1) * t767 - mrSges(4,3) * t749;
t736 = mrSges(5,1) * t749 - mrSges(5,2) * t767;
t670 = -pkin(4) * t721 - pkin(10) * t747 - t767 * t737 - t672;
t664 = -pkin(5) * t680 - qJ(6) * t730 + t705 * t732 + qJDD(6) + t670;
t799 = m(7) * t664 + t681 * mrSges(7,2) + t732 * t706;
t794 = -m(6) * t670 - t681 * mrSges(6,2) - t732 * t707 - t799;
t788 = -m(5) * t672 + t754 * mrSges(5,3) - t767 * t736 - t794;
t651 = t788 + (-mrSges(4,3) - mrSges(5,1)) * t721 + (-t726 - t727) * t748 + m(4) * t676 - t754 * mrSges(4,2) + t767 * t734 + t832;
t795 = -t645 * t779 + t831 * t651;
t636 = m(3) * t724 - mrSges(3,2) * t772 - mrSges(3,3) * t762 - t755 * t773 + t759 * t797 + t795;
t639 = t831 * t645 + t779 * t651;
t742 = -t776 * t757 - t829;
t756 = -mrSges(3,2) * t773 + mrSges(3,3) * t797;
t638 = m(3) * t742 + t762 * mrSges(3,1) + t761 * mrSges(3,2) + (t755 * t780 - t756 * t783) * t805 + t639;
t674 = t721 * pkin(3) + t786;
t814 = -t778 * t653 + t782 * t655;
t793 = -m(5) * t674 + t721 * mrSges(5,2) + t748 * t735 - t814;
t787 = -m(4) * t694 - t721 * mrSges(4,1) - t748 * t733 + (-t734 + t736) * t749 + (-mrSges(4,2) + mrSges(5,3)) * t722 + t793;
t643 = m(3) * t723 + t772 * mrSges(3,1) - t761 * mrSges(3,3) + t773 * t756 - t759 * t798 + t787;
t627 = t636 * t817 - t638 * t776 + t643 * t816;
t625 = m(2) * t769 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t785 + t627;
t632 = t783 * t636 - t643 * t780;
t631 = m(2) * t770 - mrSges(2,1) * t785 - qJDD(1) * mrSges(2,2) + t632;
t815 = t784 * t625 + t781 * t631;
t813 = -t731 * t822 - t732 * t825 - t746 * t834;
t812 = t731 * t836 + t732 * t827 + t746 * t822;
t811 = -t731 * t827 - t732 * t838 - t746 * t825;
t809 = t748 * t824 - t749 * t826 + t767 * t835;
t808 = t748 * t837 + t749 * t823 + t767 * t824;
t807 = -t823 * t748 - t749 * t839 + t826 * t767;
t626 = t636 * t819 + t777 * t638 + t643 * t818;
t796 = -t625 * t781 + t784 * t631;
t640 = -mrSges(6,1) * t670 + mrSges(6,3) * t662 - mrSges(7,1) * t664 + mrSges(7,3) * t660 - pkin(5) * (t799 - t821) + qJ(6) * t800 + (-qJ(6) * t706 - t811) * t746 + t813 * t732 + (-mrSges(7,2) * qJ(6) + t822) * t718 + t827 * t681 + (mrSges(7,1) * pkin(5) + t836) * t680;
t646 = -t722 * mrSges(5,3) - t749 * t736 - t793;
t656 = -t681 * mrSges(7,3) - t732 * t698 + t801;
t647 = mrSges(6,2) * t670 + mrSges(7,2) * t664 - mrSges(6,3) * t661 - mrSges(7,3) * t658 - qJ(6) * t656 + t827 * t680 + t681 * t838 + t825 * t718 - t813 * t731 - t812 * t746;
t623 = -mrSges(4,1) * t694 + mrSges(4,3) * t676 - mrSges(5,1) * t672 + mrSges(5,2) * t674 - t778 * t647 - t782 * t640 - pkin(4) * (t794 - t832) - pkin(10) * t814 - pkin(3) * t646 + t807 * t767 + t824 * t754 + t809 * t749 - t823 * t722 - t837 * t721;
t628 = mrSges(5,1) * t673 + mrSges(6,1) * t661 + mrSges(7,1) * t658 + mrSges(4,2) * t694 - mrSges(6,2) * t662 - mrSges(7,2) * t660 - mrSges(4,3) * t675 - mrSges(5,3) * t674 + pkin(4) * t648 + pkin(5) * t656 - qJ(4) * t646 - t808 * t767 + t826 * t754 + t809 * t748 + t812 * t732 + t811 * t731 + t839 * t722 + t823 * t721 + t834 * t718 + t825 * t681 + t822 * t680;
t739 = Ifges(3,3) * t773 + (Ifges(3,5) * t780 + Ifges(3,6) * t783) * t805;
t740 = Ifges(3,6) * t773 + (Ifges(3,4) * t780 + Ifges(3,2) * t783) * t805;
t621 = mrSges(3,2) * t742 - mrSges(3,3) * t723 + Ifges(3,1) * t761 - Ifges(3,4) * t762 + Ifges(3,5) * t772 - pkin(9) * t639 - t779 * t623 + t628 * t831 + t739 * t797 - t773 * t740;
t741 = Ifges(3,5) * t773 + (Ifges(3,1) * t780 + Ifges(3,4) * t783) * t805;
t622 = t778 * t640 - pkin(3) * (t767 * t735 + t789) - t782 * t647 + pkin(10) * t648 - pkin(2) * t639 + mrSges(5,3) * t672 - mrSges(5,2) * t673 - mrSges(4,1) * t675 + mrSges(4,2) * t676 + Ifges(3,6) * t772 + t773 * t741 + Ifges(3,4) * t761 - Ifges(3,2) * t762 + mrSges(3,3) * t724 - qJ(4) * (-t731 * t704 + t788 - t821 - t833) - mrSges(3,1) * t742 - t739 * t798 + (mrSges(5,2) * pkin(3) - t835) * t754 + t808 * t749 + (qJ(4) * t727 + t807) * t748 - t826 * t722 + (mrSges(5,1) * qJ(4) + t824) * t721;
t790 = pkin(8) * t632 + t621 * t780 + t622 * t783;
t620 = Ifges(3,5) * t761 - Ifges(3,6) * t762 + Ifges(3,3) * t772 + mrSges(3,1) * t723 - mrSges(3,2) * t724 + t779 * t628 + t831 * t623 + pkin(2) * t787 + pkin(9) * t795 + (t740 * t780 - t741 * t783) * t805;
t619 = -mrSges(2,2) * g(3) - mrSges(2,3) * t769 + Ifges(2,5) * qJDD(1) - t785 * Ifges(2,6) + t783 * t621 - t780 * t622 + (-t626 * t776 - t627 * t777) * pkin(8);
t618 = mrSges(2,1) * g(3) + mrSges(2,3) * t770 + t785 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t626 - t776 * t620 + t777 * t790;
t1 = [-m(1) * g(1) + t796; -m(1) * g(2) + t815; (-m(1) - m(2)) * g(3) + t626; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t815 - t781 * t618 + t784 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t796 + t784 * t618 + t781 * t619; -mrSges(1,1) * g(2) + mrSges(2,1) * t769 + mrSges(1,2) * g(1) - mrSges(2,2) * t770 + Ifges(2,3) * qJDD(1) + pkin(1) * t627 + t777 * t620 + t776 * t790;];
tauB  = t1;

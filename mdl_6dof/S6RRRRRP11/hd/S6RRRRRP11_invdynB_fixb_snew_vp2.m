% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 06:45:50
% EndTime: 2019-05-08 06:47:17
% DurationCPUTime: 66.58s
% Computational Cost: add. (1096963->395), mult. (2706003->513), div. (0->0), fcn. (2276972->14), ass. (0->172)
t842 = Ifges(6,1) + Ifges(7,1);
t833 = Ifges(6,4) + Ifges(7,4);
t832 = Ifges(6,5) + Ifges(7,5);
t841 = Ifges(6,2) + Ifges(7,2);
t840 = Ifges(6,6) + Ifges(7,6);
t839 = Ifges(6,3) + Ifges(7,3);
t783 = cos(pkin(6));
t778 = t783 * qJD(1) + qJD(2);
t780 = sin(pkin(7));
t782 = cos(pkin(7));
t781 = sin(pkin(6));
t792 = cos(qJ(2));
t814 = qJD(1) * t792;
t808 = t781 * t814;
t763 = (t778 * t780 + t782 * t808) * pkin(10);
t787 = sin(qJ(2));
t816 = qJD(1) * t781;
t837 = pkin(10) * t780;
t767 = (-pkin(2) * t792 - t787 * t837) * t816;
t813 = qJD(1) * qJD(2);
t773 = (qJDD(1) * t787 + t792 * t813) * t781;
t777 = t783 * qJDD(1) + qJDD(2);
t788 = sin(qJ(1));
t793 = cos(qJ(1));
t775 = t788 * g(1) - t793 * g(2);
t794 = qJD(1) ^ 2;
t838 = pkin(9) * t781;
t770 = qJDD(1) * pkin(1) + t794 * t838 + t775;
t776 = -t793 * g(1) - t788 * g(2);
t771 = -t794 * pkin(1) + qJDD(1) * t838 + t776;
t823 = t783 * t792;
t804 = t770 * t823 - t787 * t771;
t815 = qJD(1) * t787;
t836 = pkin(10) * t782;
t723 = -t773 * t836 + t777 * pkin(2) + t778 * t763 + (-g(3) * t792 - t767 * t815) * t781 + t804;
t809 = t781 * t815;
t766 = t778 * pkin(2) - t809 * t836;
t774 = (qJDD(1) * t792 - t787 * t813) * t781;
t801 = t774 * t782 + t777 * t780;
t824 = t783 * t787;
t817 = t770 * t824 + t792 * t771;
t724 = -t778 * t766 + (-g(3) * t787 + t767 * t814) * t781 + t801 * pkin(10) + t817;
t835 = t783 * g(3);
t729 = -t773 * t837 - t774 * pkin(2) - t835 + (-t770 + (-t763 * t792 + t766 * t787) * qJD(1)) * t781;
t786 = sin(qJ(3));
t791 = cos(qJ(3));
t690 = -t786 * t724 + (t723 * t782 + t729 * t780) * t791;
t825 = t782 * t792;
t830 = t780 * t786;
t754 = t778 * t830 + (t786 * t825 + t787 * t791) * t816;
t740 = -t754 * qJD(3) - t786 * t773 + t801 * t791;
t829 = t780 * t791;
t753 = (-t786 * t787 + t791 * t825) * t816 + t778 * t829;
t834 = -mrSges(6,2) - mrSges(7,2);
t828 = t781 * t787;
t827 = t781 * t792;
t826 = t782 * t786;
t691 = t723 * t826 + t791 * t724 + t729 * t830;
t742 = -t753 * mrSges(4,1) + t754 * mrSges(4,2);
t764 = t782 * t778 - t780 * t808 + qJD(3);
t748 = t764 * mrSges(4,1) - t754 * mrSges(4,3);
t755 = -t780 * t774 + t782 * t777 + qJDD(3);
t743 = -t753 * pkin(3) - t754 * pkin(11);
t762 = t764 ^ 2;
t682 = -t762 * pkin(3) + t755 * pkin(11) + t753 * t743 + t691;
t701 = -t780 * t723 + t782 * t729;
t741 = t753 * qJD(3) + t791 * t773 + t801 * t786;
t684 = (-t753 * t764 - t741) * pkin(11) + (t754 * t764 - t740) * pkin(3) + t701;
t785 = sin(qJ(4));
t790 = cos(qJ(4));
t678 = t790 * t682 + t785 * t684;
t746 = t790 * t754 + t785 * t764;
t708 = -t746 * qJD(4) - t785 * t741 + t790 * t755;
t745 = -t785 * t754 + t790 * t764;
t725 = -t745 * mrSges(5,1) + t746 * mrSges(5,2);
t752 = qJD(4) - t753;
t735 = t752 * mrSges(5,1) - t746 * mrSges(5,3);
t739 = qJDD(4) - t740;
t726 = -t745 * pkin(4) - t746 * pkin(12);
t751 = t752 ^ 2;
t673 = -t751 * pkin(4) + t739 * pkin(12) + t745 * t726 + t678;
t681 = -t755 * pkin(3) - t762 * pkin(11) + t754 * t743 - t690;
t709 = t745 * qJD(4) + t790 * t741 + t785 * t755;
t676 = (-t745 * t752 - t709) * pkin(12) + (t746 * t752 - t708) * pkin(4) + t681;
t784 = sin(qJ(5));
t789 = cos(qJ(5));
t668 = -t784 * t673 + t789 * t676;
t732 = -t784 * t746 + t789 * t752;
t689 = t732 * qJD(5) + t789 * t709 + t784 * t739;
t733 = t789 * t746 + t784 * t752;
t703 = -t732 * mrSges(7,1) + t733 * mrSges(7,2);
t704 = -t732 * mrSges(6,1) + t733 * mrSges(6,2);
t707 = qJDD(5) - t708;
t744 = qJD(5) - t745;
t712 = -t744 * mrSges(6,2) + t732 * mrSges(6,3);
t665 = -0.2e1 * qJD(6) * t733 + (t732 * t744 - t689) * qJ(6) + (t732 * t733 + t707) * pkin(5) + t668;
t711 = -t744 * mrSges(7,2) + t732 * mrSges(7,3);
t811 = m(7) * t665 + t707 * mrSges(7,1) + t744 * t711;
t658 = m(6) * t668 + t707 * mrSges(6,1) + t744 * t712 + (-t703 - t704) * t733 + (-mrSges(6,3) - mrSges(7,3)) * t689 + t811;
t669 = t789 * t673 + t784 * t676;
t688 = -t733 * qJD(5) - t784 * t709 + t789 * t739;
t713 = t744 * pkin(5) - t733 * qJ(6);
t731 = t732 ^ 2;
t667 = -t731 * pkin(5) + t688 * qJ(6) + 0.2e1 * qJD(6) * t732 - t744 * t713 + t669;
t810 = m(7) * t667 + t688 * mrSges(7,3) + t732 * t703;
t714 = t744 * mrSges(7,1) - t733 * mrSges(7,3);
t818 = -t744 * mrSges(6,1) + t733 * mrSges(6,3) - t714;
t660 = m(6) * t669 + t688 * mrSges(6,3) + t732 * t704 + t834 * t707 + t818 * t744 + t810;
t805 = -t784 * t658 + t789 * t660;
t655 = m(5) * t678 - t739 * mrSges(5,2) + t708 * mrSges(5,3) + t745 * t725 - t752 * t735 + t805;
t677 = -t785 * t682 + t790 * t684;
t734 = -t752 * mrSges(5,2) + t745 * mrSges(5,3);
t672 = -t739 * pkin(4) - t751 * pkin(12) + t746 * t726 - t677;
t670 = -t688 * pkin(5) - t731 * qJ(6) + t733 * t713 + qJDD(6) + t672;
t803 = m(7) * t670 - t688 * mrSges(7,1) - t732 * t711;
t795 = -m(6) * t672 + t688 * mrSges(6,1) + t834 * t689 + t732 * t712 + t818 * t733 - t803;
t662 = m(5) * t677 + t739 * mrSges(5,1) - t709 * mrSges(5,3) - t746 * t725 + t752 * t734 + t795;
t806 = t790 * t655 - t785 * t662;
t645 = m(4) * t691 - t755 * mrSges(4,2) + t740 * mrSges(4,3) + t753 * t742 - t764 * t748 + t806;
t648 = t785 * t655 + t790 * t662;
t747 = -t764 * mrSges(4,2) + t753 * mrSges(4,3);
t647 = m(4) * t701 - t740 * mrSges(4,1) + t741 * mrSges(4,2) - t753 * t747 + t754 * t748 + t648;
t657 = t789 * t658 + t784 * t660;
t796 = -m(5) * t681 + t708 * mrSges(5,1) - t709 * mrSges(5,2) + t745 * t734 - t746 * t735 - t657;
t652 = m(4) * t690 + t755 * mrSges(4,1) - t741 * mrSges(4,3) - t754 * t742 + t764 * t747 + t796;
t634 = t782 * t791 * t652 + t645 * t826 - t780 * t647;
t749 = -g(3) * t827 + t804;
t769 = -t778 * mrSges(3,2) + mrSges(3,3) * t808;
t772 = (-mrSges(3,1) * t792 + mrSges(3,2) * t787) * t816;
t630 = m(3) * t749 + t777 * mrSges(3,1) - t773 * mrSges(3,3) + t778 * t769 - t772 * t809 + t634;
t633 = t645 * t830 + t782 * t647 + t652 * t829;
t759 = -t781 * t770 - t835;
t768 = t778 * mrSges(3,1) - mrSges(3,3) * t809;
t632 = m(3) * t759 - t774 * mrSges(3,1) + t773 * mrSges(3,2) + (t768 * t787 - t769 * t792) * t816 + t633;
t640 = t791 * t645 - t786 * t652;
t750 = -g(3) * t828 + t817;
t639 = m(3) * t750 - t777 * mrSges(3,2) + t774 * mrSges(3,3) - t778 * t768 + t772 * t808 + t640;
t620 = t630 * t823 - t781 * t632 + t639 * t824;
t618 = m(2) * t775 + qJDD(1) * mrSges(2,1) - t794 * mrSges(2,2) + t620;
t626 = -t787 * t630 + t792 * t639;
t625 = m(2) * t776 - t794 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t626;
t822 = t793 * t618 + t788 * t625;
t821 = t840 * t732 + t832 * t733 + t839 * t744;
t820 = -t841 * t732 - t833 * t733 - t840 * t744;
t819 = t833 * t732 + t842 * t733 + t832 * t744;
t619 = t630 * t827 + t783 * t632 + t639 * t828;
t807 = -t788 * t618 + t793 * t625;
t649 = -mrSges(6,1) * t672 + mrSges(6,3) * t669 - mrSges(7,1) * t670 + mrSges(7,3) * t667 - pkin(5) * t803 + qJ(6) * t810 + (-qJ(6) * t714 + t819) * t744 + (-pkin(5) * t714 - t821) * t733 + (-qJ(6) * mrSges(7,2) + t840) * t707 + (-pkin(5) * mrSges(7,2) + t833) * t689 + t841 * t688;
t663 = -t689 * mrSges(7,3) - t733 * t703 + t811;
t656 = mrSges(6,2) * t672 + mrSges(7,2) * t670 - mrSges(6,3) * t668 - mrSges(7,3) * t665 - qJ(6) * t663 + t833 * t688 + t842 * t689 + t832 * t707 + t821 * t732 + t820 * t744;
t716 = Ifges(5,5) * t746 + Ifges(5,6) * t745 + Ifges(5,3) * t752;
t717 = Ifges(5,4) * t746 + Ifges(5,2) * t745 + Ifges(5,6) * t752;
t635 = mrSges(5,2) * t681 - mrSges(5,3) * t677 + Ifges(5,1) * t709 + Ifges(5,4) * t708 + Ifges(5,5) * t739 - pkin(12) * t657 - t784 * t649 + t789 * t656 + t745 * t716 - t752 * t717;
t718 = Ifges(5,1) * t746 + Ifges(5,4) * t745 + Ifges(5,5) * t752;
t641 = -mrSges(5,1) * t681 - mrSges(6,1) * t668 - mrSges(7,1) * t665 + mrSges(6,2) * t669 + mrSges(7,2) * t667 + mrSges(5,3) * t678 + Ifges(5,4) * t709 + Ifges(5,2) * t708 + Ifges(5,6) * t739 - pkin(4) * t657 - pkin(5) * t663 - t746 * t716 + t752 * t718 + t820 * t733 + t819 * t732 - t839 * t707 - t832 * t689 - t840 * t688;
t737 = Ifges(4,4) * t754 + Ifges(4,2) * t753 + Ifges(4,6) * t764;
t738 = Ifges(4,1) * t754 + Ifges(4,4) * t753 + Ifges(4,5) * t764;
t621 = mrSges(4,1) * t690 - mrSges(4,2) * t691 + Ifges(4,5) * t741 + Ifges(4,6) * t740 + Ifges(4,3) * t755 + pkin(3) * t796 + pkin(11) * t806 + t785 * t635 + t790 * t641 + t754 * t737 - t753 * t738;
t756 = Ifges(3,3) * t778 + (Ifges(3,5) * t787 + Ifges(3,6) * t792) * t816;
t758 = Ifges(3,5) * t778 + (Ifges(3,1) * t787 + Ifges(3,4) * t792) * t816;
t736 = Ifges(4,5) * t754 + Ifges(4,6) * t753 + Ifges(4,3) * t764;
t622 = mrSges(4,2) * t701 - mrSges(4,3) * t690 + Ifges(4,1) * t741 + Ifges(4,4) * t740 + Ifges(4,5) * t755 - pkin(11) * t648 + t790 * t635 - t785 * t641 + t753 * t736 - t764 * t737;
t627 = Ifges(4,4) * t741 + Ifges(4,2) * t740 + Ifges(4,6) * t755 - t754 * t736 + t764 * t738 - mrSges(4,1) * t701 + mrSges(4,3) * t691 - Ifges(5,5) * t709 - Ifges(5,6) * t708 - Ifges(5,3) * t739 - t746 * t717 + t745 * t718 - mrSges(5,1) * t677 + mrSges(5,2) * t678 - t784 * t656 - t789 * t649 - pkin(4) * t795 - pkin(12) * t805 - pkin(3) * t648;
t797 = pkin(10) * t640 + t622 * t786 + t627 * t791;
t615 = -mrSges(3,1) * t759 + mrSges(3,3) * t750 + Ifges(3,4) * t773 + Ifges(3,2) * t774 + Ifges(3,6) * t777 - pkin(2) * t633 - t780 * t621 - t756 * t809 + t778 * t758 + t797 * t782;
t757 = Ifges(3,6) * t778 + (Ifges(3,4) * t787 + Ifges(3,2) * t792) * t816;
t616 = t756 * t808 + mrSges(3,2) * t759 - mrSges(3,3) * t749 + Ifges(3,1) * t773 + Ifges(3,4) * t774 + Ifges(3,5) * t777 + t791 * t622 - t786 * t627 - t778 * t757 + (-t633 * t780 - t634 * t782) * pkin(10);
t798 = pkin(9) * t626 + t615 * t792 + t616 * t787;
t614 = mrSges(3,1) * t749 - mrSges(3,2) * t750 + Ifges(3,5) * t773 + Ifges(3,6) * t774 + Ifges(3,3) * t777 + pkin(2) * t634 + t782 * t621 + (t757 * t787 - t758 * t792) * t816 + t797 * t780;
t613 = -mrSges(2,2) * g(3) - mrSges(2,3) * t775 + Ifges(2,5) * qJDD(1) - t794 * Ifges(2,6) - t787 * t615 + t792 * t616 + (-t619 * t781 - t620 * t783) * pkin(9);
t612 = mrSges(2,1) * g(3) + mrSges(2,3) * t776 + t794 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t619 - t781 * t614 + t798 * t783;
t1 = [-m(1) * g(1) + t807; -m(1) * g(2) + t822; (-m(1) - m(2)) * g(3) + t619; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t822 - t788 * t612 + t793 * t613; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t807 + t793 * t612 + t788 * t613; -mrSges(1,1) * g(2) + mrSges(2,1) * t775 + mrSges(1,2) * g(1) - mrSges(2,2) * t776 + Ifges(2,3) * qJDD(1) + pkin(1) * t620 + t783 * t614 + t798 * t781;];
tauB  = t1;

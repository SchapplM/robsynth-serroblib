% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP12
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
% Datum: 2019-05-08 07:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:22:24
% EndTime: 2019-05-08 07:23:52
% DurationCPUTime: 64.71s
% Computational Cost: add. (1070686->393), mult. (2639644->513), div. (0->0), fcn. (2217859->14), ass. (0->171)
t839 = Ifges(6,1) + Ifges(7,1);
t829 = Ifges(6,4) - Ifges(7,5);
t838 = -Ifges(6,5) - Ifges(7,4);
t837 = Ifges(6,2) + Ifges(7,3);
t827 = Ifges(6,6) - Ifges(7,6);
t836 = -Ifges(6,3) - Ifges(7,2);
t781 = cos(pkin(6));
t776 = qJD(1) * t781 + qJD(2);
t778 = sin(pkin(7));
t780 = cos(pkin(7));
t779 = sin(pkin(6));
t789 = cos(qJ(2));
t810 = qJD(1) * t789;
t805 = t779 * t810;
t761 = (t776 * t778 + t780 * t805) * pkin(10);
t785 = sin(qJ(2));
t812 = qJD(1) * t779;
t833 = pkin(10) * t778;
t765 = (-pkin(2) * t789 - t785 * t833) * t812;
t809 = qJD(1) * qJD(2);
t771 = (qJDD(1) * t785 + t789 * t809) * t779;
t775 = qJDD(1) * t781 + qJDD(2);
t786 = sin(qJ(1));
t790 = cos(qJ(1));
t773 = t786 * g(1) - g(2) * t790;
t791 = qJD(1) ^ 2;
t834 = pkin(9) * t779;
t768 = qJDD(1) * pkin(1) + t791 * t834 + t773;
t774 = -g(1) * t790 - g(2) * t786;
t769 = -pkin(1) * t791 + qJDD(1) * t834 + t774;
t819 = t781 * t789;
t801 = t768 * t819 - t785 * t769;
t811 = qJD(1) * t785;
t832 = pkin(10) * t780;
t719 = -t771 * t832 + t775 * pkin(2) + t776 * t761 + (-g(3) * t789 - t765 * t811) * t779 + t801;
t806 = t779 * t811;
t764 = pkin(2) * t776 - t806 * t832;
t772 = (qJDD(1) * t789 - t785 * t809) * t779;
t798 = t772 * t780 + t775 * t778;
t820 = t781 * t785;
t813 = t768 * t820 + t789 * t769;
t720 = -t776 * t764 + (-g(3) * t785 + t765 * t810) * t779 + t798 * pkin(10) + t813;
t831 = t781 * g(3);
t725 = -t771 * t833 - t772 * pkin(2) - t831 + (-t768 + (-t761 * t789 + t764 * t785) * qJD(1)) * t779;
t784 = sin(qJ(3));
t788 = cos(qJ(3));
t685 = -t784 * t720 + (t719 * t780 + t725 * t778) * t788;
t821 = t780 * t789;
t826 = t778 * t784;
t752 = t776 * t826 + (t784 * t821 + t785 * t788) * t812;
t736 = -t752 * qJD(3) - t784 * t771 + t788 * t798;
t825 = t778 * t788;
t751 = (-t784 * t785 + t788 * t821) * t812 + t776 * t825;
t835 = cos(qJ(5));
t830 = -mrSges(6,3) - mrSges(7,2);
t824 = t779 * t785;
t823 = t779 * t789;
t822 = t780 * t784;
t686 = t719 * t822 + t788 * t720 + t725 * t826;
t738 = -mrSges(4,1) * t751 + mrSges(4,2) * t752;
t762 = t776 * t780 - t778 * t805 + qJD(3);
t745 = mrSges(4,1) * t762 - mrSges(4,3) * t752;
t753 = -t772 * t778 + t775 * t780 + qJDD(3);
t739 = -pkin(3) * t751 - pkin(11) * t752;
t760 = t762 ^ 2;
t678 = -pkin(3) * t760 + pkin(11) * t753 + t739 * t751 + t686;
t695 = -t778 * t719 + t780 * t725;
t737 = t751 * qJD(3) + t788 * t771 + t784 * t798;
t680 = (-t751 * t762 - t737) * pkin(11) + (t752 * t762 - t736) * pkin(3) + t695;
t783 = sin(qJ(4));
t787 = cos(qJ(4));
t674 = t787 * t678 + t783 * t680;
t743 = t752 * t787 + t762 * t783;
t705 = -qJD(4) * t743 - t737 * t783 + t753 * t787;
t742 = -t752 * t783 + t762 * t787;
t721 = -mrSges(5,1) * t742 + mrSges(5,2) * t743;
t750 = qJD(4) - t751;
t730 = mrSges(5,1) * t750 - mrSges(5,3) * t743;
t735 = qJDD(4) - t736;
t722 = -pkin(4) * t742 - pkin(12) * t743;
t749 = t750 ^ 2;
t670 = -pkin(4) * t749 + pkin(12) * t735 + t722 * t742 + t674;
t677 = -t753 * pkin(3) - t760 * pkin(11) + t752 * t739 - t685;
t706 = qJD(4) * t742 + t737 * t787 + t753 * t783;
t672 = (-t742 * t750 - t706) * pkin(12) + (t743 * t750 - t705) * pkin(4) + t677;
t782 = sin(qJ(5));
t667 = t835 * t670 + t782 * t672;
t728 = t743 * t835 + t782 * t750;
t683 = t728 * qJD(5) + t782 * t706 - t735 * t835;
t704 = qJDD(5) - t705;
t741 = qJD(5) - t742;
t710 = mrSges(6,1) * t741 - mrSges(6,3) * t728;
t727 = t782 * t743 - t750 * t835;
t698 = pkin(5) * t727 - qJ(6) * t728;
t740 = t741 ^ 2;
t663 = -pkin(5) * t740 + qJ(6) * t704 + 0.2e1 * qJD(6) * t741 - t698 * t727 + t667;
t711 = -mrSges(7,1) * t741 + mrSges(7,2) * t728;
t807 = m(7) * t663 + t704 * mrSges(7,3) + t741 * t711;
t699 = mrSges(7,1) * t727 - mrSges(7,3) * t728;
t814 = -mrSges(6,1) * t727 - mrSges(6,2) * t728 - t699;
t659 = m(6) * t667 - t704 * mrSges(6,2) + t683 * t830 - t741 * t710 + t727 * t814 + t807;
t666 = -t782 * t670 + t672 * t835;
t684 = -t727 * qJD(5) + t706 * t835 + t782 * t735;
t709 = -mrSges(6,2) * t741 - mrSges(6,3) * t727;
t664 = -t704 * pkin(5) - t740 * qJ(6) + t728 * t698 + qJDD(6) - t666;
t708 = -mrSges(7,2) * t727 + mrSges(7,3) * t741;
t800 = -m(7) * t664 + t704 * mrSges(7,1) + t741 * t708;
t660 = m(6) * t666 + t704 * mrSges(6,1) + t684 * t830 + t741 * t709 + t728 * t814 + t800;
t802 = t835 * t659 - t660 * t782;
t652 = m(5) * t674 - mrSges(5,2) * t735 + mrSges(5,3) * t705 + t721 * t742 - t730 * t750 + t802;
t673 = -t783 * t678 + t787 * t680;
t729 = -mrSges(5,2) * t750 + mrSges(5,3) * t742;
t669 = -t735 * pkin(4) - t749 * pkin(12) + t743 * t722 - t673;
t665 = -0.2e1 * qJD(6) * t728 + (t727 * t741 - t684) * qJ(6) + (t728 * t741 + t683) * pkin(5) + t669;
t661 = m(7) * t665 + mrSges(7,1) * t683 - t684 * mrSges(7,3) + t708 * t727 - t728 * t711;
t792 = -m(6) * t669 - t683 * mrSges(6,1) - mrSges(6,2) * t684 - t727 * t709 - t710 * t728 - t661;
t657 = m(5) * t673 + mrSges(5,1) * t735 - mrSges(5,3) * t706 - t721 * t743 + t729 * t750 + t792;
t803 = t787 * t652 - t657 * t783;
t643 = m(4) * t686 - mrSges(4,2) * t753 + mrSges(4,3) * t736 + t738 * t751 - t745 * t762 + t803;
t646 = t783 * t652 + t787 * t657;
t744 = -mrSges(4,2) * t762 + mrSges(4,3) * t751;
t645 = m(4) * t695 - mrSges(4,1) * t736 + mrSges(4,2) * t737 - t744 * t751 + t745 * t752 + t646;
t655 = t782 * t659 + t660 * t835;
t793 = -m(5) * t677 + t705 * mrSges(5,1) - t706 * mrSges(5,2) + t742 * t729 - t743 * t730 - t655;
t649 = m(4) * t685 + t753 * mrSges(4,1) - t737 * mrSges(4,3) - t752 * t738 + t762 * t744 + t793;
t632 = t780 * t788 * t649 + t643 * t822 - t645 * t778;
t746 = -g(3) * t823 + t801;
t767 = -mrSges(3,2) * t776 + mrSges(3,3) * t805;
t770 = (-mrSges(3,1) * t789 + mrSges(3,2) * t785) * t812;
t628 = m(3) * t746 + mrSges(3,1) * t775 - mrSges(3,3) * t771 + t767 * t776 - t770 * t806 + t632;
t631 = t643 * t826 + t780 * t645 + t649 * t825;
t757 = -t779 * t768 - t831;
t766 = mrSges(3,1) * t776 - mrSges(3,3) * t806;
t630 = m(3) * t757 - t772 * mrSges(3,1) + t771 * mrSges(3,2) + (t766 * t785 - t767 * t789) * t812 + t631;
t638 = t788 * t643 - t649 * t784;
t747 = -g(3) * t824 + t813;
t637 = m(3) * t747 - mrSges(3,2) * t775 + mrSges(3,3) * t772 - t766 * t776 + t770 * t805 + t638;
t618 = t628 * t819 - t630 * t779 + t637 * t820;
t616 = m(2) * t773 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t791 + t618;
t624 = -t628 * t785 + t789 * t637;
t623 = m(2) * t774 - mrSges(2,1) * t791 - qJDD(1) * mrSges(2,2) + t624;
t818 = t790 * t616 + t786 * t623;
t817 = t727 * t837 - t728 * t829 - t741 * t827;
t816 = t727 * t827 + t728 * t838 + t741 * t836;
t815 = -t829 * t727 + t728 * t839 - t838 * t741;
t617 = t628 * t823 + t781 * t630 + t637 * t824;
t804 = -t616 * t786 + t790 * t623;
t653 = -mrSges(6,1) * t669 - mrSges(7,1) * t665 + mrSges(7,2) * t663 + mrSges(6,3) * t667 - pkin(5) * t661 - t683 * t837 + t829 * t684 + t827 * t704 + t816 * t728 + t815 * t741;
t654 = mrSges(6,2) * t669 + mrSges(7,2) * t664 - mrSges(6,3) * t666 - mrSges(7,3) * t665 - qJ(6) * t661 - t829 * t683 + t684 * t839 - t704 * t838 + t816 * t727 + t817 * t741;
t712 = Ifges(5,5) * t743 + Ifges(5,6) * t742 + Ifges(5,3) * t750;
t713 = Ifges(5,4) * t743 + Ifges(5,2) * t742 + Ifges(5,6) * t750;
t633 = mrSges(5,2) * t677 - mrSges(5,3) * t673 + Ifges(5,1) * t706 + Ifges(5,4) * t705 + Ifges(5,5) * t735 - pkin(12) * t655 - t782 * t653 + t654 * t835 + t742 * t712 - t750 * t713;
t714 = Ifges(5,1) * t743 + Ifges(5,4) * t742 + Ifges(5,5) * t750;
t639 = Ifges(5,4) * t706 + Ifges(5,2) * t705 + Ifges(5,6) * t735 - t743 * t712 + t750 * t714 - mrSges(5,1) * t677 + mrSges(5,3) * t674 - mrSges(6,1) * t666 + mrSges(6,2) * t667 + mrSges(7,1) * t664 - mrSges(7,3) * t663 - pkin(5) * t800 - qJ(6) * t807 - pkin(4) * t655 + (pkin(5) * t699 + t817) * t728 + (qJ(6) * t699 - t815) * t727 + t836 * t704 + (mrSges(7,2) * pkin(5) + t838) * t684 + (mrSges(7,2) * qJ(6) + t827) * t683;
t732 = Ifges(4,4) * t752 + Ifges(4,2) * t751 + Ifges(4,6) * t762;
t733 = Ifges(4,1) * t752 + Ifges(4,4) * t751 + Ifges(4,5) * t762;
t619 = mrSges(4,1) * t685 - mrSges(4,2) * t686 + Ifges(4,5) * t737 + Ifges(4,6) * t736 + Ifges(4,3) * t753 + pkin(3) * t793 + pkin(11) * t803 + t783 * t633 + t787 * t639 + t752 * t732 - t751 * t733;
t754 = Ifges(3,3) * t776 + (Ifges(3,5) * t785 + Ifges(3,6) * t789) * t812;
t756 = Ifges(3,5) * t776 + (Ifges(3,1) * t785 + Ifges(3,4) * t789) * t812;
t731 = Ifges(4,5) * t752 + Ifges(4,6) * t751 + Ifges(4,3) * t762;
t620 = mrSges(4,2) * t695 - mrSges(4,3) * t685 + Ifges(4,1) * t737 + Ifges(4,4) * t736 + Ifges(4,5) * t753 - pkin(11) * t646 + t633 * t787 - t639 * t783 + t731 * t751 - t732 * t762;
t625 = Ifges(4,4) * t737 + Ifges(4,2) * t736 + Ifges(4,6) * t753 - t752 * t731 + t762 * t733 - mrSges(4,1) * t695 + mrSges(4,3) * t686 - Ifges(5,5) * t706 - Ifges(5,6) * t705 - Ifges(5,3) * t735 - t743 * t713 + t742 * t714 - mrSges(5,1) * t673 + mrSges(5,2) * t674 - t782 * t654 - t835 * t653 - pkin(4) * t792 - pkin(12) * t802 - pkin(3) * t646;
t794 = pkin(10) * t638 + t620 * t784 + t625 * t788;
t613 = -mrSges(3,1) * t757 + mrSges(3,3) * t747 + Ifges(3,4) * t771 + Ifges(3,2) * t772 + Ifges(3,6) * t775 - pkin(2) * t631 - t778 * t619 - t754 * t806 + t776 * t756 + t780 * t794;
t755 = Ifges(3,6) * t776 + (Ifges(3,4) * t785 + Ifges(3,2) * t789) * t812;
t614 = t754 * t805 + mrSges(3,2) * t757 - mrSges(3,3) * t746 + Ifges(3,1) * t771 + Ifges(3,4) * t772 + Ifges(3,5) * t775 + t788 * t620 - t784 * t625 - t776 * t755 + (-t631 * t778 - t632 * t780) * pkin(10);
t795 = pkin(9) * t624 + t613 * t789 + t614 * t785;
t612 = mrSges(3,1) * t746 - mrSges(3,2) * t747 + Ifges(3,5) * t771 + Ifges(3,6) * t772 + Ifges(3,3) * t775 + pkin(2) * t632 + t780 * t619 + (t755 * t785 - t756 * t789) * t812 + t794 * t778;
t611 = -mrSges(2,2) * g(3) - mrSges(2,3) * t773 + Ifges(2,5) * qJDD(1) - t791 * Ifges(2,6) - t785 * t613 + t789 * t614 + (-t617 * t779 - t618 * t781) * pkin(9);
t610 = mrSges(2,1) * g(3) + mrSges(2,3) * t774 + t791 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t617 - t779 * t612 + t781 * t795;
t1 = [-m(1) * g(1) + t804; -m(1) * g(2) + t818; (-m(1) - m(2)) * g(3) + t617; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t818 - t786 * t610 + t790 * t611; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t804 + t790 * t610 + t786 * t611; -mrSges(1,1) * g(2) + mrSges(2,1) * t773 + mrSges(1,2) * g(1) - mrSges(2,2) * t774 + Ifges(2,3) * qJDD(1) + pkin(1) * t618 + t781 * t612 + t779 * t795;];
tauB  = t1;

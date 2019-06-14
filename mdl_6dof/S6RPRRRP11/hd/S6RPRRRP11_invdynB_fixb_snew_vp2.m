% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:03:43
% EndTime: 2019-05-06 02:04:41
% DurationCPUTime: 56.33s
% Computational Cost: add. (810666->381), mult. (2517882->508), div. (0->0), fcn. (2144242->14), ass. (0->175)
t792 = sin(pkin(12));
t794 = sin(pkin(6));
t795 = cos(pkin(12));
t797 = cos(pkin(6));
t800 = sin(qJ(3));
t796 = cos(pkin(7));
t804 = cos(qJ(3));
t843 = t796 * t804;
t793 = sin(pkin(7));
t848 = t793 * t804;
t809 = (-t792 * t800 + t795 * t843) * t794 + t797 * t848;
t767 = t809 * qJD(1);
t844 = t796 * t800;
t849 = t793 * t800;
t811 = t797 * t849 + (t792 * t804 + t795 * t844) * t794;
t768 = t811 * qJD(1);
t756 = -t768 * qJD(3) + t809 * qJDD(1);
t863 = Ifges(6,1) + Ifges(7,1);
t857 = Ifges(6,4) + Ifges(7,4);
t856 = Ifges(6,5) + Ifges(7,5);
t862 = Ifges(6,2) + Ifges(7,2);
t861 = Ifges(6,6) + Ifges(7,6);
t860 = Ifges(6,3) + Ifges(7,3);
t846 = t794 * t796;
t779 = (t793 * t797 + t795 * t846) * qJD(1) * pkin(9);
t801 = sin(qJ(1));
t805 = cos(qJ(1));
t790 = -t805 * g(1) - t801 * g(2);
t806 = qJD(1) ^ 2;
t852 = qJ(2) * t794;
t783 = -t806 * pkin(1) + qJDD(1) * t852 + t790;
t859 = pkin(9) * t792;
t821 = -pkin(2) * t795 - t793 * t859;
t837 = qJD(1) * t794;
t853 = pkin(9) * qJDD(1);
t816 = qJD(1) * t821 * t837 + t796 * t853;
t789 = t801 * g(1) - t805 * g(2);
t782 = qJDD(1) * pkin(1) + t806 * t852 + t789;
t831 = qJD(2) * t837;
t845 = t795 * t797;
t847 = t794 * t795;
t822 = -g(3) * t847 + t782 * t845 - 0.2e1 * t792 * t831;
t737 = (pkin(2) * qJDD(1) + qJD(1) * t779) * t797 + (-t816 * t794 - t783) * t792 + t822;
t784 = (pkin(2) * t797 - t846 * t859) * qJD(1);
t850 = t792 * t797;
t832 = t782 * t850 + (t783 + 0.2e1 * t831) * t795;
t738 = (-qJD(1) * t784 + t793 * t853) * t797 + (-g(3) * t792 + t816 * t795) * t794 + t832;
t830 = -t797 * g(3) + qJDD(2);
t747 = (-t782 + t821 * qJDD(1) + (-t779 * t795 + t784 * t792) * qJD(1)) * t794 + t830;
t699 = -t800 * t738 + (t737 * t796 + t747 * t793) * t804;
t858 = -mrSges(6,2) - mrSges(7,2);
t854 = Ifges(3,3) * t797;
t851 = t792 * t794;
t700 = t737 * t844 + t804 * t738 + t747 * t849;
t754 = -t767 * mrSges(4,1) + t768 * mrSges(4,2);
t817 = -t793 * t847 + t796 * t797;
t780 = t817 * qJD(1) + qJD(3);
t764 = t780 * mrSges(4,1) - t768 * mrSges(4,3);
t777 = t817 * qJDD(1) + qJDD(3);
t755 = -t767 * pkin(3) - t768 * pkin(10);
t776 = t780 ^ 2;
t696 = -t776 * pkin(3) + t777 * pkin(10) + t767 * t755 + t700;
t715 = -t793 * t737 + t796 * t747;
t757 = t767 * qJD(3) + t811 * qJDD(1);
t698 = (-t767 * t780 - t757) * pkin(10) + (t768 * t780 - t756) * pkin(3) + t715;
t799 = sin(qJ(4));
t803 = cos(qJ(4));
t692 = t803 * t696 + t799 * t698;
t762 = t803 * t768 + t799 * t780;
t732 = -t762 * qJD(4) - t799 * t757 + t803 * t777;
t761 = -t799 * t768 + t803 * t780;
t739 = -t761 * mrSges(5,1) + t762 * mrSges(5,2);
t766 = qJD(4) - t767;
t749 = t766 * mrSges(5,1) - t762 * mrSges(5,3);
t753 = qJDD(4) - t756;
t740 = -t761 * pkin(4) - t762 * pkin(11);
t765 = t766 ^ 2;
t687 = -t765 * pkin(4) + t753 * pkin(11) + t761 * t740 + t692;
t695 = -t777 * pkin(3) - t776 * pkin(10) + t768 * t755 - t699;
t733 = t761 * qJD(4) + t803 * t757 + t799 * t777;
t690 = (-t761 * t766 - t733) * pkin(11) + (t762 * t766 - t732) * pkin(4) + t695;
t798 = sin(qJ(5));
t802 = cos(qJ(5));
t682 = -t798 * t687 + t802 * t690;
t745 = -t798 * t762 + t802 * t766;
t705 = t745 * qJD(5) + t802 * t733 + t798 * t753;
t746 = t802 * t762 + t798 * t766;
t717 = -t745 * mrSges(7,1) + t746 * mrSges(7,2);
t718 = -t745 * mrSges(6,1) + t746 * mrSges(6,2);
t758 = qJD(5) - t761;
t721 = -t758 * mrSges(6,2) + t745 * mrSges(6,3);
t730 = qJDD(5) - t732;
t679 = -0.2e1 * qJD(6) * t746 + (t745 * t758 - t705) * qJ(6) + (t745 * t746 + t730) * pkin(5) + t682;
t720 = -t758 * mrSges(7,2) + t745 * mrSges(7,3);
t834 = m(7) * t679 + t730 * mrSges(7,1) + t758 * t720;
t672 = m(6) * t682 + t730 * mrSges(6,1) + t758 * t721 + (-t717 - t718) * t746 + (-mrSges(6,3) - mrSges(7,3)) * t705 + t834;
t683 = t802 * t687 + t798 * t690;
t704 = -t746 * qJD(5) - t798 * t733 + t802 * t753;
t722 = t758 * pkin(5) - t746 * qJ(6);
t744 = t745 ^ 2;
t681 = -t744 * pkin(5) + t704 * qJ(6) + 0.2e1 * qJD(6) * t745 - t758 * t722 + t683;
t833 = m(7) * t681 + t704 * mrSges(7,3) + t745 * t717;
t723 = t758 * mrSges(7,1) - t746 * mrSges(7,3);
t838 = -t758 * mrSges(6,1) + t746 * mrSges(6,3) - t723;
t674 = m(6) * t683 + t704 * mrSges(6,3) + t745 * t718 + t858 * t730 + t838 * t758 + t833;
t827 = -t798 * t672 + t802 * t674;
t669 = m(5) * t692 - t753 * mrSges(5,2) + t732 * mrSges(5,3) + t761 * t739 - t766 * t749 + t827;
t691 = -t799 * t696 + t803 * t698;
t748 = -t766 * mrSges(5,2) + t761 * mrSges(5,3);
t686 = -t753 * pkin(4) - t765 * pkin(11) + t762 * t740 - t691;
t684 = -t704 * pkin(5) - t744 * qJ(6) + t746 * t722 + qJDD(6) + t686;
t826 = m(7) * t684 - t704 * mrSges(7,1) - t745 * t720;
t807 = -m(6) * t686 + t704 * mrSges(6,1) + t858 * t705 + t745 * t721 + t838 * t746 - t826;
t676 = m(5) * t691 + t753 * mrSges(5,1) - t733 * mrSges(5,3) - t762 * t739 + t766 * t748 + t807;
t828 = t803 * t669 - t799 * t676;
t659 = m(4) * t700 - t777 * mrSges(4,2) + t756 * mrSges(4,3) + t767 * t754 - t780 * t764 + t828;
t662 = t799 * t669 + t803 * t676;
t763 = -t780 * mrSges(4,2) + t767 * mrSges(4,3);
t661 = m(4) * t715 - t756 * mrSges(4,1) + t757 * mrSges(4,2) - t767 * t763 + t768 * t764 + t662;
t671 = t802 * t672 + t798 * t674;
t808 = -m(5) * t695 + t732 * mrSges(5,1) - t733 * mrSges(5,2) + t761 * t748 - t762 * t749 - t671;
t666 = m(4) * t699 + t777 * mrSges(4,1) - t757 * mrSges(4,3) - t768 * t754 + t780 * t763 + t808;
t648 = t659 * t844 - t793 * t661 + t666 * t843;
t759 = -t792 * t783 + t822;
t825 = -mrSges(3,1) * t795 + mrSges(3,2) * t792;
t781 = t825 * t837;
t819 = -mrSges(3,2) * t797 + mrSges(3,3) * t847;
t786 = t819 * qJD(1);
t820 = mrSges(3,1) * t797 - mrSges(3,3) * t851;
t644 = m(3) * t759 + t820 * qJDD(1) + (-t781 * t851 + t786 * t797) * qJD(1) + t648;
t647 = t659 * t849 + t796 * t661 + t666 * t848;
t769 = -t794 * t782 + t830;
t785 = t820 * qJD(1);
t646 = m(3) * t769 + (t825 * qJDD(1) + (t785 * t792 - t786 * t795) * qJD(1)) * t794 + t647;
t654 = t804 * t659 - t800 * t666;
t760 = -g(3) * t851 + t832;
t653 = m(3) * t760 + t819 * qJDD(1) + (t781 * t847 - t785 * t797) * qJD(1) + t654;
t634 = t644 * t845 - t794 * t646 + t653 * t850;
t632 = m(2) * t789 + qJDD(1) * mrSges(2,1) - t806 * mrSges(2,2) + t634;
t640 = -t792 * t644 + t795 * t653;
t639 = m(2) * t790 - t806 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t640;
t842 = t805 * t632 + t801 * t639;
t841 = t861 * t745 + t856 * t746 + t860 * t758;
t840 = -t862 * t745 - t857 * t746 - t861 * t758;
t839 = t857 * t745 + t863 * t746 + t856 * t758;
t633 = t644 * t847 + t797 * t646 + t653 * t851;
t829 = -t801 * t632 + t805 * t639;
t824 = Ifges(3,5) * t792 + Ifges(3,6) * t795;
t663 = -mrSges(6,1) * t686 + mrSges(6,3) * t683 - mrSges(7,1) * t684 + mrSges(7,3) * t681 - pkin(5) * t826 + qJ(6) * t833 + (-qJ(6) * t723 + t839) * t758 + (-pkin(5) * t723 - t841) * t746 + (-qJ(6) * mrSges(7,2) + t861) * t730 + (-pkin(5) * mrSges(7,2) + t857) * t705 + t862 * t704;
t677 = -t705 * mrSges(7,3) - t746 * t717 + t834;
t670 = mrSges(6,2) * t686 + mrSges(7,2) * t684 - mrSges(6,3) * t682 - mrSges(7,3) * t679 - qJ(6) * t677 + t857 * t704 + t863 * t705 + t856 * t730 + t841 * t745 + t840 * t758;
t726 = Ifges(5,5) * t762 + Ifges(5,6) * t761 + Ifges(5,3) * t766;
t727 = Ifges(5,4) * t762 + Ifges(5,2) * t761 + Ifges(5,6) * t766;
t649 = mrSges(5,2) * t695 - mrSges(5,3) * t691 + Ifges(5,1) * t733 + Ifges(5,4) * t732 + Ifges(5,5) * t753 - pkin(11) * t671 - t798 * t663 + t802 * t670 + t761 * t726 - t766 * t727;
t728 = Ifges(5,1) * t762 + Ifges(5,4) * t761 + Ifges(5,5) * t766;
t655 = -mrSges(5,1) * t695 - mrSges(6,1) * t682 - mrSges(7,1) * t679 + mrSges(6,2) * t683 + mrSges(7,2) * t681 + mrSges(5,3) * t692 + Ifges(5,4) * t733 + Ifges(5,2) * t732 + Ifges(5,6) * t753 - pkin(4) * t671 - pkin(5) * t677 - t762 * t726 + t766 * t728 + t840 * t746 + t839 * t745 - t860 * t730 - t856 * t705 - t861 * t704;
t750 = Ifges(4,5) * t768 + Ifges(4,6) * t767 + Ifges(4,3) * t780;
t751 = Ifges(4,4) * t768 + Ifges(4,2) * t767 + Ifges(4,6) * t780;
t636 = mrSges(4,2) * t715 - mrSges(4,3) * t699 + Ifges(4,1) * t757 + Ifges(4,4) * t756 + Ifges(4,5) * t777 - pkin(10) * t662 + t803 * t649 - t799 * t655 + t767 * t750 - t780 * t751;
t752 = Ifges(4,1) * t768 + Ifges(4,4) * t767 + Ifges(4,5) * t780;
t641 = Ifges(4,4) * t757 + Ifges(4,2) * t756 + Ifges(4,6) * t777 - t768 * t750 + t780 * t752 - mrSges(4,1) * t715 + mrSges(4,3) * t700 - Ifges(5,5) * t733 - Ifges(5,6) * t732 - Ifges(5,3) * t753 - t762 * t727 + t761 * t728 - mrSges(5,1) * t691 + mrSges(5,2) * t692 - t798 * t670 - t802 * t663 - pkin(4) * t807 - pkin(11) * t827 - pkin(3) * t662;
t815 = pkin(9) * t654 + t636 * t800 + t641 * t804;
t635 = mrSges(4,1) * t699 - mrSges(4,2) * t700 + Ifges(4,5) * t757 + Ifges(4,6) * t756 + Ifges(4,3) * t777 + pkin(3) * t808 + pkin(10) * t828 + t799 * t649 + t803 * t655 + t768 * t751 - t767 * t752;
t772 = (t824 * t794 + t854) * qJD(1);
t813 = Ifges(3,5) * t797 + (Ifges(3,1) * t792 + Ifges(3,4) * t795) * t794;
t774 = t813 * qJD(1);
t812 = Ifges(3,6) * t797 + (Ifges(3,4) * t792 + Ifges(3,2) * t795) * t794;
t629 = -mrSges(3,1) * t769 + mrSges(3,3) * t760 - pkin(2) * t647 - t793 * t635 + (-t772 * t851 + t774 * t797) * qJD(1) + t815 * t796 + t812 * qJDD(1);
t773 = t812 * qJD(1);
t630 = mrSges(3,2) * t769 - mrSges(3,3) * t759 + t804 * t636 - t800 * t641 + (t772 * t847 - t773 * t797) * qJD(1) + (-t647 * t793 - t648 * t796) * pkin(9) + t813 * qJDD(1);
t814 = qJ(2) * t640 + t629 * t795 + t630 * t792;
t628 = qJDD(1) * t854 + mrSges(3,1) * t759 - mrSges(3,2) * t760 + pkin(2) * t648 + t796 * t635 + t815 * t793 + (t824 * qJDD(1) + (t773 * t792 - t774 * t795) * qJD(1)) * t794;
t627 = -mrSges(2,2) * g(3) - mrSges(2,3) * t789 + Ifges(2,5) * qJDD(1) - t806 * Ifges(2,6) - t792 * t629 + t795 * t630 + (-t633 * t794 - t634 * t797) * qJ(2);
t626 = mrSges(2,1) * g(3) + mrSges(2,3) * t790 + t806 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t633 - t794 * t628 + t814 * t797;
t1 = [-m(1) * g(1) + t829; -m(1) * g(2) + t842; (-m(1) - m(2)) * g(3) + t633; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t842 - t801 * t626 + t805 * t627; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t829 + t805 * t626 + t801 * t627; -mrSges(1,1) * g(2) + mrSges(2,1) * t789 + mrSges(1,2) * g(1) - mrSges(2,2) * t790 + Ifges(2,3) * qJDD(1) + pkin(1) * t634 + t797 * t628 + t814 * t794;];
tauB  = t1;

% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:40:08
% EndTime: 2019-12-31 22:40:25
% DurationCPUTime: 17.06s
% Computational Cost: add. (277604->328), mult. (593617->428), div. (0->0), fcn. (465021->12), ass. (0->140)
t792 = sin(pkin(5));
t797 = sin(qJ(2));
t802 = cos(qJ(2));
t817 = qJD(1) * qJD(2);
t779 = (-qJDD(1) * t802 + t797 * t817) * t792;
t827 = pkin(7) * t792;
t793 = cos(pkin(5));
t826 = t793 * g(3);
t825 = t792 * t797;
t824 = t792 * t802;
t823 = t793 * t797;
t822 = t793 * t802;
t798 = sin(qJ(1));
t803 = cos(qJ(1));
t784 = t798 * g(1) - t803 * g(2);
t804 = qJD(1) ^ 2;
t774 = qJDD(1) * pkin(1) + t804 * t827 + t784;
t785 = -t803 * g(1) - t798 * g(2);
t775 = -t804 * pkin(1) + qJDD(1) * t827 + t785;
t820 = t774 * t823 + t802 * t775;
t749 = -g(3) * t825 + t820;
t788 = t793 * qJD(1) + qJD(2);
t819 = qJD(1) * t792;
t816 = t797 * t819;
t772 = t788 * mrSges(3,1) - mrSges(3,3) * t816;
t776 = (-mrSges(3,1) * t802 + mrSges(3,2) * t797) * t819;
t787 = t793 * qJDD(1) + qJDD(2);
t777 = (-pkin(2) * t802 - pkin(8) * t797) * t819;
t786 = t788 ^ 2;
t818 = qJD(1) * t802;
t729 = -t786 * pkin(2) + t787 * pkin(8) + (-g(3) * t797 + t777 * t818) * t792 + t820;
t778 = (qJDD(1) * t797 + t802 * t817) * t792;
t730 = t779 * pkin(2) - t778 * pkin(8) - t826 + (-t774 + (pkin(2) * t797 - pkin(8) * t802) * t788 * qJD(1)) * t792;
t796 = sin(qJ(3));
t801 = cos(qJ(3));
t710 = t801 * t729 + t796 * t730;
t766 = t801 * t788 - t796 * t816;
t767 = t796 * t788 + t801 * t816;
t751 = -t766 * pkin(3) - t767 * pkin(9);
t771 = qJDD(3) + t779;
t815 = t792 * t818;
t783 = qJD(3) - t815;
t781 = t783 ^ 2;
t704 = -t781 * pkin(3) + t771 * pkin(9) + t766 * t751 + t710;
t748 = -g(3) * t824 + t774 * t822 - t797 * t775;
t728 = -t787 * pkin(2) - t786 * pkin(8) + t777 * t816 - t748;
t746 = -t767 * qJD(3) - t796 * t778 + t801 * t787;
t747 = t766 * qJD(3) + t801 * t778 + t796 * t787;
t708 = (-t766 * t783 - t747) * pkin(9) + (t767 * t783 - t746) * pkin(3) + t728;
t795 = sin(qJ(4));
t800 = cos(qJ(4));
t694 = -t795 * t704 + t800 * t708;
t753 = -t795 * t767 + t800 * t783;
t718 = t753 * qJD(4) + t800 * t747 + t795 * t771;
t744 = qJDD(4) - t746;
t754 = t800 * t767 + t795 * t783;
t765 = qJD(4) - t766;
t692 = (t753 * t765 - t718) * pkin(10) + (t753 * t754 + t744) * pkin(4) + t694;
t695 = t800 * t704 + t795 * t708;
t717 = -t754 * qJD(4) - t795 * t747 + t800 * t771;
t737 = t765 * pkin(4) - t754 * pkin(10);
t752 = t753 ^ 2;
t693 = -t752 * pkin(4) + t717 * pkin(10) - t765 * t737 + t695;
t794 = sin(qJ(5));
t799 = cos(qJ(5));
t690 = t799 * t692 - t794 * t693;
t731 = t799 * t753 - t794 * t754;
t701 = t731 * qJD(5) + t794 * t717 + t799 * t718;
t732 = t794 * t753 + t799 * t754;
t715 = -t731 * mrSges(6,1) + t732 * mrSges(6,2);
t763 = qJD(5) + t765;
t719 = -t763 * mrSges(6,2) + t731 * mrSges(6,3);
t739 = qJDD(5) + t744;
t686 = m(6) * t690 + t739 * mrSges(6,1) - t701 * mrSges(6,3) - t732 * t715 + t763 * t719;
t691 = t794 * t692 + t799 * t693;
t700 = -t732 * qJD(5) + t799 * t717 - t794 * t718;
t720 = t763 * mrSges(6,1) - t732 * mrSges(6,3);
t687 = m(6) * t691 - t739 * mrSges(6,2) + t700 * mrSges(6,3) + t731 * t715 - t763 * t720;
t678 = t799 * t686 + t794 * t687;
t733 = -t753 * mrSges(5,1) + t754 * mrSges(5,2);
t735 = -t765 * mrSges(5,2) + t753 * mrSges(5,3);
t676 = m(5) * t694 + t744 * mrSges(5,1) - t718 * mrSges(5,3) - t754 * t733 + t765 * t735 + t678;
t736 = t765 * mrSges(5,1) - t754 * mrSges(5,3);
t812 = -t794 * t686 + t799 * t687;
t677 = m(5) * t695 - t744 * mrSges(5,2) + t717 * mrSges(5,3) + t753 * t733 - t765 * t736 + t812;
t674 = -t795 * t676 + t800 * t677;
t750 = -t766 * mrSges(4,1) + t767 * mrSges(4,2);
t756 = t783 * mrSges(4,1) - t767 * mrSges(4,3);
t672 = m(4) * t710 - t771 * mrSges(4,2) + t746 * mrSges(4,3) + t766 * t750 - t783 * t756 + t674;
t709 = -t796 * t729 + t801 * t730;
t703 = -t771 * pkin(3) - t781 * pkin(9) + t767 * t751 - t709;
t696 = -t717 * pkin(4) - t752 * pkin(10) + t754 * t737 + t703;
t810 = m(6) * t696 - t700 * mrSges(6,1) + t701 * mrSges(6,2) - t731 * t719 + t732 * t720;
t688 = -m(5) * t703 + t717 * mrSges(5,1) - t718 * mrSges(5,2) + t753 * t735 - t754 * t736 - t810;
t755 = -t783 * mrSges(4,2) + t766 * mrSges(4,3);
t682 = m(4) * t709 + t771 * mrSges(4,1) - t747 * mrSges(4,3) - t767 * t750 + t783 * t755 + t688;
t813 = t801 * t672 - t796 * t682;
t661 = m(3) * t749 - t787 * mrSges(3,2) - t779 * mrSges(3,3) - t788 * t772 + t776 * t815 + t813;
t664 = t796 * t672 + t801 * t682;
t760 = -t792 * t774 - t826;
t773 = -t788 * mrSges(3,2) + mrSges(3,3) * t815;
t663 = m(3) * t760 + t779 * mrSges(3,1) + t778 * mrSges(3,2) + (t772 * t797 - t773 * t802) * t819 + t664;
t673 = t800 * t676 + t795 * t677;
t807 = -m(4) * t728 + t746 * mrSges(4,1) - t747 * mrSges(4,2) + t766 * t755 - t767 * t756 - t673;
t669 = m(3) * t748 + t787 * mrSges(3,1) - t778 * mrSges(3,3) + t788 * t773 - t776 * t816 + t807;
t650 = t661 * t823 - t792 * t663 + t669 * t822;
t647 = m(2) * t784 + qJDD(1) * mrSges(2,1) - t804 * mrSges(2,2) + t650;
t656 = t802 * t661 - t797 * t669;
t654 = m(2) * t785 - t804 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t656;
t821 = t803 * t647 + t798 * t654;
t649 = t661 * t825 + t793 * t663 + t669 * t824;
t814 = -t798 * t647 + t803 * t654;
t711 = Ifges(6,5) * t732 + Ifges(6,6) * t731 + Ifges(6,3) * t763;
t713 = Ifges(6,1) * t732 + Ifges(6,4) * t731 + Ifges(6,5) * t763;
t679 = -mrSges(6,1) * t696 + mrSges(6,3) * t691 + Ifges(6,4) * t701 + Ifges(6,2) * t700 + Ifges(6,6) * t739 - t732 * t711 + t763 * t713;
t712 = Ifges(6,4) * t732 + Ifges(6,2) * t731 + Ifges(6,6) * t763;
t680 = mrSges(6,2) * t696 - mrSges(6,3) * t690 + Ifges(6,1) * t701 + Ifges(6,4) * t700 + Ifges(6,5) * t739 + t731 * t711 - t763 * t712;
t721 = Ifges(5,5) * t754 + Ifges(5,6) * t753 + Ifges(5,3) * t765;
t723 = Ifges(5,1) * t754 + Ifges(5,4) * t753 + Ifges(5,5) * t765;
t665 = -mrSges(5,1) * t703 + mrSges(5,3) * t695 + Ifges(5,4) * t718 + Ifges(5,2) * t717 + Ifges(5,6) * t744 - pkin(4) * t810 + pkin(10) * t812 + t799 * t679 + t794 * t680 - t754 * t721 + t765 * t723;
t722 = Ifges(5,4) * t754 + Ifges(5,2) * t753 + Ifges(5,6) * t765;
t666 = mrSges(5,2) * t703 - mrSges(5,3) * t694 + Ifges(5,1) * t718 + Ifges(5,4) * t717 + Ifges(5,5) * t744 - pkin(10) * t678 - t794 * t679 + t799 * t680 + t753 * t721 - t765 * t722;
t740 = Ifges(4,5) * t767 + Ifges(4,6) * t766 + Ifges(4,3) * t783;
t741 = Ifges(4,4) * t767 + Ifges(4,2) * t766 + Ifges(4,6) * t783;
t651 = mrSges(4,2) * t728 - mrSges(4,3) * t709 + Ifges(4,1) * t747 + Ifges(4,4) * t746 + Ifges(4,5) * t771 - pkin(9) * t673 - t795 * t665 + t800 * t666 + t766 * t740 - t783 * t741;
t742 = Ifges(4,1) * t767 + Ifges(4,4) * t766 + Ifges(4,5) * t783;
t808 = -mrSges(6,1) * t690 + mrSges(6,2) * t691 - Ifges(6,5) * t701 - Ifges(6,6) * t700 - Ifges(6,3) * t739 - t732 * t712 + t731 * t713;
t805 = mrSges(5,1) * t694 - mrSges(5,2) * t695 + Ifges(5,5) * t718 + Ifges(5,6) * t717 + Ifges(5,3) * t744 + pkin(4) * t678 + t754 * t722 - t753 * t723 - t808;
t657 = -mrSges(4,1) * t728 + mrSges(4,3) * t710 + Ifges(4,4) * t747 + Ifges(4,2) * t746 + Ifges(4,6) * t771 - pkin(3) * t673 - t767 * t740 + t783 * t742 - t805;
t758 = Ifges(3,6) * t788 + (Ifges(3,4) * t797 + Ifges(3,2) * t802) * t819;
t759 = Ifges(3,5) * t788 + (Ifges(3,1) * t797 + Ifges(3,4) * t802) * t819;
t641 = Ifges(3,5) * t778 - Ifges(3,6) * t779 + Ifges(3,3) * t787 + mrSges(3,1) * t748 - mrSges(3,2) * t749 + t796 * t651 + t801 * t657 + pkin(2) * t807 + pkin(8) * t813 + (t758 * t797 - t759 * t802) * t819;
t757 = Ifges(3,3) * t788 + (Ifges(3,5) * t797 + Ifges(3,6) * t802) * t819;
t643 = mrSges(3,2) * t760 - mrSges(3,3) * t748 + Ifges(3,1) * t778 - Ifges(3,4) * t779 + Ifges(3,5) * t787 - pkin(8) * t664 + t801 * t651 - t796 * t657 + t757 * t815 - t788 * t758;
t806 = mrSges(4,1) * t709 - mrSges(4,2) * t710 + Ifges(4,5) * t747 + Ifges(4,6) * t746 + Ifges(4,3) * t771 + pkin(3) * t688 + pkin(9) * t674 + t800 * t665 + t795 * t666 + t767 * t741 - t766 * t742;
t645 = -mrSges(3,1) * t760 + mrSges(3,3) * t749 + Ifges(3,4) * t778 - Ifges(3,2) * t779 + Ifges(3,6) * t787 - pkin(2) * t664 - t757 * t816 + t788 * t759 - t806;
t809 = mrSges(2,1) * t784 - mrSges(2,2) * t785 + Ifges(2,3) * qJDD(1) + pkin(1) * t650 + t793 * t641 + t643 * t825 + t645 * t824 + t656 * t827;
t639 = -mrSges(2,2) * g(3) - mrSges(2,3) * t784 + Ifges(2,5) * qJDD(1) - t804 * Ifges(2,6) + t802 * t643 - t797 * t645 + (-t649 * t792 - t650 * t793) * pkin(7);
t638 = mrSges(2,1) * g(3) + mrSges(2,3) * t785 + t804 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t649 - t792 * t641 + (pkin(7) * t656 + t643 * t797 + t645 * t802) * t793;
t1 = [-m(1) * g(1) + t814; -m(1) * g(2) + t821; (-m(1) - m(2)) * g(3) + t649; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t821 - t798 * t638 + t803 * t639; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t814 + t803 * t638 + t798 * t639; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t809; t809; t641; t806; t805; -t808;];
tauJB = t1;

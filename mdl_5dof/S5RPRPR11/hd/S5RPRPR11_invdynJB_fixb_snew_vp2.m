% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:29
% EndTime: 2019-12-31 18:27:33
% DurationCPUTime: 3.83s
% Computational Cost: add. (34299->272), mult. (83010->327), div. (0->0), fcn. (56426->8), ass. (0->122)
t777 = sin(pkin(8));
t770 = t777 ^ 2;
t778 = cos(pkin(8));
t771 = t778 ^ 2;
t785 = qJD(1) ^ 2;
t781 = sin(qJ(1));
t783 = cos(qJ(1));
t756 = t781 * g(1) - t783 * g(2);
t805 = -qJDD(2) + t756;
t827 = pkin(2) * t778;
t735 = -(qJ(2) + (t770 + t771) * pkin(6)) * t785 - (pkin(1) + t827) * qJDD(1) - t805;
t780 = sin(qJ(3));
t828 = cos(qJ(3));
t793 = t828 * t777 + t778 * t780;
t806 = t778 * t828;
t811 = t777 * qJD(1);
t750 = -qJD(1) * t806 + t780 * t811;
t813 = t750 * qJD(3);
t737 = t793 * qJDD(1) - t813;
t837 = t735 + (-t737 + t813) * qJ(4);
t836 = Ifges(4,1) + Ifges(5,1);
t825 = Ifges(4,4) - Ifges(5,5);
t824 = Ifges(4,5) + Ifges(5,4);
t835 = -Ifges(4,2) - Ifges(5,3);
t823 = Ifges(4,6) - Ifges(5,6);
t834 = Ifges(4,3) + Ifges(5,2);
t751 = t793 * qJD(1);
t730 = t750 * mrSges(5,1) - t751 * mrSges(5,3);
t757 = -t783 * g(1) - t781 * g(2);
t752 = -t785 * pkin(1) + qJDD(1) * qJ(2) + t757;
t810 = qJD(1) * qJD(2);
t804 = -t778 * g(3) - 0.2e1 * t777 * t810;
t716 = (-pkin(6) * qJDD(1) + t785 * t827 - t752) * t777 + t804;
t739 = -t777 * g(3) + (t752 + 0.2e1 * t810) * t778;
t808 = qJDD(1) * t778;
t820 = t771 * t785;
t717 = -pkin(2) * t820 + pkin(6) * t808 + t739;
t701 = t828 * t716 - t780 * t717;
t729 = t750 * pkin(3) - t751 * qJ(4);
t784 = qJD(3) ^ 2;
t693 = -qJDD(3) * pkin(3) - t784 * qJ(4) + t751 * t729 + qJDD(4) - t701;
t687 = (-t737 - t813) * pkin(7) + (t750 * t751 - qJDD(3)) * pkin(4) + t693;
t702 = t780 * t716 + t828 * t717;
t829 = 2 * qJD(4);
t692 = -t784 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t829 - t750 * t729 + t702;
t809 = qJDD(1) * t777;
t812 = t751 * qJD(3);
t736 = -qJDD(1) * t806 + t780 * t809 + t812;
t746 = -qJD(3) * pkin(4) - t751 * pkin(7);
t749 = t750 ^ 2;
t688 = -t749 * pkin(4) + t736 * pkin(7) + qJD(3) * t746 + t692;
t779 = sin(qJ(5));
t782 = cos(qJ(5));
t683 = t782 * t687 - t779 * t688;
t724 = t782 * t750 - t779 * t751;
t700 = t724 * qJD(5) + t779 * t736 + t782 * t737;
t725 = t779 * t750 + t782 * t751;
t708 = -t724 * mrSges(6,1) + t725 * mrSges(6,2);
t772 = -qJD(3) + qJD(5);
t712 = -t772 * mrSges(6,2) + t724 * mrSges(6,3);
t769 = -qJDD(3) + qJDD(5);
t680 = m(6) * t683 + t769 * mrSges(6,1) - t700 * mrSges(6,3) - t725 * t708 + t772 * t712;
t684 = t779 * t687 + t782 * t688;
t699 = -t725 * qJD(5) + t782 * t736 - t779 * t737;
t713 = t772 * mrSges(6,1) - t725 * mrSges(6,3);
t681 = m(6) * t684 - t769 * mrSges(6,2) + t699 * mrSges(6,3) + t724 * t708 - t772 * t713;
t671 = t782 * t680 + t779 * t681;
t745 = -t750 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t790 = -m(5) * t693 + qJDD(3) * mrSges(5,1) + qJD(3) * t745 - t671;
t670 = t737 * mrSges(5,2) + t751 * t730 - t790;
t704 = Ifges(6,4) * t725 + Ifges(6,2) * t724 + Ifges(6,6) * t772;
t705 = Ifges(6,1) * t725 + Ifges(6,4) * t724 + Ifges(6,5) * t772;
t789 = -mrSges(6,1) * t683 + mrSges(6,2) * t684 - Ifges(6,5) * t700 - Ifges(6,6) * t699 - Ifges(6,3) * t769 - t725 * t704 + t724 * t705;
t744 = -qJD(3) * mrSges(5,1) + t751 * mrSges(5,2);
t800 = -t779 * t680 + t782 * t681;
t792 = m(5) * t692 + qJDD(3) * mrSges(5,3) + qJD(3) * t744 + t800;
t816 = t824 * qJD(3) - t825 * t750 + t751 * t836;
t817 = qJD(3) * t823 + t750 * t835 + t751 * t825;
t830 = t834 * qJDD(3) - t823 * t736 + t824 * t737 + t816 * t750 + t817 * t751 + mrSges(4,1) * t701 - mrSges(5,1) * t693 - mrSges(4,2) * t702 + mrSges(5,3) * t692 - pkin(3) * t670 - pkin(4) * t671 + qJ(4) * (-t736 * mrSges(5,2) - t750 * t730 + t792) + t789;
t826 = -mrSges(4,3) - mrSges(5,2);
t821 = mrSges(3,2) * t777;
t743 = qJD(3) * mrSges(4,1) - t751 * mrSges(4,3);
t815 = -t750 * mrSges(4,1) - t751 * mrSges(4,2) - t730;
t667 = m(4) * t702 - qJDD(3) * mrSges(4,2) - qJD(3) * t743 + t826 * t736 + t815 * t750 + t792;
t742 = -qJD(3) * mrSges(4,2) - t750 * mrSges(4,3);
t668 = m(4) * t701 + qJDD(3) * mrSges(4,1) + qJD(3) * t742 + t826 * t737 + t815 * t751 + t790;
t663 = t780 * t667 + t828 * t668;
t738 = -t777 * t752 + t804;
t794 = mrSges(3,3) * qJDD(1) + t785 * (-mrSges(3,1) * t778 + t821);
t661 = m(3) * t738 - t794 * t777 + t663;
t801 = t828 * t667 - t780 * t668;
t662 = m(3) * t739 + t794 * t778 + t801;
t802 = -t777 * t661 + t778 * t662;
t653 = m(2) * t757 - t785 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t802;
t748 = -qJDD(1) * pkin(1) - t785 * qJ(2) - t805;
t690 = -0.2e1 * qJD(4) * t751 + (t736 + t812) * pkin(3) + t837;
t686 = -t749 * pkin(7) + (-pkin(3) - pkin(4)) * t736 + (-pkin(3) * qJD(3) + t746 + t829) * t751 - t837;
t795 = -m(6) * t686 + t699 * mrSges(6,1) - t700 * mrSges(6,2) + t724 * t712 - t725 * t713;
t678 = m(5) * t690 + t736 * mrSges(5,1) - t737 * mrSges(5,3) - t751 * t744 + t750 * t745 + t795;
t788 = m(4) * t735 + t736 * mrSges(4,1) + t737 * mrSges(4,2) + t750 * t742 + t751 * t743 + t678;
t787 = -m(3) * t748 + mrSges(3,1) * t808 - t788 + (t770 * t785 + t820) * mrSges(3,3);
t673 = t787 + (mrSges(2,1) - t821) * qJDD(1) - t785 * mrSges(2,2) + m(2) * t756;
t819 = t781 * t653 + t783 * t673;
t655 = t778 * t661 + t777 * t662;
t818 = -qJD(3) * t834 + t750 * t823 - t751 * t824;
t803 = t783 * t653 - t781 * t673;
t798 = Ifges(3,1) * t777 + Ifges(3,4) * t778;
t797 = Ifges(3,4) * t777 + Ifges(3,2) * t778;
t796 = Ifges(3,5) * t777 + Ifges(3,6) * t778;
t703 = Ifges(6,5) * t725 + Ifges(6,6) * t724 + Ifges(6,3) * t772;
t674 = -mrSges(6,1) * t686 + mrSges(6,3) * t684 + Ifges(6,4) * t700 + Ifges(6,2) * t699 + Ifges(6,6) * t769 - t725 * t703 + t772 * t705;
t675 = mrSges(6,2) * t686 - mrSges(6,3) * t683 + Ifges(6,1) * t700 + Ifges(6,4) * t699 + Ifges(6,5) * t769 + t724 * t703 - t772 * t704;
t656 = -mrSges(4,1) * t735 - mrSges(5,1) * t690 + mrSges(5,2) * t692 + mrSges(4,3) * t702 - pkin(3) * t678 - pkin(4) * t795 - pkin(7) * t800 + t816 * qJD(3) + t823 * qJDD(3) - t782 * t674 - t779 * t675 + t736 * t835 + t825 * t737 + t818 * t751;
t657 = mrSges(4,2) * t735 + mrSges(5,2) * t693 - mrSges(4,3) * t701 - mrSges(5,3) * t690 - pkin(7) * t671 - qJ(4) * t678 - t817 * qJD(3) + t824 * qJDD(3) - t779 * t674 + t782 * t675 - t825 * t736 + t737 * t836 + t818 * t750;
t754 = t796 * qJD(1);
t648 = -mrSges(3,1) * t748 + mrSges(3,3) * t739 - pkin(2) * t788 + pkin(6) * t801 + t797 * qJDD(1) + t828 * t656 + t780 * t657 - t754 * t811;
t650 = t778 * qJD(1) * t754 + mrSges(3,2) * t748 - mrSges(3,3) * t738 - pkin(6) * t663 + t798 * qJDD(1) - t780 * t656 + t828 * t657;
t677 = mrSges(3,2) * t809 - t787;
t791 = mrSges(2,1) * t756 - mrSges(2,2) * t757 + Ifges(2,3) * qJDD(1) - pkin(1) * t677 + qJ(2) * t802 + t778 * t648 + t777 * t650;
t646 = mrSges(2,1) * g(3) + (Ifges(2,6) - t796) * qJDD(1) + mrSges(2,3) * t757 - mrSges(3,1) * t738 + mrSges(3,2) * t739 - pkin(2) * t663 - pkin(1) * t655 + (-t777 * t797 + t778 * t798 + Ifges(2,5)) * t785 - t830;
t645 = -mrSges(2,2) * g(3) - mrSges(2,3) * t756 + Ifges(2,5) * qJDD(1) - t785 * Ifges(2,6) - qJ(2) * t655 - t777 * t648 + t778 * t650;
t1 = [-m(1) * g(1) + t803; -m(1) * g(2) + t819; (-m(1) - m(2)) * g(3) + t655; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t819 + t783 * t645 - t781 * t646; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t803 + t781 * t645 + t783 * t646; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t791; t791; t677; t830; t670; -t789;];
tauJB = t1;

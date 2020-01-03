% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:16
% EndTime: 2019-12-31 19:29:20
% DurationCPUTime: 4.07s
% Computational Cost: add. (39029->293), mult. (90803->356), div. (0->0), fcn. (58062->8), ass. (0->119)
t829 = -2 * qJD(3);
t828 = Ifges(4,1) + Ifges(5,1);
t819 = Ifges(4,4) - Ifges(5,5);
t818 = Ifges(4,5) + Ifges(5,4);
t827 = Ifges(4,2) + Ifges(5,3);
t826 = -Ifges(5,2) - Ifges(4,3);
t817 = Ifges(4,6) - Ifges(5,6);
t780 = sin(qJ(2));
t783 = cos(qJ(2));
t805 = qJD(1) * qJD(2);
t759 = qJDD(1) * t783 - t780 * t805;
t809 = qJD(1) * t780;
t761 = qJD(2) * pkin(2) - qJ(3) * t809;
t777 = t783 ^ 2;
t786 = qJD(1) ^ 2;
t781 = sin(qJ(1));
t784 = cos(qJ(1));
t764 = t781 * g(1) - t784 * g(2);
t797 = qJDD(1) * pkin(1) + t764;
t711 = -pkin(2) * t759 + t761 * t809 - (qJ(3) * t777 + pkin(6)) * t786 + qJDD(3) - t797;
t758 = qJDD(1) * t780 + t783 * t805;
t778 = sin(pkin(8));
t816 = cos(pkin(8));
t732 = t816 * t758 + t778 * t759;
t808 = qJD(1) * t783;
t746 = t778 * t809 - t816 * t808;
t807 = qJD(2) * t746;
t825 = t711 + (-t732 + t807) * qJ(4);
t765 = -g(1) * t784 - g(2) * t781;
t753 = -pkin(1) * t786 + qJDD(1) * pkin(6) + t765;
t815 = t753 * t780;
t821 = pkin(2) * t786;
t707 = qJDD(2) * pkin(2) - qJ(3) * t758 - t815 + (qJ(3) * t805 + t780 * t821 - g(3)) * t783;
t736 = -g(3) * t780 + t783 * t753;
t708 = qJ(3) * t759 - qJD(2) * t761 - t777 * t821 + t736;
t747 = (t778 * t783 + t816 * t780) * qJD(1);
t689 = t816 * t707 - t778 * t708 + t747 * t829;
t690 = t778 * t707 + t816 * t708 + t746 * t829;
t731 = t758 * t778 - t816 * t759;
t738 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t747;
t725 = pkin(3) * t746 - qJ(4) * t747;
t785 = qJD(2) ^ 2;
t822 = 2 * qJD(4);
t685 = -pkin(3) * t785 + qJDD(2) * qJ(4) + qJD(2) * t822 - t746 * t725 + t690;
t739 = -qJD(2) * mrSges(5,1) + mrSges(5,2) * t747;
t686 = -qJDD(2) * pkin(3) - t785 * qJ(4) + t747 * t725 + qJDD(4) - t689;
t680 = (-t732 - t807) * pkin(7) + (t746 * t747 - qJDD(2)) * pkin(4) + t686;
t741 = -qJD(2) * pkin(4) - pkin(7) * t747;
t745 = t746 ^ 2;
t681 = -pkin(4) * t745 + pkin(7) * t731 + qJD(2) * t741 + t685;
t779 = sin(qJ(5));
t782 = cos(qJ(5));
t678 = t680 * t782 - t681 * t779;
t720 = t746 * t782 - t747 * t779;
t697 = qJD(5) * t720 + t731 * t779 + t732 * t782;
t721 = t746 * t779 + t747 * t782;
t703 = -mrSges(6,1) * t720 + mrSges(6,2) * t721;
t772 = -qJD(2) + qJD(5);
t712 = -mrSges(6,2) * t772 + mrSges(6,3) * t720;
t771 = -qJDD(2) + qJDD(5);
t675 = m(6) * t678 + mrSges(6,1) * t771 - t697 * mrSges(6,3) - t703 * t721 + t712 * t772;
t679 = t680 * t779 + t681 * t782;
t696 = -qJD(5) * t721 + t731 * t782 - t732 * t779;
t713 = mrSges(6,1) * t772 - mrSges(6,3) * t721;
t676 = m(6) * t679 - mrSges(6,2) * t771 + t696 * mrSges(6,3) + t703 * t720 - t713 * t772;
t799 = -t675 * t779 + t782 * t676;
t793 = m(5) * t685 + qJDD(2) * mrSges(5,3) + qJD(2) * t739 + t799;
t726 = mrSges(5,1) * t746 - mrSges(5,3) * t747;
t810 = -mrSges(4,1) * t746 - mrSges(4,2) * t747 - t726;
t820 = -mrSges(4,3) - mrSges(5,2);
t662 = m(4) * t690 - qJDD(2) * mrSges(4,2) - qJD(2) * t738 + t820 * t731 + t810 * t746 + t793;
t737 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t746;
t666 = t675 * t782 + t676 * t779;
t740 = -mrSges(5,2) * t746 + qJD(2) * mrSges(5,3);
t790 = -m(5) * t686 + qJDD(2) * mrSges(5,1) + qJD(2) * t740 - t666;
t663 = m(4) * t689 + qJDD(2) * mrSges(4,1) + qJD(2) * t737 + t820 * t732 + t810 * t747 + t790;
t658 = t778 * t662 + t816 * t663;
t665 = mrSges(5,2) * t732 + t726 * t747 - t790;
t735 = -g(3) * t783 - t815;
t749 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t780 + Ifges(3,2) * t783) * qJD(1);
t750 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t780 + Ifges(3,4) * t783) * qJD(1);
t699 = Ifges(6,4) * t721 + Ifges(6,2) * t720 + Ifges(6,6) * t772;
t700 = Ifges(6,1) * t721 + Ifges(6,4) * t720 + Ifges(6,5) * t772;
t789 = -mrSges(6,1) * t678 + mrSges(6,2) * t679 - Ifges(6,5) * t697 - Ifges(6,6) * t696 - Ifges(6,3) * t771 - t721 * t699 + t720 * t700;
t811 = t818 * qJD(2) - t819 * t746 + t828 * t747;
t813 = -t817 * qJD(2) + t827 * t746 - t819 * t747;
t823 = (t749 * t780 - t750 * t783) * qJD(1) + (Ifges(3,3) - t826) * qJDD(2) - t817 * t731 + t818 * t732 + t811 * t746 - t813 * t747 + mrSges(3,1) * t735 + mrSges(4,1) * t689 - mrSges(5,1) * t686 - mrSges(3,2) * t736 - mrSges(4,2) * t690 + mrSges(5,3) * t685 + Ifges(3,5) * t758 + Ifges(3,6) * t759 + pkin(2) * t658 - pkin(3) * t665 - pkin(4) * t666 + qJ(4) * (-mrSges(5,2) * t731 - t726 * t746 + t793) + t789;
t757 = (-mrSges(3,1) * t783 + mrSges(3,2) * t780) * qJD(1);
t763 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t808;
t656 = m(3) * t735 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t758 + qJD(2) * t763 - t757 * t809 + t658;
t762 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t809;
t800 = t816 * t662 - t778 * t663;
t657 = m(3) * t736 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t759 - qJD(2) * t762 + t757 * t808 + t800;
t801 = -t656 * t780 + t783 * t657;
t648 = m(2) * t765 - mrSges(2,1) * t786 - qJDD(1) * mrSges(2,2) + t801;
t688 = -0.2e1 * qJD(4) * t747 + (qJD(2) * t747 + t731) * pkin(3) + t825;
t683 = -pkin(7) * t745 + (-pkin(3) - pkin(4)) * t731 + (-pkin(3) * qJD(2) + t741 + t822) * t747 - t825;
t795 = -m(6) * t683 + t696 * mrSges(6,1) - t697 * mrSges(6,2) + t720 * t712 - t721 * t713;
t673 = m(5) * t688 + t731 * mrSges(5,1) - t732 * mrSges(5,3) - t747 * t739 + t746 * t740 + t795;
t672 = m(4) * t711 + t731 * mrSges(4,1) + t732 * mrSges(4,2) + t746 * t737 + t747 * t738 + t673;
t752 = -pkin(6) * t786 - t797;
t788 = -m(3) * t752 + t759 * mrSges(3,1) - t758 * mrSges(3,2) - t762 * t809 + t763 * t808 - t672;
t670 = m(2) * t764 + qJDD(1) * mrSges(2,1) - t786 * mrSges(2,2) + t788;
t814 = t781 * t648 + t784 * t670;
t650 = t783 * t656 + t780 * t657;
t812 = t826 * qJD(2) + t817 * t746 - t818 * t747;
t802 = t784 * t648 - t670 * t781;
t698 = Ifges(6,5) * t721 + Ifges(6,6) * t720 + Ifges(6,3) * t772;
t667 = -mrSges(6,1) * t683 + mrSges(6,3) * t679 + Ifges(6,4) * t697 + Ifges(6,2) * t696 + Ifges(6,6) * t771 - t698 * t721 + t700 * t772;
t668 = mrSges(6,2) * t683 - mrSges(6,3) * t678 + Ifges(6,1) * t697 + Ifges(6,4) * t696 + Ifges(6,5) * t771 + t698 * t720 - t699 * t772;
t651 = -mrSges(4,1) * t711 - mrSges(5,1) * t688 + mrSges(5,2) * t685 + mrSges(4,3) * t690 - pkin(3) * t673 - pkin(4) * t795 - pkin(7) * t799 + t811 * qJD(2) + t817 * qJDD(2) - t782 * t667 - t779 * t668 - t827 * t731 + t819 * t732 + t812 * t747;
t652 = mrSges(4,2) * t711 + mrSges(5,2) * t686 - mrSges(4,3) * t689 - mrSges(5,3) * t688 - pkin(7) * t666 - qJ(4) * t673 + t813 * qJD(2) + t818 * qJDD(2) - t667 * t779 + t668 * t782 - t819 * t731 + t828 * t732 + t812 * t746;
t748 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t780 + Ifges(3,6) * t783) * qJD(1);
t643 = -mrSges(3,1) * t752 + mrSges(3,3) * t736 + Ifges(3,4) * t758 + Ifges(3,2) * t759 + Ifges(3,6) * qJDD(2) - pkin(2) * t672 + qJ(3) * t800 + qJD(2) * t750 + t816 * t651 + t778 * t652 - t748 * t809;
t645 = mrSges(3,2) * t752 - mrSges(3,3) * t735 + Ifges(3,1) * t758 + Ifges(3,4) * t759 + Ifges(3,5) * qJDD(2) - qJ(3) * t658 - qJD(2) * t749 - t778 * t651 + t816 * t652 + t748 * t808;
t792 = mrSges(2,1) * t764 - mrSges(2,2) * t765 + Ifges(2,3) * qJDD(1) + pkin(1) * t788 + pkin(6) * t801 + t783 * t643 + t780 * t645;
t641 = mrSges(2,1) * g(3) + mrSges(2,3) * t765 + t786 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t650 - t823;
t640 = -mrSges(2,2) * g(3) - mrSges(2,3) * t764 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t786 - pkin(6) * t650 - t643 * t780 + t645 * t783;
t1 = [-m(1) * g(1) + t802; -m(1) * g(2) + t814; (-m(1) - m(2)) * g(3) + t650; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t814 + t784 * t640 - t781 * t641; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t802 + t781 * t640 + t784 * t641; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t792; t792; t823; t672; t665; -t789;];
tauJB = t1;

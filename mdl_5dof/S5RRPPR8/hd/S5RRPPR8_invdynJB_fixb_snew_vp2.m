% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:30
% EndTime: 2019-12-31 19:38:34
% DurationCPUTime: 3.92s
% Computational Cost: add. (37939->290), mult. (85763->357), div. (0->0), fcn. (50501->8), ass. (0->116)
t837 = Ifges(3,1) + Ifges(4,1);
t829 = Ifges(3,4) - Ifges(4,5);
t828 = Ifges(3,5) + Ifges(4,4);
t836 = Ifges(3,2) + Ifges(4,3);
t827 = Ifges(3,6) - Ifges(4,6);
t835 = Ifges(3,3) + Ifges(4,2);
t797 = sin(qJ(2));
t800 = cos(qJ(2));
t763 = (-mrSges(4,1) * t800 - mrSges(4,3) * t797) * qJD(1);
t819 = qJD(1) * qJD(2);
t816 = t800 * t819;
t765 = qJDD(1) * t797 + t816;
t798 = sin(qJ(1));
t801 = cos(qJ(1));
t776 = -g(1) * t801 - g(2) * t798;
t803 = qJD(1) ^ 2;
t753 = -pkin(1) * t803 + qJDD(1) * pkin(6) + t776;
t732 = -g(3) * t797 + t800 * t753;
t762 = (-pkin(2) * t800 - qJ(3) * t797) * qJD(1);
t802 = qJD(2) ^ 2;
t820 = qJD(1) * t800;
t831 = 2 * qJD(3);
t715 = -pkin(2) * t802 + qJDD(2) * qJ(3) + qJD(2) * t831 + t762 * t820 + t732;
t817 = t797 * t819;
t766 = qJDD(1) * t800 - t817;
t821 = qJD(1) * t797;
t770 = -qJD(2) * pkin(3) - qJ(4) * t821;
t826 = t800 ^ 2 * t803;
t711 = -pkin(3) * t826 - qJ(4) * t766 + qJD(2) * t770 + t715;
t731 = -t800 * g(3) - t797 * t753;
t716 = -qJDD(2) * pkin(2) - t802 * qJ(3) + t762 * t821 + qJDD(3) - t731;
t712 = (-t765 + t816) * qJ(4) + (-t797 * t800 * t803 - qJDD(2)) * pkin(3) + t716;
t793 = sin(pkin(8));
t794 = cos(pkin(8));
t745 = (-t793 * t800 + t794 * t797) * qJD(1);
t689 = -0.2e1 * qJD(4) * t745 - t793 * t711 + t712 * t794;
t730 = t765 * t794 - t766 * t793;
t744 = (-t793 * t797 - t794 * t800) * qJD(1);
t687 = (-qJD(2) * t744 - t730) * pkin(7) + (t744 * t745 - qJDD(2)) * pkin(4) + t689;
t690 = 0.2e1 * qJD(4) * t744 + t711 * t794 + t712 * t793;
t729 = -t765 * t793 - t766 * t794;
t735 = -qJD(2) * pkin(4) - pkin(7) * t745;
t743 = t744 ^ 2;
t688 = -pkin(4) * t743 + pkin(7) * t729 + qJD(2) * t735 + t690;
t796 = sin(qJ(5));
t799 = cos(qJ(5));
t685 = t687 * t799 - t688 * t796;
t722 = t744 * t799 - t745 * t796;
t699 = qJD(5) * t722 + t729 * t796 + t730 * t799;
t723 = t744 * t796 + t745 * t799;
t710 = -mrSges(6,1) * t722 + mrSges(6,2) * t723;
t786 = -qJD(2) + qJD(5);
t717 = -mrSges(6,2) * t786 + mrSges(6,3) * t722;
t785 = -qJDD(2) + qJDD(5);
t681 = m(6) * t685 + mrSges(6,1) * t785 - t699 * mrSges(6,3) - t710 * t723 + t717 * t786;
t686 = t687 * t796 + t688 * t799;
t698 = -qJD(5) * t723 + t729 * t799 - t730 * t796;
t718 = mrSges(6,1) * t786 - mrSges(6,3) * t723;
t682 = m(6) * t686 - mrSges(6,2) * t785 + t698 * mrSges(6,3) + t710 * t722 - t718 * t786;
t671 = t681 * t799 + t682 * t796;
t726 = -mrSges(5,1) * t744 + mrSges(5,2) * t745;
t733 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t744;
t669 = m(5) * t689 - qJDD(2) * mrSges(5,1) - mrSges(5,3) * t730 - qJD(2) * t733 - t726 * t745 + t671;
t734 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t745;
t812 = -t681 * t796 + t682 * t799;
t670 = m(5) * t690 + qJDD(2) * mrSges(5,2) + mrSges(5,3) * t729 + qJD(2) * t734 + t726 * t744 + t812;
t667 = t794 * t669 + t793 * t670;
t774 = mrSges(4,2) * t820 + qJD(2) * mrSges(4,3);
t807 = -m(4) * t716 + qJDD(2) * mrSges(4,1) + qJD(2) * t774 - t667;
t666 = t765 * mrSges(4,2) + t763 * t821 - t807;
t720 = Ifges(5,4) * t745 + Ifges(5,2) * t744 - Ifges(5,6) * qJD(2);
t721 = Ifges(5,1) * t745 + Ifges(5,4) * t744 - Ifges(5,5) * qJD(2);
t703 = Ifges(6,4) * t723 + Ifges(6,2) * t722 + Ifges(6,6) * t786;
t704 = Ifges(6,1) * t723 + Ifges(6,4) * t722 + Ifges(6,5) * t786;
t806 = -mrSges(6,1) * t685 + mrSges(6,2) * t686 - Ifges(6,5) * t699 - Ifges(6,6) * t698 - Ifges(6,3) * t785 - t703 * t723 + t722 * t704;
t772 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t821;
t813 = -t793 * t669 + t670 * t794;
t809 = m(4) * t715 + qJDD(2) * mrSges(4,3) + qJD(2) * t772 + t763 * t820 + t813;
t822 = t828 * qJD(2) + (t797 * t837 + t800 * t829) * qJD(1);
t823 = -t827 * qJD(2) + (-t797 * t829 - t800 * t836) * qJD(1);
t834 = -(t797 * t823 + t800 * t822) * qJD(1) + (Ifges(5,3) + t835) * qJDD(2) + t828 * t765 + t827 * t766 + mrSges(3,1) * t731 - mrSges(4,1) * t716 - mrSges(5,1) * t689 - mrSges(3,2) * t732 + mrSges(5,2) * t690 + mrSges(4,3) * t715 - Ifges(5,5) * t730 - Ifges(5,6) * t729 - pkin(2) * t666 - pkin(3) * t667 - pkin(4) * t671 + qJ(3) * (mrSges(4,2) * t766 + t809) - t745 * t720 + t744 * t721 + t806;
t830 = mrSges(3,3) + mrSges(4,2);
t764 = (-mrSges(3,1) * t800 + mrSges(3,2) * t797) * qJD(1);
t771 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t821;
t663 = m(3) * t732 - qJDD(2) * mrSges(3,2) - qJD(2) * t771 + t764 * t820 + t766 * t830 + t809;
t773 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t820;
t664 = m(3) * t731 + qJDD(2) * mrSges(3,1) + qJD(2) * t773 - t830 * t765 + (-t763 - t764) * t821 + t807;
t814 = t663 * t800 - t664 * t797;
t655 = m(2) * t776 - mrSges(2,1) * t803 - qJDD(1) * mrSges(2,2) + t814;
t775 = t798 * g(1) - t801 * g(2);
t752 = -qJDD(1) * pkin(1) - t803 * pkin(6) - t775;
t810 = -t766 * pkin(2) + t752 + (-t765 - t816) * qJ(3);
t701 = -pkin(2) * t817 + t766 * pkin(3) - qJ(4) * t826 + qJDD(4) - t810 + (t770 + t831) * t821;
t692 = -pkin(4) * t729 - pkin(7) * t743 + t735 * t745 + t701;
t811 = m(6) * t692 - t698 * mrSges(6,1) + t699 * mrSges(6,2) - t717 * t722 + t718 * t723;
t683 = m(5) * t701 - t729 * mrSges(5,1) + t730 * mrSges(5,2) - t744 * t733 + t745 * t734 + t811;
t713 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t821 + t810;
t677 = m(4) * t713 - mrSges(4,1) * t766 - t765 * mrSges(4,3) - t772 * t821 - t774 * t820 - t683;
t805 = -m(3) * t752 + t766 * mrSges(3,1) - mrSges(3,2) * t765 - t771 * t821 + t773 * t820 - t677;
t675 = m(2) * t775 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t803 + t805;
t825 = t655 * t798 + t675 * t801;
t657 = t663 * t797 + t664 * t800;
t824 = t835 * qJD(2) + (t797 * t828 + t800 * t827) * qJD(1);
t815 = t655 * t801 - t675 * t798;
t702 = Ifges(6,5) * t723 + Ifges(6,6) * t722 + Ifges(6,3) * t786;
t672 = -mrSges(6,1) * t692 + mrSges(6,3) * t686 + Ifges(6,4) * t699 + Ifges(6,2) * t698 + Ifges(6,6) * t785 - t702 * t723 + t704 * t786;
t673 = mrSges(6,2) * t692 - mrSges(6,3) * t685 + Ifges(6,1) * t699 + Ifges(6,4) * t698 + Ifges(6,5) * t785 + t702 * t722 - t703 * t786;
t719 = Ifges(5,5) * t745 + Ifges(5,6) * t744 - Ifges(5,3) * qJD(2);
t658 = -mrSges(5,1) * t701 + mrSges(5,3) * t690 + Ifges(5,4) * t730 + Ifges(5,2) * t729 - Ifges(5,6) * qJDD(2) - pkin(4) * t811 + pkin(7) * t812 - qJD(2) * t721 + t799 * t672 + t796 * t673 - t745 * t719;
t659 = mrSges(5,2) * t701 - mrSges(5,3) * t689 + Ifges(5,1) * t730 + Ifges(5,4) * t729 - Ifges(5,5) * qJDD(2) - pkin(7) * t671 + qJD(2) * t720 - t672 * t796 + t673 * t799 + t719 * t744;
t650 = -mrSges(3,1) * t752 - mrSges(4,1) * t713 + mrSges(4,2) * t715 + mrSges(3,3) * t732 - pkin(2) * t677 + pkin(3) * t683 - qJ(4) * t813 + t822 * qJD(2) + t827 * qJDD(2) - t794 * t658 - t793 * t659 + t829 * t765 + t766 * t836 - t824 * t821;
t652 = mrSges(3,2) * t752 + mrSges(4,2) * t716 - mrSges(3,3) * t731 - mrSges(4,3) * t713 - qJ(3) * t677 - qJ(4) * t667 + t823 * qJD(2) + t828 * qJDD(2) - t793 * t658 + t794 * t659 + t765 * t837 + t829 * t766 + t824 * t820;
t808 = mrSges(2,1) * t775 - mrSges(2,2) * t776 + Ifges(2,3) * qJDD(1) + pkin(1) * t805 + pkin(6) * t814 + t650 * t800 + t652 * t797;
t648 = mrSges(2,1) * g(3) + mrSges(2,3) * t776 + t803 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t657 - t834;
t647 = -mrSges(2,2) * g(3) - mrSges(2,3) * t775 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t803 - pkin(6) * t657 - t650 * t797 + t652 * t800;
t1 = [-m(1) * g(1) + t815; -m(1) * g(2) + t825; (-m(1) - m(2)) * g(3) + t657; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t825 + t647 * t801 - t648 * t798; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t815 + t798 * t647 + t801 * t648; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t808; t808; t834; t666; t683; -t806;];
tauJB = t1;

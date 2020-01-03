% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:51
% EndTime: 2019-12-31 20:57:55
% DurationCPUTime: 3.24s
% Computational Cost: add. (22048->268), mult. (46165->315), div. (0->0), fcn. (28015->6), ass. (0->105)
t806 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t836 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t825 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t782 = sin(qJ(3));
t785 = cos(qJ(2));
t808 = qJD(1) * t785;
t783 = sin(qJ(2));
t809 = qJD(1) * t783;
t818 = cos(qJ(3));
t751 = t782 * t809 - t818 * t808;
t752 = (t782 * t785 + t818 * t783) * qJD(1);
t779 = qJD(2) + qJD(3);
t835 = t836 * t751 + t806 * t752 + t825 * t779;
t834 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t827 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t739 = t779 * mrSges(6,2) + t751 * mrSges(6,3);
t778 = qJDD(2) + qJDD(3);
t807 = qJD(1) * qJD(2);
t760 = t783 * qJDD(1) + t785 * t807;
t784 = sin(qJ(1));
t786 = cos(qJ(1));
t767 = -t786 * g(1) - t784 * g(2);
t787 = qJD(1) ^ 2;
t754 = -t787 * pkin(1) + qJDD(1) * pkin(6) + t767;
t813 = t783 * t754;
t816 = pkin(2) * t787;
t697 = qJDD(2) * pkin(2) - t760 * pkin(7) - t813 + (pkin(7) * t807 + t783 * t816 - g(3)) * t785;
t738 = -t783 * g(3) + t785 * t754;
t761 = t785 * qJDD(1) - t783 * t807;
t765 = qJD(2) * pkin(2) - pkin(7) * t809;
t781 = t785 ^ 2;
t698 = t761 * pkin(7) - qJD(2) * t765 - t781 * t816 + t738;
t689 = t818 * t697 - t782 * t698;
t731 = t751 * pkin(3) - t752 * qJ(4);
t777 = t779 ^ 2;
t687 = -t778 * pkin(3) - t777 * qJ(4) + t752 * t731 + qJDD(4) - t689;
t715 = -t751 * qJD(3) + t818 * t760 + t782 * t761;
t814 = t751 * t779;
t820 = -0.2e1 * t752;
t678 = qJD(5) * t820 + (-t715 - t814) * qJ(5) + (t751 * t752 - t778) * pkin(4) + t687;
t733 = -t751 * mrSges(6,1) + t752 * mrSges(6,2);
t797 = -m(6) * t678 + t715 * mrSges(6,3) + t752 * t733;
t675 = -t778 * mrSges(6,1) - t779 * t739 - t797;
t732 = t751 * mrSges(5,1) - t752 * mrSges(5,3);
t745 = -t751 * mrSges(5,2) + t779 * mrSges(5,3);
t822 = -m(5) * t687 + t778 * mrSges(5,1) + t779 * t745;
t671 = t715 * mrSges(5,2) + t752 * t732 + t675 - t822;
t690 = t782 * t697 + t818 * t698;
t819 = 2 * qJD(4);
t686 = -t777 * pkin(3) + t778 * qJ(4) - t751 * t731 + t779 * t819 + t690;
t714 = t752 * qJD(3) + t782 * t760 - t818 * t761;
t741 = -t779 * pkin(4) - t752 * qJ(5);
t747 = t751 ^ 2;
t681 = -t747 * pkin(4) + t714 * qJ(5) + 0.2e1 * qJD(5) * t751 + t779 * t741 + t686;
t742 = -t779 * mrSges(6,1) - t752 * mrSges(6,3);
t744 = -t779 * mrSges(5,1) + t752 * mrSges(5,2);
t805 = m(6) * t681 + t714 * mrSges(6,3) + t751 * t733;
t794 = m(5) * t686 + t778 * mrSges(5,3) + t779 * t744 + t805;
t824 = -t806 * t751 + t834 * t752 + t827 * t779;
t826 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t833 = -mrSges(5,1) * t687 - mrSges(6,1) * t678 - mrSges(4,2) * t690 - pkin(3) * t671 - pkin(4) * t675 + qJ(4) * (t779 * t742 + t794) + mrSges(6,2) * t681 + mrSges(5,3) * t686 + mrSges(4,1) * t689 + (qJ(4) * mrSges(6,2) - t826) * t778 + (-qJ(4) * t732 + t824) * t751 + t827 * t715 + (-qJ(4) * mrSges(5,2) - t825) * t714 + t835 * t752;
t766 = t784 * g(1) - t786 * g(2);
t796 = qJDD(1) * pkin(1) + t766;
t716 = -t761 * pkin(2) + t765 * t809 - (pkin(7) * t781 + pkin(6)) * t787 - t796;
t829 = t716 + (-t715 + t814) * qJ(4);
t740 = -t779 * mrSges(4,2) - t751 * mrSges(4,3);
t810 = -t751 * mrSges(4,1) - t752 * mrSges(4,2) - t732;
t815 = -mrSges(4,3) - mrSges(5,2);
t665 = m(4) * t689 + (t739 + t740) * t779 + (mrSges(4,1) + mrSges(6,1)) * t778 + t810 * t752 + t815 * t715 + t797 + t822;
t743 = t779 * mrSges(4,1) - t752 * mrSges(4,3);
t668 = m(4) * t690 + (t742 - t743) * t779 + (-mrSges(4,2) + mrSges(6,2)) * t778 + t810 * t751 + t815 * t714 + t794;
t660 = t818 * t665 + t782 * t668;
t737 = -t785 * g(3) - t813;
t749 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t783 + Ifges(3,2) * t785) * qJD(1);
t750 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t783 + Ifges(3,4) * t785) * qJD(1);
t828 = mrSges(3,1) * t737 - mrSges(3,2) * t738 + Ifges(3,5) * t760 + Ifges(3,6) * t761 + Ifges(3,3) * qJDD(2) + pkin(2) * t660 + (t783 * t749 - t785 * t750) * qJD(1) + t833;
t759 = (-mrSges(3,1) * t785 + mrSges(3,2) * t783) * qJD(1);
t764 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t808;
t658 = m(3) * t737 + qJDD(2) * mrSges(3,1) - t760 * mrSges(3,3) + qJD(2) * t764 - t759 * t809 + t660;
t763 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t809;
t799 = -t782 * t665 + t818 * t668;
t659 = m(3) * t738 - qJDD(2) * mrSges(3,2) + t761 * mrSges(3,3) - qJD(2) * t763 + t759 * t808 + t799;
t800 = -t783 * t658 + t785 * t659;
t650 = m(2) * t767 - t787 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t800;
t753 = -t787 * pkin(6) - t796;
t677 = -t747 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t714 + (-pkin(3) * t779 + t741 + t819) * t752 - t829;
t674 = m(6) * t677 - t714 * mrSges(6,1) + t715 * mrSges(6,2) - t751 * t739 + t752 * t742;
t683 = qJD(4) * t820 + (t752 * t779 + t714) * pkin(3) + t829;
t669 = m(5) * t683 + t714 * mrSges(5,1) - t715 * mrSges(5,3) - t752 * t744 + t751 * t745 - t674;
t791 = m(4) * t716 + t714 * mrSges(4,1) + t715 * mrSges(4,2) + t751 * t740 + t752 * t743 + t669;
t790 = -m(3) * t753 + t761 * mrSges(3,1) - t760 * mrSges(3,2) - t763 * t809 + t764 * t808 - t791;
t662 = m(2) * t766 + qJDD(1) * mrSges(2,1) - t787 * mrSges(2,2) + t790;
t812 = t784 * t650 + t786 * t662;
t652 = t785 * t658 + t783 * t659;
t804 = t825 * t751 - t827 * t752 + t826 * t779;
t801 = t786 * t650 - t784 * t662;
t653 = -mrSges(4,1) * t716 + mrSges(4,3) * t690 - mrSges(5,1) * t683 + mrSges(5,2) * t686 + mrSges(6,1) * t677 - mrSges(6,3) * t681 + pkin(4) * t674 - qJ(5) * t805 - pkin(3) * t669 + (-qJ(5) * t742 + t824) * t779 + (-qJ(5) * mrSges(6,2) + t825) * t778 + t804 * t752 + t806 * t715 + t836 * t714;
t654 = mrSges(4,2) * t716 + mrSges(5,2) * t687 + mrSges(6,2) * t677 - mrSges(4,3) * t689 - mrSges(5,3) * t683 - mrSges(6,3) * t678 - qJ(4) * t669 - qJ(5) * t675 - t806 * t714 + t834 * t715 + t804 * t751 + t827 * t778 - t835 * t779;
t748 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t783 + Ifges(3,6) * t785) * qJD(1);
t645 = -mrSges(3,1) * t753 + mrSges(3,3) * t738 + Ifges(3,4) * t760 + Ifges(3,2) * t761 + Ifges(3,6) * qJDD(2) - pkin(2) * t791 + pkin(7) * t799 + qJD(2) * t750 + t818 * t653 + t782 * t654 - t748 * t809;
t647 = mrSges(3,2) * t753 - mrSges(3,3) * t737 + Ifges(3,1) * t760 + Ifges(3,4) * t761 + Ifges(3,5) * qJDD(2) - pkin(7) * t660 - qJD(2) * t749 - t782 * t653 + t818 * t654 + t748 * t808;
t792 = mrSges(2,1) * t766 - mrSges(2,2) * t767 + Ifges(2,3) * qJDD(1) + pkin(1) * t790 + pkin(6) * t800 + t785 * t645 + t783 * t647;
t643 = mrSges(2,1) * g(3) + mrSges(2,3) * t767 + t787 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t652 - t828;
t642 = -mrSges(2,2) * g(3) - mrSges(2,3) * t766 + Ifges(2,5) * qJDD(1) - t787 * Ifges(2,6) - pkin(6) * t652 - t783 * t645 + t785 * t647;
t1 = [-m(1) * g(1) + t801; -m(1) * g(2) + t812; (-m(1) - m(2)) * g(3) + t652; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t812 + t786 * t642 - t784 * t643; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t801 + t784 * t642 + t786 * t643; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t792; t792; t828; t833; t671; t674;];
tauJB = t1;

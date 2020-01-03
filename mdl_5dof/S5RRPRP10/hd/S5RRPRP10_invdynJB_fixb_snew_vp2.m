% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:40
% EndTime: 2019-12-31 20:09:44
% DurationCPUTime: 2.52s
% Computational Cost: add. (17343->275), mult. (35311->314), div. (0->0), fcn. (18468->6), ass. (0->111)
t818 = Ifges(5,4) + Ifges(6,4);
t835 = Ifges(5,2) + Ifges(6,2);
t829 = Ifges(5,6) + Ifges(6,6);
t834 = -2 * qJD(3);
t833 = Ifges(3,1) + Ifges(4,2);
t832 = Ifges(5,1) + Ifges(6,1);
t819 = Ifges(3,4) + Ifges(4,6);
t817 = Ifges(3,5) - Ifges(4,4);
t831 = Ifges(5,5) + Ifges(6,5);
t830 = Ifges(3,2) + Ifges(4,3);
t816 = Ifges(3,6) - Ifges(4,5);
t828 = Ifges(3,3) + Ifges(4,1);
t827 = Ifges(5,3) + Ifges(6,3);
t777 = sin(qJ(4));
t780 = cos(qJ(4));
t781 = cos(qJ(2));
t806 = qJD(1) * t781;
t749 = -t777 * qJD(2) - t780 * t806;
t750 = t780 * qJD(2) - t777 * t806;
t778 = sin(qJ(2));
t805 = t778 * qJD(1);
t767 = qJD(4) + t805;
t826 = t835 * t749 + t818 * t750 + t829 * t767;
t804 = qJD(1) * qJD(2);
t799 = t778 * t804;
t755 = t781 * qJDD(1) - t799;
t714 = -t750 * qJD(4) - t777 * qJDD(2) - t780 * t755;
t719 = -t767 * mrSges(6,2) + t749 * mrSges(6,3);
t779 = sin(qJ(1));
t782 = cos(qJ(1));
t765 = -t782 * g(1) - t779 * g(2);
t784 = qJD(1) ^ 2;
t738 = -t784 * pkin(1) + qJDD(1) * pkin(6) + t765;
t725 = -t778 * g(3) + t781 * t738;
t751 = (-pkin(2) * t781 - qJ(3) * t778) * qJD(1);
t783 = qJD(2) ^ 2;
t696 = t783 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t834 - t751 * t806 - t725;
t763 = pkin(3) * t805 - qJD(2) * pkin(7);
t776 = t781 ^ 2;
t691 = -t776 * t784 * pkin(7) + t755 * pkin(3) + qJD(2) * t763 - t696;
t721 = t767 * pkin(4) - t750 * qJ(5);
t747 = t749 ^ 2;
t685 = -t714 * pkin(4) - t747 * qJ(5) + t750 * t721 + qJDD(5) + t691;
t715 = t749 * qJD(4) + t780 * qJDD(2) - t777 * t755;
t722 = t767 * mrSges(6,1) - t750 * mrSges(6,3);
t800 = m(6) * t685 + t715 * mrSges(6,2) + t750 * t722;
t675 = -t714 * mrSges(6,1) - t749 * t719 + t800;
t798 = t781 * t804;
t754 = t778 * qJDD(1) + t798;
t764 = t779 * g(1) - t782 * g(2);
t795 = -qJDD(1) * pkin(1) - t764;
t789 = pkin(2) * t799 + t805 * t834 + (-t754 - t798) * qJ(3) + t795;
t687 = -t763 * t805 + (-pkin(3) * t776 - pkin(6)) * t784 + (-pkin(2) - pkin(7)) * t755 + t789;
t724 = -t781 * g(3) - t778 * t738;
t697 = -qJDD(2) * pkin(2) - t783 * qJ(3) + t751 * t805 + qJDD(3) - t724;
t692 = (-t778 * t781 * t784 - qJDD(2)) * pkin(7) + (t754 - t798) * pkin(3) + t697;
t683 = t780 * t687 + t777 * t692;
t680 = -t747 * pkin(4) + t714 * qJ(5) + 0.2e1 * qJD(5) * t749 - t767 * t721 + t683;
t748 = qJDD(4) + t754;
t717 = -t749 * mrSges(6,1) + t750 * mrSges(6,2);
t801 = m(6) * t680 + t714 * mrSges(6,3) + t749 * t717;
t811 = -t818 * t749 - t832 * t750 - t831 * t767;
t812 = -t829 * t749 - t831 * t750 - t827 * t767;
t655 = -mrSges(5,1) * t691 + mrSges(5,3) * t683 - mrSges(6,1) * t685 + mrSges(6,3) * t680 - pkin(4) * t675 + qJ(5) * t801 + (-qJ(5) * t722 - t811) * t767 + t812 * t750 + (-qJ(5) * mrSges(6,2) + t829) * t748 + t818 * t715 + t835 * t714;
t752 = (mrSges(4,2) * t781 - mrSges(4,3) * t778) * qJD(1);
t761 = -mrSges(4,1) * t806 - qJD(2) * mrSges(4,3);
t682 = -t777 * t687 + t780 * t692;
t718 = -t749 * mrSges(5,1) + t750 * mrSges(5,2);
t720 = -t767 * mrSges(5,2) + t749 * mrSges(5,3);
t678 = -0.2e1 * qJD(5) * t750 + (t749 * t767 - t715) * qJ(5) + (t749 * t750 + t748) * pkin(4) + t682;
t802 = m(6) * t678 + t748 * mrSges(6,1) + t767 * t719;
t666 = m(5) * t682 + t748 * mrSges(5,1) + t767 * t720 + (-t717 - t718) * t750 + (-mrSges(5,3) - mrSges(6,3)) * t715 + t802;
t723 = t767 * mrSges(5,1) - t750 * mrSges(5,3);
t671 = m(5) * t683 + t714 * mrSges(5,3) + t749 * t718 + (-t722 - t723) * t767 + (-mrSges(5,2) - mrSges(6,2)) * t748 + t801;
t664 = t780 * t666 + t777 * t671;
t791 = -m(4) * t697 - t754 * mrSges(4,1) - t664;
t662 = qJDD(2) * mrSges(4,2) + qJD(2) * t761 + t752 * t805 - t791;
t674 = -t715 * mrSges(6,3) - t750 * t717 + t802;
t663 = mrSges(5,2) * t691 + mrSges(6,2) * t685 - mrSges(5,3) * t682 - mrSges(6,3) * t678 - qJ(5) * t674 + t818 * t714 + t832 * t715 + t831 * t748 - t812 * t749 - t826 * t767;
t762 = mrSges(4,1) * t805 + qJD(2) * mrSges(4,2);
t823 = -m(5) * t691 - t715 * mrSges(5,2) + (mrSges(5,1) + mrSges(6,1)) * t714 - t750 * t723 + (t719 + t720) * t749 - t800;
t788 = -m(4) * t696 + qJDD(2) * mrSges(4,3) + qJD(2) * t762 + t752 * t806 - t823;
t807 = t817 * qJD(2) + (t833 * t778 + t819 * t781) * qJD(1);
t808 = t816 * qJD(2) + (t819 * t778 + t830 * t781) * qJD(1);
t825 = (t808 * t778 - t807 * t781) * qJD(1) + t828 * qJDD(2) + t817 * t754 + t816 * t755 + mrSges(3,1) * t724 - mrSges(3,2) * t725 + mrSges(4,2) * t697 - mrSges(4,3) * t696 - pkin(2) * t662 - pkin(7) * t664 + qJ(3) * (t755 * mrSges(4,1) + t788) - t777 * t655 + t780 * t663;
t821 = t784 * pkin(6);
t753 = (-mrSges(3,1) * t781 + mrSges(3,2) * t778) * qJD(1);
t760 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t806;
t660 = m(3) * t724 - t754 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t760 - t761) * qJD(2) + (-t752 - t753) * t805 + t791;
t759 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t805;
t670 = t753 * t806 + t788 + (mrSges(3,3) + mrSges(4,1)) * t755 - qJDD(2) * mrSges(3,2) - qJD(2) * t759 + m(3) * t725;
t796 = -t778 * t660 + t781 * t670;
t652 = m(2) * t765 - t784 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t796;
t737 = t795 - t821;
t693 = -t755 * pkin(2) + t789 - t821;
t813 = -t777 * t666 + t780 * t671;
t793 = -m(4) * t693 - t755 * mrSges(4,2) + t762 * t805 - t813;
t786 = -m(3) * t737 + t760 * t806 + t755 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t754 + (-t759 * t778 - t761 * t781) * qJD(1) + t793;
t657 = m(2) * t764 + qJDD(1) * mrSges(2,1) - t784 * mrSges(2,2) + t786;
t814 = t779 * t652 + t782 * t657;
t654 = t781 * t660 + t778 * t670;
t809 = t828 * qJD(2) + (t817 * t778 + t816 * t781) * qJD(1);
t797 = t782 * t652 - t779 * t657;
t661 = -t754 * mrSges(4,3) + t761 * t806 - t793;
t647 = -mrSges(3,1) * t737 - mrSges(4,1) * t696 + mrSges(4,2) * t693 + mrSges(3,3) * t725 - pkin(2) * t661 - pkin(3) * t823 - pkin(7) * t813 + t807 * qJD(2) + t816 * qJDD(2) - t780 * t655 - t777 * t663 + t819 * t754 + t830 * t755 - t809 * t805;
t787 = mrSges(5,1) * t682 + mrSges(6,1) * t678 - mrSges(5,2) * t683 - mrSges(6,2) * t680 + pkin(4) * t674 + t829 * t714 + t831 * t715 + t827 * t748 + t811 * t749 + t826 * t750;
t649 = mrSges(4,1) * t697 + mrSges(3,2) * t737 - mrSges(3,3) * t724 - mrSges(4,3) * t693 + pkin(3) * t664 - qJ(3) * t661 - t808 * qJD(2) + t817 * qJDD(2) + t833 * t754 + t819 * t755 + t809 * t806 + t787;
t790 = mrSges(2,1) * t764 - mrSges(2,2) * t765 + Ifges(2,3) * qJDD(1) + pkin(1) * t786 + pkin(6) * t796 + t781 * t647 + t778 * t649;
t645 = mrSges(2,1) * g(3) + mrSges(2,3) * t765 + t784 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t654 - t825;
t644 = -mrSges(2,2) * g(3) - mrSges(2,3) * t764 + Ifges(2,5) * qJDD(1) - t784 * Ifges(2,6) - pkin(6) * t654 - t778 * t647 + t781 * t649;
t1 = [-m(1) * g(1) + t797; -m(1) * g(2) + t814; (-m(1) - m(2)) * g(3) + t654; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t814 + t782 * t644 - t779 * t645; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t797 + t779 * t644 + t782 * t645; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t790; t790; t825; t662; t787; t675;];
tauJB = t1;

% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR16
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR16_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:45:19
% EndTime: 2019-12-31 20:45:26
% DurationCPUTime: 7.09s
% Computational Cost: add. (77808->311), mult. (176131->390), div. (0->0), fcn. (123008->10), ass. (0->137)
t883 = -2 * qJD(3);
t882 = Ifges(3,1) + Ifges(4,2);
t874 = Ifges(3,4) + Ifges(4,6);
t873 = Ifges(3,5) - Ifges(4,4);
t881 = Ifges(3,2) + Ifges(4,3);
t872 = Ifges(3,6) - Ifges(4,5);
t880 = Ifges(3,3) + Ifges(4,1);
t831 = cos(pkin(5));
t824 = t831 * qJD(1) + qJD(2);
t834 = sin(qJ(2));
t830 = sin(pkin(5));
t861 = qJD(1) * t830;
t854 = t834 * t861;
t879 = (pkin(2) * t824 + t883) * t854;
t835 = sin(qJ(1));
t839 = cos(qJ(1));
t819 = t835 * g(1) - t839 * g(2);
t840 = qJD(1) ^ 2;
t877 = pkin(7) * t830;
t803 = qJDD(1) * pkin(1) + t840 * t877 + t819;
t820 = -t839 * g(1) - t835 * g(2);
t858 = qJDD(1) * t830;
t804 = -t840 * pkin(1) + pkin(7) * t858 + t820;
t838 = cos(qJ(2));
t868 = t831 * t834;
t870 = t830 * t834;
t769 = -g(3) * t870 + t803 * t868 + t838 * t804;
t805 = (-pkin(2) * t838 - qJ(3) * t834) * t861;
t822 = t824 ^ 2;
t823 = t831 * qJDD(1) + qJDD(2);
t860 = qJD(1) * t838;
t855 = t830 * t860;
t747 = t822 * pkin(2) - t823 * qJ(3) - t805 * t855 + t824 * t883 - t769;
t878 = -pkin(2) - pkin(8);
t876 = t831 * g(3);
t875 = mrSges(3,1) - mrSges(4,2);
t871 = t830 ^ 2 * t840;
t869 = t830 * t838;
t867 = t831 * t838;
t782 = -t830 * t803 - t876;
t799 = t824 * mrSges(3,1) - mrSges(3,3) * t854;
t800 = -t824 * mrSges(3,2) + mrSges(3,3) * t855;
t802 = mrSges(4,1) * t854 + t824 * mrSges(4,2);
t809 = (qJD(2) * t860 + qJDD(1) * t834) * t830;
t810 = -qJD(2) * t854 + t838 * t858;
t748 = -t810 * pkin(2) + (-t824 * t855 - t809) * qJ(3) + t782 + t879;
t801 = -mrSges(4,1) * t855 - t824 * mrSges(4,3);
t808 = pkin(3) * t854 - t824 * pkin(8);
t857 = t838 ^ 2 * t871;
t740 = -pkin(3) * t857 - t876 - t809 * qJ(3) + t878 * t810 + (-t803 + (-qJ(3) * t824 * t838 - t808 * t834) * qJD(1)) * t830 + t879;
t862 = g(3) * t869 + t834 * t804;
t849 = -t822 * qJ(3) + t805 * t854 + qJDD(3) + t862;
t742 = t809 * pkin(3) + t878 * t823 + (-pkin(3) * t824 * t861 - pkin(8) * t834 * t871 - t803 * t831) * t838 + t849;
t833 = sin(qJ(4));
t837 = cos(qJ(4));
t736 = t837 * t740 + t833 * t742;
t792 = t837 * t824 - t833 * t855;
t766 = -t792 * qJD(4) - t837 * t810 - t833 * t823;
t791 = -t833 * t824 - t837 * t855;
t770 = -t791 * mrSges(5,1) + t792 * mrSges(5,2);
t815 = qJD(4) + t854;
t775 = t815 * mrSges(5,1) - t792 * mrSges(5,3);
t798 = qJDD(4) + t809;
t771 = -t791 * pkin(4) - t792 * pkin(9);
t813 = t815 ^ 2;
t732 = -t813 * pkin(4) + t798 * pkin(9) + t791 * t771 + t736;
t739 = t810 * pkin(3) - pkin(8) * t857 + t824 * t808 - t747;
t767 = t791 * qJD(4) - t833 * t810 + t837 * t823;
t733 = (-t791 * t815 - t767) * pkin(9) + (t792 * t815 - t766) * pkin(4) + t739;
t832 = sin(qJ(5));
t836 = cos(qJ(5));
t729 = -t832 * t732 + t836 * t733;
t772 = -t832 * t792 + t836 * t815;
t745 = t772 * qJD(5) + t836 * t767 + t832 * t798;
t773 = t836 * t792 + t832 * t815;
t754 = -t772 * mrSges(6,1) + t773 * mrSges(6,2);
t790 = qJD(5) - t791;
t756 = -t790 * mrSges(6,2) + t772 * mrSges(6,3);
t764 = qJDD(5) - t766;
t726 = m(6) * t729 + t764 * mrSges(6,1) - t745 * mrSges(6,3) - t773 * t754 + t790 * t756;
t730 = t836 * t732 + t832 * t733;
t744 = -t773 * qJD(5) - t832 * t767 + t836 * t798;
t757 = t790 * mrSges(6,1) - t773 * mrSges(6,3);
t727 = m(6) * t730 - t764 * mrSges(6,2) + t744 * mrSges(6,3) + t772 * t754 - t790 * t757;
t851 = -t832 * t726 + t836 * t727;
t715 = m(5) * t736 - t798 * mrSges(5,2) + t766 * mrSges(5,3) + t791 * t770 - t815 * t775 + t851;
t735 = -t833 * t740 + t837 * t742;
t774 = -t815 * mrSges(5,2) + t791 * mrSges(5,3);
t731 = -t798 * pkin(4) - t813 * pkin(9) + t792 * t771 - t735;
t847 = -m(6) * t731 + t744 * mrSges(6,1) - t745 * mrSges(6,2) + t772 * t756 - t773 * t757;
t722 = m(5) * t735 + t798 * mrSges(5,1) - t767 * mrSges(5,3) - t792 * t770 + t815 * t774 + t847;
t852 = t837 * t715 - t833 * t722;
t850 = m(4) * t748 - t809 * mrSges(4,3) + t801 * t855 + t852;
t705 = m(3) * t782 + t809 * mrSges(3,2) - t875 * t810 + (-t800 * t838 + (t799 - t802) * t834) * t861 + t850;
t856 = t803 * t867;
t768 = t856 - t862;
t806 = (mrSges(4,2) * t838 - mrSges(4,3) * t834) * t861;
t807 = (-mrSges(3,1) * t838 + mrSges(3,2) * t834) * t861;
t709 = t833 * t715 + t837 * t722;
t753 = -t823 * pkin(2) + t849 - t856;
t848 = -m(4) * t753 - t809 * mrSges(4,1) - t709;
t706 = m(3) * t768 - t809 * mrSges(3,3) + (t800 - t801) * t824 + t875 * t823 + (-t806 - t807) * t854 + t848;
t717 = t836 * t726 + t832 * t727;
t844 = -m(5) * t739 + t766 * mrSges(5,1) - t767 * mrSges(5,2) + t791 * t774 - t792 * t775 - t717;
t842 = -m(4) * t747 + t823 * mrSges(4,3) + t824 * t802 + t806 * t855 - t844;
t713 = (mrSges(3,3) + mrSges(4,1)) * t810 + t842 - t824 * t799 - t823 * mrSges(3,2) + m(3) * t769 + t807 * t855;
t694 = -t830 * t705 + t706 * t867 + t713 * t868;
t691 = m(2) * t819 + qJDD(1) * mrSges(2,1) - t840 * mrSges(2,2) + t694;
t699 = -t834 * t706 + t838 * t713;
t697 = m(2) * t820 - t840 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t699;
t866 = t839 * t691 + t835 * t697;
t865 = (t873 * t834 + t872 * t838) * t861 + t880 * t824;
t864 = (t874 * t834 + t881 * t838) * t861 + t872 * t824;
t863 = (t882 * t834 + t874 * t838) * t861 + t873 * t824;
t693 = t831 * t705 + t706 * t869 + t713 * t870;
t853 = -t835 * t691 + t839 * t697;
t749 = Ifges(6,5) * t773 + Ifges(6,6) * t772 + Ifges(6,3) * t790;
t751 = Ifges(6,1) * t773 + Ifges(6,4) * t772 + Ifges(6,5) * t790;
t720 = -mrSges(6,1) * t731 + mrSges(6,3) * t730 + Ifges(6,4) * t745 + Ifges(6,2) * t744 + Ifges(6,6) * t764 - t773 * t749 + t790 * t751;
t750 = Ifges(6,4) * t773 + Ifges(6,2) * t772 + Ifges(6,6) * t790;
t721 = mrSges(6,2) * t731 - mrSges(6,3) * t729 + Ifges(6,1) * t745 + Ifges(6,4) * t744 + Ifges(6,5) * t764 + t772 * t749 - t790 * t750;
t758 = Ifges(5,5) * t792 + Ifges(5,6) * t791 + Ifges(5,3) * t815;
t759 = Ifges(5,4) * t792 + Ifges(5,2) * t791 + Ifges(5,6) * t815;
t700 = mrSges(5,2) * t739 - mrSges(5,3) * t735 + Ifges(5,1) * t767 + Ifges(5,4) * t766 + Ifges(5,5) * t798 - pkin(9) * t717 - t832 * t720 + t836 * t721 + t791 * t758 - t815 * t759;
t760 = Ifges(5,1) * t792 + Ifges(5,4) * t791 + Ifges(5,5) * t815;
t841 = mrSges(6,1) * t729 - mrSges(6,2) * t730 + Ifges(6,5) * t745 + Ifges(6,6) * t744 + Ifges(6,3) * t764 + t773 * t750 - t772 * t751;
t701 = -mrSges(5,1) * t739 + mrSges(5,3) * t736 + Ifges(5,4) * t767 + Ifges(5,2) * t766 + Ifges(5,6) * t798 - pkin(4) * t717 - t792 * t758 + t815 * t760 - t841;
t708 = t823 * mrSges(4,2) + t824 * t801 + t806 * t854 - t848;
t685 = mrSges(3,1) * t768 - mrSges(3,2) * t769 + mrSges(4,2) * t753 - mrSges(4,3) * t747 + t837 * t700 - t833 * t701 - pkin(8) * t709 - pkin(2) * t708 + qJ(3) * t842 + t880 * t823 + (qJ(3) * mrSges(4,1) + t872) * t810 + t873 * t809 + (t864 * t834 - t863 * t838) * t861;
t707 = t810 * mrSges(4,2) - t802 * t854 + t850;
t687 = -mrSges(3,1) * t782 - mrSges(4,1) * t747 + mrSges(4,2) * t748 + mrSges(3,3) * t769 - pkin(2) * t707 - pkin(3) * t844 - pkin(8) * t852 - t833 * t700 - t837 * t701 + t874 * t809 + t881 * t810 + t872 * t823 + t863 * t824 - t865 * t854;
t843 = mrSges(5,1) * t735 - mrSges(5,2) * t736 + Ifges(5,5) * t767 + Ifges(5,6) * t766 + Ifges(5,3) * t798 + pkin(4) * t847 + pkin(9) * t851 + t836 * t720 + t832 * t721 + t792 * t759 - t791 * t760;
t689 = mrSges(4,1) * t753 + mrSges(3,2) * t782 - mrSges(3,3) * t768 - mrSges(4,3) * t748 + pkin(3) * t709 - qJ(3) * t707 + t882 * t809 + t874 * t810 + t873 * t823 - t864 * t824 + t865 * t855 + t843;
t846 = mrSges(2,1) * t819 - mrSges(2,2) * t820 + Ifges(2,3) * qJDD(1) + pkin(1) * t694 + t831 * t685 + t687 * t869 + t689 * t870 + t699 * t877;
t683 = -mrSges(2,2) * g(3) - mrSges(2,3) * t819 + Ifges(2,5) * qJDD(1) - t840 * Ifges(2,6) - t834 * t687 + t838 * t689 + (-t693 * t830 - t694 * t831) * pkin(7);
t682 = mrSges(2,1) * g(3) + mrSges(2,3) * t820 + t840 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t693 - t830 * t685 + (pkin(7) * t699 + t687 * t838 + t689 * t834) * t831;
t1 = [-m(1) * g(1) + t853; -m(1) * g(2) + t866; (-m(1) - m(2)) * g(3) + t693; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t866 - t835 * t682 + t839 * t683; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t853 + t839 * t682 + t835 * t683; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t846; t846; t685; t708; t843; t841;];
tauJB = t1;

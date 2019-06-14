% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:29:32
% EndTime: 2019-05-05 04:29:54
% DurationCPUTime: 21.21s
% Computational Cost: add. (355004->342), mult. (768302->442), div. (0->0), fcn. (559393->14), ass. (0->145)
t874 = sin(pkin(11));
t877 = cos(pkin(11));
t863 = g(1) * t874 - g(2) * t877;
t864 = -g(1) * t877 - g(2) * t874;
t872 = -g(3) + qJDD(1);
t886 = cos(qJ(2));
t878 = cos(pkin(6));
t882 = sin(qJ(2));
t910 = t878 * t882;
t875 = sin(pkin(6));
t911 = t875 * t882;
t826 = t863 * t910 + t886 * t864 + t872 * t911;
t888 = qJD(2) ^ 2;
t821 = -pkin(2) * t888 + qJDD(2) * pkin(8) + t826;
t843 = -t863 * t875 + t872 * t878;
t881 = sin(qJ(3));
t885 = cos(qJ(3));
t803 = -t881 * t821 + t885 * t843;
t905 = qJD(2) * qJD(3);
t904 = t885 * t905;
t861 = qJDD(2) * t881 + t904;
t798 = (-t861 + t904) * qJ(4) + (t881 * t885 * t888 + qJDD(3)) * pkin(3) + t803;
t804 = t885 * t821 + t881 * t843;
t862 = qJDD(2) * t885 - t881 * t905;
t908 = qJD(2) * t881;
t865 = qJD(3) * pkin(3) - qJ(4) * t908;
t871 = t885 ^ 2;
t799 = -pkin(3) * t871 * t888 + qJ(4) * t862 - qJD(3) * t865 + t804;
t873 = sin(pkin(12));
t876 = cos(pkin(12));
t849 = (t873 * t885 + t876 * t881) * qJD(2);
t779 = -0.2e1 * qJD(4) * t849 + t798 * t876 - t873 * t799;
t825 = -t882 * t864 + (t863 * t878 + t872 * t875) * t886;
t907 = qJD(2) * t885;
t848 = -t873 * t908 + t876 * t907;
t780 = 0.2e1 * qJD(4) * t848 + t873 * t798 + t876 * t799;
t830 = -pkin(4) * t848 - pkin(9) * t849;
t887 = qJD(3) ^ 2;
t778 = -pkin(4) * t887 + qJDD(3) * pkin(9) + t830 * t848 + t780;
t892 = -qJDD(2) * pkin(2) - t825;
t800 = -t862 * pkin(3) + qJDD(4) + t865 * t908 + (-qJ(4) * t871 - pkin(8)) * t888 + t892;
t834 = -t873 * t861 + t862 * t876;
t835 = t861 * t876 + t862 * t873;
t789 = (-qJD(3) * t848 - t835) * pkin(9) + (qJD(3) * t849 - t834) * pkin(4) + t800;
t880 = sin(qJ(5));
t884 = cos(qJ(5));
t773 = -t880 * t778 + t884 * t789;
t837 = qJD(3) * t884 - t849 * t880;
t811 = qJD(5) * t837 + qJDD(3) * t880 + t835 * t884;
t833 = qJDD(5) - t834;
t838 = qJD(3) * t880 + t849 * t884;
t847 = qJD(5) - t848;
t771 = (t837 * t847 - t811) * pkin(10) + (t837 * t838 + t833) * pkin(5) + t773;
t774 = t884 * t778 + t880 * t789;
t810 = -qJD(5) * t838 + qJDD(3) * t884 - t835 * t880;
t819 = pkin(5) * t847 - pkin(10) * t838;
t836 = t837 ^ 2;
t772 = -pkin(5) * t836 + pkin(10) * t810 - t819 * t847 + t774;
t879 = sin(qJ(6));
t883 = cos(qJ(6));
t769 = t771 * t883 - t772 * t879;
t812 = t837 * t883 - t838 * t879;
t785 = qJD(6) * t812 + t810 * t879 + t811 * t883;
t813 = t837 * t879 + t838 * t883;
t794 = -mrSges(7,1) * t812 + mrSges(7,2) * t813;
t844 = qJD(6) + t847;
t801 = -mrSges(7,2) * t844 + mrSges(7,3) * t812;
t831 = qJDD(6) + t833;
t765 = m(7) * t769 + mrSges(7,1) * t831 - t785 * mrSges(7,3) - t794 * t813 + t801 * t844;
t770 = t771 * t879 + t772 * t883;
t784 = -qJD(6) * t813 + t810 * t883 - t811 * t879;
t802 = mrSges(7,1) * t844 - mrSges(7,3) * t813;
t766 = m(7) * t770 - mrSges(7,2) * t831 + t784 * mrSges(7,3) + t794 * t812 - t802 * t844;
t757 = t883 * t765 + t879 * t766;
t814 = -mrSges(6,1) * t837 + mrSges(6,2) * t838;
t817 = -mrSges(6,2) * t847 + mrSges(6,3) * t837;
t755 = m(6) * t773 + mrSges(6,1) * t833 - mrSges(6,3) * t811 - t814 * t838 + t817 * t847 + t757;
t818 = mrSges(6,1) * t847 - mrSges(6,3) * t838;
t900 = -t765 * t879 + t883 * t766;
t756 = m(6) * t774 - mrSges(6,2) * t833 + mrSges(6,3) * t810 + t814 * t837 - t818 * t847 + t900;
t751 = -t755 * t880 + t884 * t756;
t828 = -mrSges(5,1) * t848 + mrSges(5,2) * t849;
t842 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t849;
t748 = m(5) * t780 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t834 - qJD(3) * t842 + t828 * t848 + t751;
t777 = -qJDD(3) * pkin(4) - pkin(9) * t887 + t849 * t830 - t779;
t775 = -pkin(5) * t810 - pkin(10) * t836 + t819 * t838 + t777;
t894 = m(7) * t775 - t784 * mrSges(7,1) + t785 * mrSges(7,2) - t812 * t801 + t802 * t813;
t767 = -m(6) * t777 + t810 * mrSges(6,1) - mrSges(6,2) * t811 + t837 * t817 - t818 * t838 - t894;
t841 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t848;
t761 = m(5) * t779 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t835 + qJD(3) * t841 - t828 * t849 + t767;
t740 = t873 * t748 + t876 * t761;
t790 = Ifges(7,5) * t813 + Ifges(7,6) * t812 + Ifges(7,3) * t844;
t792 = Ifges(7,1) * t813 + Ifges(7,4) * t812 + Ifges(7,5) * t844;
t758 = -mrSges(7,1) * t775 + mrSges(7,3) * t770 + Ifges(7,4) * t785 + Ifges(7,2) * t784 + Ifges(7,6) * t831 - t790 * t813 + t792 * t844;
t791 = Ifges(7,4) * t813 + Ifges(7,2) * t812 + Ifges(7,6) * t844;
t759 = mrSges(7,2) * t775 - mrSges(7,3) * t769 + Ifges(7,1) * t785 + Ifges(7,4) * t784 + Ifges(7,5) * t831 + t790 * t812 - t791 * t844;
t805 = Ifges(6,5) * t838 + Ifges(6,6) * t837 + Ifges(6,3) * t847;
t807 = Ifges(6,1) * t838 + Ifges(6,4) * t837 + Ifges(6,5) * t847;
t741 = -mrSges(6,1) * t777 + mrSges(6,3) * t774 + Ifges(6,4) * t811 + Ifges(6,2) * t810 + Ifges(6,6) * t833 - pkin(5) * t894 + pkin(10) * t900 + t883 * t758 + t879 * t759 - t838 * t805 + t847 * t807;
t806 = Ifges(6,4) * t838 + Ifges(6,2) * t837 + Ifges(6,6) * t847;
t742 = mrSges(6,2) * t777 - mrSges(6,3) * t773 + Ifges(6,1) * t811 + Ifges(6,4) * t810 + Ifges(6,5) * t833 - pkin(10) * t757 - t758 * t879 + t759 * t883 + t805 * t837 - t806 * t847;
t823 = Ifges(5,4) * t849 + Ifges(5,2) * t848 + Ifges(5,6) * qJD(3);
t824 = Ifges(5,1) * t849 + Ifges(5,4) * t848 + Ifges(5,5) * qJD(3);
t852 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t881 + Ifges(4,2) * t885) * qJD(2);
t853 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t881 + Ifges(4,4) * t885) * qJD(2);
t914 = (t852 * t881 - t853 * t885) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t803 + mrSges(5,1) * t779 - mrSges(4,2) * t804 - mrSges(5,2) * t780 + Ifges(4,5) * t861 + Ifges(5,5) * t835 + Ifges(4,6) * t862 + Ifges(5,6) * t834 + pkin(3) * t740 + pkin(4) * t767 + pkin(9) * t751 + t884 * t741 + t880 * t742 + t849 * t823 - t848 * t824;
t750 = t884 * t755 + t880 * t756;
t749 = m(5) * t800 - t834 * mrSges(5,1) + mrSges(5,2) * t835 - t848 * t841 + t842 * t849 + t750;
t820 = -t888 * pkin(8) + t892;
t866 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t908;
t867 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t907;
t891 = -m(4) * t820 + t862 * mrSges(4,1) - mrSges(4,2) * t861 - t866 * t908 + t867 * t907 - t749;
t745 = m(3) * t825 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t888 + t891;
t912 = t745 * t886;
t860 = (-mrSges(4,1) * t885 + mrSges(4,2) * t881) * qJD(2);
t738 = m(4) * t803 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t861 + qJD(3) * t867 - t860 * t908 + t740;
t901 = t876 * t748 - t761 * t873;
t739 = m(4) * t804 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t862 - qJD(3) * t866 + t860 * t907 + t901;
t902 = -t738 * t881 + t885 * t739;
t730 = m(3) * t826 - mrSges(3,1) * t888 - qJDD(2) * mrSges(3,2) + t902;
t733 = t885 * t738 + t881 * t739;
t732 = m(3) * t843 + t733;
t721 = t730 * t910 - t732 * t875 + t878 * t912;
t719 = m(2) * t863 + t721;
t725 = t886 * t730 - t745 * t882;
t724 = m(2) * t864 + t725;
t909 = t877 * t719 + t874 * t724;
t720 = t730 * t911 + t878 * t732 + t875 * t912;
t903 = -t719 * t874 + t877 * t724;
t898 = m(2) * t872 + t720;
t822 = Ifges(5,5) * t849 + Ifges(5,6) * t848 + Ifges(5,3) * qJD(3);
t726 = mrSges(5,2) * t800 - mrSges(5,3) * t779 + Ifges(5,1) * t835 + Ifges(5,4) * t834 + Ifges(5,5) * qJDD(3) - pkin(9) * t750 - qJD(3) * t823 - t741 * t880 + t742 * t884 + t822 * t848;
t893 = -mrSges(7,1) * t769 + mrSges(7,2) * t770 - Ifges(7,5) * t785 - Ifges(7,6) * t784 - Ifges(7,3) * t831 - t813 * t791 + t812 * t792;
t890 = mrSges(6,1) * t773 - mrSges(6,2) * t774 + Ifges(6,5) * t811 + Ifges(6,6) * t810 + Ifges(6,3) * t833 + pkin(5) * t757 + t838 * t806 - t837 * t807 - t893;
t734 = -mrSges(5,1) * t800 + mrSges(5,3) * t780 + Ifges(5,4) * t835 + Ifges(5,2) * t834 + Ifges(5,6) * qJDD(3) - pkin(4) * t750 + qJD(3) * t824 - t849 * t822 - t890;
t851 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t881 + Ifges(4,6) * t885) * qJD(2);
t715 = -mrSges(4,1) * t820 + mrSges(4,3) * t804 + Ifges(4,4) * t861 + Ifges(4,2) * t862 + Ifges(4,6) * qJDD(3) - pkin(3) * t749 + qJ(4) * t901 + qJD(3) * t853 + t873 * t726 + t876 * t734 - t851 * t908;
t717 = mrSges(4,2) * t820 - mrSges(4,3) * t803 + Ifges(4,1) * t861 + Ifges(4,4) * t862 + Ifges(4,5) * qJDD(3) - qJ(4) * t740 - qJD(3) * t852 + t726 * t876 - t734 * t873 + t851 * t907;
t714 = mrSges(3,2) * t843 - mrSges(3,3) * t825 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t888 - pkin(8) * t733 - t715 * t881 + t717 * t885;
t716 = -mrSges(3,1) * t843 + mrSges(3,3) * t826 + t888 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t733 - t914;
t895 = pkin(7) * t725 + t714 * t882 + t716 * t886;
t713 = mrSges(3,1) * t825 - mrSges(3,2) * t826 + Ifges(3,3) * qJDD(2) + pkin(2) * t891 + pkin(8) * t902 + t885 * t715 + t881 * t717;
t712 = mrSges(2,2) * t872 - mrSges(2,3) * t863 + t886 * t714 - t882 * t716 + (-t720 * t875 - t721 * t878) * pkin(7);
t711 = -mrSges(2,1) * t872 + mrSges(2,3) * t864 - pkin(1) * t720 - t875 * t713 + t878 * t895;
t1 = [-m(1) * g(1) + t903; -m(1) * g(2) + t909; -m(1) * g(3) + t898; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t909 - t874 * t711 + t877 * t712; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t903 + t877 * t711 + t874 * t712; -mrSges(1,1) * g(2) + mrSges(2,1) * t863 + mrSges(1,2) * g(1) - mrSges(2,2) * t864 + pkin(1) * t721 + t878 * t713 + t875 * t895; t898; t713; t914; t749; t890; -t893;];
tauJB  = t1;

% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:14:03
% EndTime: 2019-05-05 07:14:26
% DurationCPUTime: 22.42s
% Computational Cost: add. (398284->342), mult. (802046->442), div. (0->0), fcn. (585342->14), ass. (0->145)
t873 = sin(pkin(11));
t876 = cos(pkin(11));
t858 = g(1) * t873 - g(2) * t876;
t859 = -g(1) * t876 - g(2) * t873;
t871 = -g(3) + qJDD(1);
t874 = sin(pkin(6));
t877 = cos(pkin(6));
t881 = sin(qJ(2));
t884 = cos(qJ(2));
t827 = -t881 * t859 + (t858 * t877 + t871 * t874) * t884;
t908 = t877 * t881;
t909 = t874 * t881;
t828 = t858 * t908 + t884 * t859 + t871 * t909;
t885 = qJD(2) ^ 2;
t823 = -pkin(2) * t885 + qJDD(2) * pkin(8) + t828;
t841 = -t858 * t874 + t871 * t877;
t880 = sin(qJ(3));
t883 = cos(qJ(3));
t803 = -t880 * t823 + t883 * t841;
t904 = qJD(2) * qJD(3);
t903 = t883 * t904;
t856 = qJDD(2) * t880 + t903;
t794 = (-t856 + t903) * pkin(9) + (t880 * t883 * t885 + qJDD(3)) * pkin(3) + t803;
t804 = t883 * t823 + t880 * t841;
t857 = qJDD(2) * t883 - t880 * t904;
t906 = qJD(2) * t880;
t863 = qJD(3) * pkin(3) - pkin(9) * t906;
t870 = t883 ^ 2;
t795 = -pkin(3) * t870 * t885 + pkin(9) * t857 - qJD(3) * t863 + t804;
t879 = sin(qJ(4));
t911 = cos(qJ(4));
t782 = t879 * t794 + t911 * t795;
t848 = (t879 * t883 + t911 * t880) * qJD(2);
t818 = qJD(4) * t848 + t856 * t879 - t911 * t857;
t905 = qJD(2) * t883;
t847 = t879 * t906 - t911 * t905;
t831 = mrSges(5,1) * t847 + mrSges(5,2) * t848;
t869 = qJD(3) + qJD(4);
t840 = mrSges(5,1) * t869 - mrSges(5,3) * t848;
t868 = qJDD(3) + qJDD(4);
t830 = pkin(4) * t847 - qJ(5) * t848;
t867 = t869 ^ 2;
t776 = -pkin(4) * t867 + qJ(5) * t868 - t830 * t847 + t782;
t891 = -qJDD(2) * pkin(2) - t827;
t798 = -t857 * pkin(3) + t863 * t906 + (-pkin(9) * t870 - pkin(8)) * t885 + t891;
t819 = -t847 * qJD(4) + t911 * t856 + t879 * t857;
t779 = (t847 * t869 - t819) * qJ(5) + (t848 * t869 + t818) * pkin(4) + t798;
t872 = sin(pkin(12));
t875 = cos(pkin(12));
t836 = t848 * t875 + t869 * t872;
t771 = -0.2e1 * qJD(5) * t836 - t872 * t776 + t875 * t779;
t809 = t819 * t875 + t868 * t872;
t835 = -t848 * t872 + t869 * t875;
t769 = (t835 * t847 - t809) * pkin(10) + (t835 * t836 + t818) * pkin(5) + t771;
t772 = 0.2e1 * qJD(5) * t835 + t875 * t776 + t872 * t779;
t808 = -t819 * t872 + t868 * t875;
t821 = pkin(5) * t847 - pkin(10) * t836;
t834 = t835 ^ 2;
t770 = -pkin(5) * t834 + pkin(10) * t808 - t821 * t847 + t772;
t878 = sin(qJ(6));
t882 = cos(qJ(6));
t767 = t769 * t882 - t770 * t878;
t806 = t835 * t882 - t836 * t878;
t785 = qJD(6) * t806 + t808 * t878 + t809 * t882;
t807 = t835 * t878 + t836 * t882;
t790 = -mrSges(7,1) * t806 + mrSges(7,2) * t807;
t842 = qJD(6) + t847;
t796 = -mrSges(7,2) * t842 + mrSges(7,3) * t806;
t816 = qJDD(6) + t818;
t763 = m(7) * t767 + mrSges(7,1) * t816 - mrSges(7,3) * t785 - t790 * t807 + t796 * t842;
t768 = t769 * t878 + t770 * t882;
t784 = -qJD(6) * t807 + t808 * t882 - t809 * t878;
t797 = mrSges(7,1) * t842 - mrSges(7,3) * t807;
t764 = m(7) * t768 - mrSges(7,2) * t816 + mrSges(7,3) * t784 + t790 * t806 - t797 * t842;
t755 = t882 * t763 + t878 * t764;
t810 = -mrSges(6,1) * t835 + mrSges(6,2) * t836;
t896 = -mrSges(6,2) * t847 + mrSges(6,3) * t835;
t753 = m(6) * t771 + t818 * mrSges(6,1) - t809 * mrSges(6,3) - t836 * t810 + t847 * t896 + t755;
t820 = mrSges(6,1) * t847 - mrSges(6,3) * t836;
t898 = -t763 * t878 + t882 * t764;
t754 = m(6) * t772 - mrSges(6,2) * t818 + mrSges(6,3) * t808 + t810 * t835 - t820 * t847 + t898;
t899 = -t753 * t872 + t875 * t754;
t746 = m(5) * t782 - mrSges(5,2) * t868 - mrSges(5,3) * t818 - t831 * t847 - t840 * t869 + t899;
t781 = t911 * t794 - t879 * t795;
t775 = -t868 * pkin(4) - t867 * qJ(5) + t848 * t830 + qJDD(5) - t781;
t773 = -t808 * pkin(5) - t834 * pkin(10) + t836 * t821 + t775;
t892 = m(7) * t773 - t784 * mrSges(7,1) + mrSges(7,2) * t785 - t806 * t796 + t797 * t807;
t766 = m(6) * t775 - t808 * mrSges(6,1) + mrSges(6,2) * t809 + t820 * t836 - t835 * t896 + t892;
t839 = -mrSges(5,2) * t869 - mrSges(5,3) * t847;
t759 = m(5) * t781 + mrSges(5,1) * t868 - mrSges(5,3) * t819 - t831 * t848 + t839 * t869 - t766;
t736 = t879 * t746 + t911 * t759;
t845 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t880 + Ifges(4,2) * t883) * qJD(2);
t846 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t880 + Ifges(4,4) * t883) * qJD(2);
t786 = Ifges(7,5) * t807 + Ifges(7,6) * t806 + Ifges(7,3) * t842;
t788 = Ifges(7,1) * t807 + Ifges(7,4) * t806 + Ifges(7,5) * t842;
t756 = -mrSges(7,1) * t773 + mrSges(7,3) * t768 + Ifges(7,4) * t785 + Ifges(7,2) * t784 + Ifges(7,6) * t816 - t786 * t807 + t788 * t842;
t787 = Ifges(7,4) * t807 + Ifges(7,2) * t806 + Ifges(7,6) * t842;
t757 = mrSges(7,2) * t773 - mrSges(7,3) * t767 + Ifges(7,1) * t785 + Ifges(7,4) * t784 + Ifges(7,5) * t816 + t786 * t806 - t787 * t842;
t799 = Ifges(6,5) * t836 + Ifges(6,6) * t835 + Ifges(6,3) * t847;
t801 = Ifges(6,1) * t836 + Ifges(6,4) * t835 + Ifges(6,5) * t847;
t738 = -mrSges(6,1) * t775 + mrSges(6,3) * t772 + Ifges(6,4) * t809 + Ifges(6,2) * t808 + Ifges(6,6) * t818 - pkin(5) * t892 + pkin(10) * t898 + t882 * t756 + t878 * t757 - t836 * t799 + t847 * t801;
t800 = Ifges(6,4) * t836 + Ifges(6,2) * t835 + Ifges(6,6) * t847;
t740 = mrSges(6,2) * t775 - mrSges(6,3) * t771 + Ifges(6,1) * t809 + Ifges(6,4) * t808 + Ifges(6,5) * t818 - pkin(10) * t755 - t756 * t878 + t757 * t882 + t799 * t835 - t800 * t847;
t825 = Ifges(5,4) * t848 - Ifges(5,2) * t847 + Ifges(5,6) * t869;
t826 = Ifges(5,1) * t848 - Ifges(5,4) * t847 + Ifges(5,5) * t869;
t889 = -mrSges(5,1) * t781 + mrSges(5,2) * t782 - Ifges(5,5) * t819 + Ifges(5,6) * t818 - Ifges(5,3) * t868 + pkin(4) * t766 - qJ(5) * t899 - t875 * t738 - t872 * t740 - t848 * t825 - t847 * t826;
t912 = mrSges(4,1) * t803 - mrSges(4,2) * t804 + Ifges(4,5) * t856 + Ifges(4,6) * t857 + Ifges(4,3) * qJDD(3) + pkin(3) * t736 + (t845 * t880 - t846 * t883) * qJD(2) - t889;
t822 = -t885 * pkin(8) + t891;
t860 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t906;
t861 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t905;
t748 = t875 * t753 + t872 * t754;
t890 = m(5) * t798 + t818 * mrSges(5,1) + mrSges(5,2) * t819 + t847 * t839 + t840 * t848 + t748;
t887 = -m(4) * t822 + t857 * mrSges(4,1) - mrSges(4,2) * t856 - t860 * t906 + t861 * t905 - t890;
t743 = m(3) * t827 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t885 + t887;
t910 = t743 * t884;
t855 = (-mrSges(4,1) * t883 + mrSges(4,2) * t880) * qJD(2);
t734 = m(4) * t803 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t856 + qJD(3) * t861 - t855 * t906 + t736;
t900 = t911 * t746 - t759 * t879;
t735 = m(4) * t804 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t857 - qJD(3) * t860 + t855 * t905 + t900;
t901 = -t734 * t880 + t883 * t735;
t726 = m(3) * t828 - mrSges(3,1) * t885 - qJDD(2) * mrSges(3,2) + t901;
t729 = t883 * t734 + t880 * t735;
t728 = m(3) * t841 + t729;
t717 = t726 * t908 - t728 * t874 + t877 * t910;
t715 = m(2) * t858 + t717;
t721 = t884 * t726 - t743 * t881;
t720 = m(2) * t859 + t721;
t907 = t876 * t715 + t873 * t720;
t716 = t726 * t909 + t877 * t728 + t874 * t910;
t902 = -t715 * t873 + t876 * t720;
t897 = m(2) * t871 + t716;
t824 = Ifges(5,5) * t848 - Ifges(5,6) * t847 + Ifges(5,3) * t869;
t722 = mrSges(5,2) * t798 - mrSges(5,3) * t781 + Ifges(5,1) * t819 - Ifges(5,4) * t818 + Ifges(5,5) * t868 - qJ(5) * t748 - t738 * t872 + t740 * t875 - t824 * t847 - t825 * t869;
t888 = mrSges(7,1) * t767 - mrSges(7,2) * t768 + Ifges(7,5) * t785 + Ifges(7,6) * t784 + Ifges(7,3) * t816 + t807 * t787 - t806 * t788;
t730 = -t888 + (-Ifges(5,2) - Ifges(6,3)) * t818 + Ifges(5,6) * t868 + t869 * t826 - t848 * t824 + t835 * t801 - t836 * t800 + Ifges(5,4) * t819 - Ifges(6,6) * t808 - Ifges(6,5) * t809 - mrSges(5,1) * t798 + mrSges(5,3) * t782 + mrSges(6,2) * t772 - mrSges(6,1) * t771 - pkin(5) * t755 - pkin(4) * t748;
t844 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t880 + Ifges(4,6) * t883) * qJD(2);
t711 = -mrSges(4,1) * t822 + mrSges(4,3) * t804 + Ifges(4,4) * t856 + Ifges(4,2) * t857 + Ifges(4,6) * qJDD(3) - pkin(3) * t890 + pkin(9) * t900 + qJD(3) * t846 + t879 * t722 + t911 * t730 - t844 * t906;
t713 = mrSges(4,2) * t822 - mrSges(4,3) * t803 + Ifges(4,1) * t856 + Ifges(4,4) * t857 + Ifges(4,5) * qJDD(3) - pkin(9) * t736 - qJD(3) * t845 + t911 * t722 - t879 * t730 + t844 * t905;
t710 = mrSges(3,2) * t841 - mrSges(3,3) * t827 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t885 - pkin(8) * t729 - t711 * t880 + t713 * t883;
t712 = -mrSges(3,1) * t841 + mrSges(3,3) * t828 + t885 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t729 - t912;
t893 = pkin(7) * t721 + t710 * t881 + t712 * t884;
t709 = mrSges(3,1) * t827 - mrSges(3,2) * t828 + Ifges(3,3) * qJDD(2) + pkin(2) * t887 + pkin(8) * t901 + t883 * t711 + t880 * t713;
t708 = mrSges(2,2) * t871 - mrSges(2,3) * t858 + t884 * t710 - t881 * t712 + (-t716 * t874 - t717 * t877) * pkin(7);
t707 = -mrSges(2,1) * t871 + mrSges(2,3) * t859 - pkin(1) * t716 - t874 * t709 + t893 * t877;
t1 = [-m(1) * g(1) + t902; -m(1) * g(2) + t907; -m(1) * g(3) + t897; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t907 - t873 * t707 + t876 * t708; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t902 + t876 * t707 + t873 * t708; -mrSges(1,1) * g(2) + mrSges(2,1) * t858 + mrSges(1,2) * g(1) - mrSges(2,2) * t859 + pkin(1) * t717 + t877 * t709 + t893 * t874; t897; t709; t912; -t889; t766; t888;];
tauJB  = t1;

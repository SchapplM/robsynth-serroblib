% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:14:24
% EndTime: 2019-05-05 16:14:33
% DurationCPUTime: 8.28s
% Computational Cost: add. (121320->316), mult. (275086->384), div. (0->0), fcn. (195942->10), ass. (0->139)
t870 = sin(qJ(1));
t874 = cos(qJ(1));
t842 = t870 * g(1) - t874 * g(2);
t876 = qJD(1) ^ 2;
t887 = -t876 * qJ(2) + qJDD(2) - t842;
t909 = -pkin(1) - qJ(3);
t917 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t909 + t887;
t843 = -t874 * g(1) - t870 * g(2);
t916 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t843;
t831 = t876 * pkin(1) - t916;
t915 = -m(3) * t831 + t876 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t913 = pkin(3) * t876;
t912 = mrSges(2,1) - mrSges(3,2);
t911 = -Ifges(3,4) + Ifges(2,5);
t910 = -Ifges(2,6) + Ifges(3,5);
t865 = sin(pkin(10));
t866 = cos(pkin(10));
t816 = t865 * g(3) + t866 * t917;
t798 = (-pkin(7) * qJDD(1) - t865 * t913) * t866 + t816;
t817 = -g(3) * t866 + t865 * t917;
t856 = t865 ^ 2;
t902 = qJDD(1) * t865;
t799 = -pkin(7) * t902 - t856 * t913 + t817;
t869 = sin(qJ(4));
t873 = cos(qJ(4));
t780 = t869 * t798 + t873 * t799;
t905 = qJD(1) * t866;
t906 = qJD(1) * t865;
t837 = -t869 * t905 - t873 * t906;
t890 = -t865 * t869 + t866 * t873;
t838 = t890 * qJD(1);
t811 = -mrSges(5,1) * t837 + mrSges(5,2) * t838;
t834 = t838 * qJD(4);
t901 = qJDD(1) * t866;
t819 = -t869 * t901 - t873 * t902 - t834;
t829 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t838;
t818 = -pkin(4) * t837 - pkin(8) * t838;
t875 = qJD(4) ^ 2;
t769 = -pkin(4) * t875 + qJDD(4) * pkin(8) + t818 * t837 + t780;
t886 = qJDD(3) + t916;
t907 = -t866 ^ 2 - t856;
t804 = pkin(3) * t902 + (pkin(7) * t907 + t909) * t876 + t886;
t904 = t837 * qJD(4);
t820 = qJDD(1) * t890 + t904;
t772 = (-t820 - t904) * pkin(8) + (-t819 + t834) * pkin(4) + t804;
t868 = sin(qJ(5));
t872 = cos(qJ(5));
t759 = -t868 * t769 + t872 * t772;
t823 = qJD(4) * t872 - t838 * t868;
t789 = qJD(5) * t823 + qJDD(4) * t868 + t820 * t872;
t815 = qJDD(5) - t819;
t824 = qJD(4) * t868 + t838 * t872;
t835 = qJD(5) - t837;
t756 = (t823 * t835 - t789) * pkin(9) + (t823 * t824 + t815) * pkin(5) + t759;
t760 = t872 * t769 + t868 * t772;
t788 = -qJD(5) * t824 + qJDD(4) * t872 - t820 * t868;
t802 = pkin(5) * t835 - pkin(9) * t824;
t822 = t823 ^ 2;
t757 = -pkin(5) * t822 + pkin(9) * t788 - t802 * t835 + t760;
t867 = sin(qJ(6));
t871 = cos(qJ(6));
t754 = t756 * t871 - t757 * t867;
t790 = t823 * t871 - t824 * t867;
t765 = qJD(6) * t790 + t788 * t867 + t789 * t871;
t791 = t823 * t867 + t824 * t871;
t777 = -mrSges(7,1) * t790 + mrSges(7,2) * t791;
t832 = qJD(6) + t835;
t781 = -mrSges(7,2) * t832 + mrSges(7,3) * t790;
t810 = qJDD(6) + t815;
t750 = m(7) * t754 + mrSges(7,1) * t810 - t765 * mrSges(7,3) - t777 * t791 + t781 * t832;
t755 = t756 * t867 + t757 * t871;
t764 = -qJD(6) * t791 + t788 * t871 - t789 * t867;
t782 = mrSges(7,1) * t832 - mrSges(7,3) * t791;
t751 = m(7) * t755 - mrSges(7,2) * t810 + t764 * mrSges(7,3) + t777 * t790 - t782 * t832;
t742 = t871 * t750 + t867 * t751;
t794 = -mrSges(6,1) * t823 + mrSges(6,2) * t824;
t800 = -mrSges(6,2) * t835 + mrSges(6,3) * t823;
t740 = m(6) * t759 + mrSges(6,1) * t815 - mrSges(6,3) * t789 - t794 * t824 + t800 * t835 + t742;
t801 = mrSges(6,1) * t835 - mrSges(6,3) * t824;
t894 = -t750 * t867 + t871 * t751;
t741 = m(6) * t760 - mrSges(6,2) * t815 + mrSges(6,3) * t788 + t794 * t823 - t801 * t835 + t894;
t895 = -t740 * t868 + t872 * t741;
t734 = m(5) * t780 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t819 - qJD(4) * t829 + t811 * t837 + t895;
t779 = t798 * t873 - t869 * t799;
t828 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t837;
t768 = -qJDD(4) * pkin(4) - pkin(8) * t875 + t838 * t818 - t779;
t758 = -pkin(5) * t788 - pkin(9) * t822 + t802 * t824 + t768;
t884 = m(7) * t758 - t764 * mrSges(7,1) + t765 * mrSges(7,2) - t790 * t781 + t782 * t791;
t878 = -m(6) * t768 + t788 * mrSges(6,1) - mrSges(6,2) * t789 + t823 * t800 - t801 * t824 - t884;
t746 = m(5) * t779 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t820 + qJD(4) * t828 - t811 * t838 + t878;
t724 = t869 * t734 + t873 * t746;
t889 = -qJDD(1) * mrSges(4,3) - t876 * (mrSges(4,1) * t865 + mrSges(4,2) * t866);
t720 = m(4) * t816 + t866 * t889 + t724;
t896 = t873 * t734 - t869 * t746;
t721 = m(4) * t817 + t865 * t889 + t896;
t717 = t866 * t720 + t865 * t721;
t836 = -qJDD(1) * pkin(1) + t887;
t885 = -m(3) * t836 + t876 * mrSges(3,3) - t717;
t713 = m(2) * t842 - t876 * mrSges(2,2) + qJDD(1) * t912 + t885;
t827 = t876 * t909 + t886;
t736 = t872 * t740 + t868 * t741;
t883 = m(5) * t804 - t819 * mrSges(5,1) + t820 * mrSges(5,2) - t837 * t828 + t838 * t829 + t736;
t881 = m(4) * t827 + mrSges(4,1) * t902 + mrSges(4,2) * t901 + t883;
t899 = t907 * mrSges(4,3);
t729 = t881 + (-mrSges(2,1) + t899) * t876 - qJDD(1) * mrSges(2,2) + m(2) * t843 + t915;
t908 = t874 * t713 + t870 * t729;
t898 = -t713 * t870 + t874 * t729;
t897 = -t865 * t720 + t866 * t721;
t893 = Ifges(4,1) * t866 - Ifges(4,4) * t865;
t892 = Ifges(4,4) * t866 - Ifges(4,2) * t865;
t891 = Ifges(4,5) * t866 - Ifges(4,6) * t865;
t774 = Ifges(7,4) * t791 + Ifges(7,2) * t790 + Ifges(7,6) * t832;
t775 = Ifges(7,1) * t791 + Ifges(7,4) * t790 + Ifges(7,5) * t832;
t882 = -mrSges(7,1) * t754 + mrSges(7,2) * t755 - Ifges(7,5) * t765 - Ifges(7,6) * t764 - Ifges(7,3) * t810 - t791 * t774 + t790 * t775;
t773 = Ifges(7,5) * t791 + Ifges(7,6) * t790 + Ifges(7,3) * t832;
t743 = -mrSges(7,1) * t758 + mrSges(7,3) * t755 + Ifges(7,4) * t765 + Ifges(7,2) * t764 + Ifges(7,6) * t810 - t773 * t791 + t775 * t832;
t744 = mrSges(7,2) * t758 - mrSges(7,3) * t754 + Ifges(7,1) * t765 + Ifges(7,4) * t764 + Ifges(7,5) * t810 + t773 * t790 - t774 * t832;
t783 = Ifges(6,5) * t824 + Ifges(6,6) * t823 + Ifges(6,3) * t835;
t785 = Ifges(6,1) * t824 + Ifges(6,4) * t823 + Ifges(6,5) * t835;
t723 = -mrSges(6,1) * t768 + mrSges(6,3) * t760 + Ifges(6,4) * t789 + Ifges(6,2) * t788 + Ifges(6,6) * t815 - pkin(5) * t884 + pkin(9) * t894 + t871 * t743 + t867 * t744 - t824 * t783 + t835 * t785;
t784 = Ifges(6,4) * t824 + Ifges(6,2) * t823 + Ifges(6,6) * t835;
t726 = mrSges(6,2) * t768 - mrSges(6,3) * t759 + Ifges(6,1) * t789 + Ifges(6,4) * t788 + Ifges(6,5) * t815 - pkin(9) * t742 - t743 * t867 + t744 * t871 + t783 * t823 - t784 * t835;
t806 = Ifges(5,4) * t838 + Ifges(5,2) * t837 + Ifges(5,6) * qJD(4);
t807 = Ifges(5,1) * t838 + Ifges(5,4) * t837 + Ifges(5,5) * qJD(4);
t880 = mrSges(5,1) * t779 - mrSges(5,2) * t780 + Ifges(5,5) * t820 + Ifges(5,6) * t819 + Ifges(5,3) * qJDD(4) + pkin(4) * t878 + pkin(8) * t895 + t872 * t723 + t868 * t726 + t838 * t806 - t837 * t807;
t805 = Ifges(5,5) * t838 + Ifges(5,6) * t837 + Ifges(5,3) * qJD(4);
t711 = mrSges(5,2) * t804 - mrSges(5,3) * t779 + Ifges(5,1) * t820 + Ifges(5,4) * t819 + Ifges(5,5) * qJDD(4) - pkin(8) * t736 - qJD(4) * t806 - t723 * t868 + t726 * t872 + t805 * t837;
t877 = mrSges(6,1) * t759 - mrSges(6,2) * t760 + Ifges(6,5) * t789 + Ifges(6,6) * t788 + Ifges(6,3) * t815 + pkin(5) * t742 + t824 * t784 - t823 * t785 - t882;
t718 = -mrSges(5,1) * t804 + mrSges(5,3) * t780 + Ifges(5,4) * t820 + Ifges(5,2) * t819 + Ifges(5,6) * qJDD(4) - pkin(4) * t736 + qJD(4) * t807 - t838 * t805 - t877;
t840 = t891 * qJD(1);
t708 = -mrSges(4,1) * t827 + mrSges(4,3) * t817 - pkin(3) * t883 + pkin(7) * t896 + qJDD(1) * t892 + t869 * t711 + t873 * t718 - t840 * t905;
t710 = mrSges(4,2) * t827 - mrSges(4,3) * t816 - pkin(7) * t724 + qJDD(1) * t893 + t873 * t711 - t869 * t718 - t840 * t906;
t715 = qJDD(1) * mrSges(3,2) - t885;
t731 = t876 * t899 + t881;
t879 = -mrSges(2,2) * t843 - mrSges(3,3) * t831 - pkin(1) * t715 - qJ(3) * t717 - t708 * t865 + t866 * t710 + qJ(2) * (t731 + t915) + mrSges(3,2) * t836 + mrSges(2,1) * t842 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t716 = -m(3) * g(3) + t897;
t707 = t880 + pkin(3) * t724 + pkin(2) * t717 - qJ(2) * t716 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t891 + t911) * qJDD(1) - mrSges(2,3) * t842 + mrSges(3,1) * t836 + mrSges(4,1) * t816 - mrSges(4,2) * t817 + (t865 * t893 + t866 * t892 + t910) * t876;
t706 = -mrSges(3,1) * t831 + mrSges(2,3) * t843 - pkin(1) * t716 + pkin(2) * t731 + g(3) * t912 - qJ(3) * t897 - qJDD(1) * t910 - t866 * t708 - t865 * t710 + t876 * t911;
t1 = [-m(1) * g(1) + t898; -m(1) * g(2) + t908; (-m(1) - m(2) - m(3)) * g(3) + t897; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t908 - t870 * t706 + t874 * t707; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t898 + t874 * t706 + t870 * t707; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t879; t879; t715; t731; t880; t877; -t882;];
tauJB  = t1;

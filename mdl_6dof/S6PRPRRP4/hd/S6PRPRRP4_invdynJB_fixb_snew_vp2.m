% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:47:57
% EndTime: 2019-05-04 23:48:07
% DurationCPUTime: 9.31s
% Computational Cost: add. (145400->296), mult. (319185->367), div. (0->0), fcn. (236323->12), ass. (0->139)
t912 = Ifges(6,1) + Ifges(7,1);
t905 = Ifges(6,4) - Ifges(7,5);
t904 = -Ifges(6,5) - Ifges(7,4);
t911 = Ifges(6,2) + Ifges(7,3);
t903 = Ifges(6,6) - Ifges(7,6);
t910 = -Ifges(6,3) - Ifges(7,2);
t866 = qJD(2) ^ 2;
t854 = sin(pkin(11));
t857 = cos(pkin(11));
t861 = sin(qJ(4));
t863 = cos(qJ(4));
t874 = t854 * t861 - t857 * t863;
t833 = t874 * qJD(2);
t855 = sin(pkin(10));
t858 = cos(pkin(10));
t840 = g(1) * t855 - g(2) * t858;
t841 = -g(1) * t858 - g(2) * t855;
t853 = -g(3) + qJDD(1);
t856 = sin(pkin(6));
t859 = cos(pkin(6));
t862 = sin(qJ(2));
t864 = cos(qJ(2));
t813 = -t862 * t841 + (t840 * t859 + t853 * t856) * t864;
t875 = t854 * t863 + t857 * t861;
t834 = t875 * qJD(2);
t888 = t834 * qJD(4);
t822 = -t874 * qJDD(2) - t888;
t889 = t833 * qJD(4);
t823 = t875 * qJDD(2) - t889;
t860 = sin(qJ(5));
t908 = cos(qJ(5));
t824 = -t908 * qJD(4) + t860 * t834;
t795 = -t824 * qJD(5) + t860 * qJDD(4) + t908 * t823;
t825 = t860 * qJD(4) + t908 * t834;
t800 = mrSges(7,1) * t824 - mrSges(7,3) * t825;
t897 = t859 * t862;
t898 = t856 * t862;
t814 = t840 * t897 + t864 * t841 + t853 * t898;
t809 = -pkin(2) * t866 + qJDD(2) * qJ(3) + t814;
t830 = -t840 * t856 + t853 * t859;
t887 = qJD(2) * qJD(3);
t891 = t857 * t830 - 0.2e1 * t854 * t887;
t907 = pkin(3) * t857;
t780 = (-pkin(8) * qJDD(2) + t866 * t907 - t809) * t854 + t891;
t783 = t854 * t830 + (t809 + 0.2e1 * t887) * t857;
t886 = qJDD(2) * t857;
t852 = t857 ^ 2;
t899 = t852 * t866;
t781 = -pkin(3) * t899 + pkin(8) * t886 + t783;
t774 = t861 * t780 + t863 * t781;
t821 = pkin(4) * t833 - pkin(9) * t834;
t865 = qJD(4) ^ 2;
t772 = -pkin(4) * t865 + qJDD(4) * pkin(9) - t821 * t833 + t774;
t851 = t854 ^ 2;
t871 = qJDD(3) - t813;
t796 = (-pkin(2) - t907) * qJDD(2) + (-qJ(3) + (-t851 - t852) * pkin(8)) * t866 + t871;
t776 = (-t823 + t889) * pkin(9) + (-t822 + t888) * pkin(4) + t796;
t768 = -t860 * t772 + t908 * t776;
t799 = pkin(5) * t824 - qJ(6) * t825;
t820 = qJDD(5) - t822;
t832 = qJD(5) + t833;
t831 = t832 ^ 2;
t766 = -t820 * pkin(5) - t831 * qJ(6) + t825 * t799 + qJDD(6) - t768;
t805 = -mrSges(7,2) * t824 + mrSges(7,3) * t832;
t880 = -m(7) * t766 + t820 * mrSges(7,1) + t832 * t805;
t762 = t795 * mrSges(7,2) + t825 * t800 - t880;
t769 = t908 * t772 + t860 * t776;
t765 = -pkin(5) * t831 + qJ(6) * t820 + 0.2e1 * qJD(6) * t832 - t799 * t824 + t769;
t794 = t825 * qJD(5) - t908 * qJDD(4) + t860 * t823;
t808 = -mrSges(7,1) * t832 + mrSges(7,2) * t825;
t885 = m(7) * t765 + t820 * mrSges(7,3) + t832 * t808;
t893 = t905 * t824 - t825 * t912 + t904 * t832;
t894 = t824 * t911 - t825 * t905 - t832 * t903;
t909 = -t903 * t794 - t904 * t795 - t910 * t820 - t893 * t824 - t894 * t825 + mrSges(6,1) * t768 - mrSges(7,1) * t766 - mrSges(6,2) * t769 + mrSges(7,3) * t765 - pkin(5) * t762 + qJ(6) * (-t794 * mrSges(7,2) - t824 * t800 + t885);
t906 = -mrSges(6,3) - mrSges(7,2);
t901 = mrSges(4,2) * t854;
t804 = -qJDD(2) * pkin(2) - t866 * qJ(3) + t871;
t807 = mrSges(6,1) * t832 - mrSges(6,3) * t825;
t892 = -mrSges(6,1) * t824 - mrSges(6,2) * t825 - t800;
t757 = m(6) * t769 - t820 * mrSges(6,2) + t906 * t794 - t832 * t807 + t892 * t824 + t885;
t806 = -mrSges(6,2) * t832 - mrSges(6,3) * t824;
t759 = m(6) * t768 + t820 * mrSges(6,1) + t906 * t795 + t832 * t806 + t892 * t825 + t880;
t751 = t860 * t757 + t908 * t759;
t828 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t833;
t829 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t834;
t870 = m(5) * t796 - t822 * mrSges(5,1) + t823 * mrSges(5,2) + t833 * t828 + t834 * t829 + t751;
t868 = -m(4) * t804 + mrSges(4,1) * t886 - t870 + (t851 * t866 + t899) * mrSges(4,3);
t744 = (mrSges(3,1) - t901) * qJDD(2) + t868 - t866 * mrSges(3,2) + m(3) * t813;
t900 = t744 * t864;
t752 = t908 * t757 - t759 * t860;
t818 = mrSges(5,1) * t833 + mrSges(5,2) * t834;
t747 = m(5) * t774 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t822 - qJD(4) * t829 - t818 * t833 + t752;
t773 = t863 * t780 - t861 * t781;
t771 = -qJDD(4) * pkin(4) - t865 * pkin(9) + t834 * t821 - t773;
t767 = -0.2e1 * qJD(6) * t825 + (t824 * t832 - t795) * qJ(6) + (t825 * t832 + t794) * pkin(5) + t771;
t763 = m(7) * t767 + t794 * mrSges(7,1) - t795 * mrSges(7,3) + t805 * t824 - t825 * t808;
t760 = -m(6) * t771 - t794 * mrSges(6,1) - t795 * mrSges(6,2) - t824 * t806 - t807 * t825 - t763;
t754 = m(5) * t773 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t823 + qJD(4) * t828 - t818 * t834 + t760;
t741 = t861 * t747 + t863 * t754;
t782 = -t809 * t854 + t891;
t873 = mrSges(4,3) * qJDD(2) + t866 * (-mrSges(4,1) * t857 + t901);
t739 = m(4) * t782 - t873 * t854 + t741;
t882 = t863 * t747 - t861 * t754;
t740 = m(4) * t783 + t873 * t857 + t882;
t883 = -t739 * t854 + t857 * t740;
t730 = m(3) * t814 - mrSges(3,1) * t866 - qJDD(2) * mrSges(3,2) + t883;
t733 = t857 * t739 + t854 * t740;
t732 = m(3) * t830 + t733;
t721 = t730 * t897 - t732 * t856 + t859 * t900;
t719 = m(2) * t840 + t721;
t726 = t864 * t730 - t744 * t862;
t725 = m(2) * t841 + t726;
t896 = t858 * t719 + t855 * t725;
t895 = t824 * t903 + t825 * t904 + t832 * t910;
t877 = Ifges(4,5) * t854 + Ifges(4,6) * t857;
t890 = t866 * t877;
t720 = t730 * t898 + t859 * t732 + t856 * t900;
t884 = -t719 * t855 + t858 * t725;
t881 = m(2) * t853 + t720;
t879 = Ifges(4,1) * t854 + Ifges(4,4) * t857;
t878 = Ifges(4,4) * t854 + Ifges(4,2) * t857;
t749 = -mrSges(6,1) * t771 - mrSges(7,1) * t767 + mrSges(7,2) * t765 + mrSges(6,3) * t769 - pkin(5) * t763 - t794 * t911 + t905 * t795 + t903 * t820 + t895 * t825 - t893 * t832;
t750 = mrSges(6,2) * t771 + mrSges(7,2) * t766 - mrSges(6,3) * t768 - mrSges(7,3) * t767 - qJ(6) * t763 - t905 * t794 + t795 * t912 - t904 * t820 + t895 * t824 + t894 * t832;
t810 = Ifges(5,5) * t834 - Ifges(5,6) * t833 + Ifges(5,3) * qJD(4);
t811 = Ifges(5,4) * t834 - Ifges(5,2) * t833 + Ifges(5,6) * qJD(4);
t734 = mrSges(5,2) * t796 - mrSges(5,3) * t773 + Ifges(5,1) * t823 + Ifges(5,4) * t822 + Ifges(5,5) * qJDD(4) - pkin(9) * t751 - qJD(4) * t811 - t860 * t749 + t908 * t750 - t833 * t810;
t812 = Ifges(5,1) * t834 - Ifges(5,4) * t833 + Ifges(5,5) * qJD(4);
t735 = -mrSges(5,1) * t796 + mrSges(5,3) * t774 + Ifges(5,4) * t823 + Ifges(5,2) * t822 + Ifges(5,6) * qJDD(4) - pkin(4) * t751 + qJD(4) * t812 - t834 * t810 - t909;
t717 = -mrSges(4,1) * t804 + mrSges(4,3) * t783 - pkin(3) * t870 + pkin(8) * t882 + t878 * qJDD(2) + t861 * t734 + t863 * t735 - t854 * t890;
t722 = mrSges(4,2) * t804 - mrSges(4,3) * t782 - pkin(8) * t741 + t879 * qJDD(2) + t863 * t734 - t861 * t735 + t857 * t890;
t715 = mrSges(3,2) * t830 - mrSges(3,3) * t813 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t866 - qJ(3) * t733 - t717 * t854 + t722 * t857;
t867 = mrSges(5,1) * t773 - mrSges(5,2) * t774 + Ifges(5,5) * t823 + Ifges(5,6) * t822 + Ifges(5,3) * qJDD(4) + pkin(4) * t760 + pkin(9) * t752 + t908 * t749 + t860 * t750 + t834 * t811 + t833 * t812;
t716 = -pkin(2) * t733 - t867 + (Ifges(3,6) - t877) * qJDD(2) - mrSges(3,1) * t830 + mrSges(3,3) * t814 - mrSges(4,1) * t782 + mrSges(4,2) * t783 - pkin(3) * t741 + (-t854 * t878 + t857 * t879 + Ifges(3,5)) * t866;
t872 = pkin(7) * t726 + t715 * t862 + t716 * t864;
t748 = qJDD(2) * t901 - t868;
t714 = mrSges(3,1) * t813 - mrSges(3,2) * t814 + Ifges(3,3) * qJDD(2) - pkin(2) * t748 + qJ(3) * t883 + t857 * t717 + t854 * t722;
t713 = mrSges(2,2) * t853 - mrSges(2,3) * t840 + t864 * t715 - t862 * t716 + (-t720 * t856 - t721 * t859) * pkin(7);
t712 = -mrSges(2,1) * t853 + mrSges(2,3) * t841 - pkin(1) * t720 - t856 * t714 + t872 * t859;
t1 = [-m(1) * g(1) + t884; -m(1) * g(2) + t896; -m(1) * g(3) + t881; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t896 - t855 * t712 + t858 * t713; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t884 + t858 * t712 + t855 * t713; -mrSges(1,1) * g(2) + mrSges(2,1) * t840 + mrSges(1,2) * g(1) - mrSges(2,2) * t841 + pkin(1) * t721 + t859 * t714 + t872 * t856; t881; t714; t748; t867; t909; t762;];
tauJB  = t1;

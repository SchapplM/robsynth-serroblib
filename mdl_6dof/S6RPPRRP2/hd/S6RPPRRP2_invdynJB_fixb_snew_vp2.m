% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:47:24
% EndTime: 2019-05-05 14:47:32
% DurationCPUTime: 7.71s
% Computational Cost: add. (90378->298), mult. (199250->358), div. (0->0), fcn. (136185->10), ass. (0->132)
t881 = Ifges(6,1) + Ifges(7,1);
t874 = Ifges(6,4) - Ifges(7,5);
t873 = -Ifges(6,5) - Ifges(7,4);
t880 = Ifges(6,2) + Ifges(7,3);
t872 = Ifges(6,6) - Ifges(7,6);
t879 = -Ifges(6,3) - Ifges(7,2);
t839 = qJD(1) ^ 2;
t829 = sin(pkin(10));
t831 = cos(pkin(10));
t834 = sin(qJ(4));
t836 = cos(qJ(4));
t846 = t829 * t834 - t831 * t836;
t801 = t846 * qJD(1);
t847 = t829 * t836 + t831 * t834;
t802 = t847 * qJD(1);
t860 = t802 * qJD(4);
t790 = -qJDD(1) * t846 - t860;
t861 = t801 * qJD(4);
t791 = qJDD(1) * t847 - t861;
t833 = sin(qJ(5));
t877 = cos(qJ(5));
t795 = -qJD(4) * t877 + t833 * t802;
t762 = -t795 * qJD(5) + t833 * qJDD(4) + t791 * t877;
t796 = t833 * qJD(4) + t802 * t877;
t770 = mrSges(7,1) * t795 - mrSges(7,3) * t796;
t835 = sin(qJ(1));
t837 = cos(qJ(1));
t811 = t835 * g(1) - g(2) * t837;
t808 = qJDD(1) * pkin(1) + t811;
t812 = -g(1) * t837 - g(2) * t835;
t809 = -pkin(1) * t839 + t812;
t830 = sin(pkin(9));
t832 = cos(pkin(9));
t794 = t830 * t808 + t832 * t809;
t782 = -pkin(2) * t839 + qJDD(1) * qJ(3) + t794;
t828 = -g(3) + qJDD(2);
t859 = qJD(1) * qJD(3);
t863 = t831 * t828 - 0.2e1 * t829 * t859;
t876 = pkin(3) * t831;
t760 = (-pkin(7) * qJDD(1) + t839 * t876 - t782) * t829 + t863;
t768 = t829 * t828 + (t782 + 0.2e1 * t859) * t831;
t858 = qJDD(1) * t831;
t825 = t831 ^ 2;
t869 = t825 * t839;
t763 = -pkin(3) * t869 + pkin(7) * t858 + t768;
t746 = t834 * t760 + t836 * t763;
t789 = pkin(4) * t801 - pkin(8) * t802;
t838 = qJD(4) ^ 2;
t742 = -pkin(4) * t838 + qJDD(4) * pkin(8) - t789 * t801 + t746;
t824 = t829 ^ 2;
t793 = t832 * t808 - t830 * t809;
t848 = qJDD(3) - t793;
t766 = (-pkin(2) - t876) * qJDD(1) + (-qJ(3) + (-t824 - t825) * pkin(7)) * t839 + t848;
t744 = (-t791 + t861) * pkin(8) + (-t790 + t860) * pkin(4) + t766;
t738 = -t833 * t742 + t744 * t877;
t769 = pkin(5) * t795 - qJ(6) * t796;
t788 = qJDD(5) - t790;
t800 = qJD(5) + t801;
t799 = t800 ^ 2;
t736 = -t788 * pkin(5) - t799 * qJ(6) + t796 * t769 + qJDD(6) - t738;
t773 = -mrSges(7,2) * t795 + mrSges(7,3) * t800;
t852 = -m(7) * t736 + t788 * mrSges(7,1) + t800 * t773;
t732 = t762 * mrSges(7,2) + t796 * t770 - t852;
t739 = t877 * t742 + t833 * t744;
t735 = -pkin(5) * t799 + qJ(6) * t788 + 0.2e1 * qJD(6) * t800 - t769 * t795 + t739;
t761 = t796 * qJD(5) - qJDD(4) * t877 + t833 * t791;
t776 = -mrSges(7,1) * t800 + mrSges(7,2) * t796;
t857 = m(7) * t735 + t788 * mrSges(7,3) + t800 * t776;
t865 = t874 * t795 - t796 * t881 + t873 * t800;
t866 = t795 * t880 - t796 * t874 - t800 * t872;
t878 = -t761 * t872 - t762 * t873 - t879 * t788 - t795 * t865 - t796 * t866 + mrSges(6,1) * t738 - mrSges(7,1) * t736 - mrSges(6,2) * t739 + mrSges(7,3) * t735 - pkin(5) * t732 + qJ(6) * (-t761 * mrSges(7,2) - t795 * t770 + t857);
t875 = -mrSges(6,3) - mrSges(7,2);
t870 = mrSges(4,2) * t829;
t775 = mrSges(6,1) * t800 - mrSges(6,3) * t796;
t864 = -mrSges(6,1) * t795 - mrSges(6,2) * t796 - t770;
t727 = m(6) * t739 - t788 * mrSges(6,2) + t761 * t875 - t800 * t775 + t795 * t864 + t857;
t774 = -mrSges(6,2) * t800 - mrSges(6,3) * t795;
t729 = m(6) * t738 + t788 * mrSges(6,1) + t762 * t875 + t800 * t774 + t796 * t864 + t852;
t722 = t877 * t727 - t729 * t833;
t786 = mrSges(5,1) * t801 + mrSges(5,2) * t802;
t798 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t802;
t718 = m(5) * t746 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t790 - qJD(4) * t798 - t786 * t801 + t722;
t745 = t836 * t760 - t834 * t763;
t741 = -qJDD(4) * pkin(4) - t838 * pkin(8) + t802 * t789 - t745;
t737 = -0.2e1 * qJD(6) * t796 + (t795 * t800 - t762) * qJ(6) + (t796 * t800 + t761) * pkin(5) + t741;
t733 = m(7) * t737 + mrSges(7,1) * t761 - t762 * mrSges(7,3) + t773 * t795 - t796 * t776;
t730 = -m(6) * t741 - t761 * mrSges(6,1) - mrSges(6,2) * t762 - t795 * t774 - t775 * t796 - t733;
t797 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t801;
t724 = m(5) * t745 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t791 + qJD(4) * t797 - t786 * t802 + t730;
t711 = t834 * t718 + t836 * t724;
t767 = -t782 * t829 + t863;
t845 = mrSges(4,3) * qJDD(1) + t839 * (-mrSges(4,1) * t831 + t870);
t709 = m(4) * t767 - t829 * t845 + t711;
t853 = t836 * t718 - t834 * t724;
t710 = m(4) * t768 + t831 * t845 + t853;
t854 = -t709 * t829 + t831 * t710;
t700 = m(3) * t794 - mrSges(3,1) * t839 - qJDD(1) * mrSges(3,2) + t854;
t778 = -qJDD(1) * pkin(2) - t839 * qJ(3) + t848;
t721 = t833 * t727 + t877 * t729;
t844 = m(5) * t766 - t790 * mrSges(5,1) + t791 * mrSges(5,2) + t801 * t797 + t802 * t798 + t721;
t841 = -m(4) * t778 + mrSges(4,1) * t858 - t844 + (t824 * t839 + t869) * mrSges(4,3);
t713 = (mrSges(3,1) - t870) * qJDD(1) + t841 + m(3) * t793 - t839 * mrSges(3,2);
t697 = t830 * t700 + t832 * t713;
t694 = m(2) * t811 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t839 + t697;
t855 = t832 * t700 - t713 * t830;
t695 = m(2) * t812 - mrSges(2,1) * t839 - qJDD(1) * mrSges(2,2) + t855;
t868 = t837 * t694 + t835 * t695;
t703 = t831 * t709 + t829 * t710;
t867 = t795 * t872 + t796 * t873 + t800 * t879;
t849 = Ifges(4,5) * t829 + Ifges(4,6) * t831;
t862 = t839 * t849;
t701 = m(3) * t828 + t703;
t856 = -t694 * t835 + t837 * t695;
t851 = Ifges(4,1) * t829 + Ifges(4,4) * t831;
t850 = Ifges(4,4) * t829 + Ifges(4,2) * t831;
t719 = -mrSges(6,1) * t741 - mrSges(7,1) * t737 + mrSges(7,2) * t735 + mrSges(6,3) * t739 - pkin(5) * t733 - t761 * t880 + t874 * t762 + t872 * t788 + t867 * t796 - t865 * t800;
t720 = mrSges(6,2) * t741 + mrSges(7,2) * t736 - mrSges(6,3) * t738 - mrSges(7,3) * t737 - qJ(6) * t733 - t874 * t761 + t762 * t881 - t873 * t788 + t867 * t795 + t866 * t800;
t779 = Ifges(5,5) * t802 - Ifges(5,6) * t801 + Ifges(5,3) * qJD(4);
t780 = Ifges(5,4) * t802 - Ifges(5,2) * t801 + Ifges(5,6) * qJD(4);
t704 = mrSges(5,2) * t766 - mrSges(5,3) * t745 + Ifges(5,1) * t791 + Ifges(5,4) * t790 + Ifges(5,5) * qJDD(4) - pkin(8) * t721 - qJD(4) * t780 - t833 * t719 + t720 * t877 - t801 * t779;
t781 = Ifges(5,1) * t802 - Ifges(5,4) * t801 + Ifges(5,5) * qJD(4);
t705 = -mrSges(5,1) * t766 + mrSges(5,3) * t746 + Ifges(5,4) * t791 + Ifges(5,2) * t790 + Ifges(5,6) * qJDD(4) - pkin(4) * t721 + qJD(4) * t781 - t802 * t779 - t878;
t688 = -mrSges(4,1) * t778 + mrSges(4,3) * t768 - pkin(3) * t844 + pkin(7) * t853 + qJDD(1) * t850 + t834 * t704 + t836 * t705 - t829 * t862;
t690 = mrSges(4,2) * t778 - mrSges(4,3) * t767 - pkin(7) * t711 + qJDD(1) * t851 + t836 * t704 - t834 * t705 + t831 * t862;
t715 = qJDD(1) * t870 - t841;
t843 = mrSges(2,1) * t811 + mrSges(3,1) * t793 - mrSges(2,2) * t812 - mrSges(3,2) * t794 + pkin(1) * t697 - pkin(2) * t715 + qJ(3) * t854 + t831 * t688 + t829 * t690 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t840 = mrSges(5,1) * t745 - mrSges(5,2) * t746 + Ifges(5,5) * t791 + Ifges(5,6) * t790 + Ifges(5,3) * qJDD(4) + pkin(4) * t730 + pkin(8) * t722 + t719 * t877 + t833 * t720 + t802 * t780 + t801 * t781;
t686 = -pkin(3) * t711 - t840 - pkin(2) * t703 + (Ifges(3,6) - t849) * qJDD(1) - mrSges(4,1) * t767 + mrSges(4,2) * t768 + mrSges(3,3) * t794 - mrSges(3,1) * t828 + (-t829 * t850 + t831 * t851 + Ifges(3,5)) * t839;
t685 = mrSges(3,2) * t828 - mrSges(3,3) * t793 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t839 - qJ(3) * t703 - t688 * t829 + t690 * t831;
t684 = -mrSges(2,2) * g(3) - mrSges(2,3) * t811 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t839 - qJ(2) * t697 + t685 * t832 - t686 * t830;
t683 = mrSges(2,1) * g(3) + mrSges(2,3) * t812 + t839 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t701 + qJ(2) * t855 + t830 * t685 + t832 * t686;
t1 = [-m(1) * g(1) + t856; -m(1) * g(2) + t868; (-m(1) - m(2)) * g(3) + t701; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t868 - t835 * t683 + t837 * t684; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t856 + t837 * t683 + t835 * t684; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t843; t843; t701; t715; t840; t878; t732;];
tauJB  = t1;

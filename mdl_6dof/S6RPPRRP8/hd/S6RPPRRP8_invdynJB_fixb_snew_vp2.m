% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:06:37
% EndTime: 2019-05-05 15:06:43
% DurationCPUTime: 5.07s
% Computational Cost: add. (52496->290), mult. (117245->339), div. (0->0), fcn. (79353->8), ass. (0->128)
t898 = Ifges(6,1) + Ifges(7,1);
t885 = Ifges(6,4) - Ifges(7,5);
t884 = -Ifges(6,5) - Ifges(7,4);
t897 = Ifges(6,2) + Ifges(7,3);
t882 = Ifges(6,6) - Ifges(7,6);
t896 = -Ifges(6,3) - Ifges(7,2);
t839 = sin(qJ(1));
t841 = cos(qJ(1));
t814 = t839 * g(1) - g(2) * t841;
t843 = qJD(1) ^ 2;
t852 = -t843 * qJ(2) + qJDD(2) - t814;
t880 = -pkin(1) - qJ(3);
t895 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t880 + t852;
t815 = -t841 * g(1) - t839 * g(2);
t894 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t815;
t835 = sin(pkin(9));
t836 = cos(pkin(9));
t838 = sin(qJ(4));
t840 = cos(qJ(4));
t856 = t835 * t840 + t836 * t838;
t809 = t856 * qJD(1);
t804 = t843 * pkin(1) - t894;
t893 = -m(3) * t804 + mrSges(3,2) * t843 + qJDD(1) * mrSges(3,3);
t855 = -t835 * t838 + t836 * t840;
t810 = t855 * qJD(1);
t870 = t810 * qJD(4);
t793 = -qJDD(1) * t856 - t870;
t871 = t809 * qJD(4);
t794 = qJDD(1) * t855 - t871;
t837 = sin(qJ(5));
t891 = cos(qJ(5));
t796 = -qJD(4) * t891 + t810 * t837;
t760 = -qJD(5) * t796 + qJDD(4) * t837 + t794 * t891;
t797 = qJD(4) * t837 + t810 * t891;
t766 = mrSges(7,1) * t796 - mrSges(7,3) * t797;
t790 = t835 * g(3) + t836 * t895;
t889 = pkin(3) * t843;
t771 = (-pkin(7) * qJDD(1) - t835 * t889) * t836 + t790;
t791 = -t836 * g(3) + t835 * t895;
t826 = t835 ^ 2;
t868 = qJDD(1) * t835;
t772 = -pkin(7) * t868 - t826 * t889 + t791;
t748 = t771 * t838 + t772 * t840;
t792 = pkin(4) * t809 - pkin(8) * t810;
t842 = qJD(4) ^ 2;
t743 = -pkin(4) * t842 + qJDD(4) * pkin(8) - t792 * t809 + t748;
t851 = qJDD(3) + t894;
t873 = -t836 ^ 2 - t826;
t778 = pkin(3) * t868 + (pkin(7) * t873 + t880) * t843 + t851;
t745 = (-t794 + t871) * pkin(8) + (-t793 + t870) * pkin(4) + t778;
t739 = -t743 * t837 + t745 * t891;
t765 = pkin(5) * t796 - qJ(6) * t797;
t789 = qJDD(5) - t793;
t807 = qJD(5) + t809;
t806 = t807 ^ 2;
t737 = -pkin(5) * t789 - qJ(6) * t806 + t765 * t797 + qJDD(6) - t739;
t773 = -mrSges(7,2) * t796 + mrSges(7,3) * t807;
t860 = -m(7) * t737 + mrSges(7,1) * t789 + t773 * t807;
t733 = t760 * mrSges(7,2) + t797 * t766 - t860;
t740 = t743 * t891 + t745 * t837;
t736 = -pkin(5) * t806 + qJ(6) * t789 + 0.2e1 * qJD(6) * t807 - t765 * t796 + t740;
t759 = qJD(5) * t797 - qJDD(4) * t891 + t794 * t837;
t776 = -mrSges(7,1) * t807 + mrSges(7,2) * t797;
t866 = m(7) * t736 + mrSges(7,3) * t789 + t776 * t807;
t875 = t796 * t885 - t797 * t898 + t807 * t884;
t876 = t796 * t897 - t797 * t885 - t807 * t882;
t892 = -t882 * t759 - t884 * t760 - t896 * t789 - t875 * t796 - t876 * t797 + mrSges(6,1) * t739 - mrSges(7,1) * t737 - mrSges(6,2) * t740 + mrSges(7,3) * t736 - pkin(5) * t733 + qJ(6) * (-t759 * mrSges(7,2) - t796 * t766 + t866);
t888 = mrSges(2,1) - mrSges(3,2);
t887 = -mrSges(6,3) - mrSges(7,2);
t886 = -Ifges(3,4) + Ifges(2,5);
t883 = -Ifges(2,6) + Ifges(3,5);
t879 = mrSges(4,2) * t836;
t785 = mrSges(5,1) * t809 + mrSges(5,2) * t810;
t802 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t810;
t775 = mrSges(6,1) * t807 - mrSges(6,3) * t797;
t874 = -mrSges(6,1) * t796 - mrSges(6,2) * t797 - t766;
t728 = m(6) * t740 - t789 * mrSges(6,2) + t759 * t887 - t807 * t775 + t796 * t874 + t866;
t774 = -mrSges(6,2) * t807 - mrSges(6,3) * t796;
t730 = m(6) * t739 + t789 * mrSges(6,1) + t760 * t887 + t807 * t774 + t797 * t874 + t860;
t861 = t728 * t891 - t730 * t837;
t720 = m(5) * t748 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t793 - qJD(4) * t802 - t785 * t809 + t861;
t747 = t840 * t771 - t772 * t838;
t801 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t809;
t742 = -qJDD(4) * pkin(4) - t842 * pkin(8) + t792 * t810 - t747;
t738 = -0.2e1 * qJD(6) * t797 + (t796 * t807 - t760) * qJ(6) + (t797 * t807 + t759) * pkin(5) + t742;
t734 = m(7) * t738 + mrSges(7,1) * t759 - mrSges(7,3) * t760 + t773 * t796 - t776 * t797;
t844 = -m(6) * t742 - mrSges(6,1) * t759 - mrSges(6,2) * t760 - t774 * t796 - t775 * t797 - t734;
t725 = m(5) * t747 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t794 + qJD(4) * t801 - t785 * t810 + t844;
t709 = t720 * t838 + t725 * t840;
t854 = -qJDD(1) * mrSges(4,3) - t843 * (mrSges(4,1) * t835 + t879);
t707 = m(4) * t790 + t836 * t854 + t709;
t862 = t720 * t840 - t838 * t725;
t708 = m(4) * t791 + t835 * t854 + t862;
t704 = t836 * t707 + t835 * t708;
t808 = -qJDD(1) * pkin(1) + t852;
t850 = -m(3) * t808 + mrSges(3,3) * t843 - t704;
t699 = m(2) * t814 - t843 * mrSges(2,2) + qJDD(1) * t888 + t850;
t800 = t843 * t880 + t851;
t723 = t728 * t837 + t730 * t891;
t849 = m(5) * t778 - t793 * mrSges(5,1) + mrSges(5,2) * t794 + t809 * t801 + t802 * t810 + t723;
t848 = m(4) * t800 + mrSges(4,1) * t868 + qJDD(1) * t879 + t849;
t865 = t873 * mrSges(4,3);
t712 = -qJDD(1) * mrSges(2,2) + m(2) * t815 + t848 + (-mrSges(2,1) + t865) * t843 + t893;
t878 = t699 * t841 + t712 * t839;
t877 = t796 * t882 + t797 * t884 + t807 * t896;
t857 = Ifges(4,5) * t836 - Ifges(4,6) * t835;
t872 = t843 * t857;
t864 = -t699 * t839 + t712 * t841;
t863 = -t835 * t707 + t708 * t836;
t859 = Ifges(4,1) * t836 - Ifges(4,4) * t835;
t858 = Ifges(4,4) * t836 - Ifges(4,2) * t835;
t716 = -mrSges(6,1) * t742 - mrSges(7,1) * t738 + mrSges(7,2) * t736 + mrSges(6,3) * t740 - pkin(5) * t734 - t759 * t897 + t760 * t885 + t789 * t882 + t797 * t877 - t807 * t875;
t721 = mrSges(6,2) * t742 + mrSges(7,2) * t737 - mrSges(6,3) * t739 - mrSges(7,3) * t738 - qJ(6) * t734 - t759 * t885 + t760 * t898 - t789 * t884 + t796 * t877 + t807 * t876;
t780 = Ifges(5,4) * t810 - Ifges(5,2) * t809 + Ifges(5,6) * qJD(4);
t781 = Ifges(5,1) * t810 - Ifges(5,4) * t809 + Ifges(5,5) * qJD(4);
t847 = mrSges(5,1) * t747 - mrSges(5,2) * t748 + Ifges(5,5) * t794 + Ifges(5,6) * t793 + Ifges(5,3) * qJDD(4) + pkin(4) * t844 + pkin(8) * t861 + t716 * t891 + t721 * t837 + t780 * t810 + t809 * t781;
t779 = Ifges(5,5) * t810 - Ifges(5,6) * t809 + Ifges(5,3) * qJD(4);
t700 = mrSges(5,2) * t778 - mrSges(5,3) * t747 + Ifges(5,1) * t794 + Ifges(5,4) * t793 + Ifges(5,5) * qJDD(4) - pkin(8) * t723 - qJD(4) * t780 - t716 * t837 + t721 * t891 - t779 * t809;
t705 = -mrSges(5,1) * t778 + mrSges(5,3) * t748 + Ifges(5,4) * t794 + Ifges(5,2) * t793 + Ifges(5,6) * qJDD(4) - pkin(4) * t723 + qJD(4) * t781 - t810 * t779 - t892;
t695 = -mrSges(4,1) * t800 + mrSges(4,3) * t791 - pkin(3) * t849 + pkin(7) * t862 + qJDD(1) * t858 + t838 * t700 + t840 * t705 - t836 * t872;
t697 = mrSges(4,2) * t800 - mrSges(4,3) * t790 - pkin(7) * t709 + qJDD(1) * t859 + t840 * t700 - t838 * t705 - t835 * t872;
t702 = qJDD(1) * mrSges(3,2) - t850;
t714 = t843 * t865 + t848;
t845 = -mrSges(2,2) * t815 - mrSges(3,3) * t804 - pkin(1) * t702 - qJ(3) * t704 - t835 * t695 + t836 * t697 + qJ(2) * (t714 + t893) + mrSges(3,2) * t808 + mrSges(2,1) * t814 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t703 = -m(3) * g(3) + t863;
t694 = (mrSges(3,3) - mrSges(2,2)) * g(3) - mrSges(2,3) * t814 + mrSges(3,1) * t808 + mrSges(4,1) * t790 - mrSges(4,2) * t791 + t847 + pkin(3) * t709 + pkin(2) * t704 - qJ(2) * t703 + (t857 + t886) * qJDD(1) + (t835 * t859 + t836 * t858 + t883) * t843;
t693 = -mrSges(3,1) * t804 + mrSges(2,3) * t815 - pkin(1) * t703 + pkin(2) * t714 + g(3) * t888 - qJ(3) * t863 - qJDD(1) * t883 - t836 * t695 - t835 * t697 + t843 * t886;
t1 = [-m(1) * g(1) + t864; -m(1) * g(2) + t878; (-m(1) - m(2) - m(3)) * g(3) + t863; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t878 - t693 * t839 + t694 * t841; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t864 + t841 * t693 + t839 * t694; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t845; t845; t702; t714; t847; t892; t733;];
tauJB  = t1;

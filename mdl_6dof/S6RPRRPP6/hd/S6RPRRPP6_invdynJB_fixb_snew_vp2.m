% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-05-05 21:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:41:54
% EndTime: 2019-05-05 21:42:04
% DurationCPUTime: 6.02s
% Computational Cost: add. (69554->314), mult. (138549->370), div. (0->0), fcn. (87249->8), ass. (0->128)
t870 = Ifges(6,1) + Ifges(7,1);
t860 = Ifges(6,4) - Ifges(7,5);
t858 = Ifges(6,5) + Ifges(7,4);
t869 = -Ifges(6,2) - Ifges(7,3);
t868 = -Ifges(7,2) - Ifges(6,3);
t857 = Ifges(6,6) - Ifges(7,6);
t822 = sin(qJ(1));
t825 = cos(qJ(1));
t807 = -t825 * g(1) - t822 * g(2);
t867 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t807;
t827 = qJD(1) ^ 2;
t864 = (-pkin(1) - pkin(7));
t776 = (t864 * t827) - t867;
t821 = sin(qJ(3));
t824 = cos(qJ(3));
t848 = qJD(1) * qJD(3);
t843 = t824 * t848;
t801 = -qJDD(1) * t821 - t843;
t844 = t821 * t848;
t802 = qJDD(1) * t824 - t844;
t745 = (-t802 + t844) * pkin(8) + (-t801 + t843) * pkin(3) + t776;
t806 = t822 * g(1) - t825 * g(2);
t835 = -t827 * qJ(2) + qJDD(2) - t806;
t777 = t864 * qJDD(1) + t835;
t770 = -g(3) * t824 + t821 * t777;
t800 = (pkin(3) * t821 - pkin(8) * t824) * qJD(1);
t826 = qJD(3) ^ 2;
t850 = qJD(1) * t821;
t750 = -pkin(3) * t826 + qJDD(3) * pkin(8) - t800 * t850 + t770;
t820 = sin(qJ(4));
t823 = cos(qJ(4));
t724 = t823 * t745 - t820 * t750;
t849 = qJD(1) * t824;
t797 = qJD(3) * t823 - t820 * t849;
t764 = qJD(4) * t797 + qJDD(3) * t820 + t802 * t823;
t796 = qJDD(4) - t801;
t798 = qJD(3) * t820 + t823 * t849;
t809 = qJD(4) + t850;
t720 = (t797 * t809 - t764) * qJ(5) + (t797 * t798 + t796) * pkin(4) + t724;
t725 = t820 * t745 + t823 * t750;
t763 = -qJD(4) * t798 + qJDD(3) * t823 - t802 * t820;
t772 = pkin(4) * t809 - qJ(5) * t798;
t795 = t797 ^ 2;
t722 = -pkin(4) * t795 + qJ(5) * t763 - t772 * t809 + t725;
t819 = sin(pkin(9));
t856 = cos(pkin(9));
t765 = -t856 * t797 + t819 * t798;
t865 = -2 * qJD(5);
t716 = t819 * t720 + t856 * t722 + t765 * t865;
t735 = -t856 * t763 + t819 * t764;
t766 = t819 * t797 + t856 * t798;
t753 = mrSges(6,1) * t809 - mrSges(6,3) * t766;
t740 = pkin(5) * t765 - qJ(6) * t766;
t808 = t809 ^ 2;
t713 = -pkin(5) * t808 + qJ(6) * t796 + 0.2e1 * qJD(6) * t809 - t740 * t765 + t716;
t754 = -mrSges(7,1) * t809 + mrSges(7,2) * t766;
t845 = m(7) * t713 + t796 * mrSges(7,3) + t809 * t754;
t741 = mrSges(7,1) * t765 - mrSges(7,3) * t766;
t851 = -mrSges(6,1) * t765 - mrSges(6,2) * t766 - t741;
t862 = -mrSges(6,3) - mrSges(7,2);
t703 = m(6) * t716 - t796 * mrSges(6,2) + t862 * t735 - t809 * t753 + t851 * t765 + t845;
t836 = t856 * t720 - t819 * t722;
t715 = t766 * t865 + t836;
t736 = t819 * t763 + t856 * t764;
t751 = -mrSges(6,2) * t809 - mrSges(6,3) * t765;
t714 = -t796 * pkin(5) - t808 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t740) * t766 - t836;
t752 = -mrSges(7,2) * t765 + mrSges(7,3) * t809;
t838 = -m(7) * t714 + t796 * mrSges(7,1) + t809 * t752;
t705 = m(6) * t715 + t796 * mrSges(6,1) + t862 * t736 + t809 * t751 + t851 * t766 + t838;
t698 = t819 * t703 + t856 * t705;
t710 = t736 * mrSges(7,2) + t766 * t741 - t838;
t757 = Ifges(5,4) * t798 + Ifges(5,2) * t797 + Ifges(5,6) * t809;
t758 = Ifges(5,1) * t798 + Ifges(5,4) * t797 + Ifges(5,5) * t809;
t852 = -t860 * t765 + t870 * t766 + t858 * t809;
t853 = t869 * t765 + t860 * t766 + t857 * t809;
t866 = -t857 * t735 + t858 * t736 + t852 * t765 + t853 * t766 + (Ifges(5,3) - t868) * t796 + mrSges(5,1) * t724 + mrSges(6,1) * t715 - mrSges(7,1) * t714 - mrSges(5,2) * t725 - mrSges(6,2) * t716 + mrSges(7,3) * t713 + Ifges(5,5) * t764 + Ifges(5,6) * t763 + pkin(4) * t698 - pkin(5) * t710 + qJ(6) * (-t735 * mrSges(7,2) - t765 * t741 + t845) + t798 * t757 - t797 * t758;
t863 = mrSges(2,1) - mrSges(3,2);
t861 = -Ifges(3,4) + Ifges(2,5);
t859 = (Ifges(3,5) - Ifges(2,6));
t799 = (mrSges(4,1) * t821 + mrSges(4,2) * t824) * qJD(1);
t805 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t849;
t768 = -mrSges(5,1) * t797 + mrSges(5,2) * t798;
t771 = -mrSges(5,2) * t809 + mrSges(5,3) * t797;
t696 = m(5) * t724 + mrSges(5,1) * t796 - mrSges(5,3) * t764 - t768 * t798 + t771 * t809 + t698;
t773 = mrSges(5,1) * t809 - mrSges(5,3) * t798;
t839 = t856 * t703 - t705 * t819;
t697 = m(5) * t725 - mrSges(5,2) * t796 + mrSges(5,3) * t763 + t768 * t797 - t773 * t809 + t839;
t840 = -t696 * t820 + t823 * t697;
t690 = m(4) * t770 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t801 - qJD(3) * t805 - t799 * t850 + t840;
t769 = t821 * g(3) + t824 * t777;
t804 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t850;
t749 = -qJDD(3) * pkin(3) - t826 * pkin(8) + t800 * t849 - t769;
t723 = -t763 * pkin(4) - t795 * qJ(5) + t798 * t772 + qJDD(5) + t749;
t718 = -0.2e1 * qJD(6) * t766 + (t765 * t809 - t736) * qJ(6) + (t766 * t809 + t735) * pkin(5) + t723;
t711 = m(7) * t718 + t735 * mrSges(7,1) - t736 * mrSges(7,3) + t765 * t752 - t766 * t754;
t708 = m(6) * t723 + t735 * mrSges(6,1) + t736 * mrSges(6,2) + t765 * t751 + t766 * t753 + t711;
t829 = -m(5) * t749 + t763 * mrSges(5,1) - t764 * mrSges(5,2) + t797 * t771 - t798 * t773 - t708;
t706 = m(4) * t769 + qJDD(3) * mrSges(4,1) - t802 * mrSges(4,3) + qJD(3) * t804 - t799 * t849 + t829;
t684 = t821 * t690 + t824 * t706;
t782 = -qJDD(1) * pkin(1) + t835;
t834 = -m(3) * t782 + (t827 * mrSges(3,3)) - t684;
t680 = m(2) * t806 - (t827 * mrSges(2,2)) + t863 * qJDD(1) + t834;
t780 = t827 * pkin(1) + t867;
t692 = t823 * t696 + t820 * t697;
t833 = -m(4) * t776 + mrSges(4,1) * t801 - t802 * mrSges(4,2) - t804 * t850 - t805 * t849 - t692;
t831 = -m(3) * t780 + (t827 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t833;
t687 = m(2) * t807 - (mrSges(2,1) * t827) - qJDD(1) * mrSges(2,2) + t831;
t855 = t825 * t680 + t822 * t687;
t854 = t857 * t765 - t858 * t766 + t868 * t809;
t842 = -t680 * t822 + t825 * t687;
t841 = t824 * t690 - t821 * t706;
t699 = -mrSges(6,1) * t723 - mrSges(7,1) * t718 + mrSges(7,2) * t713 + mrSges(6,3) * t716 - pkin(5) * t711 + t869 * t735 + t860 * t736 + t854 * t766 + t857 * t796 + t852 * t809;
t700 = mrSges(6,2) * t723 + mrSges(7,2) * t714 - mrSges(6,3) * t715 - mrSges(7,3) * t718 - qJ(6) * t711 - t860 * t735 + t870 * t736 + t854 * t765 + t858 * t796 - t853 * t809;
t756 = Ifges(5,5) * t798 + Ifges(5,6) * t797 + Ifges(5,3) * t809;
t676 = -mrSges(5,1) * t749 + mrSges(5,3) * t725 + Ifges(5,4) * t764 + Ifges(5,2) * t763 + Ifges(5,6) * t796 - pkin(4) * t708 + qJ(5) * t839 + t856 * t699 + t819 * t700 - t798 * t756 + t809 * t758;
t678 = mrSges(5,2) * t749 - mrSges(5,3) * t724 + Ifges(5,1) * t764 + Ifges(5,4) * t763 + Ifges(5,5) * t796 - qJ(5) * t698 - t819 * t699 + t856 * t700 + t797 * t756 - t809 * t757;
t786 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t824 - Ifges(4,2) * t821) * qJD(1);
t787 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t824 - Ifges(4,4) * t821) * qJD(1);
t832 = mrSges(4,1) * t769 - mrSges(4,2) * t770 + Ifges(4,5) * t802 + Ifges(4,6) * t801 + Ifges(4,3) * qJDD(3) + pkin(3) * t829 + pkin(8) * t840 + t823 * t676 + t820 * t678 + t786 * t849 + t787 * t850;
t785 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t824 - Ifges(4,6) * t821) * qJD(1);
t673 = mrSges(4,2) * t776 - mrSges(4,3) * t769 + Ifges(4,1) * t802 + Ifges(4,4) * t801 + Ifges(4,5) * qJDD(3) - pkin(8) * t692 - qJD(3) * t786 - t676 * t820 + t678 * t823 - t785 * t850;
t674 = -mrSges(4,1) * t776 + mrSges(4,3) * t770 + Ifges(4,4) * t802 + Ifges(4,2) * t801 + Ifges(4,6) * qJDD(3) - pkin(3) * t692 + qJD(3) * t787 - t785 * t849 - t866;
t682 = qJDD(1) * mrSges(3,2) - t834;
t830 = mrSges(2,1) * t806 - mrSges(2,2) * t807 + mrSges(3,2) * t782 - mrSges(3,3) * t780 - pkin(1) * t682 - pkin(7) * t684 + qJ(2) * t831 + t824 * t673 - t674 * t821 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t683 = -m(3) * g(3) + t841;
t671 = -mrSges(2,3) * t806 + mrSges(3,1) * t782 + pkin(2) * t684 - qJ(2) * t683 + t832 + (t859 * t827) + t861 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t670 = -mrSges(3,1) * t780 + mrSges(2,3) * t807 - pkin(1) * t683 - pkin(2) * t833 - pkin(7) * t841 + t863 * g(3) - t859 * qJDD(1) - t821 * t673 - t824 * t674 + t861 * t827;
t1 = [-m(1) * g(1) + t842; -m(1) * g(2) + t855; (-m(1) - m(2) - m(3)) * g(3) + t841; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t855 - t822 * t670 + t825 * t671; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t842 + t825 * t670 + t822 * t671; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t830; t830; t682; t832; t866; t708; t710;];
tauJB  = t1;

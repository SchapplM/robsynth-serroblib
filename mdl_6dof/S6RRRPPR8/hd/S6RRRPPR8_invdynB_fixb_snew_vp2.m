% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 06:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:10:01
% EndTime: 2019-05-07 06:10:18
% DurationCPUTime: 8.14s
% Computational Cost: add. (110152->360), mult. (237866->434), div. (0->0), fcn. (177188->10), ass. (0->150)
t824 = cos(pkin(6));
t820 = qJD(1) * t824 + qJD(2);
t826 = sin(qJ(3));
t827 = sin(qJ(2));
t823 = sin(pkin(6));
t859 = qJD(1) * t823;
t850 = t827 * t859;
t872 = cos(qJ(3));
t790 = t826 * t820 + t872 * t850;
t830 = cos(qJ(2));
t858 = qJD(1) * t830;
t849 = t823 * t858;
t811 = -qJD(3) + t849;
t772 = pkin(4) * t811 - qJ(5) * t790;
t882 = (2 * qJD(4)) + t772;
t881 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t880 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t855 = Ifges(5,5) - Ifges(6,4) - Ifges(4,4);
t854 = Ifges(5,6) - Ifges(6,5) - Ifges(4,6);
t879 = -Ifges(5,3) - Ifges(6,1) - Ifges(4,2);
t878 = Ifges(6,3) + Ifges(4,3) + Ifges(5,2);
t877 = t811 ^ 2;
t876 = -2 * qJD(4);
t874 = -2 * qJD(5);
t873 = -pkin(4) - pkin(10);
t871 = pkin(3) * t811;
t870 = t824 * g(3);
t869 = -mrSges(4,3) - mrSges(5,2);
t789 = -t872 * t820 + t826 * t850;
t868 = t789 * t811;
t867 = t823 * t827;
t866 = t823 * t830;
t865 = t824 * t827;
t864 = t824 * t830;
t828 = sin(qJ(1));
t831 = cos(qJ(1));
t813 = t828 * g(1) - g(2) * t831;
t832 = qJD(1) ^ 2;
t800 = pkin(8) * t823 * t832 + qJDD(1) * pkin(1) + t813;
t814 = -g(1) * t831 - g(2) * t828;
t857 = qJDD(1) * t823;
t801 = -pkin(1) * t832 + pkin(8) * t857 + t814;
t860 = t800 * t865 + t830 * t801;
t758 = -g(3) * t867 + t860;
t798 = mrSges(3,1) * t820 - mrSges(3,3) * t850;
t802 = (-mrSges(3,1) * t830 + mrSges(3,2) * t827) * t859;
t805 = -qJD(2) * t850 + t830 * t857;
t819 = qJDD(1) * t824 + qJDD(2);
t803 = (-pkin(2) * t830 - pkin(9) * t827) * t859;
t818 = t820 ^ 2;
t725 = -t818 * pkin(2) + t819 * pkin(9) + (-g(3) * t827 + t803 * t858) * t823 + t860;
t804 = (qJD(2) * t858 + qJDD(1) * t827) * t823;
t726 = -t805 * pkin(2) - t804 * pkin(9) - t870 + (-t800 + (pkin(2) * t827 - pkin(9) * t830) * t820 * qJD(1)) * t823;
t712 = -t826 * t725 + t872 * t726;
t756 = -t789 * qJD(3) + t872 * t804 + t826 * t819;
t759 = mrSges(6,1) * t790 + mrSges(6,2) * t789;
t771 = mrSges(4,2) * t811 - mrSges(4,3) * t789;
t797 = -qJDD(3) + t805;
t760 = pkin(3) * t789 - qJ(4) * t790;
t840 = -qJ(4) * t877 + t790 * t760 + qJDD(4) - t712;
t710 = t797 * pkin(3) + t840;
t776 = -mrSges(5,2) * t789 - mrSges(5,3) * t811;
t834 = (-t756 + t868) * qJ(5) + t840 + (t789 * pkin(4) + t874) * t790;
t704 = (pkin(3) + pkin(4)) * t797 + t834;
t770 = mrSges(6,1) * t811 - mrSges(6,3) * t789;
t755 = qJD(3) * t790 + t804 * t826 - t872 * t819;
t788 = t789 ^ 2;
t757 = -g(3) * t866 + t800 * t864 - t827 * t801;
t724 = -t819 * pkin(2) - t818 * pkin(9) + t803 * t850 - t757;
t841 = t755 * pkin(3) + t724 + (-t756 - t868) * qJ(4);
t836 = -qJ(5) * t788 + t790 * t882 + qJDD(5) - t841;
t700 = t836 + pkin(5) * t756 + t873 * t755 + (pkin(5) * t789 + (pkin(3) + pkin(10)) * t790) * t811;
t763 = pkin(5) * t790 - pkin(10) * t789;
t701 = -t877 * pkin(5) - t790 * t763 + (pkin(3) - t873) * t797 + t834;
t825 = sin(qJ(6));
t829 = cos(qJ(6));
t698 = t700 * t829 - t701 * t825;
t768 = -t789 * t825 + t811 * t829;
t716 = qJD(6) * t768 + t755 * t829 + t797 * t825;
t769 = t789 * t829 + t811 * t825;
t727 = -mrSges(7,1) * t768 + mrSges(7,2) * t769;
t787 = qJD(6) + t790;
t731 = -mrSges(7,2) * t787 + mrSges(7,3) * t768;
t748 = qJDD(6) + t756;
t696 = m(7) * t698 + mrSges(7,1) * t748 - mrSges(7,3) * t716 - t727 * t769 + t731 * t787;
t699 = t700 * t825 + t701 * t829;
t715 = -qJD(6) * t769 - t755 * t825 + t797 * t829;
t732 = mrSges(7,1) * t787 - mrSges(7,3) * t769;
t697 = m(7) * t699 - mrSges(7,2) * t748 + mrSges(7,3) * t715 + t727 * t768 - t732 * t787;
t862 = -t825 * t696 + t829 * t697;
t846 = m(6) * t704 - t797 * mrSges(6,2) - t811 * t770 + t862;
t839 = -m(5) * t710 - t797 * mrSges(5,1) - t811 * t776 - t846;
t761 = mrSges(5,1) * t789 - mrSges(5,3) * t790;
t861 = -mrSges(4,1) * t789 - mrSges(4,2) * t790 - t761;
t683 = m(4) * t712 - mrSges(4,1) * t797 - t771 * t811 + (t759 + t861) * t790 + (mrSges(6,3) + t869) * t756 + t839;
t713 = t872 * t725 + t826 * t726;
t773 = -mrSges(4,1) * t811 - mrSges(4,3) * t790;
t775 = -mrSges(6,2) * t811 - mrSges(6,3) * t790;
t806 = t811 * t876;
t845 = pkin(3) * t877 + t797 * qJ(4) + t789 * t760 - t713;
t709 = t806 - t845;
t774 = mrSges(5,1) * t811 + mrSges(5,2) * t790;
t838 = pkin(4) * t788 - qJ(5) * t755 + t845;
t705 = t789 * t874 + t811 * t882 + t838;
t703 = -pkin(5) * t797 - pkin(10) * t877 - t772 * t811 + t806 + ((2 * qJD(5)) + t763) * t789 - t838;
t842 = -m(7) * t703 + mrSges(7,1) * t715 - t716 * mrSges(7,2) + t731 * t768 - t769 * t732;
t837 = -m(6) * t705 + t755 * mrSges(6,3) + t789 * t759 - t842;
t835 = m(5) * t709 - t797 * mrSges(5,3) - t811 * t774 + t837;
t689 = t835 + (t773 - t775) * t811 + (mrSges(4,2) - mrSges(6,1)) * t797 + t861 * t789 + t869 * t755 + m(4) * t713;
t847 = -t683 * t826 + t872 * t689;
t675 = m(3) * t758 - mrSges(3,2) * t819 + mrSges(3,3) * t805 - t798 * t820 + t802 * t849 + t847;
t678 = t872 * t683 + t826 * t689;
t781 = -t823 * t800 - t870;
t799 = -mrSges(3,2) * t820 + mrSges(3,3) * t849;
t677 = m(3) * t781 - t805 * mrSges(3,1) + t804 * mrSges(3,2) + (t798 * t827 - t799 * t830) * t859 + t678;
t711 = (t876 - t871) * t790 + t841;
t686 = t829 * t696 + t825 * t697;
t707 = -pkin(4) * t755 + t790 * t871 + t836;
t844 = -m(6) * t707 - t756 * mrSges(6,1) - t755 * mrSges(6,2) - t789 * t770 - t790 * t775 - t686;
t684 = m(5) * t711 + t755 * mrSges(5,1) - t756 * mrSges(5,3) - t790 * t774 + t789 * t776 + t844;
t833 = -m(4) * t724 - t755 * mrSges(4,1) - t756 * mrSges(4,2) - t789 * t771 - t790 * t773 - t684;
t681 = m(3) * t757 + t819 * mrSges(3,1) - t804 * mrSges(3,3) + t820 * t799 - t802 * t850 + t833;
t666 = t675 * t865 - t677 * t823 + t681 * t864;
t664 = m(2) * t813 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t832 + t666;
t671 = t830 * t675 - t681 * t827;
t670 = m(2) * t814 - mrSges(2,1) * t832 - qJDD(1) * mrSges(2,2) + t671;
t863 = t831 * t664 + t828 * t670;
t665 = t675 * t867 + t824 * t677 + t681 * t866;
t853 = -t789 * t854 - t790 * t880 + t811 * t878;
t852 = -t789 * t855 - t790 * t881 + t811 * t880;
t851 = t789 * t879 - t790 * t855 + t811 * t854;
t848 = -t664 * t828 + t831 * t670;
t717 = Ifges(7,5) * t769 + Ifges(7,6) * t768 + Ifges(7,3) * t787;
t719 = Ifges(7,1) * t769 + Ifges(7,4) * t768 + Ifges(7,5) * t787;
t690 = -mrSges(7,1) * t703 + mrSges(7,3) * t699 + Ifges(7,4) * t716 + Ifges(7,2) * t715 + Ifges(7,6) * t748 - t717 * t769 + t719 * t787;
t718 = Ifges(7,4) * t769 + Ifges(7,2) * t768 + Ifges(7,6) * t787;
t691 = mrSges(7,2) * t703 - mrSges(7,3) * t698 + Ifges(7,1) * t716 + Ifges(7,4) * t715 + Ifges(7,5) * t748 + t717 * t768 - t718 * t787;
t662 = t825 * t690 - t829 * t691 - qJ(5) * t837 - pkin(4) * t844 - mrSges(4,1) * t724 + mrSges(5,2) * t709 - mrSges(5,1) * t711 + mrSges(4,3) * t713 + mrSges(6,3) * t705 - mrSges(6,2) * t707 + pkin(10) * t686 - pkin(3) * t684 + (qJ(5) * t775 + t852) * t811 + (mrSges(6,1) * qJ(5) + t854) * t797 + t853 * t790 - t855 * t756 + t879 * t755;
t685 = -mrSges(6,3) * t756 - t759 * t790 + t846;
t667 = t769 * t718 - t768 * t719 + Ifges(7,3) * t748 + Ifges(7,6) * t715 + Ifges(7,5) * t716 + mrSges(4,2) * t724 + mrSges(5,2) * t710 - mrSges(5,3) * t711 - mrSges(4,3) * t712 - mrSges(6,3) * t704 + mrSges(6,1) * t707 - mrSges(7,2) * t699 + mrSges(7,1) * t698 + pkin(5) * t686 - qJ(5) * t685 - qJ(4) * t684 + t851 * t811 - t880 * t797 + t853 * t789 + t881 * t756 + t855 * t755;
t778 = Ifges(3,3) * t820 + (Ifges(3,5) * t827 + Ifges(3,6) * t830) * t859;
t779 = Ifges(3,6) * t820 + (Ifges(3,4) * t827 + Ifges(3,2) * t830) * t859;
t660 = mrSges(3,2) * t781 - mrSges(3,3) * t757 + Ifges(3,1) * t804 + Ifges(3,4) * t805 + Ifges(3,5) * t819 - pkin(9) * t678 - t826 * t662 + t872 * t667 + t778 * t849 - t820 * t779;
t780 = Ifges(3,5) * t820 + (Ifges(3,1) * t827 + Ifges(3,4) * t830) * t859;
t661 = pkin(10) * t862 + pkin(5) * t842 - qJ(4) * (-t775 * t811 + t835) - pkin(3) * t839 + (-pkin(3) * (-mrSges(5,2) + mrSges(6,3)) - t880) * t756 + (qJ(4) * t761 + t852) * t789 + (mrSges(5,2) * qJ(4) - t854) * t755 - t778 * t850 + (-pkin(3) * (t759 - t761) - t851) * t790 + (mrSges(6,1) * qJ(4) + t878) * t797 + t825 * t691 + t829 * t690 + Ifges(3,6) * t819 + t820 * t780 + Ifges(3,4) * t804 + Ifges(3,2) * t805 - mrSges(3,1) * t781 + mrSges(3,3) * t758 - mrSges(5,3) * t709 + mrSges(5,1) * t710 - mrSges(4,1) * t712 + mrSges(4,2) * t713 - mrSges(6,2) * t704 + mrSges(6,1) * t705 + pkin(4) * t685 - pkin(2) * t678;
t843 = pkin(8) * t671 + t660 * t827 + t661 * t830;
t659 = Ifges(3,5) * t804 + Ifges(3,6) * t805 + Ifges(3,3) * t819 + mrSges(3,1) * t757 - mrSges(3,2) * t758 + t826 * t667 + t872 * t662 + pkin(2) * t833 + pkin(9) * t847 + (t779 * t827 - t780 * t830) * t859;
t658 = -mrSges(2,2) * g(3) - mrSges(2,3) * t813 + Ifges(2,5) * qJDD(1) - t832 * Ifges(2,6) + t830 * t660 - t827 * t661 + (-t665 * t823 - t666 * t824) * pkin(8);
t657 = mrSges(2,1) * g(3) + mrSges(2,3) * t814 + Ifges(2,5) * t832 + Ifges(2,6) * qJDD(1) - pkin(1) * t665 - t659 * t823 + t843 * t824;
t1 = [-m(1) * g(1) + t848; -m(1) * g(2) + t863; (-m(1) - m(2)) * g(3) + t665; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t863 - t828 * t657 + t831 * t658; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t848 + t831 * t657 + t828 * t658; -mrSges(1,1) * g(2) + mrSges(2,1) * t813 + mrSges(1,2) * g(1) - mrSges(2,2) * t814 + Ifges(2,3) * qJDD(1) + pkin(1) * t666 + t659 * t824 + t843 * t823;];
tauB  = t1;

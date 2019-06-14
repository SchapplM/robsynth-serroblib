% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:20:25
% EndTime: 2019-05-05 19:20:34
% DurationCPUTime: 8.41s
% Computational Cost: add. (133177->339), mult. (292003->415), div. (0->0), fcn. (196636->10), ass. (0->135)
t865 = sin(qJ(1));
t869 = cos(qJ(1));
t843 = t865 * g(1) - t869 * g(2);
t871 = qJD(1) ^ 2;
t881 = -t871 * qJ(2) + qJDD(2) - t843;
t900 = -pkin(1) - pkin(7);
t814 = qJDD(1) * t900 + t881;
t864 = sin(qJ(3));
t868 = cos(qJ(3));
t803 = t864 * g(3) + t868 * t814;
t892 = qJD(1) * qJD(3);
t890 = t864 * t892;
t838 = qJDD(1) * t868 - t890;
t776 = (-t838 - t890) * qJ(4) + (-t864 * t868 * t871 + qJDD(3)) * pkin(3) + t803;
t804 = -g(3) * t868 + t864 * t814;
t837 = -qJDD(1) * t864 - t868 * t892;
t894 = qJD(1) * t868;
t841 = qJD(3) * pkin(3) - qJ(4) * t894;
t857 = t864 ^ 2;
t777 = -pkin(3) * t857 * t871 + qJ(4) * t837 - qJD(3) * t841 + t804;
t860 = sin(pkin(10));
t861 = cos(pkin(10));
t895 = qJD(1) * t864;
t825 = -t860 * t895 + t861 * t894;
t757 = -0.2e1 * qJD(4) * t825 + t776 * t861 - t860 * t777;
t844 = -t869 * g(1) - t865 * g(2);
t882 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t844;
t899 = mrSges(2,1) - mrSges(3,2);
t898 = Ifges(2,5) - Ifges(3,4);
t897 = -Ifges(2,6) + Ifges(3,5);
t824 = -t860 * t894 - t861 * t895;
t758 = 0.2e1 * qJD(4) * t824 + t860 * t776 + t861 * t777;
t792 = -mrSges(5,1) * t824 + mrSges(5,2) * t825;
t800 = t837 * t861 - t860 * t838;
t813 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t825;
t794 = -pkin(4) * t824 - pkin(8) * t825;
t870 = qJD(3) ^ 2;
t746 = -pkin(4) * t870 + qJDD(3) * pkin(8) + t794 * t824 + t758;
t779 = -t837 * pkin(3) + qJDD(4) + t841 * t894 + (-qJ(4) * t857 + t900) * t871 + t882;
t801 = t837 * t860 + t838 * t861;
t755 = (-qJD(3) * t824 - t801) * pkin(8) + (qJD(3) * t825 - t800) * pkin(4) + t779;
t863 = sin(qJ(5));
t867 = cos(qJ(5));
t741 = -t863 * t746 + t867 * t755;
t806 = qJD(3) * t867 - t825 * t863;
t772 = qJD(5) * t806 + qJDD(3) * t863 + t801 * t867;
t799 = qJDD(5) - t800;
t807 = qJD(3) * t863 + t825 * t867;
t822 = qJD(5) - t824;
t739 = (t806 * t822 - t772) * pkin(9) + (t806 * t807 + t799) * pkin(5) + t741;
t742 = t867 * t746 + t863 * t755;
t771 = -qJD(5) * t807 + qJDD(3) * t867 - t801 * t863;
t787 = pkin(5) * t822 - pkin(9) * t807;
t805 = t806 ^ 2;
t740 = -pkin(5) * t805 + t771 * pkin(9) - t787 * t822 + t742;
t862 = sin(qJ(6));
t866 = cos(qJ(6));
t737 = t739 * t866 - t740 * t862;
t780 = t806 * t866 - t807 * t862;
t751 = qJD(6) * t780 + t771 * t862 + t772 * t866;
t781 = t806 * t862 + t807 * t866;
t763 = -mrSges(7,1) * t780 + mrSges(7,2) * t781;
t818 = qJD(6) + t822;
t764 = -mrSges(7,2) * t818 + mrSges(7,3) * t780;
t795 = qJDD(6) + t799;
t733 = m(7) * t737 + mrSges(7,1) * t795 - t751 * mrSges(7,3) - t763 * t781 + t764 * t818;
t738 = t739 * t862 + t740 * t866;
t750 = -qJD(6) * t781 + t771 * t866 - t772 * t862;
t765 = mrSges(7,1) * t818 - mrSges(7,3) * t781;
t734 = m(7) * t738 - mrSges(7,2) * t795 + t750 * mrSges(7,3) + t763 * t780 - t765 * t818;
t725 = t866 * t733 + t862 * t734;
t783 = -mrSges(6,1) * t806 + mrSges(6,2) * t807;
t785 = -mrSges(6,2) * t822 + mrSges(6,3) * t806;
t723 = m(6) * t741 + mrSges(6,1) * t799 - t772 * mrSges(6,3) - t783 * t807 + t785 * t822 + t725;
t786 = mrSges(6,1) * t822 - mrSges(6,3) * t807;
t885 = -t733 * t862 + t866 * t734;
t724 = m(6) * t742 - mrSges(6,2) * t799 + t771 * mrSges(6,3) + t783 * t806 - t786 * t822 + t885;
t886 = -t723 * t863 + t867 * t724;
t716 = m(5) * t758 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t800 - qJD(3) * t813 + t792 * t824 + t886;
t812 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t824;
t745 = -qJDD(3) * pkin(4) - pkin(8) * t870 + t825 * t794 - t757;
t743 = -t771 * pkin(5) - pkin(9) * t805 + t787 * t807 + t745;
t879 = m(7) * t743 - t750 * mrSges(7,1) + t751 * mrSges(7,2) - t780 * t764 + t765 * t781;
t875 = -m(6) * t745 + t771 * mrSges(6,1) - t772 * mrSges(6,2) + t806 * t785 - t786 * t807 - t879;
t729 = m(5) * t757 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t801 + qJD(3) * t812 - t792 * t825 + t875;
t705 = t860 * t716 + t861 * t729;
t836 = (mrSges(4,1) * t864 + mrSges(4,2) * t868) * qJD(1);
t840 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t895;
t702 = m(4) * t803 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t838 + qJD(3) * t840 - t836 * t894 + t705;
t842 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t894;
t887 = t861 * t716 - t729 * t860;
t703 = m(4) * t804 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t837 - qJD(3) * t842 - t836 * t895 + t887;
t699 = t868 * t702 + t864 * t703;
t823 = -qJDD(1) * pkin(1) + t881;
t880 = -m(3) * t823 + t871 * mrSges(3,3) - t699;
t695 = m(2) * t843 - t871 * mrSges(2,2) + qJDD(1) * t899 + t880;
t817 = t871 * pkin(1) - t882;
t719 = t867 * t723 + t863 * t724;
t717 = m(5) * t779 - mrSges(5,1) * t800 + t801 * mrSges(5,2) - t812 * t824 + t825 * t813 + t719;
t811 = t871 * t900 + t882;
t877 = -m(4) * t811 + mrSges(4,1) * t837 - t838 * mrSges(4,2) - t840 * t895 - t842 * t894 - t717;
t874 = -m(3) * t817 + t871 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t877;
t712 = m(2) * t844 - mrSges(2,1) * t871 - qJDD(1) * mrSges(2,2) + t874;
t896 = t869 * t695 + t865 * t712;
t889 = -t695 * t865 + t869 * t712;
t888 = -t864 * t702 + t868 * t703;
t760 = Ifges(7,4) * t781 + Ifges(7,2) * t780 + Ifges(7,6) * t818;
t761 = Ifges(7,1) * t781 + Ifges(7,4) * t780 + Ifges(7,5) * t818;
t878 = -mrSges(7,1) * t737 + mrSges(7,2) * t738 - Ifges(7,5) * t751 - Ifges(7,6) * t750 - Ifges(7,3) * t795 - t781 * t760 + t780 * t761;
t759 = Ifges(7,5) * t781 + Ifges(7,6) * t780 + Ifges(7,3) * t818;
t726 = -mrSges(7,1) * t743 + mrSges(7,3) * t738 + Ifges(7,4) * t751 + Ifges(7,2) * t750 + Ifges(7,6) * t795 - t759 * t781 + t761 * t818;
t727 = mrSges(7,2) * t743 - mrSges(7,3) * t737 + Ifges(7,1) * t751 + Ifges(7,4) * t750 + Ifges(7,5) * t795 + t759 * t780 - t760 * t818;
t766 = Ifges(6,5) * t807 + Ifges(6,6) * t806 + Ifges(6,3) * t822;
t768 = Ifges(6,1) * t807 + Ifges(6,4) * t806 + Ifges(6,5) * t822;
t707 = -mrSges(6,1) * t745 + mrSges(6,3) * t742 + Ifges(6,4) * t772 + Ifges(6,2) * t771 + Ifges(6,6) * t799 - pkin(5) * t879 + pkin(9) * t885 + t866 * t726 + t862 * t727 - t807 * t766 + t822 * t768;
t767 = Ifges(6,4) * t807 + Ifges(6,2) * t806 + Ifges(6,6) * t822;
t709 = mrSges(6,2) * t745 - mrSges(6,3) * t741 + Ifges(6,1) * t772 + Ifges(6,4) * t771 + Ifges(6,5) * t799 - pkin(9) * t725 - t726 * t862 + t727 * t866 + t766 * t806 - t767 * t822;
t788 = Ifges(5,5) * t825 + Ifges(5,6) * t824 + Ifges(5,3) * qJD(3);
t789 = Ifges(5,4) * t825 + Ifges(5,2) * t824 + Ifges(5,6) * qJD(3);
t693 = mrSges(5,2) * t779 - mrSges(5,3) * t757 + Ifges(5,1) * t801 + Ifges(5,4) * t800 + Ifges(5,5) * qJDD(3) - pkin(8) * t719 - qJD(3) * t789 - t707 * t863 + t709 * t867 + t788 * t824;
t790 = Ifges(5,1) * t825 + Ifges(5,4) * t824 + Ifges(5,5) * qJD(3);
t872 = mrSges(6,1) * t741 - mrSges(6,2) * t742 + Ifges(6,5) * t772 + Ifges(6,6) * t771 + Ifges(6,3) * t799 + pkin(5) * t725 + t807 * t767 - t806 * t768 - t878;
t700 = -mrSges(5,1) * t779 + mrSges(5,3) * t758 + Ifges(5,4) * t801 + Ifges(5,2) * t800 + Ifges(5,6) * qJDD(3) - pkin(4) * t719 + qJD(3) * t790 - t825 * t788 - t872;
t826 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t868 - Ifges(4,6) * t864) * qJD(1);
t828 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t868 - Ifges(4,4) * t864) * qJD(1);
t690 = -mrSges(4,1) * t811 + mrSges(4,3) * t804 + Ifges(4,4) * t838 + Ifges(4,2) * t837 + Ifges(4,6) * qJDD(3) - pkin(3) * t717 + qJ(4) * t887 + qJD(3) * t828 + t860 * t693 + t861 * t700 - t826 * t894;
t827 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t868 - Ifges(4,2) * t864) * qJD(1);
t692 = mrSges(4,2) * t811 - mrSges(4,3) * t803 + Ifges(4,1) * t838 + Ifges(4,4) * t837 + Ifges(4,5) * qJDD(3) - qJ(4) * t705 - qJD(3) * t827 + t693 * t861 - t700 * t860 - t826 * t895;
t697 = qJDD(1) * mrSges(3,2) - t880;
t876 = mrSges(2,1) * t843 - mrSges(2,2) * t844 + mrSges(3,2) * t823 - mrSges(3,3) * t817 - pkin(1) * t697 - pkin(7) * t699 + qJ(2) * t874 - t690 * t864 + t868 * t692 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t873 = mrSges(5,1) * t757 - mrSges(4,2) * t804 - mrSges(5,2) * t758 + Ifges(5,5) * t801 + Ifges(5,6) * t800 + pkin(3) * t705 + pkin(4) * t875 + pkin(8) * t886 + t867 * t707 + t863 * t709 + t825 * t789 - t824 * t790 + mrSges(4,1) * t803 + t828 * t895 + t827 * t894 + Ifges(4,6) * t837 + Ifges(4,5) * t838 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t698 = -m(3) * g(3) + t888;
t689 = t873 + pkin(2) * t699 - qJ(2) * t698 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t898 * qJDD(1) + t897 * t871 - mrSges(2,3) * t843 + mrSges(3,1) * t823;
t688 = -mrSges(3,1) * t817 + mrSges(2,3) * t844 - pkin(1) * t698 - pkin(2) * t877 - pkin(7) * t888 + g(3) * t899 - qJDD(1) * t897 - t868 * t690 - t864 * t692 + t871 * t898;
t1 = [-m(1) * g(1) + t889; -m(1) * g(2) + t896; (-m(1) - m(2) - m(3)) * g(3) + t888; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t896 - t865 * t688 + t869 * t689; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t889 + t869 * t688 + t865 * t689; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t876; t876; t697; t873; t717; t872; -t878;];
tauJB  = t1;

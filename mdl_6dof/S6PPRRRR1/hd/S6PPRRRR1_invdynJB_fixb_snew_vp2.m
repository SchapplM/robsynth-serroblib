% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:40:04
% EndTime: 2019-05-04 20:40:25
% DurationCPUTime: 20.93s
% Computational Cost: add. (398919->296), mult. (735688->391), div. (0->0), fcn. (578085->16), ass. (0->138)
t786 = sin(pkin(12));
t790 = cos(pkin(12));
t774 = -g(1) * t790 - g(2) * t786;
t785 = sin(pkin(13));
t789 = cos(pkin(13));
t773 = g(1) * t786 - g(2) * t790;
t784 = -g(3) + qJDD(1);
t788 = sin(pkin(6));
t792 = cos(pkin(6));
t812 = t773 * t792 + t784 * t788;
t747 = -t785 * t774 + t789 * t812;
t748 = t789 * t774 + t785 * t812;
t760 = -t773 * t788 + t784 * t792 + qJDD(2);
t800 = cos(qJ(3));
t791 = cos(pkin(7));
t796 = sin(qJ(3));
t825 = t791 * t796;
t787 = sin(pkin(7));
t826 = t787 * t796;
t723 = t747 * t825 + t800 * t748 + t760 * t826;
t801 = qJD(3) ^ 2;
t721 = -pkin(3) * t801 + qJDD(3) * pkin(9) + t723;
t733 = -t747 * t787 + t760 * t791;
t795 = sin(qJ(4));
t799 = cos(qJ(4));
t716 = -t795 * t721 + t799 * t733;
t821 = qJD(3) * qJD(4);
t820 = t799 * t821;
t771 = qJDD(3) * t795 + t820;
t714 = (-t771 + t820) * pkin(10) + (t795 * t799 * t801 + qJDD(4)) * pkin(4) + t716;
t717 = t799 * t721 + t795 * t733;
t772 = qJDD(3) * t799 - t795 * t821;
t823 = qJD(3) * t795;
t777 = qJD(4) * pkin(4) - pkin(10) * t823;
t783 = t799 ^ 2;
t715 = -pkin(4) * t783 * t801 + pkin(10) * t772 - qJD(4) * t777 + t717;
t794 = sin(qJ(5));
t798 = cos(qJ(5));
t710 = t794 * t714 + t798 * t715;
t766 = (t794 * t799 + t795 * t798) * qJD(3);
t740 = -qJD(5) * t766 - t771 * t794 + t772 * t798;
t765 = (t794 * t795 - t798 * t799) * qJD(3);
t753 = mrSges(6,1) * t765 + mrSges(6,2) * t766;
t782 = qJD(4) + qJD(5);
t759 = mrSges(6,1) * t782 - mrSges(6,3) * t766;
t781 = qJDD(4) + qJDD(5);
t754 = pkin(5) * t765 - pkin(11) * t766;
t780 = t782 ^ 2;
t707 = -pkin(5) * t780 + pkin(11) * t781 - t754 * t765 + t710;
t722 = -t796 * t748 + (t747 * t791 + t760 * t787) * t800;
t807 = -qJDD(3) * pkin(3) - t722;
t718 = -t772 * pkin(4) + t777 * t823 + (-pkin(10) * t783 - pkin(9)) * t801 + t807;
t741 = -qJD(5) * t765 + t771 * t798 + t772 * t794;
t711 = (t765 * t782 - t741) * pkin(11) + (t766 * t782 - t740) * pkin(5) + t718;
t793 = sin(qJ(6));
t797 = cos(qJ(6));
t704 = -t707 * t793 + t711 * t797;
t756 = -t766 * t793 + t782 * t797;
t726 = qJD(6) * t756 + t741 * t797 + t781 * t793;
t757 = t766 * t797 + t782 * t793;
t734 = -mrSges(7,1) * t756 + mrSges(7,2) * t757;
t739 = qJDD(6) - t740;
t761 = qJD(6) + t765;
t742 = -mrSges(7,2) * t761 + mrSges(7,3) * t756;
t700 = m(7) * t704 + mrSges(7,1) * t739 - mrSges(7,3) * t726 - t734 * t757 + t742 * t761;
t705 = t707 * t797 + t711 * t793;
t725 = -qJD(6) * t757 - t741 * t793 + t781 * t797;
t743 = mrSges(7,1) * t761 - mrSges(7,3) * t757;
t701 = m(7) * t705 - mrSges(7,2) * t739 + mrSges(7,3) * t725 + t734 * t756 - t743 * t761;
t816 = -t700 * t793 + t797 * t701;
t687 = m(6) * t710 - mrSges(6,2) * t781 + mrSges(6,3) * t740 - t753 * t765 - t759 * t782 + t816;
t709 = t714 * t798 - t715 * t794;
t758 = -mrSges(6,2) * t782 - mrSges(6,3) * t765;
t706 = -pkin(5) * t781 - pkin(11) * t780 + t754 * t766 - t709;
t808 = -m(7) * t706 + t725 * mrSges(7,1) - mrSges(7,2) * t726 + t756 * t742 - t743 * t757;
t696 = m(6) * t709 + mrSges(6,1) * t781 - mrSges(6,3) * t741 - t753 * t766 + t758 * t782 + t808;
t681 = t794 * t687 + t798 * t696;
t770 = (-mrSges(5,1) * t799 + mrSges(5,2) * t795) * qJD(3);
t822 = qJD(3) * t799;
t776 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t822;
t679 = m(5) * t716 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t771 + qJD(4) * t776 - t770 * t823 + t681;
t775 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t823;
t817 = t798 * t687 - t696 * t794;
t680 = m(5) * t717 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t772 - qJD(4) * t775 + t770 * t822 + t817;
t818 = -t679 * t795 + t799 * t680;
t670 = m(4) * t723 - mrSges(4,1) * t801 - qJDD(3) * mrSges(4,2) + t818;
t673 = t799 * t679 + t795 * t680;
t672 = m(4) * t733 + t673;
t720 = -t801 * pkin(9) + t807;
t689 = t797 * t700 + t793 * t701;
t806 = m(6) * t718 - t740 * mrSges(6,1) + mrSges(6,2) * t741 + t765 * t758 + t759 * t766 + t689;
t803 = -m(5) * t720 + t772 * mrSges(5,1) - mrSges(5,2) * t771 - t775 * t823 + t776 * t822 - t806;
t684 = m(4) * t722 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t801 + t803;
t827 = t684 * t800;
t659 = t670 * t825 - t672 * t787 + t791 * t827;
t655 = m(3) * t747 + t659;
t666 = t800 * t670 - t684 * t796;
t665 = m(3) * t748 + t666;
t831 = t655 * t789 + t665 * t785;
t763 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t795 + Ifges(5,2) * t799) * qJD(3);
t764 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t795 + Ifges(5,4) * t799) * qJD(3);
t727 = Ifges(7,5) * t757 + Ifges(7,6) * t756 + Ifges(7,3) * t761;
t729 = Ifges(7,1) * t757 + Ifges(7,4) * t756 + Ifges(7,5) * t761;
t693 = -mrSges(7,1) * t706 + mrSges(7,3) * t705 + Ifges(7,4) * t726 + Ifges(7,2) * t725 + Ifges(7,6) * t739 - t727 * t757 + t729 * t761;
t728 = Ifges(7,4) * t757 + Ifges(7,2) * t756 + Ifges(7,6) * t761;
t694 = mrSges(7,2) * t706 - mrSges(7,3) * t704 + Ifges(7,1) * t726 + Ifges(7,4) * t725 + Ifges(7,5) * t739 + t727 * t756 - t728 * t761;
t750 = Ifges(6,4) * t766 - Ifges(6,2) * t765 + Ifges(6,6) * t782;
t751 = Ifges(6,1) * t766 - Ifges(6,4) * t765 + Ifges(6,5) * t782;
t805 = -mrSges(6,1) * t709 + mrSges(6,2) * t710 - Ifges(6,5) * t741 - Ifges(6,6) * t740 - Ifges(6,3) * t781 - pkin(5) * t808 - pkin(11) * t816 - t797 * t693 - t793 * t694 - t766 * t750 - t765 * t751;
t830 = mrSges(5,1) * t716 - mrSges(5,2) * t717 + Ifges(5,5) * t771 + Ifges(5,6) * t772 + Ifges(5,3) * qJDD(4) + pkin(4) * t681 + (t763 * t795 - t764 * t799) * qJD(3) - t805;
t658 = t670 * t826 + t791 * t672 + t787 * t827;
t657 = m(3) * t760 + t658;
t645 = -t657 * t788 + t792 * t831;
t643 = m(2) * t773 + t645;
t651 = -t655 * t785 + t789 * t665;
t650 = m(2) * t774 + t651;
t824 = t790 * t643 + t786 * t650;
t644 = t792 * t657 + t788 * t831;
t819 = -t643 * t786 + t790 * t650;
t815 = m(2) * t784 + t644;
t749 = Ifges(6,5) * t766 - Ifges(6,6) * t765 + Ifges(6,3) * t782;
t674 = mrSges(6,2) * t718 - mrSges(6,3) * t709 + Ifges(6,1) * t741 + Ifges(6,4) * t740 + Ifges(6,5) * t781 - pkin(11) * t689 - t693 * t793 + t694 * t797 - t749 * t765 - t750 * t782;
t804 = mrSges(7,1) * t704 - mrSges(7,2) * t705 + Ifges(7,5) * t726 + Ifges(7,6) * t725 + Ifges(7,3) * t739 + t728 * t757 - t729 * t756;
t675 = -mrSges(6,1) * t718 + mrSges(6,3) * t710 + Ifges(6,4) * t741 + Ifges(6,2) * t740 + Ifges(6,6) * t781 - pkin(5) * t689 - t749 * t766 + t751 * t782 - t804;
t762 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t795 + Ifges(5,6) * t799) * qJD(3);
t660 = -mrSges(5,1) * t720 + mrSges(5,3) * t717 + Ifges(5,4) * t771 + Ifges(5,2) * t772 + Ifges(5,6) * qJDD(4) - pkin(4) * t806 + pkin(10) * t817 + qJD(4) * t764 + t794 * t674 + t798 * t675 - t762 * t823;
t661 = mrSges(5,2) * t720 - mrSges(5,3) * t716 + Ifges(5,1) * t771 + Ifges(5,4) * t772 + Ifges(5,5) * qJDD(4) - pkin(10) * t681 - qJD(4) * t763 + t674 * t798 - t675 * t794 + t762 * t822;
t647 = mrSges(4,2) * t733 - mrSges(4,3) * t722 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t801 - pkin(9) * t673 - t660 * t795 + t661 * t799;
t652 = -mrSges(4,1) * t733 + mrSges(4,3) * t723 + t801 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t673 - t830;
t810 = pkin(8) * t666 + t647 * t796 + t652 * t800;
t646 = mrSges(4,1) * t722 - mrSges(4,2) * t723 + Ifges(4,3) * qJDD(3) + pkin(3) * t803 + pkin(9) * t818 + t799 * t660 + t795 * t661;
t640 = -mrSges(3,1) * t760 + mrSges(3,3) * t748 - pkin(2) * t658 - t787 * t646 + t791 * t810;
t641 = mrSges(3,2) * t760 - mrSges(3,3) * t747 + t800 * t647 - t796 * t652 + (-t658 * t787 - t659 * t791) * pkin(8);
t809 = qJ(2) * t651 + t640 * t789 + t641 * t785;
t639 = mrSges(3,1) * t747 - mrSges(3,2) * t748 + pkin(2) * t659 + t791 * t646 + t787 * t810;
t638 = mrSges(2,2) * t784 - mrSges(2,3) * t773 - t785 * t640 + t789 * t641 + (-t644 * t788 - t645 * t792) * qJ(2);
t637 = -mrSges(2,1) * t784 + mrSges(2,3) * t774 - pkin(1) * t644 - t788 * t639 + t792 * t809;
t1 = [-m(1) * g(1) + t819; -m(1) * g(2) + t824; -m(1) * g(3) + t815; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t824 - t786 * t637 + t790 * t638; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t819 + t790 * t637 + t786 * t638; -mrSges(1,1) * g(2) + mrSges(2,1) * t773 + mrSges(1,2) * g(1) - mrSges(2,2) * t774 + pkin(1) * t645 + t792 * t639 + t788 * t809; t815; t657; t646; t830; -t805; t804;];
tauJB  = t1;

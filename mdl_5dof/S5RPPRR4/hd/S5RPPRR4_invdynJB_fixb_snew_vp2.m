% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:43
% EndTime: 2022-01-23 09:16:50
% DurationCPUTime: 7.30s
% Computational Cost: add. (73722->278), mult. (201875->377), div. (0->0), fcn. (138541->10), ass. (0->134)
t787 = sin(qJ(1));
t790 = cos(qJ(1));
t766 = -t790 * g(1) - t787 * g(2);
t791 = qJD(1) ^ 2;
t835 = -t791 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t766;
t765 = t787 * g(1) - t790 * g(2);
t800 = -t791 * qJ(2) + qJDD(2) - t765;
t782 = sin(pkin(8));
t784 = cos(pkin(8));
t808 = -pkin(2) * t784 - qJ(3) * t782;
t826 = t782 * qJD(1);
t834 = (-pkin(1) + t808) * qJDD(1) + t800 - 0.2e1 * qJD(3) * t826;
t738 = -t784 * g(3) - t835 * t782;
t833 = mrSges(3,2) * t782;
t779 = t782 ^ 2;
t832 = t779 * t791;
t781 = sin(pkin(9));
t831 = t781 * t782;
t783 = cos(pkin(9));
t830 = t782 * t783;
t739 = -t782 * g(3) + t835 * t784;
t758 = (-mrSges(3,1) * t784 + t833) * qJD(1);
t757 = t808 * qJD(1);
t825 = t784 * qJD(1);
t728 = t757 * t825 + t739;
t805 = -pkin(3) * t784 - pkin(6) * t830;
t828 = t834 * t783;
t707 = t805 * qJDD(1) + (-t728 + (-pkin(3) * t779 * t783 + pkin(6) * t782 * t784) * t791) * t781 + t828;
t716 = t783 * t728 + t834 * t781;
t756 = t805 * qJD(1);
t823 = qJDD(1) * t782;
t819 = t781 * t823;
t821 = t781 ^ 2 * t832;
t708 = -pkin(3) * t821 - pkin(6) * t819 + t756 * t825 + t716;
t786 = sin(qJ(4));
t789 = cos(qJ(4));
t693 = t789 * t707 - t786 * t708;
t802 = (-t781 * t789 - t783 * t786) * t782;
t746 = qJD(1) * t802;
t801 = (-t781 * t786 + t783 * t789) * t782;
t732 = t746 * qJD(4) + qJDD(1) * t801;
t747 = qJD(1) * t801;
t822 = t784 * qJDD(1);
t769 = qJDD(4) - t822;
t770 = qJD(4) - t825;
t691 = (t746 * t770 - t732) * pkin(7) + (t746 * t747 + t769) * pkin(4) + t693;
t694 = t786 * t707 + t789 * t708;
t731 = -t747 * qJD(4) + qJDD(1) * t802;
t737 = t770 * pkin(4) - t747 * pkin(7);
t745 = t746 ^ 2;
t692 = -t745 * pkin(4) + t731 * pkin(7) - t770 * t737 + t694;
t785 = sin(qJ(5));
t788 = cos(qJ(5));
t689 = t788 * t691 - t785 * t692;
t725 = t788 * t746 - t785 * t747;
t703 = t725 * qJD(5) + t785 * t731 + t788 * t732;
t726 = t785 * t746 + t788 * t747;
t714 = -t725 * mrSges(6,1) + t726 * mrSges(6,2);
t768 = qJD(5) + t770;
t718 = -t768 * mrSges(6,2) + t725 * mrSges(6,3);
t764 = qJDD(5) + t769;
t685 = m(6) * t689 + t764 * mrSges(6,1) - t703 * mrSges(6,3) - t726 * t714 + t768 * t718;
t690 = t785 * t691 + t788 * t692;
t702 = -t726 * qJD(5) + t788 * t731 - t785 * t732;
t719 = t768 * mrSges(6,1) - t726 * mrSges(6,3);
t686 = m(6) * t690 - t764 * mrSges(6,2) + t702 * mrSges(6,3) + t725 * t714 - t768 * t719;
t677 = t788 * t685 + t785 * t686;
t729 = -t746 * mrSges(5,1) + t747 * mrSges(5,2);
t733 = -t770 * mrSges(5,2) + t746 * mrSges(5,3);
t675 = m(5) * t693 + t769 * mrSges(5,1) - t732 * mrSges(5,3) - t747 * t729 + t770 * t733 + t677;
t734 = t770 * mrSges(5,1) - t747 * mrSges(5,3);
t814 = -t785 * t685 + t788 * t686;
t676 = m(5) * t694 - t769 * mrSges(5,2) + t731 * mrSges(5,3) + t746 * t729 - t770 * t734 + t814;
t671 = t789 * t675 + t786 * t676;
t715 = -t781 * t728 + t828;
t812 = mrSges(4,1) * t781 + mrSges(4,2) * t783;
t750 = t812 * t826;
t803 = mrSges(4,2) * t784 - mrSges(4,3) * t831;
t753 = t803 * qJD(1);
t804 = -mrSges(4,1) * t784 - mrSges(4,3) * t830;
t669 = m(4) * t715 + t804 * qJDD(1) + (-t750 * t830 - t753 * t784) * qJD(1) + t671;
t754 = t804 * qJD(1);
t815 = -t786 * t675 + t789 * t676;
t670 = m(4) * t716 + t803 * qJDD(1) + (-t750 * t831 + t754 * t784) * qJD(1) + t815;
t816 = -t781 * t669 + t783 * t670;
t662 = m(3) * t739 + (qJDD(1) * mrSges(3,3) + qJD(1) * t758) * t784 + t816;
t727 = t757 * t826 + qJDD(3) - t738;
t717 = t783 * t756 * t826 + pkin(3) * t819 - pkin(6) * t821 + t727;
t696 = -t731 * pkin(4) - t745 * pkin(7) + t747 * t737 + t717;
t807 = m(6) * t696 - t702 * mrSges(6,1) + t703 * mrSges(6,2) - t725 * t718 + t726 * t719;
t794 = m(5) * t717 - t731 * mrSges(5,1) + t732 * mrSges(5,2) - t746 * t733 + t747 * t734 + t807;
t793 = m(4) * t727 + t794;
t806 = t753 * t781 + t754 * t783;
t681 = -t793 + m(3) * t738 + ((-mrSges(3,3) - t812) * qJDD(1) + (-t758 - t806) * qJD(1)) * t782;
t817 = t784 * t662 - t782 * t681;
t655 = m(2) * t766 - t791 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t817;
t665 = t783 * t669 + t781 * t670;
t752 = -qJDD(1) * pkin(1) + t800;
t795 = -m(3) * t752 + mrSges(3,1) * t822 - t665 + (t784 ^ 2 * t791 + t832) * mrSges(3,3);
t659 = m(2) * t765 - t791 * mrSges(2,2) + (mrSges(2,1) - t833) * qJDD(1) + t795;
t829 = t787 * t655 + t790 * t659;
t657 = t782 * t662 + t784 * t681;
t818 = t790 * t655 - t787 * t659;
t811 = Ifges(3,1) * t782 + Ifges(3,4) * t784;
t810 = Ifges(3,5) * t782 + Ifges(3,6) * t784;
t809 = Ifges(4,5) * t783 - Ifges(4,6) * t781;
t709 = Ifges(6,5) * t726 + Ifges(6,6) * t725 + Ifges(6,3) * t768;
t711 = Ifges(6,1) * t726 + Ifges(6,4) * t725 + Ifges(6,5) * t768;
t678 = -mrSges(6,1) * t696 + mrSges(6,3) * t690 + Ifges(6,4) * t703 + Ifges(6,2) * t702 + Ifges(6,6) * t764 - t726 * t709 + t768 * t711;
t710 = Ifges(6,4) * t726 + Ifges(6,2) * t725 + Ifges(6,6) * t768;
t679 = mrSges(6,2) * t696 - mrSges(6,3) * t689 + Ifges(6,1) * t703 + Ifges(6,4) * t702 + Ifges(6,5) * t764 + t725 * t709 - t768 * t710;
t720 = Ifges(5,5) * t747 + Ifges(5,6) * t746 + Ifges(5,3) * t770;
t722 = Ifges(5,1) * t747 + Ifges(5,4) * t746 + Ifges(5,5) * t770;
t666 = -mrSges(5,1) * t717 + mrSges(5,3) * t694 + Ifges(5,4) * t732 + Ifges(5,2) * t731 + Ifges(5,6) * t769 - pkin(4) * t807 + pkin(7) * t814 + t788 * t678 + t785 * t679 - t747 * t720 + t770 * t722;
t721 = Ifges(5,4) * t747 + Ifges(5,2) * t746 + Ifges(5,6) * t770;
t667 = mrSges(5,2) * t717 - mrSges(5,3) * t693 + Ifges(5,1) * t732 + Ifges(5,4) * t731 + Ifges(5,5) * t769 - pkin(7) * t677 - t785 * t678 + t788 * t679 + t746 * t720 - t770 * t721;
t741 = (-Ifges(4,3) * t784 + t809 * t782) * qJD(1);
t798 = -Ifges(4,5) * t784 + (Ifges(4,1) * t783 - Ifges(4,4) * t781) * t782;
t743 = t798 * qJD(1);
t797 = -Ifges(4,6) * t784 + (Ifges(4,4) * t783 - Ifges(4,2) * t781) * t782;
t651 = -mrSges(4,1) * t727 + mrSges(4,3) * t716 + t786 * t667 + t789 * t666 - pkin(3) * t794 + pkin(6) * t815 + (-t741 * t830 - t784 * t743) * qJD(1) + t797 * qJDD(1);
t742 = t797 * qJD(1);
t652 = mrSges(4,2) * t727 - mrSges(4,3) * t715 - pkin(6) * t671 - t786 * t666 + t789 * t667 + (-t741 * t831 + t742 * t784) * qJD(1) + t798 * qJDD(1);
t759 = t810 * qJD(1);
t648 = mrSges(3,2) * t752 - mrSges(3,3) * t738 - qJ(3) * t665 + t811 * qJDD(1) - t781 * t651 + t783 * t652 + t759 * t825;
t796 = -mrSges(6,1) * t689 + mrSges(6,2) * t690 - Ifges(6,5) * t703 - Ifges(6,6) * t702 - Ifges(6,3) * t764 - t726 * t710 + t725 * t711;
t792 = mrSges(5,1) * t693 - mrSges(5,2) * t694 + Ifges(5,5) * t732 + Ifges(5,6) * t731 + Ifges(5,3) * t769 + pkin(4) * t677 + t747 * t721 - t746 * t722 - t796;
t650 = (Ifges(3,2) + Ifges(4,3)) * t822 - t792 - pkin(3) * t671 + ((Ifges(3,4) - t809) * qJDD(1) + (-t742 * t783 - t743 * t781 - t759) * qJD(1)) * t782 - mrSges(4,1) * t715 + mrSges(4,2) * t716 - pkin(2) * t665 + mrSges(3,3) * t739 - mrSges(3,1) * t752;
t664 = mrSges(3,2) * t823 - t795;
t799 = mrSges(2,1) * t765 - mrSges(2,2) * t766 + Ifges(2,3) * qJDD(1) - pkin(1) * t664 + qJ(2) * t817 + t782 * t648 + t784 * t650;
t687 = (t806 * qJD(1) + t812 * qJDD(1)) * t782 + t793;
t646 = mrSges(2,1) * g(3) + mrSges(2,3) * t766 - mrSges(3,1) * t738 + mrSges(3,2) * t739 - t781 * t652 - t783 * t651 + pkin(2) * t687 - qJ(3) * t816 - pkin(1) * t657 + (Ifges(2,6) - t810) * qJDD(1) + (Ifges(2,5) - t782 * (Ifges(3,4) * t782 + Ifges(3,2) * t784) + t784 * t811) * t791;
t645 = -mrSges(2,2) * g(3) - mrSges(2,3) * t765 + Ifges(2,5) * qJDD(1) - t791 * Ifges(2,6) - qJ(2) * t657 + t784 * t648 - t782 * t650;
t1 = [-m(1) * g(1) + t818; -m(1) * g(2) + t829; (-m(1) - m(2)) * g(3) + t657; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t829 + t790 * t645 - t787 * t646; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t818 + t787 * t645 + t790 * t646; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t799; t799; t664; t687; t792; -t796;];
tauJB = t1;

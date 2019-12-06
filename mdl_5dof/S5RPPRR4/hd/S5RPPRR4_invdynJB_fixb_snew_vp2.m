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
% m_mdh [6x1]
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:44:45
% EndTime: 2019-12-05 17:44:53
% DurationCPUTime: 7.93s
% Computational Cost: add. (73722->278), mult. (201875->377), div. (0->0), fcn. (138541->10), ass. (0->134)
t783 = sin(qJ(1));
t786 = cos(qJ(1));
t761 = t783 * g(2) - t786 * g(3);
t787 = qJD(1) ^ 2;
t831 = -t787 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t761;
t762 = t786 * g(2) + t783 * g(3);
t796 = -t787 * qJ(2) + qJDD(2) - t762;
t778 = sin(pkin(8));
t780 = cos(pkin(8));
t805 = -pkin(2) * t780 - qJ(3) * t778;
t823 = t778 * qJD(1);
t830 = (-pkin(1) + t805) * qJDD(1) + t796 - 0.2e1 * qJD(3) * t823;
t734 = -t780 * g(1) - t831 * t778;
t829 = mrSges(3,2) * t778;
t775 = t778 ^ 2;
t828 = t775 * t787;
t777 = sin(pkin(9));
t827 = t777 * t778;
t779 = cos(pkin(9));
t826 = t778 * t779;
t735 = -t778 * g(1) + t831 * t780;
t754 = (-mrSges(3,1) * t780 + t829) * qJD(1);
t753 = t805 * qJD(1);
t822 = t780 * qJD(1);
t724 = t753 * t822 + t735;
t801 = -pkin(3) * t780 - pkin(6) * t826;
t825 = t830 * t779;
t703 = t801 * qJDD(1) + (-t724 + (-pkin(3) * t775 * t779 + pkin(6) * t778 * t780) * t787) * t777 + t825;
t712 = t779 * t724 + t830 * t777;
t752 = t801 * qJD(1);
t820 = qJDD(1) * t778;
t816 = t777 * t820;
t818 = t777 ^ 2 * t828;
t704 = -pkin(3) * t818 - pkin(6) * t816 + t752 * t822 + t712;
t782 = sin(qJ(4));
t785 = cos(qJ(4));
t689 = t785 * t703 - t782 * t704;
t798 = (-t777 * t785 - t779 * t782) * t778;
t742 = qJD(1) * t798;
t797 = (-t777 * t782 + t779 * t785) * t778;
t728 = t742 * qJD(4) + qJDD(1) * t797;
t743 = qJD(1) * t797;
t819 = t780 * qJDD(1);
t765 = qJDD(4) - t819;
t766 = qJD(4) - t822;
t687 = (t742 * t766 - t728) * pkin(7) + (t742 * t743 + t765) * pkin(4) + t689;
t690 = t782 * t703 + t785 * t704;
t727 = -t743 * qJD(4) + qJDD(1) * t798;
t733 = t766 * pkin(4) - t743 * pkin(7);
t741 = t742 ^ 2;
t688 = -t741 * pkin(4) + t727 * pkin(7) - t766 * t733 + t690;
t781 = sin(qJ(5));
t784 = cos(qJ(5));
t685 = t784 * t687 - t781 * t688;
t721 = t784 * t742 - t781 * t743;
t699 = t721 * qJD(5) + t781 * t727 + t784 * t728;
t722 = t781 * t742 + t784 * t743;
t710 = -t721 * mrSges(6,1) + t722 * mrSges(6,2);
t764 = qJD(5) + t766;
t714 = -t764 * mrSges(6,2) + t721 * mrSges(6,3);
t760 = qJDD(5) + t765;
t681 = m(6) * t685 + t760 * mrSges(6,1) - t699 * mrSges(6,3) - t722 * t710 + t764 * t714;
t686 = t781 * t687 + t784 * t688;
t698 = -t722 * qJD(5) + t784 * t727 - t781 * t728;
t715 = t764 * mrSges(6,1) - t722 * mrSges(6,3);
t682 = m(6) * t686 - t760 * mrSges(6,2) + t698 * mrSges(6,3) + t721 * t710 - t764 * t715;
t673 = t784 * t681 + t781 * t682;
t725 = -t742 * mrSges(5,1) + t743 * mrSges(5,2);
t729 = -t766 * mrSges(5,2) + t742 * mrSges(5,3);
t671 = m(5) * t689 + t765 * mrSges(5,1) - t728 * mrSges(5,3) - t743 * t725 + t766 * t729 + t673;
t730 = t766 * mrSges(5,1) - t743 * mrSges(5,3);
t811 = -t781 * t681 + t784 * t682;
t672 = m(5) * t690 - t765 * mrSges(5,2) + t727 * mrSges(5,3) + t742 * t725 - t766 * t730 + t811;
t667 = t785 * t671 + t782 * t672;
t711 = -t777 * t724 + t825;
t809 = mrSges(4,1) * t777 + mrSges(4,2) * t779;
t746 = t809 * t823;
t799 = mrSges(4,2) * t780 - mrSges(4,3) * t827;
t749 = t799 * qJD(1);
t800 = -mrSges(4,1) * t780 - mrSges(4,3) * t826;
t665 = m(4) * t711 + t800 * qJDD(1) + (-t746 * t826 - t749 * t780) * qJD(1) + t667;
t750 = t800 * qJD(1);
t812 = -t782 * t671 + t785 * t672;
t666 = m(4) * t712 + t799 * qJDD(1) + (-t746 * t827 + t750 * t780) * qJD(1) + t812;
t813 = -t777 * t665 + t779 * t666;
t658 = m(3) * t735 + (qJDD(1) * mrSges(3,3) + qJD(1) * t754) * t780 + t813;
t723 = t753 * t823 + qJDD(3) - t734;
t713 = t779 * t752 * t823 + pkin(3) * t816 - pkin(6) * t818 + t723;
t692 = -t727 * pkin(4) - t741 * pkin(7) + t743 * t733 + t713;
t804 = m(6) * t692 - t698 * mrSges(6,1) + t699 * mrSges(6,2) - t721 * t714 + t722 * t715;
t790 = m(5) * t713 - t727 * mrSges(5,1) + t728 * mrSges(5,2) - t742 * t729 + t743 * t730 + t804;
t789 = m(4) * t723 + t790;
t802 = t749 * t777 + t750 * t779;
t677 = -t789 + ((-mrSges(3,3) - t809) * qJDD(1) + (-t754 - t802) * qJD(1)) * t778 + m(3) * t734;
t654 = t778 * t658 + t780 * t677;
t814 = t780 * t658 - t778 * t677;
t652 = m(2) * t761 - t787 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t814;
t661 = t779 * t665 + t777 * t666;
t748 = -qJDD(1) * pkin(1) + t796;
t791 = -m(3) * t748 + mrSges(3,1) * t819 - t661 + (t780 ^ 2 * t787 + t828) * mrSges(3,3);
t655 = m(2) * t762 - t787 * mrSges(2,2) + (mrSges(2,1) - t829) * qJDD(1) + t791;
t815 = t786 * t652 - t783 * t655;
t808 = Ifges(3,1) * t778 + Ifges(3,4) * t780;
t807 = Ifges(3,5) * t778 + Ifges(3,6) * t780;
t806 = Ifges(4,5) * t779 - Ifges(4,6) * t777;
t803 = -t783 * t652 - t786 * t655;
t705 = Ifges(6,5) * t722 + Ifges(6,6) * t721 + Ifges(6,3) * t764;
t707 = Ifges(6,1) * t722 + Ifges(6,4) * t721 + Ifges(6,5) * t764;
t674 = -mrSges(6,1) * t692 + mrSges(6,3) * t686 + Ifges(6,4) * t699 + Ifges(6,2) * t698 + Ifges(6,6) * t760 - t722 * t705 + t764 * t707;
t706 = Ifges(6,4) * t722 + Ifges(6,2) * t721 + Ifges(6,6) * t764;
t675 = mrSges(6,2) * t692 - mrSges(6,3) * t685 + Ifges(6,1) * t699 + Ifges(6,4) * t698 + Ifges(6,5) * t760 + t721 * t705 - t764 * t706;
t716 = Ifges(5,5) * t743 + Ifges(5,6) * t742 + Ifges(5,3) * t766;
t718 = Ifges(5,1) * t743 + Ifges(5,4) * t742 + Ifges(5,5) * t766;
t662 = -mrSges(5,1) * t713 + mrSges(5,3) * t690 + Ifges(5,4) * t728 + Ifges(5,2) * t727 + Ifges(5,6) * t765 - pkin(4) * t804 + pkin(7) * t811 + t784 * t674 + t781 * t675 - t743 * t716 + t766 * t718;
t717 = Ifges(5,4) * t743 + Ifges(5,2) * t742 + Ifges(5,6) * t766;
t663 = mrSges(5,2) * t713 - mrSges(5,3) * t689 + Ifges(5,1) * t728 + Ifges(5,4) * t727 + Ifges(5,5) * t765 - pkin(7) * t673 - t781 * t674 + t784 * t675 + t742 * t716 - t766 * t717;
t737 = (-Ifges(4,3) * t780 + t806 * t778) * qJD(1);
t794 = -Ifges(4,5) * t780 + (Ifges(4,1) * t779 - Ifges(4,4) * t777) * t778;
t739 = t794 * qJD(1);
t793 = -Ifges(4,6) * t780 + (Ifges(4,4) * t779 - Ifges(4,2) * t777) * t778;
t649 = -mrSges(4,1) * t723 + mrSges(4,3) * t712 + t782 * t663 + t785 * t662 - pkin(3) * t790 + pkin(6) * t812 + (-t737 * t826 - t780 * t739) * qJD(1) + t793 * qJDD(1);
t738 = t793 * qJD(1);
t650 = mrSges(4,2) * t723 - mrSges(4,3) * t711 - pkin(6) * t667 - t782 * t662 + t785 * t663 + (-t737 * t827 + t738 * t780) * qJD(1) + t794 * qJDD(1);
t755 = t807 * qJD(1);
t646 = mrSges(3,2) * t748 - mrSges(3,3) * t734 - qJ(3) * t661 + t808 * qJDD(1) - t777 * t649 + t779 * t650 + t755 * t822;
t792 = -mrSges(6,1) * t685 + mrSges(6,2) * t686 - Ifges(6,5) * t699 - Ifges(6,6) * t698 - Ifges(6,3) * t760 - t722 * t706 + t721 * t707;
t788 = mrSges(5,1) * t689 - mrSges(5,2) * t690 + Ifges(5,5) * t728 + Ifges(5,6) * t727 + Ifges(5,3) * t765 + pkin(4) * t673 + t743 * t717 - t742 * t718 - t792;
t648 = -t788 + ((Ifges(3,4) - t806) * qJDD(1) + (-t738 * t779 - t739 * t777 - t755) * qJD(1)) * t778 + mrSges(4,2) * t712 - mrSges(3,1) * t748 - mrSges(4,1) * t711 + (Ifges(3,2) + Ifges(4,3)) * t819 - pkin(3) * t667 - pkin(2) * t661 + mrSges(3,3) * t735;
t660 = mrSges(3,2) * t820 - t791;
t795 = mrSges(2,1) * t762 - mrSges(2,2) * t761 + Ifges(2,3) * qJDD(1) - pkin(1) * t660 + qJ(2) * t814 + t778 * t646 + t780 * t648;
t683 = (t802 * qJD(1) + t809 * qJDD(1)) * t778 + t789;
t644 = mrSges(2,1) * g(1) + mrSges(2,3) * t761 - mrSges(3,1) * t734 + mrSges(3,2) * t735 - t777 * t650 - t779 * t649 + pkin(2) * t683 - qJ(3) * t813 - pkin(1) * t654 + (Ifges(2,6) - t807) * qJDD(1) + (Ifges(2,5) - t778 * (Ifges(3,4) * t778 + Ifges(3,2) * t780) + t780 * t808) * t787;
t643 = -mrSges(2,2) * g(1) - mrSges(2,3) * t762 + Ifges(2,5) * qJDD(1) - t787 * Ifges(2,6) - qJ(2) * t654 + t780 * t646 - t778 * t648;
t1 = [(-m(1) - m(2)) * g(1) + t654; -m(1) * g(2) + t803; -m(1) * g(3) + t815; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t795; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t815 - t783 * t643 - t786 * t644; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t803 + t786 * t643 - t783 * t644; t795; t660; t683; t788; -t792;];
tauJB = t1;

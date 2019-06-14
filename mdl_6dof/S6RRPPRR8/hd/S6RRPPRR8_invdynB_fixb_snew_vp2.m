% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-05-06 11:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:18:27
% EndTime: 2019-05-06 11:18:41
% DurationCPUTime: 9.46s
% Computational Cost: add. (130897->363), mult. (288516->439), div. (0->0), fcn. (195085->10), ass. (0->141)
t807 = -2 * qJD(3);
t806 = -2 * qJD(4);
t805 = Ifges(4,1) + Ifges(5,1);
t799 = Ifges(4,4) - Ifges(5,5);
t804 = Ifges(4,5) + Ifges(5,4);
t803 = -Ifges(4,2) - Ifges(5,3);
t802 = Ifges(5,2) + Ifges(4,3);
t801 = Ifges(4,6) - Ifges(5,6);
t764 = sin(qJ(1));
t768 = cos(qJ(1));
t747 = -g(1) * t768 - g(2) * t764;
t770 = qJD(1) ^ 2;
t729 = -pkin(1) * t770 + qJDD(1) * pkin(7) + t747;
t763 = sin(qJ(2));
t767 = cos(qJ(2));
t705 = -t767 * g(3) - t763 * t729;
t740 = (-pkin(2) * t767 - qJ(3) * t763) * qJD(1);
t769 = qJD(2) ^ 2;
t789 = qJD(1) * t763;
t683 = -qJDD(2) * pkin(2) - t769 * qJ(3) + t740 * t789 + qJDD(3) - t705;
t786 = qJD(1) * qJD(2);
t784 = t767 * t786;
t742 = qJDD(1) * t763 + t784;
t760 = sin(pkin(10));
t796 = cos(pkin(10));
t715 = -qJDD(2) * t796 + t742 * t760;
t716 = t760 * qJDD(2) + t742 * t796;
t735 = t760 * qJD(2) + t789 * t796;
t754 = t767 * qJD(1);
t734 = -qJD(2) * t796 + t760 * t789;
t785 = t734 * t754;
t662 = t683 + (-t716 - t785) * qJ(4) + (-t735 * t754 + t715) * pkin(3) + t735 * t806;
t746 = t764 * g(1) - t768 * g(2);
t728 = -qJDD(1) * pkin(1) - t770 * pkin(7) - t746;
t751 = t763 * t786;
t743 = qJDD(1) * t767 - t751;
t680 = (-t742 - t784) * qJ(3) + (-t743 + t751) * pkin(2) + t728;
t706 = -g(3) * t763 + t767 * t729;
t684 = -pkin(2) * t769 + qJDD(2) * qJ(3) + t740 * t754 + t706;
t663 = t680 * t796 - t760 * t684 + t735 * t807;
t800 = -mrSges(4,3) - mrSges(5,2);
t795 = t767 ^ 2 * t770;
t741 = (-mrSges(3,1) * t767 + mrSges(3,2) * t763) * qJD(1);
t744 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t789;
t664 = t760 * t680 + t796 * t684 + t734 * t807;
t713 = -mrSges(4,1) * t754 - mrSges(4,3) * t735;
t714 = mrSges(5,1) * t754 + mrSges(5,2) * t735;
t702 = pkin(3) * t734 - qJ(4) * t735;
t654 = -pkin(3) * t795 - t743 * qJ(4) - t734 * t702 + t754 * t806 + t664;
t655 = t743 * pkin(3) - qJ(4) * t795 + t735 * t702 + qJDD(4) - t663;
t644 = (-t716 + t785) * pkin(8) + (t734 * t735 + t743) * pkin(4) + t655;
t717 = pkin(4) * t754 - pkin(8) * t735;
t732 = t734 ^ 2;
t646 = -pkin(4) * t732 + pkin(8) * t715 - t717 * t754 + t654;
t762 = sin(qJ(5));
t766 = cos(qJ(5));
t638 = t766 * t644 - t762 * t646;
t700 = t734 * t766 - t735 * t762;
t670 = qJD(5) * t700 + t715 * t762 + t716 * t766;
t701 = t734 * t762 + t735 * t766;
t739 = qJDD(5) + t743;
t749 = t754 + qJD(5);
t636 = (t700 * t749 - t670) * pkin(9) + (t700 * t701 + t739) * pkin(5) + t638;
t639 = t762 * t644 + t766 * t646;
t669 = -qJD(5) * t701 + t715 * t766 - t716 * t762;
t689 = pkin(5) * t749 - pkin(9) * t701;
t699 = t700 ^ 2;
t637 = -pkin(5) * t699 + pkin(9) * t669 - t689 * t749 + t639;
t761 = sin(qJ(6));
t765 = cos(qJ(6));
t634 = t636 * t765 - t637 * t761;
t676 = t700 * t765 - t701 * t761;
t650 = qJD(6) * t676 + t669 * t761 + t670 * t765;
t677 = t700 * t761 + t701 * t765;
t661 = -mrSges(7,1) * t676 + mrSges(7,2) * t677;
t748 = qJD(6) + t749;
t665 = -mrSges(7,2) * t748 + mrSges(7,3) * t676;
t730 = qJDD(6) + t739;
t632 = m(7) * t634 + mrSges(7,1) * t730 - mrSges(7,3) * t650 - t661 * t677 + t665 * t748;
t635 = t636 * t761 + t637 * t765;
t649 = -qJD(6) * t677 + t669 * t765 - t670 * t761;
t666 = mrSges(7,1) * t748 - mrSges(7,3) * t677;
t633 = m(7) * t635 - mrSges(7,2) * t730 + mrSges(7,3) * t649 + t661 * t676 - t666 * t748;
t623 = t765 * t632 + t761 * t633;
t678 = -mrSges(6,1) * t700 + mrSges(6,2) * t701;
t685 = -mrSges(6,2) * t749 + mrSges(6,3) * t700;
t621 = m(6) * t638 + mrSges(6,1) * t739 - mrSges(6,3) * t670 - t678 * t701 + t685 * t749 + t623;
t686 = mrSges(6,1) * t749 - mrSges(6,3) * t701;
t779 = -t632 * t761 + t765 * t633;
t622 = m(6) * t639 - mrSges(6,2) * t739 + mrSges(6,3) * t669 + t678 * t700 - t686 * t749 + t779;
t780 = -t762 * t621 + t766 * t622;
t776 = m(5) * t654 - t743 * mrSges(5,3) + t780;
t703 = mrSges(5,1) * t734 - mrSges(5,3) * t735;
t790 = -mrSges(4,1) * t734 - mrSges(4,2) * t735 - t703;
t617 = m(4) * t664 + t743 * mrSges(4,2) + t790 * t734 + t800 * t715 + (t713 - t714) * t754 + t776;
t711 = -mrSges(5,2) * t734 - mrSges(5,3) * t754;
t712 = mrSges(4,2) * t754 - mrSges(4,3) * t734;
t619 = t766 * t621 + t762 * t622;
t775 = -m(5) * t655 - t743 * mrSges(5,1) - t619;
t618 = m(4) * t663 - t743 * mrSges(4,1) + t790 * t735 + t800 * t716 + (-t711 - t712) * t754 + t775;
t781 = t796 * t617 - t618 * t760;
t612 = m(3) * t706 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t743 - qJD(2) * t744 + t741 * t754 + t781;
t745 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t754;
t652 = -pkin(4) * t715 - pkin(8) * t732 + t735 * t717 - t662;
t641 = -pkin(5) * t669 - pkin(9) * t699 + t689 * t701 + t652;
t777 = m(7) * t641 - t649 * mrSges(7,1) + t650 * mrSges(7,2) - t676 * t665 + t677 * t666;
t774 = -m(6) * t652 + t669 * mrSges(6,1) - t670 * mrSges(6,2) + t700 * t685 - t701 * t686 - t777;
t628 = m(5) * t662 + t715 * mrSges(5,1) - t716 * mrSges(5,3) + t734 * t711 - t735 * t714 + t774;
t771 = -m(4) * t683 - t715 * mrSges(4,1) - t716 * mrSges(4,2) - t734 * t712 - t735 * t713 - t628;
t627 = m(3) * t705 + qJDD(2) * mrSges(3,1) - t742 * mrSges(3,3) + qJD(2) * t745 - t741 * t789 + t771;
t782 = t767 * t612 - t627 * t763;
t606 = m(2) * t747 - mrSges(2,1) * t770 - qJDD(1) * mrSges(2,2) + t782;
t613 = t760 * t617 + t618 * t796;
t772 = -m(3) * t728 + t743 * mrSges(3,1) - t742 * mrSges(3,2) - t744 * t789 + t745 * t754 - t613;
t609 = m(2) * t746 + qJDD(1) * mrSges(2,1) - t770 * mrSges(2,2) + t772;
t794 = t764 * t606 + t768 * t609;
t607 = t763 * t612 + t767 * t627;
t793 = t803 * t734 + t799 * t735 - t801 * t754;
t792 = t801 * t734 - t804 * t735 + t802 * t754;
t791 = t799 * t734 - t805 * t735 + t804 * t754;
t783 = t768 * t606 - t609 * t764;
t727 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t763 + Ifges(3,4) * t767) * qJD(1);
t726 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t763 + Ifges(3,2) * t767) * qJD(1);
t725 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t763 + Ifges(3,6) * t767) * qJD(1);
t673 = Ifges(6,1) * t701 + Ifges(6,4) * t700 + Ifges(6,5) * t749;
t672 = Ifges(6,4) * t701 + Ifges(6,2) * t700 + Ifges(6,6) * t749;
t671 = Ifges(6,5) * t701 + Ifges(6,6) * t700 + Ifges(6,3) * t749;
t658 = Ifges(7,1) * t677 + Ifges(7,4) * t676 + Ifges(7,5) * t748;
t657 = Ifges(7,4) * t677 + Ifges(7,2) * t676 + Ifges(7,6) * t748;
t656 = Ifges(7,5) * t677 + Ifges(7,6) * t676 + Ifges(7,3) * t748;
t625 = mrSges(7,2) * t641 - mrSges(7,3) * t634 + Ifges(7,1) * t650 + Ifges(7,4) * t649 + Ifges(7,5) * t730 + t656 * t676 - t657 * t748;
t624 = -mrSges(7,1) * t641 + mrSges(7,3) * t635 + Ifges(7,4) * t650 + Ifges(7,2) * t649 + Ifges(7,6) * t730 - t656 * t677 + t658 * t748;
t615 = mrSges(6,2) * t652 - mrSges(6,3) * t638 + Ifges(6,1) * t670 + Ifges(6,4) * t669 + Ifges(6,5) * t739 - pkin(9) * t623 - t624 * t761 + t625 * t765 + t671 * t700 - t672 * t749;
t614 = -mrSges(6,1) * t652 + mrSges(6,3) * t639 + Ifges(6,4) * t670 + Ifges(6,2) * t669 + Ifges(6,6) * t739 - pkin(5) * t777 + pkin(9) * t779 + t765 * t624 + t761 * t625 - t701 * t671 + t749 * t673;
t603 = mrSges(4,2) * t683 + mrSges(5,2) * t655 - mrSges(4,3) * t663 - mrSges(5,3) * t662 - pkin(8) * t619 - qJ(4) * t628 - t762 * t614 + t766 * t615 - t799 * t715 + t805 * t716 + t792 * t734 - t743 * t804 + t793 * t754;
t602 = -mrSges(4,1) * t683 - mrSges(5,1) * t662 + mrSges(5,2) * t654 + mrSges(4,3) * t664 - pkin(3) * t628 - pkin(4) * t774 - pkin(8) * t780 - t766 * t614 - t762 * t615 + t803 * t715 + t799 * t716 + t792 * t735 - t743 * t801 + t791 * t754;
t601 = -t725 * t789 + (Ifges(3,2) + t802) * t743 - qJ(4) * (-t714 * t754 + t776) - pkin(3) * (-t711 * t754 + t775) + Ifges(3,6) * qJDD(2) + (qJ(4) * t703 + t791) * t734 + (pkin(3) * t703 - t793) * t735 + (mrSges(5,2) * qJ(4) + t801) * t715 + (mrSges(5,2) * pkin(3) - t804) * t716 + Ifges(6,3) * t739 + Ifges(3,4) * t742 + Ifges(7,3) * t730 + qJD(2) * t727 - mrSges(3,1) * t728 + mrSges(3,3) * t706 - t700 * t673 + t701 * t672 - t676 * t658 + t677 * t657 + Ifges(6,6) * t669 + Ifges(6,5) * t670 - mrSges(5,3) * t654 + mrSges(5,1) * t655 - mrSges(4,1) * t663 + mrSges(4,2) * t664 + Ifges(7,6) * t649 + Ifges(7,5) * t650 - mrSges(6,2) * t639 + mrSges(6,1) * t638 + mrSges(7,1) * t634 - mrSges(7,2) * t635 + pkin(5) * t623 + pkin(4) * t619 - pkin(2) * t613;
t600 = mrSges(3,2) * t728 - mrSges(3,3) * t705 + Ifges(3,1) * t742 + Ifges(3,4) * t743 + Ifges(3,5) * qJDD(2) - qJ(3) * t613 - qJD(2) * t726 - t760 * t602 + t603 * t796 + t725 * t754;
t599 = Ifges(2,6) * qJDD(1) + t770 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t747 - Ifges(3,5) * t742 - Ifges(3,6) * t743 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t705 + mrSges(3,2) * t706 - t760 * t603 - t796 * t602 - pkin(2) * t771 - qJ(3) * t781 - pkin(1) * t607 + (-t726 * t763 + t727 * t767) * qJD(1);
t598 = -mrSges(2,2) * g(3) - mrSges(2,3) * t746 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t770 - pkin(7) * t607 + t600 * t767 - t601 * t763;
t1 = [-m(1) * g(1) + t783; -m(1) * g(2) + t794; (-m(1) - m(2)) * g(3) + t607; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t794 + t768 * t598 - t764 * t599; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t783 + t764 * t598 + t768 * t599; -mrSges(1,1) * g(2) + mrSges(2,1) * t746 + mrSges(1,2) * g(1) - mrSges(2,2) * t747 + Ifges(2,3) * qJDD(1) + pkin(1) * t772 + pkin(7) * t782 + t763 * t600 + t767 * t601;];
tauB  = t1;

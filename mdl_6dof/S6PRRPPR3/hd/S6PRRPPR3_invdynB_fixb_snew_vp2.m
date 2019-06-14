% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-05-05 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:57:51
% EndTime: 2019-05-05 02:57:56
% DurationCPUTime: 3.83s
% Computational Cost: add. (38833->305), mult. (76983->367), div. (0->0), fcn. (42152->10), ass. (0->132)
t796 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t771 = Ifges(4,4) - Ifges(5,5) + Ifges(6,4);
t770 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t795 = Ifges(4,2) + Ifges(5,3) + Ifges(6,1);
t794 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t769 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t737 = sin(pkin(10));
t739 = cos(pkin(10));
t712 = g(1) * t737 - g(2) * t739;
t713 = -g(1) * t739 - g(2) * t737;
t735 = -g(3) + qJDD(1);
t747 = cos(qJ(2));
t740 = cos(pkin(6));
t744 = sin(qJ(2));
t781 = t740 * t744;
t738 = sin(pkin(6));
t783 = t738 * t744;
t656 = t712 * t781 + t747 * t713 + t735 * t783;
t749 = qJD(2) ^ 2;
t654 = -pkin(2) * t749 + qJDD(2) * pkin(8) + t656;
t668 = -t712 * t738 + t735 * t740;
t743 = sin(qJ(3));
t746 = cos(qJ(3));
t650 = t746 * t654 + t743 * t668;
t702 = (-pkin(3) * t746 - qJ(4) * t743) * qJD(2);
t775 = qJD(2) * t746;
t793 = qJDD(3) * qJ(4) + t702 * t775 + t650;
t773 = qJD(2) * qJD(3);
t765 = t746 * t773;
t707 = qJDD(2) * t743 + t765;
t772 = qJD(2) * qJD(5);
t792 = -0.2e1 * t743 * t772 + (-t707 + t765) * qJ(5);
t791 = -2 * qJD(4);
t790 = 2 * qJD(4);
t789 = -pkin(3) - pkin(9);
t788 = pkin(4) + pkin(9);
t748 = qJD(3) ^ 2;
t787 = pkin(3) * t748;
t786 = mrSges(4,3) + mrSges(5,2);
t784 = t746 ^ 2 * t749;
t782 = t738 * t747;
t780 = t740 * t747;
t779 = t746 * t749;
t649 = -t743 * t654 + t668 * t746;
t703 = (-mrSges(5,1) * t746 - mrSges(5,3) * t743) * qJD(2);
t704 = (-mrSges(4,1) * t746 + mrSges(4,2) * t743) * qJD(2);
t718 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t775;
t719 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t775;
t776 = qJD(2) * t743;
t760 = t702 * t776 + qJDD(4) - t649;
t647 = -qJDD(3) * pkin(3) - qJ(4) * t748 + t760;
t720 = mrSges(5,2) * t775 + qJD(3) * mrSges(5,3);
t642 = (-t743 * t779 - qJDD(3)) * pkin(4) + t647 + t792;
t705 = (mrSges(6,1) * t743 - mrSges(6,2) * t746) * qJD(2);
t764 = t743 * t773;
t708 = qJDD(2) * t746 - t764;
t714 = -qJD(3) * pkin(4) - qJ(5) * t776;
t655 = t712 * t780 - t744 * t713 + t735 * t782;
t653 = -qJDD(2) * pkin(2) - t749 * pkin(8) - t655;
t757 = -t708 * pkin(3) + t653 + (-t707 - t765) * qJ(4);
t752 = -qJ(5) * t784 + qJDD(5) - t757 + (t714 + t790) * t776;
t637 = t788 * t708 + pkin(5) * t707 + (pkin(5) * t746 + t789 * t743) * t773 + t752;
t706 = (pkin(5) * t743 + pkin(9) * t746) * qJD(2);
t640 = (-pkin(5) - qJ(4)) * t748 + (-pkin(4) * t779 - qJD(2) * t706) * t743 + (-pkin(3) - t788) * qJDD(3) + t760 + t792;
t742 = sin(qJ(6));
t745 = cos(qJ(6));
t635 = t637 * t745 - t640 * t742;
t700 = -qJD(3) * t745 + t742 * t775;
t663 = qJD(6) * t700 - qJDD(3) * t742 - t708 * t745;
t701 = -qJD(3) * t742 - t745 * t775;
t664 = -mrSges(7,1) * t700 + mrSges(7,2) * t701;
t725 = qJD(6) + t776;
t666 = -mrSges(7,2) * t725 + mrSges(7,3) * t700;
t696 = qJDD(6) + t707;
t633 = m(7) * t635 + mrSges(7,1) * t696 - mrSges(7,3) * t663 - t664 * t701 + t666 * t725;
t636 = t637 * t742 + t640 * t745;
t662 = -qJD(6) * t701 - qJDD(3) * t745 + t708 * t742;
t667 = mrSges(7,1) * t725 - mrSges(7,3) * t701;
t634 = m(7) * t636 - mrSges(7,2) * t696 + mrSges(7,3) * t662 + t664 * t700 - t667 * t725;
t777 = -t742 * t633 + t745 * t634;
t761 = -m(6) * t642 + t705 * t776 - t777;
t755 = -m(5) * t647 + qJDD(3) * mrSges(5,1) + qJD(3) * t720 + t761;
t620 = m(4) * t649 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t718 + t719) * qJD(3) + (-t703 - t704) * t776 + (mrSges(6,3) - t786) * t707 + t755;
t716 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t776;
t729 = qJD(3) * t790;
t646 = t729 - t787 + t793;
t717 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t776;
t754 = pkin(4) * t784 + qJ(5) * t708 - t793;
t641 = 0.2e1 * t746 * t772 + t787 + (t791 - t714) * qJD(3) + t754;
t715 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t776;
t639 = qJDD(3) * pkin(5) + qJD(3) * t714 + t729 + t789 * t748 + (-0.2e1 * qJD(5) - t706) * t775 - t754;
t756 = -m(7) * t639 + mrSges(7,1) * t662 - t663 * mrSges(7,2) + t666 * t700 - t701 * t667;
t753 = -m(6) * t641 + qJDD(3) * mrSges(6,1) - t708 * mrSges(6,3) + qJD(3) * t715 - t756;
t751 = m(5) * t646 + qJDD(3) * mrSges(5,3) + qJD(3) * t717 + t703 * t775 + t753;
t626 = (t704 - t705) * t775 + t786 * t708 + m(4) * t650 - qJD(3) * t716 - qJDD(3) * mrSges(4,2) + t751;
t762 = -t620 * t743 + t746 * t626;
t612 = m(3) * t656 - mrSges(3,1) * t749 - qJDD(2) * mrSges(3,2) + t762;
t615 = t746 * t620 + t743 * t626;
t614 = m(3) * t668 + t615;
t648 = (pkin(3) * qJD(3) + t791) * t776 + t757;
t623 = t745 * t633 + t742 * t634;
t644 = -pkin(3) * t764 + pkin(4) * t708 + t752;
t759 = -m(6) * t644 - t707 * mrSges(6,1) + t708 * mrSges(6,2) - t715 * t776 + t718 * t775 - t623;
t621 = m(5) * t648 - mrSges(5,1) * t708 - t707 * mrSges(5,3) - t717 * t776 - t720 * t775 + t759;
t750 = -m(4) * t653 + t708 * mrSges(4,1) - mrSges(4,2) * t707 - t716 * t776 + t719 * t775 - t621;
t618 = m(3) * t655 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t749 + t750;
t603 = t612 * t781 - t614 * t738 + t618 * t780;
t601 = m(2) * t712 + t603;
t608 = t747 * t612 - t618 * t744;
t607 = m(2) * t713 + t608;
t778 = t739 * t601 + t737 * t607;
t774 = qJD(3) * t718;
t602 = t612 * t783 + t740 * t614 + t618 * t782;
t768 = t794 * qJD(3) + (-t770 * t743 - t769 * t746) * qJD(2);
t767 = -t769 * qJD(3) + (-t771 * t743 - t795 * t746) * qJD(2);
t766 = t770 * qJD(3) + (t796 * t743 + t771 * t746) * qJD(2);
t763 = -t601 * t737 + t739 * t607;
t657 = Ifges(7,5) * t701 + Ifges(7,6) * t700 + Ifges(7,3) * t725;
t659 = Ifges(7,1) * t701 + Ifges(7,4) * t700 + Ifges(7,5) * t725;
t627 = -mrSges(7,1) * t639 + mrSges(7,3) * t636 + Ifges(7,4) * t663 + Ifges(7,2) * t662 + Ifges(7,6) * t696 - t657 * t701 + t659 * t725;
t658 = Ifges(7,4) * t701 + Ifges(7,2) * t700 + Ifges(7,6) * t725;
t628 = mrSges(7,2) * t639 - mrSges(7,3) * t635 + Ifges(7,1) * t663 + Ifges(7,4) * t662 + Ifges(7,5) * t696 + t657 * t700 - t658 * t725;
t599 = -t745 * t628 - qJ(5) * t753 + t742 * t627 - pkin(4) * t759 + mrSges(6,3) * t641 - mrSges(6,2) * t644 + mrSges(5,2) * t646 - mrSges(5,1) * t648 + mrSges(4,3) * t650 - mrSges(4,1) * t653 + pkin(9) * t623 - pkin(3) * t621 + t795 * t708 + t771 * t707 + t769 * qJDD(3) + t766 * qJD(3) + (qJ(5) * t705 * t746 + t768 * t743) * qJD(2);
t622 = qJDD(3) * mrSges(6,2) - mrSges(6,3) * t707 - t761 + t774;
t604 = Ifges(7,3) * t696 - t700 * t659 + t701 * t658 + mrSges(7,1) * t635 - mrSges(7,2) * t636 + Ifges(7,6) * t662 + Ifges(7,5) * t663 - mrSges(6,3) * t642 + mrSges(6,1) * t644 + mrSges(5,2) * t647 - mrSges(5,3) * t648 - mrSges(4,3) * t649 + mrSges(4,2) * t653 + pkin(5) * t623 - qJ(4) * t621 - qJ(5) * t622 + t771 * t708 + t796 * t707 + t770 * qJDD(3) + t767 * qJD(3) - t768 * t775;
t597 = mrSges(3,2) * t668 - mrSges(3,3) * t655 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t749 - pkin(8) * t615 - t599 * t743 + t604 * t746;
t598 = t745 * t627 + t749 * Ifges(3,5) + t742 * t628 - qJ(4) * t751 - pkin(3) * (t755 - t774) - pkin(2) * t615 + pkin(5) * t756 + mrSges(6,1) * t641 + Ifges(3,6) * qJDD(2) + mrSges(3,3) * t656 - mrSges(3,1) * t668 - mrSges(6,2) * t642 - mrSges(5,3) * t646 + mrSges(5,1) * t647 - mrSges(4,1) * t649 + mrSges(4,2) * t650 + pkin(9) * t777 + pkin(4) * t622 + (-qJ(4) * mrSges(5,2) - t769) * t708 + (pkin(3) * mrSges(6,2) + t794) * qJDD(3) + (-pkin(3) * (-mrSges(5,2) + mrSges(6,3)) - t770) * t707 + ((qJ(4) * t705 + t766) * t746 + (pkin(3) * t703 + t767) * t743) * qJD(2);
t758 = pkin(7) * t608 + t597 * t744 + t598 * t747;
t596 = mrSges(3,1) * t655 - mrSges(3,2) * t656 + Ifges(3,3) * qJDD(2) + pkin(2) * t750 + pkin(8) * t762 + t746 * t599 + t743 * t604;
t595 = mrSges(2,2) * t735 - mrSges(2,3) * t712 + t597 * t747 - t598 * t744 + (-t602 * t738 - t603 * t740) * pkin(7);
t594 = -mrSges(2,1) * t735 + mrSges(2,3) * t713 - pkin(1) * t602 - t596 * t738 + t758 * t740;
t1 = [-m(1) * g(1) + t763; -m(1) * g(2) + t778; -m(1) * g(3) + m(2) * t735 + t602; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t778 - t737 * t594 + t739 * t595; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t763 + t739 * t594 + t737 * t595; -mrSges(1,1) * g(2) + mrSges(2,1) * t712 + mrSges(1,2) * g(1) - mrSges(2,2) * t713 + pkin(1) * t603 + t596 * t740 + t758 * t738;];
tauB  = t1;

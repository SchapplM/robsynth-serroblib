% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:20:36
% EndTime: 2019-05-07 04:20:48
% DurationCPUTime: 10.72s
% Computational Cost: add. (153721->370), mult. (349924->448), div. (0->0), fcn. (253102->10), ass. (0->145)
t799 = -2 * qJD(4);
t798 = Ifges(5,1) + Ifges(6,2);
t797 = -Ifges(6,1) - Ifges(5,3);
t793 = Ifges(5,4) + Ifges(6,6);
t792 = Ifges(5,5) - Ifges(6,4);
t796 = Ifges(5,2) + Ifges(6,3);
t791 = Ifges(5,6) - Ifges(6,5);
t757 = sin(qJ(2));
t761 = cos(qJ(2));
t778 = qJD(1) * qJD(2);
t741 = qJDD(1) * t757 + t761 * t778;
t758 = sin(qJ(1));
t762 = cos(qJ(1));
t747 = -g(1) * t762 - g(2) * t758;
t763 = qJD(1) ^ 2;
t736 = -pkin(1) * t763 + qJDD(1) * pkin(7) + t747;
t788 = t757 * t736;
t794 = pkin(2) * t763;
t693 = qJDD(2) * pkin(2) - t741 * pkin(8) - t788 + (pkin(8) * t778 + t757 * t794 - g(3)) * t761;
t723 = -g(3) * t757 + t761 * t736;
t742 = qJDD(1) * t761 - t757 * t778;
t781 = qJD(1) * t757;
t745 = qJD(2) * pkin(2) - pkin(8) * t781;
t753 = t761 ^ 2;
t694 = pkin(8) * t742 - qJD(2) * t745 - t753 * t794 + t723;
t756 = sin(qJ(3));
t760 = cos(qJ(3));
t664 = t760 * t693 - t756 * t694;
t733 = (-t756 * t757 + t760 * t761) * qJD(1);
t701 = qJD(3) * t733 + t741 * t760 + t742 * t756;
t734 = (t756 * t761 + t757 * t760) * qJD(1);
t751 = qJDD(2) + qJDD(3);
t752 = qJD(2) + qJD(3);
t649 = (t733 * t752 - t701) * qJ(4) + (t733 * t734 + t751) * pkin(3) + t664;
t665 = t756 * t693 + t760 * t694;
t700 = -qJD(3) * t734 - t741 * t756 + t742 * t760;
t725 = pkin(3) * t752 - qJ(4) * t734;
t729 = t733 ^ 2;
t652 = -pkin(3) * t729 + qJ(4) * t700 - t725 * t752 + t665;
t754 = sin(pkin(10));
t790 = cos(pkin(10));
t720 = t754 * t733 + t790 * t734;
t646 = t790 * t649 - t754 * t652 + t720 * t799;
t795 = -2 * qJD(5);
t719 = -t790 * t733 + t754 * t734;
t789 = t719 * t752;
t673 = t754 * t700 + t790 * t701;
t688 = mrSges(5,1) * t719 + mrSges(5,2) * t720;
t703 = -mrSges(5,2) * t752 - mrSges(5,3) * t719;
t705 = mrSges(6,1) * t719 - mrSges(6,3) * t752;
t687 = pkin(4) * t719 - qJ(5) * t720;
t750 = t752 ^ 2;
t643 = -t751 * pkin(4) - t750 * qJ(5) + t720 * t687 + qJDD(5) - t646;
t638 = (t719 * t720 - t751) * pkin(9) + (t673 + t789) * pkin(5) + t643;
t672 = -t790 * t700 + t754 * t701;
t707 = pkin(5) * t720 - pkin(9) * t752;
t716 = t719 ^ 2;
t746 = t758 * g(1) - t762 * g(2);
t772 = -qJDD(1) * pkin(1) - t746;
t702 = -t742 * pkin(2) + t745 * t781 + (-pkin(8) * t753 - pkin(7)) * t763 + t772;
t657 = -t700 * pkin(3) - t729 * qJ(4) + t734 * t725 + qJDD(4) + t702;
t765 = (-t673 + t789) * qJ(5) + t657 + (t752 * pkin(4) + t795) * t720;
t641 = (pkin(4) + pkin(9)) * t672 - t720 * t707 - t716 * pkin(5) + t765;
t755 = sin(qJ(6));
t759 = cos(qJ(6));
t636 = t638 * t759 - t641 * t755;
t698 = t719 * t759 - t752 * t755;
t655 = qJD(6) * t698 + t672 * t755 + t751 * t759;
t671 = qJDD(6) + t673;
t699 = t719 * t755 + t752 * t759;
t674 = -mrSges(7,1) * t698 + mrSges(7,2) * t699;
t712 = qJD(6) + t720;
t675 = -mrSges(7,2) * t712 + mrSges(7,3) * t698;
t634 = m(7) * t636 + mrSges(7,1) * t671 - mrSges(7,3) * t655 - t674 * t699 + t675 * t712;
t637 = t638 * t755 + t641 * t759;
t654 = -qJD(6) * t699 + t672 * t759 - t751 * t755;
t676 = mrSges(7,1) * t712 - mrSges(7,3) * t699;
t635 = m(7) * t637 - mrSges(7,2) * t671 + mrSges(7,3) * t654 + t674 * t698 - t676 * t712;
t626 = t759 * t634 + t755 * t635;
t689 = -mrSges(6,2) * t719 - mrSges(6,3) * t720;
t769 = -m(6) * t643 - t673 * mrSges(6,1) - t720 * t689 - t626;
t624 = m(5) * t646 - t673 * mrSges(5,3) - t720 * t688 + (t703 - t705) * t752 + (mrSges(5,1) - mrSges(6,2)) * t751 + t769;
t710 = t719 * t799;
t785 = t754 * t649 + t790 * t652;
t647 = t710 + t785;
t704 = mrSges(5,1) * t752 - mrSges(5,3) * t720;
t771 = t750 * pkin(4) - t751 * qJ(5) - t785;
t642 = t752 * t795 + ((2 * qJD(4)) + t687) * t719 + t771;
t706 = mrSges(6,1) * t720 + mrSges(6,2) * t752;
t640 = -t672 * pkin(5) - t716 * pkin(9) - t719 * t687 + t710 + ((2 * qJD(5)) + t707) * t752 - t771;
t770 = -m(7) * t640 + t654 * mrSges(7,1) - t655 * mrSges(7,2) + t698 * t675 - t699 * t676;
t767 = -m(6) * t642 + t751 * mrSges(6,3) + t752 * t706 - t770;
t631 = m(5) * t647 - t751 * mrSges(5,2) - t752 * t704 + (-t688 - t689) * t719 + (-mrSges(5,3) - mrSges(6,1)) * t672 + t767;
t620 = t790 * t624 + t754 * t631;
t721 = -mrSges(4,1) * t733 + mrSges(4,2) * t734;
t724 = -mrSges(4,2) * t752 + mrSges(4,3) * t733;
t618 = m(4) * t664 + mrSges(4,1) * t751 - mrSges(4,3) * t701 - t721 * t734 + t724 * t752 + t620;
t726 = mrSges(4,1) * t752 - mrSges(4,3) * t734;
t774 = -t624 * t754 + t790 * t631;
t619 = m(4) * t665 - mrSges(4,2) * t751 + mrSges(4,3) * t700 + t721 * t733 - t726 * t752 + t774;
t613 = t760 * t618 + t756 * t619;
t722 = -t761 * g(3) - t788;
t740 = (-mrSges(3,1) * t761 + mrSges(3,2) * t757) * qJD(1);
t780 = qJD(1) * t761;
t744 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t780;
t611 = m(3) * t722 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t741 + qJD(2) * t744 - t740 * t781 + t613;
t743 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t781;
t775 = -t618 * t756 + t760 * t619;
t612 = m(3) * t723 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t742 - qJD(2) * t743 + t740 * t780 + t775;
t776 = -t611 * t757 + t761 * t612;
t605 = m(2) * t747 - mrSges(2,1) * t763 - qJDD(1) * mrSges(2,2) + t776;
t735 = -t763 * pkin(7) + t772;
t645 = t672 * pkin(4) + t765;
t786 = -t755 * t634 + t759 * t635;
t625 = m(6) * t645 - t672 * mrSges(6,2) - t673 * mrSges(6,3) - t719 * t705 - t720 * t706 + t786;
t768 = m(5) * t657 + t672 * mrSges(5,1) + t673 * mrSges(5,2) + t719 * t703 + t720 * t704 + t625;
t766 = m(4) * t702 - t700 * mrSges(4,1) + t701 * mrSges(4,2) - t733 * t724 + t734 * t726 + t768;
t764 = -m(3) * t735 + t742 * mrSges(3,1) - t741 * mrSges(3,2) - t743 * t781 + t744 * t780 - t766;
t622 = m(2) * t746 + qJDD(1) * mrSges(2,1) - t763 * mrSges(2,2) + t764;
t787 = t758 * t605 + t762 * t622;
t606 = t761 * t611 + t757 * t612;
t784 = t791 * t719 - t792 * t720 + t797 * t752;
t783 = t796 * t719 - t793 * t720 - t791 * t752;
t782 = -t793 * t719 + t798 * t720 + t792 * t752;
t777 = t762 * t605 - t622 * t758;
t732 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t757 + Ifges(3,4) * t761) * qJD(1);
t731 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t757 + Ifges(3,2) * t761) * qJD(1);
t730 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t757 + Ifges(3,6) * t761) * qJD(1);
t715 = Ifges(4,1) * t734 + Ifges(4,4) * t733 + Ifges(4,5) * t752;
t714 = Ifges(4,4) * t734 + Ifges(4,2) * t733 + Ifges(4,6) * t752;
t713 = Ifges(4,5) * t734 + Ifges(4,6) * t733 + Ifges(4,3) * t752;
t660 = Ifges(7,1) * t699 + Ifges(7,4) * t698 + Ifges(7,5) * t712;
t659 = Ifges(7,4) * t699 + Ifges(7,2) * t698 + Ifges(7,6) * t712;
t658 = Ifges(7,5) * t699 + Ifges(7,6) * t698 + Ifges(7,3) * t712;
t628 = mrSges(7,2) * t640 - mrSges(7,3) * t636 + Ifges(7,1) * t655 + Ifges(7,4) * t654 + Ifges(7,5) * t671 + t658 * t698 - t659 * t712;
t627 = -mrSges(7,1) * t640 + mrSges(7,3) * t637 + Ifges(7,4) * t655 + Ifges(7,2) * t654 + Ifges(7,6) * t671 - t658 * t699 + t660 * t712;
t614 = mrSges(6,1) * t643 + mrSges(7,1) * t636 + mrSges(5,2) * t657 - mrSges(7,2) * t637 - mrSges(5,3) * t646 - mrSges(6,3) * t645 + Ifges(7,5) * t655 + Ifges(7,6) * t654 + Ifges(7,3) * t671 + pkin(5) * t626 - qJ(5) * t625 + t699 * t659 - t698 * t660 + t783 * t752 + t792 * t751 + t784 * t719 + t798 * t673 - t793 * t672;
t607 = -mrSges(5,1) * t657 - mrSges(6,1) * t642 + mrSges(6,2) * t645 + mrSges(5,3) * t647 - pkin(4) * t625 - pkin(5) * t770 - pkin(9) * t786 - t759 * t627 - t755 * t628 - t672 * t796 + t793 * t673 + t784 * t720 + t791 * t751 + t782 * t752;
t602 = mrSges(4,2) * t702 - mrSges(4,3) * t664 + Ifges(4,1) * t701 + Ifges(4,4) * t700 + Ifges(4,5) * t751 - qJ(4) * t620 - t754 * t607 + t790 * t614 + t733 * t713 - t752 * t714;
t601 = -mrSges(4,1) * t702 + mrSges(4,3) * t665 + Ifges(4,4) * t701 + Ifges(4,2) * t700 + Ifges(4,6) * t751 - pkin(3) * t768 + qJ(4) * t774 + t790 * t607 + t754 * t614 - t734 * t713 + t752 * t715;
t600 = -qJ(5) * t767 - pkin(4) * (-t752 * t705 + t769) + (mrSges(6,1) * qJ(5) + t791) * t672 - t792 * t673 - pkin(1) * t606 - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) + (-t731 * t757 + t732 * t761) * qJD(1) + t763 * Ifges(2,5) - t759 * t628 + t755 * t627 + mrSges(2,3) * t747 - Ifges(3,5) * t741 - Ifges(3,6) * t742 + t733 * t715 - t734 * t714 - mrSges(3,1) * t722 + mrSges(3,2) * t723 + Ifges(2,6) * qJDD(1) - Ifges(4,6) * t700 - Ifges(4,5) * t701 - mrSges(4,1) * t664 + mrSges(4,2) * t665 - mrSges(5,1) * t646 + mrSges(5,2) * t647 - mrSges(6,2) * t643 + mrSges(6,3) * t642 + pkin(9) * t626 - pkin(3) * t620 + (mrSges(6,2) * pkin(4) - Ifges(4,3) + t797) * t751 - pkin(2) * t613 + (qJ(5) * t689 - t782) * t719 + t783 * t720;
t599 = mrSges(3,2) * t735 - mrSges(3,3) * t722 + Ifges(3,1) * t741 + Ifges(3,4) * t742 + Ifges(3,5) * qJDD(2) - pkin(8) * t613 - qJD(2) * t731 - t601 * t756 + t602 * t760 + t730 * t780;
t598 = -mrSges(3,1) * t735 + mrSges(3,3) * t723 + Ifges(3,4) * t741 + Ifges(3,2) * t742 + Ifges(3,6) * qJDD(2) - pkin(2) * t766 + pkin(8) * t775 + qJD(2) * t732 + t760 * t601 + t756 * t602 - t730 * t781;
t597 = -mrSges(2,2) * g(3) - mrSges(2,3) * t746 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t763 - pkin(7) * t606 - t598 * t757 + t599 * t761;
t1 = [-m(1) * g(1) + t777; -m(1) * g(2) + t787; (-m(1) - m(2)) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t787 + t762 * t597 - t758 * t600; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t777 + t758 * t597 + t762 * t600; -mrSges(1,1) * g(2) + mrSges(2,1) * t746 + mrSges(1,2) * g(1) - mrSges(2,2) * t747 + Ifges(2,3) * qJDD(1) + pkin(1) * t764 + pkin(7) * t776 + t761 * t598 + t757 * t599;];
tauB  = t1;

% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 05:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:52:41
% EndTime: 2019-05-07 05:52:57
% DurationCPUTime: 9.80s
% Computational Cost: add. (144137->363), mult. (293631->441), div. (0->0), fcn. (199651->10), ass. (0->142)
t799 = Ifges(4,1) + Ifges(5,1);
t790 = Ifges(4,4) - Ifges(5,5);
t798 = Ifges(4,5) + Ifges(5,4);
t797 = -Ifges(4,2) - Ifges(5,3);
t796 = Ifges(5,2) + Ifges(4,3);
t795 = -Ifges(5,6) + Ifges(4,6);
t794 = -2 * qJD(4);
t793 = cos(qJ(3));
t762 = cos(qJ(2));
t780 = qJD(1) * t762;
t746 = -qJD(3) + t780;
t792 = pkin(3) * t746;
t791 = -mrSges(4,3) - mrSges(5,2);
t758 = sin(qJ(3));
t759 = sin(qJ(2));
t781 = qJD(1) * t759;
t733 = -qJD(2) * t793 + t758 * t781;
t787 = t733 * t746;
t760 = sin(qJ(1));
t763 = cos(qJ(1));
t743 = -g(1) * t763 - g(2) * t760;
t765 = qJD(1) ^ 2;
t724 = -pkin(1) * t765 + qJDD(1) * pkin(7) + t743;
t714 = -g(3) * t759 + t762 * t724;
t735 = (-mrSges(3,1) * t762 + mrSges(3,2) * t759) * qJD(1);
t779 = qJD(1) * qJD(2);
t748 = t759 * t779;
t738 = qJDD(1) * t762 - t748;
t739 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t781;
t742 = t760 * g(1) - t763 * g(2);
t723 = -qJDD(1) * pkin(1) - t765 * pkin(7) - t742;
t778 = t762 * t779;
t737 = qJDD(1) * t759 + t778;
t673 = (-t737 - t778) * pkin(8) + (-t738 + t748) * pkin(2) + t723;
t736 = (-pkin(2) * t762 - pkin(8) * t759) * qJD(1);
t764 = qJD(2) ^ 2;
t677 = -pkin(2) * t764 + qJDD(2) * pkin(8) + t736 * t780 + t714;
t657 = t758 * t673 + t793 * t677;
t734 = t758 * qJD(2) + t781 * t793;
t698 = qJD(3) * t734 - qJDD(2) * t793 + t737 * t758;
t710 = -mrSges(4,1) * t746 - mrSges(4,3) * t734;
t732 = -qJDD(3) + t738;
t705 = pkin(3) * t733 - qJ(4) * t734;
t745 = t746 ^ 2;
t647 = -pkin(3) * t745 - t732 * qJ(4) - t733 * t705 + t746 * t794 + t657;
t711 = mrSges(5,1) * t746 + mrSges(5,2) * t734;
t656 = t673 * t793 - t758 * t677;
t649 = t732 * pkin(3) - t745 * qJ(4) + t734 * t705 + qJDD(4) - t656;
t699 = -t733 * qJD(3) + t758 * qJDD(2) + t737 * t793;
t637 = (-t699 + t787) * qJ(5) + (t733 * t734 + t732) * pkin(4) + t649;
t709 = pkin(4) * t746 - qJ(5) * t734;
t731 = t733 ^ 2;
t641 = -pkin(4) * t731 + qJ(5) * t698 - t709 * t746 + t647;
t754 = sin(pkin(10));
t755 = cos(pkin(10));
t702 = t733 * t754 + t734 * t755;
t631 = -0.2e1 * qJD(5) * t702 + t755 * t637 - t754 * t641;
t666 = t698 * t754 + t699 * t755;
t701 = t733 * t755 - t734 * t754;
t629 = (t701 * t746 - t666) * pkin(9) + (t701 * t702 + t732) * pkin(5) + t631;
t632 = 0.2e1 * qJD(5) * t701 + t754 * t637 + t755 * t641;
t665 = t698 * t755 - t699 * t754;
t680 = pkin(5) * t746 - pkin(9) * t702;
t700 = t701 ^ 2;
t630 = -pkin(5) * t700 + pkin(9) * t665 - t680 * t746 + t632;
t757 = sin(qJ(6));
t761 = cos(qJ(6));
t627 = t629 * t761 - t630 * t757;
t669 = t701 * t761 - t702 * t757;
t645 = qJD(6) * t669 + t665 * t757 + t666 * t761;
t670 = t701 * t757 + t702 * t761;
t655 = -mrSges(7,1) * t669 + mrSges(7,2) * t670;
t744 = qJD(6) + t746;
t658 = -mrSges(7,2) * t744 + mrSges(7,3) * t669;
t728 = qJDD(6) + t732;
t625 = m(7) * t627 + mrSges(7,1) * t728 - mrSges(7,3) * t645 - t655 * t670 + t658 * t744;
t628 = t629 * t757 + t630 * t761;
t644 = -qJD(6) * t670 + t665 * t761 - t666 * t757;
t659 = mrSges(7,1) * t744 - mrSges(7,3) * t670;
t626 = m(7) * t628 - mrSges(7,2) * t728 + mrSges(7,3) * t644 + t655 * t669 - t659 * t744;
t616 = t761 * t625 + t757 * t626;
t671 = -mrSges(6,1) * t701 + mrSges(6,2) * t702;
t678 = -mrSges(6,2) * t746 + mrSges(6,3) * t701;
t614 = m(6) * t631 + mrSges(6,1) * t732 - mrSges(6,3) * t666 - t671 * t702 + t678 * t746 + t616;
t679 = mrSges(6,1) * t746 - mrSges(6,3) * t702;
t773 = -t625 * t757 + t761 * t626;
t615 = m(6) * t632 - mrSges(6,2) * t732 + mrSges(6,3) * t665 + t671 * t701 - t679 * t746 + t773;
t774 = -t754 * t614 + t755 * t615;
t771 = m(5) * t647 - t732 * mrSges(5,3) - t746 * t711 + t774;
t706 = mrSges(5,1) * t733 - mrSges(5,3) * t734;
t782 = -mrSges(4,1) * t733 - mrSges(4,2) * t734 - t706;
t610 = m(4) * t657 + t732 * mrSges(4,2) + t698 * t791 + t746 * t710 + t733 * t782 + t771;
t708 = mrSges(4,2) * t746 - mrSges(4,3) * t733;
t612 = t755 * t614 + t754 * t615;
t712 = -mrSges(5,2) * t733 - mrSges(5,3) * t746;
t769 = -m(5) * t649 - t732 * mrSges(5,1) - t746 * t712 - t612;
t611 = m(4) * t656 - t732 * mrSges(4,1) + t699 * t791 - t746 * t708 + t734 * t782 + t769;
t775 = t793 * t610 - t611 * t758;
t605 = m(3) * t714 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t738 - qJD(2) * t739 + t735 * t780 + t775;
t713 = -t762 * g(3) - t759 * t724;
t740 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t780;
t676 = -qJDD(2) * pkin(2) - t764 * pkin(8) + t736 * t781 - t713;
t770 = t698 * pkin(3) + t676 + (-t699 - t787) * qJ(4);
t648 = (t794 - t792) * t734 + t770;
t640 = -pkin(4) * t698 - qJ(5) * t731 + qJDD(5) - t770 + ((2 * qJD(4)) + t709 + t792) * t734;
t634 = -pkin(5) * t665 - pkin(9) * t700 + t680 * t702 + t640;
t772 = m(7) * t634 - t644 * mrSges(7,1) + t645 * mrSges(7,2) - t669 * t658 + t670 * t659;
t768 = -m(6) * t640 + t665 * mrSges(6,1) - t666 * mrSges(6,2) + t701 * t678 - t702 * t679 - t772;
t621 = m(5) * t648 + t698 * mrSges(5,1) - t699 * mrSges(5,3) - t734 * t711 + t733 * t712 + t768;
t766 = -m(4) * t676 - t698 * mrSges(4,1) - t699 * mrSges(4,2) - t733 * t708 - t734 * t710 - t621;
t620 = m(3) * t713 + qJDD(2) * mrSges(3,1) - t737 * mrSges(3,3) + qJD(2) * t740 - t735 * t781 + t766;
t776 = t762 * t605 - t620 * t759;
t599 = m(2) * t743 - mrSges(2,1) * t765 - qJDD(1) * mrSges(2,2) + t776;
t606 = t758 * t610 + t611 * t793;
t767 = -m(3) * t723 + t738 * mrSges(3,1) - t737 * mrSges(3,2) - t739 * t781 + t740 * t780 - t606;
t602 = m(2) * t742 + qJDD(1) * mrSges(2,1) - t765 * mrSges(2,2) + t767;
t786 = t760 * t599 + t763 * t602;
t600 = t759 * t605 + t762 * t620;
t785 = t733 * t797 + t734 * t790 - t746 * t795;
t784 = t733 * t795 - t734 * t798 + t746 * t796;
t783 = t790 * t733 - t734 * t799 + t798 * t746;
t777 = t763 * t599 - t602 * t760;
t722 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t759 + Ifges(3,4) * t762) * qJD(1);
t721 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t759 + Ifges(3,2) * t762) * qJD(1);
t720 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t759 + Ifges(3,6) * t762) * qJD(1);
t664 = Ifges(6,1) * t702 + Ifges(6,4) * t701 + Ifges(6,5) * t746;
t663 = Ifges(6,4) * t702 + Ifges(6,2) * t701 + Ifges(6,6) * t746;
t662 = Ifges(6,5) * t702 + Ifges(6,6) * t701 + Ifges(6,3) * t746;
t652 = Ifges(7,1) * t670 + Ifges(7,4) * t669 + Ifges(7,5) * t744;
t651 = Ifges(7,4) * t670 + Ifges(7,2) * t669 + Ifges(7,6) * t744;
t650 = Ifges(7,5) * t670 + Ifges(7,6) * t669 + Ifges(7,3) * t744;
t618 = mrSges(7,2) * t634 - mrSges(7,3) * t627 + Ifges(7,1) * t645 + Ifges(7,4) * t644 + Ifges(7,5) * t728 + t650 * t669 - t651 * t744;
t617 = -mrSges(7,1) * t634 + mrSges(7,3) * t628 + Ifges(7,4) * t645 + Ifges(7,2) * t644 + Ifges(7,6) * t728 - t650 * t670 + t652 * t744;
t608 = mrSges(6,2) * t640 - mrSges(6,3) * t631 + Ifges(6,1) * t666 + Ifges(6,4) * t665 + Ifges(6,5) * t732 - pkin(9) * t616 - t617 * t757 + t618 * t761 + t662 * t701 - t663 * t746;
t607 = -mrSges(6,1) * t640 + mrSges(6,3) * t632 + Ifges(6,4) * t666 + Ifges(6,2) * t665 + Ifges(6,6) * t732 - pkin(5) * t772 + pkin(9) * t773 + t761 * t617 + t757 * t618 - t702 * t662 + t746 * t664;
t596 = mrSges(4,2) * t676 + mrSges(5,2) * t649 - mrSges(4,3) * t656 - mrSges(5,3) * t648 - qJ(4) * t621 - qJ(5) * t612 - t754 * t607 + t755 * t608 - t790 * t698 + t699 * t799 - t732 * t798 + t784 * t733 + t785 * t746;
t595 = -mrSges(4,1) * t676 - mrSges(5,1) * t648 + mrSges(5,2) * t647 + mrSges(4,3) * t657 - pkin(3) * t621 - pkin(4) * t768 - qJ(5) * t774 - t755 * t607 - t754 * t608 + t698 * t797 + t790 * t699 - t732 * t795 + t784 * t734 + t783 * t746;
t594 = -t720 * t781 + (qJ(4) * t706 + t783) * t733 + (pkin(3) * t706 - t785) * t734 - qJ(4) * t771 - pkin(3) * t769 + (Ifges(6,3) + t796) * t732 + (mrSges(5,2) * qJ(4) + t795) * t698 + (mrSges(5,2) * pkin(3) - t798) * t699 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t737 + Ifges(3,2) * t738 + mrSges(3,3) * t714 + qJD(2) * t722 - mrSges(3,1) * t723 + Ifges(7,3) * t728 - t701 * t664 + t702 * t663 + Ifges(6,5) * t666 - t669 * t652 + t670 * t651 - mrSges(4,1) * t656 + mrSges(4,2) * t657 + Ifges(6,6) * t665 + Ifges(7,6) * t644 + Ifges(7,5) * t645 - mrSges(5,3) * t647 + mrSges(5,1) * t649 + mrSges(6,1) * t631 - mrSges(6,2) * t632 + mrSges(7,1) * t627 - mrSges(7,2) * t628 + pkin(5) * t616 + pkin(4) * t612 - pkin(2) * t606;
t593 = mrSges(3,2) * t723 - mrSges(3,3) * t713 + Ifges(3,1) * t737 + Ifges(3,4) * t738 + Ifges(3,5) * qJDD(2) - pkin(8) * t606 - qJD(2) * t721 - t758 * t595 + t596 * t793 + t720 * t780;
t592 = Ifges(2,6) * qJDD(1) + t765 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t743 - Ifges(3,5) * t737 - Ifges(3,6) * t738 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t713 + mrSges(3,2) * t714 - t758 * t596 - t793 * t595 - pkin(2) * t766 - pkin(8) * t775 - pkin(1) * t600 + (-t721 * t759 + t722 * t762) * qJD(1);
t591 = -mrSges(2,2) * g(3) - mrSges(2,3) * t742 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t765 - pkin(7) * t600 + t593 * t762 - t594 * t759;
t1 = [-m(1) * g(1) + t777; -m(1) * g(2) + t786; (-m(1) - m(2)) * g(3) + t600; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t786 + t763 * t591 - t760 * t592; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t777 + t760 * t591 + t763 * t592; -mrSges(1,1) * g(2) + mrSges(2,1) * t742 + mrSges(1,2) * g(1) - mrSges(2,2) * t743 + Ifges(2,3) * qJDD(1) + pkin(1) * t767 + pkin(7) * t776 + t759 * t593 + t762 * t594;];
tauB  = t1;

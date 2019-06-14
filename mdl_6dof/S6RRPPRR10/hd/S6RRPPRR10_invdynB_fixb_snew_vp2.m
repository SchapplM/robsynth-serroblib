% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 11:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:47:24
% EndTime: 2019-05-06 11:47:37
% DurationCPUTime: 9.91s
% Computational Cost: add. (146067->368), mult. (321565->447), div. (0->0), fcn. (204604->10), ass. (0->142)
t774 = -2 * qJD(3);
t773 = Ifges(3,1) + Ifges(4,2);
t769 = Ifges(3,4) + Ifges(4,6);
t768 = Ifges(3,5) - Ifges(4,4);
t772 = Ifges(3,2) + Ifges(4,3);
t767 = Ifges(3,6) - Ifges(4,5);
t771 = (Ifges(3,3) + Ifges(4,1));
t736 = sin(qJ(1));
t740 = cos(qJ(1));
t718 = -t740 * g(1) - t736 * g(2);
t742 = qJD(1) ^ 2;
t694 = -t742 * pkin(1) + qJDD(1) * pkin(7) + t718;
t735 = sin(qJ(2));
t739 = cos(qJ(2));
t673 = -t735 * g(3) + t739 * t694;
t705 = (-pkin(2) * t739 - qJ(3) * t735) * qJD(1);
t741 = qJD(2) ^ 2;
t761 = qJD(1) * t739;
t658 = t741 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t774) - t705 * t761 - t673;
t770 = t742 * pkin(7);
t672 = -t739 * g(3) - t735 * t694;
t706 = (mrSges(4,2) * t739 - mrSges(4,3) * t735) * qJD(1);
t707 = (-mrSges(3,1) * t739 + mrSges(3,2) * t735) * qJD(1);
t760 = qJD(1) * qJD(2);
t757 = t739 * t760;
t708 = t735 * qJDD(1) + t757;
t713 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t761;
t715 = -mrSges(4,1) * t761 - qJD(2) * mrSges(4,3);
t724 = t735 * qJD(1);
t758 = t735 * t760;
t709 = t739 * qJDD(1) - t758;
t714 = pkin(3) * t724 - (qJD(2) * qJ(4));
t730 = t739 ^ 2;
t717 = t736 * g(1) - t740 * g(2);
t752 = -qJDD(1) * pkin(1) - t717;
t747 = pkin(2) * t758 + t724 * t774 + (-t708 - t757) * qJ(3) + t752;
t639 = -t714 * t724 + (-pkin(3) * t730 - pkin(7)) * t742 + (-pkin(2) - qJ(4)) * t709 + t747;
t659 = -qJDD(2) * pkin(2) - t741 * qJ(3) + t705 * t724 + qJDD(3) - t672;
t652 = (-t735 * t739 * t742 - qJDD(2)) * qJ(4) + (t708 - t757) * pkin(3) + t659;
t731 = sin(pkin(10));
t732 = cos(pkin(10));
t700 = t732 * qJD(2) - t731 * t761;
t629 = -0.2e1 * qJD(4) * t700 - t731 * t639 + t732 * t652;
t678 = t732 * qJDD(2) - t731 * t709;
t699 = -t731 * qJD(2) - t732 * t761;
t623 = (t699 * t724 - t678) * pkin(8) + (t699 * t700 + t708) * pkin(4) + t629;
t630 = 0.2e1 * qJD(4) * t699 + t732 * t639 + t731 * t652;
t677 = -t731 * qJDD(2) - t732 * t709;
t679 = pkin(4) * t724 - t700 * pkin(8);
t698 = t699 ^ 2;
t625 = -t698 * pkin(4) + t677 * pkin(8) - t679 * t724 + t630;
t734 = sin(qJ(5));
t738 = cos(qJ(5));
t617 = t738 * t623 - t734 * t625;
t669 = t738 * t699 - t734 * t700;
t644 = t669 * qJD(5) + t734 * t677 + t738 * t678;
t670 = t734 * t699 + t738 * t700;
t704 = qJDD(5) + t708;
t721 = t724 + qJD(5);
t615 = (t669 * t721 - t644) * pkin(9) + (t669 * t670 + t704) * pkin(5) + t617;
t618 = t734 * t623 + t738 * t625;
t643 = -t670 * qJD(5) + t738 * t677 - t734 * t678;
t662 = t721 * pkin(5) - t670 * pkin(9);
t668 = t669 ^ 2;
t616 = -t668 * pkin(5) + t643 * pkin(9) - t721 * t662 + t618;
t733 = sin(qJ(6));
t737 = cos(qJ(6));
t613 = t737 * t615 - t733 * t616;
t654 = t737 * t669 - t733 * t670;
t628 = t654 * qJD(6) + t733 * t643 + t737 * t644;
t655 = t733 * t669 + t737 * t670;
t637 = -t654 * mrSges(7,1) + t655 * mrSges(7,2);
t719 = qJD(6) + t721;
t640 = -t719 * mrSges(7,2) + t654 * mrSges(7,3);
t695 = qJDD(6) + t704;
t608 = m(7) * t613 + t695 * mrSges(7,1) - t628 * mrSges(7,3) - t655 * t637 + t719 * t640;
t614 = t733 * t615 + t737 * t616;
t627 = -t655 * qJD(6) + t737 * t643 - t733 * t644;
t641 = t719 * mrSges(7,1) - t655 * mrSges(7,3);
t609 = m(7) * t614 - t695 * mrSges(7,2) + t627 * mrSges(7,3) + t654 * t637 - t719 * t641;
t602 = t737 * t608 + t733 * t609;
t656 = -t669 * mrSges(6,1) + t670 * mrSges(6,2);
t660 = -t721 * mrSges(6,2) + t669 * mrSges(6,3);
t600 = m(6) * t617 + t704 * mrSges(6,1) - t644 * mrSges(6,3) - t670 * t656 + t721 * t660 + t602;
t661 = t721 * mrSges(6,1) - t670 * mrSges(6,3);
t753 = -t733 * t608 + t737 * t609;
t601 = m(6) * t618 - t704 * mrSges(6,2) + t643 * mrSges(6,3) + t669 * t656 - t721 * t661 + t753;
t596 = t738 * t600 + t734 * t601;
t671 = -t699 * mrSges(5,1) + t700 * mrSges(5,2);
t675 = -mrSges(5,2) * t724 + t699 * mrSges(5,3);
t594 = m(5) * t629 + t708 * mrSges(5,1) - t678 * mrSges(5,3) - t700 * t671 + t675 * t724 + t596;
t676 = mrSges(5,1) * t724 - t700 * mrSges(5,3);
t754 = -t734 * t600 + t738 * t601;
t595 = m(5) * t630 - t708 * mrSges(5,2) + t677 * mrSges(5,3) + t699 * t671 - t676 * t724 + t754;
t589 = t732 * t594 + t731 * t595;
t748 = -m(4) * t659 - t708 * mrSges(4,1) - t589;
t587 = m(3) * t672 - t708 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t713 - t715) * qJD(2) + (-t706 - t707) * t724 + t748;
t712 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t724;
t716 = mrSges(4,1) * t724 + qJD(2) * mrSges(4,2);
t648 = -t730 * t742 * qJ(4) + t709 * pkin(3) + qJD(2) * t714 + qJDD(4) - t658;
t636 = -t677 * pkin(4) - t698 * pkin(8) + t700 * t679 + t648;
t620 = -t643 * pkin(5) - t668 * pkin(9) + t670 * t662 + t636;
t749 = m(7) * t620 - t627 * mrSges(7,1) + t628 * mrSges(7,2) - t654 * t640 + t655 * t641;
t746 = m(6) * t636 - t643 * mrSges(6,1) + t644 * mrSges(6,2) - t669 * t660 + t670 * t661 + t749;
t744 = -m(5) * t648 + t677 * mrSges(5,1) - t678 * mrSges(5,2) + t699 * t675 - t700 * t676 - t746;
t743 = -m(4) * t658 + qJDD(2) * mrSges(4,3) + qJD(2) * t716 + t706 * t761 - t744;
t612 = t707 * t761 - qJDD(2) * mrSges(3,2) + t743 + (mrSges(3,3) + mrSges(4,1)) * t709 + m(3) * t673 - qJD(2) * t712;
t755 = -t735 * t587 + t739 * t612;
t582 = m(2) * t718 - t742 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t755;
t693 = t752 - t770;
t657 = -t709 * pkin(2) + t747 - t770;
t765 = -t731 * t594 + t732 * t595;
t751 = -m(4) * t657 - t709 * mrSges(4,2) + t716 * t724 - t765;
t745 = -m(3) * t693 + t713 * t761 + t709 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t708 + (-t712 * t735 - t715 * t739) * qJD(1) + t751;
t585 = m(2) * t717 + qJDD(1) * mrSges(2,1) - t742 * mrSges(2,2) + t745;
t766 = t736 * t582 + t740 * t585;
t583 = t739 * t587 + t735 * t612;
t764 = (t771 * qJD(2)) + (t735 * t768 + t739 * t767) * qJD(1);
t763 = -t767 * qJD(2) + (-t735 * t769 - t739 * t772) * qJD(1);
t762 = t768 * qJD(2) + (t735 * t773 + t739 * t769) * qJD(1);
t756 = t740 * t582 - t736 * t585;
t665 = Ifges(5,1) * t700 + Ifges(5,4) * t699 + Ifges(5,5) * t724;
t664 = Ifges(5,4) * t700 + Ifges(5,2) * t699 + Ifges(5,6) * t724;
t663 = Ifges(5,5) * t700 + Ifges(5,6) * t699 + Ifges(5,3) * t724;
t651 = Ifges(6,1) * t670 + Ifges(6,4) * t669 + Ifges(6,5) * t721;
t650 = Ifges(6,4) * t670 + Ifges(6,2) * t669 + Ifges(6,6) * t721;
t649 = Ifges(6,5) * t670 + Ifges(6,6) * t669 + Ifges(6,3) * t721;
t633 = Ifges(7,1) * t655 + Ifges(7,4) * t654 + Ifges(7,5) * t719;
t632 = Ifges(7,4) * t655 + Ifges(7,2) * t654 + Ifges(7,6) * t719;
t631 = Ifges(7,5) * t655 + Ifges(7,6) * t654 + Ifges(7,3) * t719;
t604 = mrSges(7,2) * t620 - mrSges(7,3) * t613 + Ifges(7,1) * t628 + Ifges(7,4) * t627 + Ifges(7,5) * t695 + t654 * t631 - t719 * t632;
t603 = -mrSges(7,1) * t620 + mrSges(7,3) * t614 + Ifges(7,4) * t628 + Ifges(7,2) * t627 + Ifges(7,6) * t695 - t655 * t631 + t719 * t633;
t591 = mrSges(6,2) * t636 - mrSges(6,3) * t617 + Ifges(6,1) * t644 + Ifges(6,4) * t643 + Ifges(6,5) * t704 - pkin(9) * t602 - t733 * t603 + t737 * t604 + t669 * t649 - t721 * t650;
t590 = -mrSges(6,1) * t636 + mrSges(6,3) * t618 + Ifges(6,4) * t644 + Ifges(6,2) * t643 + Ifges(6,6) * t704 - pkin(5) * t749 + pkin(9) * t753 + t737 * t603 + t733 * t604 - t670 * t649 + t721 * t651;
t588 = -t708 * mrSges(4,3) + t715 * t761 - t751;
t579 = mrSges(5,2) * t648 - mrSges(5,3) * t629 + Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t708 - pkin(8) * t596 - t734 * t590 + t738 * t591 + t699 * t663 - t664 * t724;
t578 = -mrSges(5,1) * t648 + mrSges(5,3) * t630 + Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t708 - pkin(4) * t746 + pkin(8) * t754 + t738 * t590 + t734 * t591 - t700 * t663 + t665 * t724;
t577 = (Ifges(5,3) + t773) * t708 + pkin(3) * t589 - qJ(3) * t588 + pkin(5) * t602 + pkin(4) * t596 + t763 * qJD(2) + t764 * t761 + t768 * qJDD(2) + t769 * t709 + mrSges(7,1) * t613 - mrSges(7,2) * t614 + mrSges(6,1) * t617 - mrSges(6,2) * t618 + Ifges(7,6) * t627 + Ifges(7,5) * t628 + mrSges(5,1) * t629 - mrSges(5,2) * t630 + Ifges(6,6) * t643 + Ifges(6,5) * t644 - t654 * t633 + t655 * t632 - mrSges(4,3) * t657 + mrSges(4,1) * t659 - t669 * t651 + t670 * t650 - mrSges(3,3) * t672 + Ifges(5,6) * t677 + Ifges(5,5) * t678 + mrSges(3,2) * t693 + Ifges(7,3) * t695 - t699 * t665 + t700 * t664 + Ifges(6,3) * t704;
t576 = -mrSges(3,1) * t693 - mrSges(4,1) * t658 + mrSges(4,2) * t657 + mrSges(3,3) * t673 - pkin(2) * t588 - pkin(3) * t744 - qJ(4) * t765 + t762 * qJD(2) + t767 * qJDD(2) - t732 * t578 - t731 * t579 + t769 * t708 + t709 * t772 - t764 * t724;
t575 = -pkin(1) * t583 + mrSges(2,3) * t718 - qJ(3) * t743 - pkin(2) * (-qJD(2) * t715 + t748) - t732 * t579 + t731 * t578 + qJ(4) * t589 - mrSges(3,1) * t672 + mrSges(3,2) * t673 - mrSges(4,2) * t659 + mrSges(4,3) * t658 + mrSges(2,1) * g(3) + t742 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-qJ(3) * mrSges(4,1) - t767) * t709 - t768 * t708 + (pkin(2) * mrSges(4,2) - t771) * qJDD(2) + (t762 * t739 + (pkin(2) * t706 + t763) * t735) * qJD(1);
t574 = -mrSges(2,2) * g(3) - mrSges(2,3) * t717 + Ifges(2,5) * qJDD(1) - t742 * Ifges(2,6) - pkin(7) * t583 - t735 * t576 + t739 * t577;
t1 = [-m(1) * g(1) + t756; -m(1) * g(2) + t766; (-m(1) - m(2)) * g(3) + t583; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t766 + t740 * t574 - t736 * t575; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t756 + t736 * t574 + t740 * t575; -mrSges(1,1) * g(2) + mrSges(2,1) * t717 + mrSges(1,2) * g(1) - mrSges(2,2) * t718 + Ifges(2,3) * qJDD(1) + pkin(1) * t745 + pkin(7) * t755 + t739 * t576 + t735 * t577;];
tauB  = t1;

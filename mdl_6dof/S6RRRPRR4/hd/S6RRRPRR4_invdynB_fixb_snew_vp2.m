% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 10:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:26:13
% EndTime: 2019-05-07 10:26:45
% DurationCPUTime: 29.19s
% Computational Cost: add. (486989->387), mult. (1004746->489), div. (0->0), fcn. (740876->12), ass. (0->149)
t756 = cos(qJ(3));
t737 = qJD(1) ^ 2;
t755 = pkin(2) * t737;
t732 = sin(qJ(1));
t736 = cos(qJ(1));
t719 = -g(1) * t736 - g(2) * t732;
t707 = -pkin(1) * t737 + qJDD(1) * pkin(7) + t719;
t731 = sin(qJ(2));
t754 = t731 * t707;
t735 = cos(qJ(2));
t750 = qJD(1) * qJD(2);
t713 = qJDD(1) * t731 + t735 * t750;
t671 = qJDD(2) * pkin(2) - t713 * pkin(8) - t754 + (pkin(8) * t750 + t731 * t755 - g(3)) * t735;
t695 = -g(3) * t731 + t735 * t707;
t714 = qJDD(1) * t735 - t731 * t750;
t752 = qJD(1) * t731;
t717 = qJD(2) * pkin(2) - pkin(8) * t752;
t725 = t735 ^ 2;
t672 = pkin(8) * t714 - qJD(2) * t717 - t725 * t755 + t695;
t730 = sin(qJ(3));
t652 = t730 * t671 + t672 * t756;
t705 = (t730 * t735 + t731 * t756) * qJD(1);
t677 = qJD(3) * t705 + t713 * t730 - t714 * t756;
t751 = qJD(1) * t735;
t704 = t730 * t752 - t751 * t756;
t688 = mrSges(4,1) * t704 + mrSges(4,2) * t705;
t724 = qJD(2) + qJD(3);
t697 = mrSges(4,1) * t724 - mrSges(4,3) * t705;
t723 = qJDD(2) + qJDD(3);
t678 = -t704 * qJD(3) + t713 * t756 + t730 * t714;
t718 = t732 * g(1) - t736 * g(2);
t742 = -qJDD(1) * pkin(1) - t718;
t679 = -t714 * pkin(2) + t717 * t752 + (-pkin(8) * t725 - pkin(7)) * t737 + t742;
t638 = (t704 * t724 - t678) * qJ(4) + (t705 * t724 + t677) * pkin(3) + t679;
t687 = pkin(3) * t704 - qJ(4) * t705;
t722 = t724 ^ 2;
t641 = -pkin(3) * t722 + qJ(4) * t723 - t687 * t704 + t652;
t726 = sin(pkin(11));
t727 = cos(pkin(11));
t693 = t705 * t727 + t724 * t726;
t624 = -0.2e1 * qJD(4) * t693 + t727 * t638 - t726 * t641;
t666 = t678 * t727 + t723 * t726;
t692 = -t705 * t726 + t724 * t727;
t617 = (t692 * t704 - t666) * pkin(9) + (t692 * t693 + t677) * pkin(4) + t624;
t625 = 0.2e1 * qJD(4) * t692 + t726 * t638 + t727 * t641;
t665 = -t678 * t726 + t723 * t727;
t682 = pkin(4) * t704 - pkin(9) * t693;
t691 = t692 ^ 2;
t619 = -pkin(4) * t691 + pkin(9) * t665 - t682 * t704 + t625;
t729 = sin(qJ(5));
t734 = cos(qJ(5));
t611 = t734 * t617 - t729 * t619;
t663 = t692 * t734 - t693 * t729;
t637 = qJD(5) * t663 + t665 * t729 + t666 * t734;
t664 = t692 * t729 + t693 * t734;
t676 = qJDD(5) + t677;
t700 = qJD(5) + t704;
t609 = (t663 * t700 - t637) * pkin(10) + (t663 * t664 + t676) * pkin(5) + t611;
t612 = t729 * t617 + t734 * t619;
t636 = -qJD(5) * t664 + t665 * t734 - t666 * t729;
t655 = pkin(5) * t700 - pkin(10) * t664;
t662 = t663 ^ 2;
t610 = -pkin(5) * t662 + pkin(10) * t636 - t655 * t700 + t612;
t728 = sin(qJ(6));
t733 = cos(qJ(6));
t607 = t609 * t733 - t610 * t728;
t648 = t663 * t733 - t664 * t728;
t623 = qJD(6) * t648 + t636 * t728 + t637 * t733;
t649 = t663 * t728 + t664 * t733;
t632 = -mrSges(7,1) * t648 + mrSges(7,2) * t649;
t699 = qJD(6) + t700;
t642 = -mrSges(7,2) * t699 + mrSges(7,3) * t648;
t673 = qJDD(6) + t676;
t603 = m(7) * t607 + mrSges(7,1) * t673 - mrSges(7,3) * t623 - t632 * t649 + t642 * t699;
t608 = t609 * t728 + t610 * t733;
t622 = -qJD(6) * t649 + t636 * t733 - t637 * t728;
t643 = mrSges(7,1) * t699 - mrSges(7,3) * t649;
t604 = m(7) * t608 - mrSges(7,2) * t673 + mrSges(7,3) * t622 + t632 * t648 - t643 * t699;
t597 = t733 * t603 + t728 * t604;
t650 = -mrSges(6,1) * t663 + mrSges(6,2) * t664;
t653 = -mrSges(6,2) * t700 + mrSges(6,3) * t663;
t595 = m(6) * t611 + mrSges(6,1) * t676 - mrSges(6,3) * t637 - t650 * t664 + t653 * t700 + t597;
t654 = mrSges(6,1) * t700 - mrSges(6,3) * t664;
t744 = -t603 * t728 + t733 * t604;
t596 = m(6) * t612 - mrSges(6,2) * t676 + mrSges(6,3) * t636 + t650 * t663 - t654 * t700 + t744;
t591 = t734 * t595 + t729 * t596;
t668 = -mrSges(5,1) * t692 + mrSges(5,2) * t693;
t680 = -mrSges(5,2) * t704 + mrSges(5,3) * t692;
t589 = m(5) * t624 + mrSges(5,1) * t677 - mrSges(5,3) * t666 - t668 * t693 + t680 * t704 + t591;
t681 = mrSges(5,1) * t704 - mrSges(5,3) * t693;
t745 = -t595 * t729 + t734 * t596;
t590 = m(5) * t625 - mrSges(5,2) * t677 + mrSges(5,3) * t665 + t668 * t692 - t681 * t704 + t745;
t746 = -t589 * t726 + t727 * t590;
t582 = m(4) * t652 - mrSges(4,2) * t723 - mrSges(4,3) * t677 - t688 * t704 - t697 * t724 + t746;
t651 = t671 * t756 - t730 * t672;
t696 = -mrSges(4,2) * t724 - mrSges(4,3) * t704;
t640 = -t723 * pkin(3) - t722 * qJ(4) + t705 * t687 + qJDD(4) - t651;
t626 = -t665 * pkin(4) - t691 * pkin(9) + t693 * t682 + t640;
t614 = -t636 * pkin(5) - t662 * pkin(10) + t664 * t655 + t626;
t743 = m(7) * t614 - t622 * mrSges(7,1) + t623 * mrSges(7,2) - t648 * t642 + t649 * t643;
t740 = m(6) * t626 - t636 * mrSges(6,1) + mrSges(6,2) * t637 - t663 * t653 + t654 * t664 + t743;
t738 = -m(5) * t640 + t665 * mrSges(5,1) - mrSges(5,2) * t666 + t692 * t680 - t681 * t693 - t740;
t606 = m(4) * t651 + mrSges(4,1) * t723 - mrSges(4,3) * t678 - t688 * t705 + t696 * t724 + t738;
t577 = t730 * t582 + t606 * t756;
t694 = -t735 * g(3) - t754;
t712 = (-mrSges(3,1) * t735 + mrSges(3,2) * t731) * qJD(1);
t716 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t751;
t575 = m(3) * t694 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t713 + qJD(2) * t716 - t712 * t752 + t577;
t715 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t752;
t747 = t582 * t756 - t606 * t730;
t576 = m(3) * t695 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t714 - qJD(2) * t715 + t712 * t751 + t747;
t748 = -t575 * t731 + t735 * t576;
t567 = m(2) * t719 - mrSges(2,1) * t737 - qJDD(1) * mrSges(2,2) + t748;
t706 = -t737 * pkin(7) + t742;
t583 = t727 * t589 + t726 * t590;
t741 = m(4) * t679 + t677 * mrSges(4,1) + mrSges(4,2) * t678 + t704 * t696 + t697 * t705 + t583;
t739 = -m(3) * t706 + t714 * mrSges(3,1) - mrSges(3,2) * t713 - t715 * t752 + t716 * t751 - t741;
t579 = m(2) * t718 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t737 + t739;
t753 = t732 * t567 + t736 * t579;
t568 = t735 * t575 + t731 * t576;
t749 = t736 * t567 - t579 * t732;
t703 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t731 + Ifges(3,4) * t735) * qJD(1);
t702 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t731 + Ifges(3,2) * t735) * qJD(1);
t701 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t731 + Ifges(3,6) * t735) * qJD(1);
t685 = Ifges(4,1) * t705 - Ifges(4,4) * t704 + Ifges(4,5) * t724;
t684 = Ifges(4,4) * t705 - Ifges(4,2) * t704 + Ifges(4,6) * t724;
t683 = Ifges(4,5) * t705 - Ifges(4,6) * t704 + Ifges(4,3) * t724;
t658 = Ifges(5,1) * t693 + Ifges(5,4) * t692 + Ifges(5,5) * t704;
t657 = Ifges(5,4) * t693 + Ifges(5,2) * t692 + Ifges(5,6) * t704;
t656 = Ifges(5,5) * t693 + Ifges(5,6) * t692 + Ifges(5,3) * t704;
t646 = Ifges(6,1) * t664 + Ifges(6,4) * t663 + Ifges(6,5) * t700;
t645 = Ifges(6,4) * t664 + Ifges(6,2) * t663 + Ifges(6,6) * t700;
t644 = Ifges(6,5) * t664 + Ifges(6,6) * t663 + Ifges(6,3) * t700;
t629 = Ifges(7,1) * t649 + Ifges(7,4) * t648 + Ifges(7,5) * t699;
t628 = Ifges(7,4) * t649 + Ifges(7,2) * t648 + Ifges(7,6) * t699;
t627 = Ifges(7,5) * t649 + Ifges(7,6) * t648 + Ifges(7,3) * t699;
t599 = mrSges(7,2) * t614 - mrSges(7,3) * t607 + Ifges(7,1) * t623 + Ifges(7,4) * t622 + Ifges(7,5) * t673 + t627 * t648 - t628 * t699;
t598 = -mrSges(7,1) * t614 + mrSges(7,3) * t608 + Ifges(7,4) * t623 + Ifges(7,2) * t622 + Ifges(7,6) * t673 - t627 * t649 + t629 * t699;
t585 = mrSges(6,2) * t626 - mrSges(6,3) * t611 + Ifges(6,1) * t637 + Ifges(6,4) * t636 + Ifges(6,5) * t676 - pkin(10) * t597 - t598 * t728 + t599 * t733 + t644 * t663 - t645 * t700;
t584 = -mrSges(6,1) * t626 + mrSges(6,3) * t612 + Ifges(6,4) * t637 + Ifges(6,2) * t636 + Ifges(6,6) * t676 - pkin(5) * t743 + pkin(10) * t744 + t733 * t598 + t728 * t599 - t664 * t644 + t700 * t646;
t571 = mrSges(5,2) * t640 - mrSges(5,3) * t624 + Ifges(5,1) * t666 + Ifges(5,4) * t665 + Ifges(5,5) * t677 - pkin(9) * t591 - t584 * t729 + t585 * t734 + t656 * t692 - t657 * t704;
t570 = -mrSges(5,1) * t640 + mrSges(5,3) * t625 + Ifges(5,4) * t666 + Ifges(5,2) * t665 + Ifges(5,6) * t677 - pkin(4) * t740 + pkin(9) * t745 + t734 * t584 + t729 * t585 - t693 * t656 + t704 * t658;
t569 = (-Ifges(5,3) - Ifges(4,2)) * t677 + Ifges(4,6) * t723 + t724 * t685 - t705 * t683 - t693 * t657 + t692 * t658 - Ifges(7,3) * t673 - Ifges(6,3) * t676 + Ifges(4,4) * t678 - mrSges(4,1) * t679 - t664 * t645 - Ifges(5,6) * t665 - Ifges(5,5) * t666 + t663 * t646 + t648 * t629 - t649 * t628 + mrSges(4,3) * t652 - Ifges(6,6) * t636 - Ifges(6,5) * t637 - mrSges(5,1) * t624 + mrSges(5,2) * t625 - Ifges(7,6) * t622 - Ifges(7,5) * t623 + mrSges(6,2) * t612 - mrSges(6,1) * t611 + mrSges(7,2) * t608 - mrSges(7,1) * t607 - pkin(5) * t597 - pkin(4) * t591 - pkin(3) * t583;
t564 = mrSges(4,2) * t679 - mrSges(4,3) * t651 + Ifges(4,1) * t678 - Ifges(4,4) * t677 + Ifges(4,5) * t723 - qJ(4) * t583 - t570 * t726 + t571 * t727 - t683 * t704 - t684 * t724;
t563 = mrSges(3,2) * t706 - mrSges(3,3) * t694 + Ifges(3,1) * t713 + Ifges(3,4) * t714 + Ifges(3,5) * qJDD(2) - pkin(8) * t577 - qJD(2) * t702 + t564 * t756 - t730 * t569 + t701 * t751;
t562 = -pkin(1) * t568 + mrSges(2,3) * t719 - pkin(2) * t577 - Ifges(3,5) * t713 - Ifges(3,6) * t714 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t694 + mrSges(3,2) * t695 - t726 * t571 - t727 * t570 - pkin(3) * t738 - qJ(4) * t746 - Ifges(4,5) * t678 + Ifges(4,6) * t677 - Ifges(4,3) * t723 - mrSges(4,1) * t651 + mrSges(4,2) * t652 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + t737 * Ifges(2,5) - t705 * t684 - t704 * t685 + (-t702 * t731 + t703 * t735) * qJD(1);
t561 = -mrSges(3,1) * t706 + mrSges(3,3) * t695 + Ifges(3,4) * t713 + Ifges(3,2) * t714 + Ifges(3,6) * qJDD(2) - pkin(2) * t741 + pkin(8) * t747 + qJD(2) * t703 + t730 * t564 + t569 * t756 - t701 * t752;
t560 = -mrSges(2,2) * g(3) - mrSges(2,3) * t718 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t737 - pkin(7) * t568 - t561 * t731 + t563 * t735;
t1 = [-m(1) * g(1) + t749; -m(1) * g(2) + t753; (-m(1) - m(2)) * g(3) + t568; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t753 + t736 * t560 - t732 * t562; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t749 + t732 * t560 + t736 * t562; -mrSges(1,1) * g(2) + mrSges(2,1) * t718 + mrSges(1,2) * g(1) - mrSges(2,2) * t719 + Ifges(2,3) * qJDD(1) + pkin(1) * t739 + pkin(7) * t748 + t735 * t561 + t731 * t563;];
tauB  = t1;

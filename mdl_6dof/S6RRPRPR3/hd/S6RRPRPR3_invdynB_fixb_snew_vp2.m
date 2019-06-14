% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:15:47
% EndTime: 2019-05-06 13:16:15
% DurationCPUTime: 26.29s
% Computational Cost: add. (406411->386), mult. (945230->490), div. (0->0), fcn. (690904->12), ass. (0->147)
t729 = sin(qJ(2));
t733 = cos(qJ(2));
t750 = qJD(1) * qJD(2);
t713 = qJDD(1) * t729 + t733 * t750;
t730 = sin(qJ(1));
t734 = cos(qJ(1));
t719 = -g(1) * t734 - g(2) * t730;
t736 = qJD(1) ^ 2;
t708 = -pkin(1) * t736 + qJDD(1) * pkin(7) + t719;
t755 = t729 * t708;
t756 = pkin(2) * t736;
t670 = qJDD(2) * pkin(2) - t713 * qJ(3) - t755 + (qJ(3) * t750 + t729 * t756 - g(3)) * t733;
t693 = -g(3) * t729 + t733 * t708;
t714 = qJDD(1) * t733 - t729 * t750;
t753 = qJD(1) * t729;
t715 = qJD(2) * pkin(2) - qJ(3) * t753;
t722 = t733 ^ 2;
t671 = qJ(3) * t714 - qJD(2) * t715 - t722 * t756 + t693;
t724 = sin(pkin(10));
t726 = cos(pkin(10));
t702 = (t724 * t733 + t726 * t729) * qJD(1);
t646 = -0.2e1 * qJD(3) * t702 + t670 * t726 - t724 * t671;
t752 = qJD(1) * t733;
t701 = -t724 * t753 + t726 * t752;
t647 = 0.2e1 * qJD(3) * t701 + t724 * t670 + t726 * t671;
t682 = -mrSges(4,1) * t701 + mrSges(4,2) * t702;
t687 = -t724 * t713 + t714 * t726;
t695 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t702;
t683 = -pkin(3) * t701 - pkin(8) * t702;
t735 = qJD(2) ^ 2;
t633 = -pkin(3) * t735 + qJDD(2) * pkin(8) + t683 * t701 + t647;
t718 = t730 * g(1) - t734 * g(2);
t741 = -qJDD(1) * pkin(1) - t718;
t674 = -t714 * pkin(2) + qJDD(3) + t715 * t753 + (-qJ(3) * t722 - pkin(7)) * t736 + t741;
t688 = t713 * t726 + t714 * t724;
t636 = (-qJD(2) * t701 - t688) * pkin(8) + (qJD(2) * t702 - t687) * pkin(3) + t674;
t728 = sin(qJ(4));
t732 = cos(qJ(4));
t623 = -t728 * t633 + t732 * t636;
t690 = qJD(2) * t732 - t702 * t728;
t661 = qJD(4) * t690 + qJDD(2) * t728 + t688 * t732;
t686 = qJDD(4) - t687;
t691 = qJD(2) * t728 + t702 * t732;
t700 = qJD(4) - t701;
t615 = (t690 * t700 - t661) * qJ(5) + (t690 * t691 + t686) * pkin(4) + t623;
t624 = t732 * t633 + t728 * t636;
t660 = -qJD(4) * t691 + qJDD(2) * t732 - t688 * t728;
t676 = pkin(4) * t700 - qJ(5) * t691;
t689 = t690 ^ 2;
t617 = -pkin(4) * t689 + qJ(5) * t660 - t676 * t700 + t624;
t723 = sin(pkin(11));
t725 = cos(pkin(11));
t669 = t690 * t723 + t691 * t725;
t609 = -0.2e1 * qJD(5) * t669 + t725 * t615 - t723 * t617;
t641 = t660 * t723 + t661 * t725;
t668 = t690 * t725 - t691 * t723;
t607 = (t668 * t700 - t641) * pkin(9) + (t668 * t669 + t686) * pkin(5) + t609;
t610 = 0.2e1 * qJD(5) * t668 + t723 * t615 + t725 * t617;
t640 = t660 * t725 - t661 * t723;
t653 = pkin(5) * t700 - pkin(9) * t669;
t665 = t668 ^ 2;
t608 = -pkin(5) * t665 + pkin(9) * t640 - t653 * t700 + t610;
t727 = sin(qJ(6));
t731 = cos(qJ(6));
t605 = t607 * t731 - t608 * t727;
t648 = t668 * t731 - t669 * t727;
t621 = qJD(6) * t648 + t640 * t727 + t641 * t731;
t649 = t668 * t727 + t669 * t731;
t630 = -mrSges(7,1) * t648 + mrSges(7,2) * t649;
t696 = qJD(6) + t700;
t637 = -mrSges(7,2) * t696 + mrSges(7,3) * t648;
t684 = qJDD(6) + t686;
t601 = m(7) * t605 + mrSges(7,1) * t684 - mrSges(7,3) * t621 - t630 * t649 + t637 * t696;
t606 = t607 * t727 + t608 * t731;
t620 = -qJD(6) * t649 + t640 * t731 - t641 * t727;
t638 = mrSges(7,1) * t696 - mrSges(7,3) * t649;
t602 = m(7) * t606 - mrSges(7,2) * t684 + mrSges(7,3) * t620 + t630 * t648 - t638 * t696;
t595 = t731 * t601 + t727 * t602;
t650 = -mrSges(6,1) * t668 + mrSges(6,2) * t669;
t651 = -mrSges(6,2) * t700 + mrSges(6,3) * t668;
t593 = m(6) * t609 + mrSges(6,1) * t686 - mrSges(6,3) * t641 - t650 * t669 + t651 * t700 + t595;
t652 = mrSges(6,1) * t700 - mrSges(6,3) * t669;
t744 = -t601 * t727 + t731 * t602;
t594 = m(6) * t610 - mrSges(6,2) * t686 + mrSges(6,3) * t640 + t650 * t668 - t652 * t700 + t744;
t589 = t725 * t593 + t723 * t594;
t672 = -mrSges(5,1) * t690 + mrSges(5,2) * t691;
t675 = -mrSges(5,2) * t700 + mrSges(5,3) * t690;
t587 = m(5) * t623 + mrSges(5,1) * t686 - mrSges(5,3) * t661 - t672 * t691 + t675 * t700 + t589;
t677 = mrSges(5,1) * t700 - mrSges(5,3) * t691;
t745 = -t593 * t723 + t725 * t594;
t588 = m(5) * t624 - mrSges(5,2) * t686 + mrSges(5,3) * t660 + t672 * t690 - t677 * t700 + t745;
t746 = -t587 * t728 + t732 * t588;
t580 = m(4) * t647 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t687 - qJD(2) * t695 + t682 * t701 + t746;
t694 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t701;
t632 = -qJDD(2) * pkin(3) - pkin(8) * t735 + t702 * t683 - t646;
t622 = -pkin(4) * t660 - qJ(5) * t689 + t691 * t676 + qJDD(5) + t632;
t612 = -pkin(5) * t640 - pkin(9) * t665 + t653 * t669 + t622;
t742 = m(7) * t612 - t620 * mrSges(7,1) + t621 * mrSges(7,2) - t648 * t637 + t649 * t638;
t739 = m(6) * t622 - t640 * mrSges(6,1) + mrSges(6,2) * t641 - t668 * t651 + t652 * t669 + t742;
t737 = -m(5) * t632 + t660 * mrSges(5,1) - mrSges(5,2) * t661 + t690 * t675 - t677 * t691 - t739;
t604 = m(4) * t646 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t688 + qJD(2) * t694 - t682 * t702 + t737;
t575 = t724 * t580 + t726 * t604;
t692 = -t733 * g(3) - t755;
t712 = (-mrSges(3,1) * t733 + mrSges(3,2) * t729) * qJD(1);
t717 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t752;
t573 = m(3) * t692 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t713 + qJD(2) * t717 - t712 * t753 + t575;
t716 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t753;
t747 = t726 * t580 - t604 * t724;
t574 = m(3) * t693 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t714 - qJD(2) * t716 + t712 * t752 + t747;
t748 = -t573 * t729 + t733 * t574;
t565 = m(2) * t719 - mrSges(2,1) * t736 - qJDD(1) * mrSges(2,2) + t748;
t707 = -t736 * pkin(7) + t741;
t581 = t732 * t587 + t728 * t588;
t740 = m(4) * t674 - t687 * mrSges(4,1) + mrSges(4,2) * t688 - t701 * t694 + t695 * t702 + t581;
t738 = -m(3) * t707 + t714 * mrSges(3,1) - mrSges(3,2) * t713 - t716 * t753 + t717 * t752 - t740;
t577 = m(2) * t718 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t736 + t738;
t754 = t730 * t565 + t734 * t577;
t566 = t733 * t573 + t729 * t574;
t749 = t734 * t565 - t577 * t730;
t705 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t729 + Ifges(3,4) * t733) * qJD(1);
t704 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t729 + Ifges(3,2) * t733) * qJD(1);
t703 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t729 + Ifges(3,6) * t733) * qJD(1);
t680 = Ifges(4,1) * t702 + Ifges(4,4) * t701 + Ifges(4,5) * qJD(2);
t679 = Ifges(4,4) * t702 + Ifges(4,2) * t701 + Ifges(4,6) * qJD(2);
t678 = Ifges(4,5) * t702 + Ifges(4,6) * t701 + Ifges(4,3) * qJD(2);
t656 = Ifges(5,1) * t691 + Ifges(5,4) * t690 + Ifges(5,5) * t700;
t655 = Ifges(5,4) * t691 + Ifges(5,2) * t690 + Ifges(5,6) * t700;
t654 = Ifges(5,5) * t691 + Ifges(5,6) * t690 + Ifges(5,3) * t700;
t644 = Ifges(6,1) * t669 + Ifges(6,4) * t668 + Ifges(6,5) * t700;
t643 = Ifges(6,4) * t669 + Ifges(6,2) * t668 + Ifges(6,6) * t700;
t642 = Ifges(6,5) * t669 + Ifges(6,6) * t668 + Ifges(6,3) * t700;
t627 = Ifges(7,1) * t649 + Ifges(7,4) * t648 + Ifges(7,5) * t696;
t626 = Ifges(7,4) * t649 + Ifges(7,2) * t648 + Ifges(7,6) * t696;
t625 = Ifges(7,5) * t649 + Ifges(7,6) * t648 + Ifges(7,3) * t696;
t597 = mrSges(7,2) * t612 - mrSges(7,3) * t605 + Ifges(7,1) * t621 + Ifges(7,4) * t620 + Ifges(7,5) * t684 + t625 * t648 - t626 * t696;
t596 = -mrSges(7,1) * t612 + mrSges(7,3) * t606 + Ifges(7,4) * t621 + Ifges(7,2) * t620 + Ifges(7,6) * t684 - t625 * t649 + t627 * t696;
t583 = mrSges(6,2) * t622 - mrSges(6,3) * t609 + Ifges(6,1) * t641 + Ifges(6,4) * t640 + Ifges(6,5) * t686 - pkin(9) * t595 - t596 * t727 + t597 * t731 + t642 * t668 - t643 * t700;
t582 = -mrSges(6,1) * t622 + mrSges(6,3) * t610 + Ifges(6,4) * t641 + Ifges(6,2) * t640 + Ifges(6,6) * t686 - pkin(5) * t742 + pkin(9) * t744 + t731 * t596 + t727 * t597 - t669 * t642 + t700 * t644;
t569 = mrSges(5,2) * t632 - mrSges(5,3) * t623 + Ifges(5,1) * t661 + Ifges(5,4) * t660 + Ifges(5,5) * t686 - qJ(5) * t589 - t582 * t723 + t583 * t725 + t654 * t690 - t655 * t700;
t568 = -mrSges(5,1) * t632 + mrSges(5,3) * t624 + Ifges(5,4) * t661 + Ifges(5,2) * t660 + Ifges(5,6) * t686 - pkin(4) * t739 + qJ(5) * t745 + t725 * t582 + t723 * t583 - t691 * t654 + t700 * t656;
t567 = Ifges(4,6) * qJDD(2) + (-Ifges(5,3) - Ifges(6,3)) * t686 - t702 * t678 + t690 * t656 - t691 * t655 - Ifges(7,3) * t684 + Ifges(4,2) * t687 + Ifges(4,4) * t688 + t668 * t644 - t669 * t643 - mrSges(4,1) * t674 + qJD(2) * t680 - Ifges(5,6) * t660 - Ifges(5,5) * t661 + t648 * t627 - t649 * t626 - Ifges(6,6) * t640 - Ifges(6,5) * t641 + mrSges(4,3) * t647 - mrSges(5,1) * t623 + mrSges(5,2) * t624 - Ifges(7,6) * t620 - Ifges(7,5) * t621 - mrSges(6,1) * t609 + mrSges(6,2) * t610 + mrSges(7,2) * t606 - mrSges(7,1) * t605 - pkin(5) * t595 - pkin(4) * t589 - pkin(3) * t581;
t562 = mrSges(4,2) * t674 - mrSges(4,3) * t646 + Ifges(4,1) * t688 + Ifges(4,4) * t687 + Ifges(4,5) * qJDD(2) - pkin(8) * t581 - qJD(2) * t679 - t568 * t728 + t569 * t732 + t678 * t701;
t561 = mrSges(3,2) * t707 - mrSges(3,3) * t692 + Ifges(3,1) * t713 + Ifges(3,4) * t714 + Ifges(3,5) * qJDD(2) - qJ(3) * t575 - qJD(2) * t704 + t562 * t726 - t567 * t724 + t703 * t752;
t560 = -pkin(1) * t566 + mrSges(2,3) * t719 - pkin(2) * t575 - Ifges(3,5) * t713 - Ifges(3,6) * t714 - mrSges(3,1) * t692 + mrSges(3,2) * t693 - pkin(8) * t746 - t728 * t569 - t732 * t568 - pkin(3) * t737 - Ifges(4,5) * t688 - Ifges(4,6) * t687 - mrSges(4,1) * t646 + mrSges(4,2) * t647 + mrSges(2,1) * g(3) + t736 * Ifges(2,5) + t701 * t680 - t702 * t679 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t704 * t729 + t705 * t733) * qJD(1);
t559 = -mrSges(3,1) * t707 + mrSges(3,3) * t693 + Ifges(3,4) * t713 + Ifges(3,2) * t714 + Ifges(3,6) * qJDD(2) - pkin(2) * t740 + qJ(3) * t747 + qJD(2) * t705 + t724 * t562 + t726 * t567 - t703 * t753;
t558 = -mrSges(2,2) * g(3) - mrSges(2,3) * t718 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t736 - pkin(7) * t566 - t559 * t729 + t561 * t733;
t1 = [-m(1) * g(1) + t749; -m(1) * g(2) + t754; (-m(1) - m(2)) * g(3) + t566; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t754 + t734 * t558 - t730 * t560; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t749 + t730 * t558 + t734 * t560; -mrSges(1,1) * g(2) + mrSges(2,1) * t718 + mrSges(1,2) * g(1) - mrSges(2,2) * t719 + Ifges(2,3) * qJDD(1) + pkin(1) * t738 + pkin(7) * t748 + t733 * t559 + t729 * t561;];
tauB  = t1;

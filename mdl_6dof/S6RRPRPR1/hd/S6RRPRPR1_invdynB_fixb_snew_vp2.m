% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR1
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
% Datum: 2019-05-06 12:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:48:39
% EndTime: 2019-05-06 12:49:06
% DurationCPUTime: 26.54s
% Computational Cost: add. (395253->386), mult. (928568->491), div. (0->0), fcn. (697256->12), ass. (0->147)
t748 = cos(qJ(4));
t729 = qJD(1) ^ 2;
t747 = pkin(2) * t729;
t725 = sin(qJ(1));
t728 = cos(qJ(1));
t712 = -g(1) * t728 - g(2) * t725;
t701 = -pkin(1) * t729 + qJDD(1) * pkin(7) + t712;
t724 = sin(qJ(2));
t746 = t724 * t701;
t727 = cos(qJ(2));
t742 = qJD(1) * qJD(2);
t706 = qJDD(1) * t724 + t727 * t742;
t662 = qJDD(2) * pkin(2) - t706 * qJ(3) - t746 + (qJ(3) * t742 + t724 * t747 - g(3)) * t727;
t686 = -g(3) * t724 + t727 * t701;
t707 = qJDD(1) * t727 - t724 * t742;
t744 = qJD(1) * t724;
t708 = qJD(2) * pkin(2) - qJ(3) * t744;
t717 = t727 ^ 2;
t663 = qJ(3) * t707 - qJD(2) * t708 - t717 * t747 + t686;
t719 = sin(pkin(10));
t721 = cos(pkin(10));
t696 = (t719 * t727 + t721 * t724) * qJD(1);
t635 = -0.2e1 * qJD(3) * t696 + t721 * t662 - t719 * t663;
t684 = t706 * t721 + t707 * t719;
t695 = (-t719 * t724 + t721 * t727) * qJD(1);
t622 = (qJD(2) * t695 - t684) * pkin(8) + (t695 * t696 + qJDD(2)) * pkin(3) + t635;
t636 = 0.2e1 * qJD(3) * t695 + t719 * t662 + t721 * t663;
t683 = -t706 * t719 + t707 * t721;
t689 = qJD(2) * pkin(3) - pkin(8) * t696;
t694 = t695 ^ 2;
t626 = -pkin(3) * t694 + pkin(8) * t683 - qJD(2) * t689 + t636;
t723 = sin(qJ(4));
t613 = t723 * t622 + t748 * t626;
t677 = t723 * t695 + t748 * t696;
t644 = qJD(4) * t677 - t748 * t683 + t684 * t723;
t676 = -t748 * t695 + t696 * t723;
t658 = mrSges(5,1) * t676 + mrSges(5,2) * t677;
t716 = qJD(2) + qJD(4);
t671 = mrSges(5,1) * t716 - mrSges(5,3) * t677;
t715 = qJDD(2) + qJDD(4);
t657 = pkin(4) * t676 - qJ(5) * t677;
t714 = t716 ^ 2;
t608 = -pkin(4) * t714 + qJ(5) * t715 - t657 * t676 + t613;
t711 = t725 * g(1) - t728 * g(2);
t735 = -qJDD(1) * pkin(1) - t711;
t669 = -t707 * pkin(2) + qJDD(3) + t708 * t744 + (-qJ(3) * t717 - pkin(7)) * t729 + t735;
t634 = -t683 * pkin(3) - t694 * pkin(8) + t696 * t689 + t669;
t645 = -t676 * qJD(4) + t723 * t683 + t748 * t684;
t611 = (t676 * t716 - t645) * qJ(5) + (t677 * t716 + t644) * pkin(4) + t634;
t718 = sin(pkin(11));
t720 = cos(pkin(11));
t668 = t677 * t720 + t716 * t718;
t603 = -0.2e1 * qJD(5) * t668 - t718 * t608 + t720 * t611;
t639 = t645 * t720 + t715 * t718;
t667 = -t677 * t718 + t716 * t720;
t601 = (t667 * t676 - t639) * pkin(9) + (t667 * t668 + t644) * pkin(5) + t603;
t604 = 0.2e1 * qJD(5) * t667 + t720 * t608 + t718 * t611;
t638 = -t645 * t718 + t715 * t720;
t651 = pkin(5) * t676 - pkin(9) * t668;
t666 = t667 ^ 2;
t602 = -pkin(5) * t666 + pkin(9) * t638 - t651 * t676 + t604;
t722 = sin(qJ(6));
t726 = cos(qJ(6));
t599 = t601 * t726 - t602 * t722;
t646 = t667 * t726 - t668 * t722;
t616 = qJD(6) * t646 + t638 * t722 + t639 * t726;
t647 = t667 * t722 + t668 * t726;
t623 = -mrSges(7,1) * t646 + mrSges(7,2) * t647;
t672 = qJD(6) + t676;
t627 = -mrSges(7,2) * t672 + mrSges(7,3) * t646;
t643 = qJDD(6) + t644;
t597 = m(7) * t599 + mrSges(7,1) * t643 - mrSges(7,3) * t616 - t623 * t647 + t627 * t672;
t600 = t601 * t722 + t602 * t726;
t615 = -qJD(6) * t647 + t638 * t726 - t639 * t722;
t628 = mrSges(7,1) * t672 - mrSges(7,3) * t647;
t598 = m(7) * t600 - mrSges(7,2) * t643 + mrSges(7,3) * t615 + t623 * t646 - t628 * t672;
t589 = t726 * t597 + t722 * t598;
t648 = -mrSges(6,1) * t667 + mrSges(6,2) * t668;
t649 = -mrSges(6,2) * t676 + mrSges(6,3) * t667;
t587 = m(6) * t603 + mrSges(6,1) * t644 - mrSges(6,3) * t639 - t648 * t668 + t649 * t676 + t589;
t650 = mrSges(6,1) * t676 - mrSges(6,3) * t668;
t736 = -t597 * t722 + t726 * t598;
t588 = m(6) * t604 - mrSges(6,2) * t644 + mrSges(6,3) * t638 + t648 * t667 - t650 * t676 + t736;
t737 = -t587 * t718 + t720 * t588;
t582 = m(5) * t613 - mrSges(5,2) * t715 - mrSges(5,3) * t644 - t658 * t676 - t671 * t716 + t737;
t612 = t748 * t622 - t723 * t626;
t670 = -mrSges(5,2) * t716 - mrSges(5,3) * t676;
t607 = -t715 * pkin(4) - t714 * qJ(5) + t677 * t657 + qJDD(5) - t612;
t605 = -t638 * pkin(5) - t666 * pkin(9) + t668 * t651 + t607;
t733 = m(7) * t605 - t615 * mrSges(7,1) + mrSges(7,2) * t616 - t646 * t627 + t628 * t647;
t731 = -m(6) * t607 + t638 * mrSges(6,1) - mrSges(6,2) * t639 + t667 * t649 - t650 * t668 - t733;
t593 = m(5) * t612 + mrSges(5,1) * t715 - mrSges(5,3) * t645 - t658 * t677 + t670 * t716 + t731;
t575 = t723 * t582 + t748 * t593;
t680 = -mrSges(4,1) * t695 + mrSges(4,2) * t696;
t687 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t695;
t573 = m(4) * t635 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t684 + qJD(2) * t687 - t680 * t696 + t575;
t688 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t696;
t738 = t748 * t582 - t593 * t723;
t574 = m(4) * t636 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t683 - qJD(2) * t688 + t680 * t695 + t738;
t568 = t721 * t573 + t719 * t574;
t685 = -t727 * g(3) - t746;
t705 = (-mrSges(3,1) * t727 + mrSges(3,2) * t724) * qJD(1);
t743 = qJD(1) * t727;
t710 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t743;
t566 = m(3) * t685 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t706 + qJD(2) * t710 - t705 * t744 + t568;
t709 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t744;
t739 = -t573 * t719 + t721 * t574;
t567 = m(3) * t686 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t707 - qJD(2) * t709 + t705 * t743 + t739;
t740 = -t566 * t724 + t727 * t567;
t560 = m(2) * t712 - mrSges(2,1) * t729 - qJDD(1) * mrSges(2,2) + t740;
t700 = -t729 * pkin(7) + t735;
t583 = t720 * t587 + t718 * t588;
t734 = m(5) * t634 + t644 * mrSges(5,1) + t645 * mrSges(5,2) + t676 * t670 + t677 * t671 + t583;
t732 = m(4) * t669 - t683 * mrSges(4,1) + mrSges(4,2) * t684 - t695 * t687 + t688 * t696 + t734;
t730 = -m(3) * t700 + t707 * mrSges(3,1) - mrSges(3,2) * t706 - t709 * t744 + t710 * t743 - t732;
t579 = m(2) * t711 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t729 + t730;
t745 = t725 * t560 + t728 * t579;
t561 = t727 * t566 + t724 * t567;
t741 = t728 * t560 - t579 * t725;
t699 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t724 + Ifges(3,4) * t727) * qJD(1);
t698 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t724 + Ifges(3,2) * t727) * qJD(1);
t697 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t724 + Ifges(3,6) * t727) * qJD(1);
t675 = Ifges(4,1) * t696 + Ifges(4,4) * t695 + Ifges(4,5) * qJD(2);
t674 = Ifges(4,4) * t696 + Ifges(4,2) * t695 + Ifges(4,6) * qJD(2);
t673 = Ifges(4,5) * t696 + Ifges(4,6) * t695 + Ifges(4,3) * qJD(2);
t654 = Ifges(5,1) * t677 - Ifges(5,4) * t676 + Ifges(5,5) * t716;
t653 = Ifges(5,4) * t677 - Ifges(5,2) * t676 + Ifges(5,6) * t716;
t652 = Ifges(5,5) * t677 - Ifges(5,6) * t676 + Ifges(5,3) * t716;
t631 = Ifges(6,1) * t668 + Ifges(6,4) * t667 + Ifges(6,5) * t676;
t630 = Ifges(6,4) * t668 + Ifges(6,2) * t667 + Ifges(6,6) * t676;
t629 = Ifges(6,5) * t668 + Ifges(6,6) * t667 + Ifges(6,3) * t676;
t619 = Ifges(7,1) * t647 + Ifges(7,4) * t646 + Ifges(7,5) * t672;
t618 = Ifges(7,4) * t647 + Ifges(7,2) * t646 + Ifges(7,6) * t672;
t617 = Ifges(7,5) * t647 + Ifges(7,6) * t646 + Ifges(7,3) * t672;
t591 = mrSges(7,2) * t605 - mrSges(7,3) * t599 + Ifges(7,1) * t616 + Ifges(7,4) * t615 + Ifges(7,5) * t643 + t617 * t646 - t618 * t672;
t590 = -mrSges(7,1) * t605 + mrSges(7,3) * t600 + Ifges(7,4) * t616 + Ifges(7,2) * t615 + Ifges(7,6) * t643 - t617 * t647 + t619 * t672;
t577 = mrSges(6,2) * t607 - mrSges(6,3) * t603 + Ifges(6,1) * t639 + Ifges(6,4) * t638 + Ifges(6,5) * t644 - pkin(9) * t589 - t590 * t722 + t591 * t726 + t629 * t667 - t630 * t676;
t576 = -mrSges(6,1) * t607 + mrSges(6,3) * t604 + Ifges(6,4) * t639 + Ifges(6,2) * t638 + Ifges(6,6) * t644 - pkin(5) * t733 + pkin(9) * t736 + t726 * t590 + t722 * t591 - t668 * t629 + t676 * t631;
t569 = Ifges(5,4) * t645 + Ifges(5,6) * t715 - t677 * t652 + t716 * t654 - mrSges(5,1) * t634 + mrSges(5,3) * t613 - Ifges(6,5) * t639 - Ifges(6,6) * t638 - t668 * t630 + t667 * t631 - mrSges(6,1) * t603 + mrSges(6,2) * t604 - Ifges(7,5) * t616 - Ifges(7,6) * t615 - Ifges(7,3) * t643 - t647 * t618 + t646 * t619 - mrSges(7,1) * t599 + mrSges(7,2) * t600 - pkin(5) * t589 - pkin(4) * t583 + (-Ifges(5,2) - Ifges(6,3)) * t644;
t562 = mrSges(5,2) * t634 - mrSges(5,3) * t612 + Ifges(5,1) * t645 - Ifges(5,4) * t644 + Ifges(5,5) * t715 - qJ(5) * t583 - t576 * t718 + t577 * t720 - t652 * t676 - t653 * t716;
t557 = mrSges(4,2) * t669 - mrSges(4,3) * t635 + Ifges(4,1) * t684 + Ifges(4,4) * t683 + Ifges(4,5) * qJDD(2) - pkin(8) * t575 - qJD(2) * t674 + t748 * t562 - t723 * t569 + t695 * t673;
t556 = -mrSges(4,1) * t669 + mrSges(4,3) * t636 + Ifges(4,4) * t684 + Ifges(4,2) * t683 + Ifges(4,6) * qJDD(2) - pkin(3) * t734 + pkin(8) * t738 + qJD(2) * t675 + t723 * t562 + t748 * t569 - t696 * t673;
t555 = mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) + (-Ifges(4,3) - Ifges(3,3)) * qJDD(2) + (-t698 * t724 + t699 * t727) * qJD(1) + t729 * Ifges(2,5) - t718 * t577 - t720 * t576 + mrSges(2,3) * t712 - Ifges(5,3) * t715 - Ifges(3,5) * t706 - Ifges(3,6) * t707 + t695 * t675 - t696 * t674 - Ifges(4,6) * t683 - Ifges(4,5) * t684 - mrSges(3,1) * t685 + mrSges(3,2) * t686 - t676 * t654 - t677 * t653 + Ifges(5,6) * t644 - Ifges(5,5) * t645 - mrSges(4,1) * t635 + mrSges(4,2) * t636 - mrSges(5,1) * t612 + mrSges(5,2) * t613 - pkin(1) * t561 - qJ(5) * t737 - pkin(4) * t731 - pkin(3) * t575 - pkin(2) * t568;
t554 = mrSges(3,2) * t700 - mrSges(3,3) * t685 + Ifges(3,1) * t706 + Ifges(3,4) * t707 + Ifges(3,5) * qJDD(2) - qJ(3) * t568 - qJD(2) * t698 - t556 * t719 + t557 * t721 + t697 * t743;
t553 = -mrSges(3,1) * t700 + mrSges(3,3) * t686 + Ifges(3,4) * t706 + Ifges(3,2) * t707 + Ifges(3,6) * qJDD(2) - pkin(2) * t732 + qJ(3) * t739 + qJD(2) * t699 + t721 * t556 + t719 * t557 - t697 * t744;
t552 = -mrSges(2,2) * g(3) - mrSges(2,3) * t711 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t729 - pkin(7) * t561 - t553 * t724 + t554 * t727;
t1 = [-m(1) * g(1) + t741; -m(1) * g(2) + t745; (-m(1) - m(2)) * g(3) + t561; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t745 + t728 * t552 - t725 * t555; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t741 + t725 * t552 + t728 * t555; -mrSges(1,1) * g(2) + mrSges(2,1) * t711 + mrSges(1,2) * g(1) - mrSges(2,2) * t712 + Ifges(2,3) * qJDD(1) + pkin(1) * t730 + pkin(7) * t740 + t727 * t553 + t724 * t554;];
tauB  = t1;

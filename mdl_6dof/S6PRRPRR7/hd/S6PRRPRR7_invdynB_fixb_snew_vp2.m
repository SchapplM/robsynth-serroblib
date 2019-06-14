% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:01:06
% EndTime: 2019-05-05 06:01:17
% DurationCPUTime: 7.57s
% Computational Cost: add. (118530->326), mult. (235261->401), div. (0->0), fcn. (149012->12), ass. (0->142)
t749 = -2 * qJD(4);
t742 = Ifges(4,4) + Ifges(5,6);
t741 = Ifges(4,5) - Ifges(5,4);
t748 = Ifges(4,2) + Ifges(5,3);
t747 = Ifges(5,2) + Ifges(4,1);
t740 = Ifges(4,6) - Ifges(5,5);
t746 = (Ifges(4,3) + Ifges(5,1));
t697 = sin(pkin(11));
t699 = cos(pkin(11));
t678 = g(1) * t697 - g(2) * t699;
t679 = -g(1) * t699 - g(2) * t697;
t696 = -g(3) + qJDD(1);
t708 = cos(qJ(2));
t700 = cos(pkin(6));
t704 = sin(qJ(2));
t736 = t700 * t704;
t698 = sin(pkin(6));
t737 = t698 * t704;
t633 = t678 * t736 + t708 * t679 + t696 * t737;
t710 = qJD(2) ^ 2;
t629 = -pkin(2) * t710 + qJDD(2) * pkin(8) + t633;
t648 = -t678 * t698 + t696 * t700;
t703 = sin(qJ(3));
t707 = cos(qJ(3));
t624 = t707 * t629 + t703 * t648;
t671 = (-pkin(3) * t707 - qJ(4) * t703) * qJD(2);
t709 = qJD(3) ^ 2;
t730 = qJD(2) * t707;
t616 = pkin(3) * t709 - qJDD(3) * qJ(4) + (qJD(3) * t749) - t671 * t730 - t624;
t632 = -t704 * t679 + (t678 * t700 + t696 * t698) * t708;
t745 = -pkin(3) - pkin(9);
t744 = pkin(8) * t710;
t743 = pkin(9) * t710;
t715 = -qJDD(2) * pkin(2) - t632;
t628 = t715 - t744;
t729 = qJD(2) * qJD(3);
t727 = t707 * t729;
t674 = qJDD(2) * t703 + t727;
t726 = t703 * t729;
t675 = qJDD(2) * t707 - t726;
t691 = t703 * qJD(2);
t680 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t691;
t681 = -(qJD(3) * mrSges(4,2)) + mrSges(4,3) * t730;
t682 = -mrSges(5,1) * t730 - (qJD(3) * mrSges(5,3));
t713 = pkin(3) * t726 + t691 * t749 + (-t674 - t727) * qJ(4) + t715;
t618 = -pkin(3) * t675 + t713 - t744;
t683 = mrSges(5,1) * t691 + (qJD(3) * mrSges(5,2));
t626 = t703 * t629;
t721 = -qJ(4) * t709 + t671 * t691 + qJDD(4) + t626;
t610 = pkin(4) * t674 + t745 * qJDD(3) + (-pkin(4) * t729 - t703 * t743 - t648) * t707 + t721;
t685 = pkin(4) * t691 - qJD(3) * pkin(9);
t695 = t707 ^ 2;
t612 = -t685 * t691 + (-pkin(4) * t695 - pkin(8)) * t710 + t745 * t675 + t713;
t702 = sin(qJ(5));
t706 = cos(qJ(5));
t602 = t706 * t610 - t612 * t702;
t669 = -qJD(3) * t702 - t706 * t730;
t640 = qJD(5) * t669 + qJDD(3) * t706 - t675 * t702;
t666 = qJDD(5) + t674;
t670 = qJD(3) * t706 - t702 * t730;
t688 = t691 + qJD(5);
t600 = (t669 * t688 - t640) * pkin(10) + (t669 * t670 + t666) * pkin(5) + t602;
t603 = t702 * t610 + t706 * t612;
t639 = -qJD(5) * t670 - qJDD(3) * t702 - t675 * t706;
t647 = pkin(5) * t688 - pkin(10) * t670;
t665 = t669 ^ 2;
t601 = -pkin(5) * t665 + pkin(10) * t639 - t647 * t688 + t603;
t701 = sin(qJ(6));
t705 = cos(qJ(6));
t598 = t600 * t705 - t601 * t701;
t641 = t669 * t705 - t670 * t701;
t615 = qJD(6) * t641 + t639 * t701 + t640 * t705;
t642 = t669 * t701 + t670 * t705;
t625 = -mrSges(7,1) * t641 + mrSges(7,2) * t642;
t686 = qJD(6) + t688;
t630 = -mrSges(7,2) * t686 + mrSges(7,3) * t641;
t659 = qJDD(6) + t666;
t596 = m(7) * t598 + mrSges(7,1) * t659 - mrSges(7,3) * t615 - t625 * t642 + t630 * t686;
t599 = t600 * t701 + t601 * t705;
t614 = -qJD(6) * t642 + t639 * t705 - t640 * t701;
t631 = mrSges(7,1) * t686 - mrSges(7,3) * t642;
t597 = m(7) * t599 - mrSges(7,2) * t659 + mrSges(7,3) * t614 + t625 * t641 - t631 * t686;
t587 = t705 * t596 + t701 * t597;
t643 = -mrSges(6,1) * t669 + mrSges(6,2) * t670;
t645 = -mrSges(6,2) * t688 + mrSges(6,3) * t669;
t585 = m(6) * t602 + mrSges(6,1) * t666 - mrSges(6,3) * t640 - t643 * t670 + t645 * t688 + t587;
t646 = mrSges(6,1) * t688 - mrSges(6,3) * t670;
t723 = -t596 * t701 + t705 * t597;
t586 = m(6) * t603 - mrSges(6,2) * t666 + mrSges(6,3) * t639 + t643 * t669 - t646 * t688 + t723;
t734 = -t702 * t585 + t706 * t586;
t720 = -m(5) * t618 - t675 * mrSges(5,2) + t683 * t691 - t734;
t711 = -m(4) * t628 + t681 * t730 + t675 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t674 + (-t680 * t703 - t682 * t707) * qJD(2) + t720;
t578 = m(3) * t632 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t710 + t711;
t739 = t578 * t708;
t738 = t648 * t707;
t623 = -t626 + t738;
t672 = (mrSges(5,2) * t707 - mrSges(5,3) * t703) * qJD(2);
t673 = (-mrSges(4,1) * t707 + mrSges(4,2) * t703) * qJD(2);
t582 = t706 * t585 + t702 * t586;
t617 = -qJDD(3) * pkin(3) + t721 - t738;
t716 = -m(5) * t617 - t674 * mrSges(5,1) - t582;
t580 = m(4) * t623 - mrSges(4,3) * t674 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t681 - t682) * qJD(3) + (-t672 - t673) * t691 + t716;
t609 = pkin(4) * t675 + qJD(3) * t685 - t695 * t743 - t616;
t605 = -pkin(5) * t639 - pkin(10) * t665 + t647 * t670 + t609;
t717 = m(7) * t605 - t614 * mrSges(7,1) + t615 * mrSges(7,2) - t641 * t630 + t642 * t631;
t714 = -m(6) * t609 + t639 * mrSges(6,1) - t640 * mrSges(6,2) + t669 * t645 - t670 * t646 - t717;
t712 = -m(5) * t616 + qJDD(3) * mrSges(5,3) + qJD(3) * t683 + t672 * t730 - t714;
t592 = t712 - qJDD(3) * mrSges(4,2) + t673 * t730 + (mrSges(4,3) + mrSges(5,1)) * t675 + m(4) * t624 - qJD(3) * t680;
t724 = -t580 * t703 + t707 * t592;
t570 = m(3) * t633 - mrSges(3,1) * t710 - qJDD(2) * mrSges(3,2) + t724;
t573 = t707 * t580 + t703 * t592;
t572 = m(3) * t648 + t573;
t561 = t570 * t736 - t572 * t698 + t700 * t739;
t559 = m(2) * t678 + t561;
t566 = t708 * t570 - t578 * t704;
t565 = m(2) * t679 + t566;
t735 = t699 * t559 + t697 * t565;
t733 = (t746 * qJD(3)) + (t741 * t703 + t740 * t707) * qJD(2);
t732 = -t740 * qJD(3) + (-t742 * t703 - t748 * t707) * qJD(2);
t731 = t741 * qJD(3) + (t747 * t703 + t742 * t707) * qJD(2);
t560 = t570 * t737 + t700 * t572 + t698 * t739;
t725 = -t559 * t697 + t699 * t565;
t619 = Ifges(7,5) * t642 + Ifges(7,6) * t641 + Ifges(7,3) * t686;
t621 = Ifges(7,1) * t642 + Ifges(7,4) * t641 + Ifges(7,5) * t686;
t588 = -mrSges(7,1) * t605 + mrSges(7,3) * t599 + Ifges(7,4) * t615 + Ifges(7,2) * t614 + Ifges(7,6) * t659 - t619 * t642 + t621 * t686;
t620 = Ifges(7,4) * t642 + Ifges(7,2) * t641 + Ifges(7,6) * t686;
t589 = mrSges(7,2) * t605 - mrSges(7,3) * t598 + Ifges(7,1) * t615 + Ifges(7,4) * t614 + Ifges(7,5) * t659 + t619 * t641 - t620 * t686;
t634 = Ifges(6,5) * t670 + Ifges(6,6) * t669 + Ifges(6,3) * t688;
t636 = Ifges(6,1) * t670 + Ifges(6,4) * t669 + Ifges(6,5) * t688;
t574 = -mrSges(6,1) * t609 + mrSges(6,3) * t603 + Ifges(6,4) * t640 + Ifges(6,2) * t639 + Ifges(6,6) * t666 - pkin(5) * t717 + pkin(10) * t723 + t705 * t588 + t701 * t589 - t670 * t634 + t688 * t636;
t635 = Ifges(6,4) * t670 + Ifges(6,2) * t669 + Ifges(6,6) * t688;
t575 = mrSges(6,2) * t609 - mrSges(6,3) * t602 + Ifges(6,1) * t640 + Ifges(6,4) * t639 + Ifges(6,5) * t666 - pkin(10) * t587 - t588 * t701 + t589 * t705 + t634 * t669 - t635 * t688;
t581 = -mrSges(5,3) * t674 + t682 * t730 - t720;
t557 = -mrSges(4,1) * t628 - mrSges(5,1) * t616 + mrSges(5,2) * t618 + mrSges(4,3) * t624 - pkin(3) * t581 - pkin(4) * t714 - pkin(9) * t734 + t731 * qJD(3) + t740 * qJDD(3) - t706 * t574 - t702 * t575 + t742 * t674 + t748 * t675 - t733 * t691;
t562 = Ifges(6,3) * t666 - t669 * t636 + t670 * t635 + Ifges(7,3) * t659 + Ifges(6,6) * t639 + Ifges(6,5) * t640 - t641 * t621 + t642 * t620 + mrSges(4,2) * t628 + Ifges(7,6) * t614 + Ifges(7,5) * t615 + mrSges(5,1) * t617 - mrSges(5,3) * t618 - mrSges(4,3) * t623 + t747 * t674 + mrSges(6,1) * t602 - mrSges(6,2) * t603 + mrSges(7,1) * t598 - mrSges(7,2) * t599 + pkin(5) * t587 + pkin(4) * t582 - qJ(4) * t581 + t741 * qJDD(3) + t742 * t675 + t732 * qJD(3) + t733 * t730;
t555 = mrSges(3,2) * t648 - mrSges(3,3) * t632 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t710 - pkin(8) * t573 - t557 * t703 + t562 * t707;
t556 = Ifges(3,6) * qJDD(2) - pkin(2) * t573 + mrSges(3,3) * t633 - mrSges(3,1) * t648 - pkin(3) * (-qJD(3) * t682 + t716) - qJ(4) * t712 + t702 * t574 + pkin(9) * t582 - mrSges(4,1) * t623 + mrSges(4,2) * t624 - mrSges(5,2) * t617 + mrSges(5,3) * t616 - t706 * t575 + t710 * Ifges(3,5) + (-qJ(4) * mrSges(5,1) - t740) * t675 - t741 * t674 + (mrSges(5,2) * pkin(3) - t746) * qJDD(3) + (t731 * t707 + (pkin(3) * t672 + t732) * t703) * qJD(2);
t718 = pkin(7) * t566 + t555 * t704 + t556 * t708;
t554 = mrSges(3,1) * t632 - mrSges(3,2) * t633 + Ifges(3,3) * qJDD(2) + pkin(2) * t711 + pkin(8) * t724 + t707 * t557 + t703 * t562;
t553 = mrSges(2,2) * t696 - mrSges(2,3) * t678 + t555 * t708 - t556 * t704 + (-t560 * t698 - t561 * t700) * pkin(7);
t552 = -mrSges(2,1) * t696 + mrSges(2,3) * t679 - pkin(1) * t560 - t554 * t698 + t700 * t718;
t1 = [-m(1) * g(1) + t725; -m(1) * g(2) + t735; -m(1) * g(3) + m(2) * t696 + t560; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t735 - t697 * t552 + t699 * t553; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t725 + t699 * t552 + t697 * t553; -mrSges(1,1) * g(2) + mrSges(2,1) * t678 + mrSges(1,2) * g(1) - mrSges(2,2) * t679 + pkin(1) * t561 + t554 * t700 + t698 * t718;];
tauB  = t1;

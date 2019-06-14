% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:20:13
% EndTime: 2019-05-06 18:20:31
% DurationCPUTime: 14.33s
% Computational Cost: add. (211311->361), mult. (459533->442), div. (0->0), fcn. (330780->10), ass. (0->138)
t749 = Ifges(6,1) + Ifges(7,1);
t744 = Ifges(6,4) - Ifges(7,5);
t743 = Ifges(7,4) + Ifges(6,5);
t748 = Ifges(6,2) + Ifges(7,3);
t742 = Ifges(6,6) - Ifges(7,6);
t747 = -Ifges(6,3) - Ifges(7,2);
t746 = cos(qJ(5));
t745 = -mrSges(6,3) - mrSges(7,2);
t716 = sin(qJ(1));
t719 = cos(qJ(1));
t704 = -g(1) * t719 - g(2) * t716;
t721 = qJD(1) ^ 2;
t689 = -pkin(1) * t721 + qJDD(1) * pkin(7) + t704;
t715 = sin(qJ(2));
t718 = cos(qJ(2));
t671 = -g(3) * t715 + t718 * t689;
t698 = (-mrSges(3,1) * t718 + mrSges(3,2) * t715) * qJD(1);
t734 = qJD(1) * qJD(2);
t708 = t715 * t734;
t700 = qJDD(1) * t718 - t708;
t736 = qJD(1) * t715;
t701 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t736;
t703 = t716 * g(1) - t719 * g(2);
t688 = -qJDD(1) * pkin(1) - t721 * pkin(7) - t703;
t732 = t718 * t734;
t699 = qJDD(1) * t715 + t732;
t653 = (-t699 - t732) * qJ(3) + (-t700 + t708) * pkin(2) + t688;
t697 = (-pkin(2) * t718 - qJ(3) * t715) * qJD(1);
t720 = qJD(2) ^ 2;
t735 = qJD(1) * t718;
t656 = -pkin(2) * t720 + qJDD(2) * qJ(3) + t697 * t735 + t671;
t711 = sin(pkin(10));
t712 = cos(pkin(10));
t694 = qJD(2) * t711 + t712 * t736;
t630 = -0.2e1 * qJD(3) * t694 + t712 * t653 - t711 * t656;
t676 = qJDD(2) * t711 + t699 * t712;
t693 = qJD(2) * t712 - t711 * t736;
t614 = (-t693 * t735 - t676) * pkin(8) + (t693 * t694 - t700) * pkin(3) + t630;
t631 = 0.2e1 * qJD(3) * t693 + t711 * t653 + t712 * t656;
t675 = qJDD(2) * t712 - t699 * t711;
t677 = -pkin(3) * t735 - pkin(8) * t694;
t692 = t693 ^ 2;
t616 = -pkin(3) * t692 + pkin(8) * t675 + t677 * t735 + t631;
t714 = sin(qJ(4));
t717 = cos(qJ(4));
t602 = t717 * t614 - t714 * t616;
t667 = t693 * t717 - t694 * t714;
t642 = qJD(4) * t667 + t675 * t714 + t676 * t717;
t668 = t693 * t714 + t694 * t717;
t696 = qJDD(4) - t700;
t707 = qJD(4) - t735;
t599 = (t667 * t707 - t642) * pkin(9) + (t667 * t668 + t696) * pkin(4) + t602;
t603 = t714 * t614 + t717 * t616;
t641 = -qJD(4) * t668 + t675 * t717 - t676 * t714;
t659 = pkin(4) * t707 - pkin(9) * t668;
t666 = t667 ^ 2;
t601 = -pkin(4) * t666 + pkin(9) * t641 - t659 * t707 + t603;
t713 = sin(qJ(5));
t595 = t713 * t599 + t746 * t601;
t649 = t713 * t667 + t746 * t668;
t610 = qJD(5) * t649 - t746 * t641 + t642 * t713;
t706 = qJD(5) + t707;
t638 = mrSges(6,1) * t706 - mrSges(6,3) * t649;
t648 = -t746 * t667 + t668 * t713;
t690 = qJDD(5) + t696;
t627 = pkin(5) * t648 - qJ(6) * t649;
t705 = t706 ^ 2;
t592 = -pkin(5) * t705 + qJ(6) * t690 + 0.2e1 * qJD(6) * t706 - t627 * t648 + t595;
t639 = -mrSges(7,1) * t706 + mrSges(7,2) * t649;
t733 = m(7) * t592 + t690 * mrSges(7,3) + t706 * t639;
t628 = mrSges(7,1) * t648 - mrSges(7,3) * t649;
t737 = -mrSges(6,1) * t648 - mrSges(6,2) * t649 - t628;
t585 = m(6) * t595 - t690 * mrSges(6,2) + t745 * t610 - t706 * t638 + t737 * t648 + t733;
t594 = t746 * t599 - t713 * t601;
t611 = -t648 * qJD(5) + t713 * t641 + t746 * t642;
t637 = -mrSges(6,2) * t706 - mrSges(6,3) * t648;
t593 = -t690 * pkin(5) - t705 * qJ(6) + t649 * t627 + qJDD(6) - t594;
t636 = -mrSges(7,2) * t648 + mrSges(7,3) * t706;
t726 = -m(7) * t593 + t690 * mrSges(7,1) + t706 * t636;
t587 = m(6) * t594 + t690 * mrSges(6,1) + t745 * t611 + t706 * t637 + t737 * t649 + t726;
t582 = t713 * t585 + t746 * t587;
t650 = -mrSges(5,1) * t667 + mrSges(5,2) * t668;
t657 = -mrSges(5,2) * t707 + mrSges(5,3) * t667;
t578 = m(5) * t602 + mrSges(5,1) * t696 - mrSges(5,3) * t642 - t650 * t668 + t657 * t707 + t582;
t658 = mrSges(5,1) * t707 - mrSges(5,3) * t668;
t727 = t746 * t585 - t587 * t713;
t579 = m(5) * t603 - mrSges(5,2) * t696 + mrSges(5,3) * t641 + t650 * t667 - t658 * t707 + t727;
t574 = t717 * t578 + t714 * t579;
t669 = -mrSges(4,1) * t693 + mrSges(4,2) * t694;
t673 = mrSges(4,2) * t735 + mrSges(4,3) * t693;
t572 = m(4) * t630 - mrSges(4,1) * t700 - mrSges(4,3) * t676 - t669 * t694 - t673 * t735 + t574;
t674 = -mrSges(4,1) * t735 - mrSges(4,3) * t694;
t728 = -t578 * t714 + t717 * t579;
t573 = m(4) * t631 + mrSges(4,2) * t700 + mrSges(4,3) * t675 + t669 * t693 + t674 * t735 + t728;
t729 = -t572 * t711 + t712 * t573;
t567 = m(3) * t671 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t700 - qJD(2) * t701 + t698 * t735 + t729;
t670 = -t718 * g(3) - t715 * t689;
t702 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t735;
t655 = -qJDD(2) * pkin(2) - t720 * qJ(3) + t697 * t736 + qJDD(3) - t670;
t632 = -t675 * pkin(3) - t692 * pkin(8) + t694 * t677 + t655;
t605 = -t641 * pkin(4) - t666 * pkin(9) + t668 * t659 + t632;
t597 = t605 - 0.2e1 * qJD(6) * t649 + (t649 * t706 + t610) * pkin(5) + (t648 * t706 - t611) * qJ(6);
t590 = m(7) * t597 + t610 * mrSges(7,1) - t611 * mrSges(7,3) + t648 * t636 - t649 * t639;
t725 = m(6) * t605 + t610 * mrSges(6,1) + t611 * mrSges(6,2) + t648 * t637 + t649 * t638 + t590;
t723 = m(5) * t632 - t641 * mrSges(5,1) + t642 * mrSges(5,2) - t667 * t657 + t668 * t658 + t725;
t722 = -m(4) * t655 + t675 * mrSges(4,1) - t676 * mrSges(4,2) + t693 * t673 - t694 * t674 - t723;
t589 = m(3) * t670 + qJDD(2) * mrSges(3,1) - t699 * mrSges(3,3) + qJD(2) * t702 - t698 * t736 + t722;
t730 = t718 * t567 - t589 * t715;
t561 = m(2) * t704 - mrSges(2,1) * t721 - qJDD(1) * mrSges(2,2) + t730;
t568 = t572 * t712 + t573 * t711;
t724 = -m(3) * t688 + t700 * mrSges(3,1) - mrSges(3,2) * t699 - t701 * t736 + t702 * t735 - t568;
t564 = m(2) * t703 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t721 + t724;
t741 = t716 * t561 + t719 * t564;
t562 = t715 * t567 + t718 * t589;
t740 = t748 * t648 - t744 * t649 - t742 * t706;
t739 = t742 * t648 - t743 * t649 + t747 * t706;
t738 = -t744 * t648 + t749 * t649 + t743 * t706;
t731 = t719 * t561 - t564 * t716;
t687 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t715 + Ifges(3,4) * t718) * qJD(1);
t686 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t715 + Ifges(3,2) * t718) * qJD(1);
t685 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t715 + Ifges(3,6) * t718) * qJD(1);
t663 = Ifges(4,1) * t694 + Ifges(4,4) * t693 - Ifges(4,5) * t735;
t662 = Ifges(4,4) * t694 + Ifges(4,2) * t693 - Ifges(4,6) * t735;
t661 = Ifges(4,5) * t694 + Ifges(4,6) * t693 - Ifges(4,3) * t735;
t645 = Ifges(5,1) * t668 + Ifges(5,4) * t667 + Ifges(5,5) * t707;
t644 = Ifges(5,4) * t668 + Ifges(5,2) * t667 + Ifges(5,6) * t707;
t643 = Ifges(5,5) * t668 + Ifges(5,6) * t667 + Ifges(5,3) * t707;
t581 = mrSges(6,2) * t605 + mrSges(7,2) * t593 - mrSges(6,3) * t594 - mrSges(7,3) * t597 - qJ(6) * t590 - t744 * t610 + t749 * t611 + t739 * t648 + t743 * t690 + t740 * t706;
t580 = -mrSges(6,1) * t605 - mrSges(7,1) * t597 + mrSges(7,2) * t592 + mrSges(6,3) * t595 - pkin(5) * t590 - t748 * t610 + t744 * t611 + t739 * t649 + t742 * t690 + t738 * t706;
t570 = mrSges(5,2) * t632 - mrSges(5,3) * t602 + Ifges(5,1) * t642 + Ifges(5,4) * t641 + Ifges(5,5) * t696 - pkin(9) * t582 - t713 * t580 + t746 * t581 + t667 * t643 - t707 * t644;
t569 = -mrSges(5,1) * t632 + mrSges(5,3) * t603 + Ifges(5,4) * t642 + Ifges(5,2) * t641 + Ifges(5,6) * t696 - pkin(4) * t725 + pkin(9) * t727 + t746 * t580 + t713 * t581 - t668 * t643 + t707 * t645;
t558 = mrSges(4,2) * t655 - mrSges(4,3) * t630 + Ifges(4,1) * t676 + Ifges(4,4) * t675 - Ifges(4,5) * t700 - pkin(8) * t574 - t569 * t714 + t570 * t717 + t661 * t693 + t662 * t735;
t557 = Ifges(3,6) * qJDD(2) - qJ(6) * t733 + t747 * t690 - pkin(5) * t726 + (qJ(6) * t628 - t738) * t648 + (pkin(5) * t628 + t740) * t649 + (mrSges(7,2) * qJ(6) + t742) * t610 + (mrSges(7,2) * pkin(5) - t743) * t611 - pkin(2) * t568 - t685 * t736 + (Ifges(3,2) + Ifges(4,3)) * t700 - Ifges(5,3) * t696 + Ifges(3,4) * t699 + qJD(2) * t687 - mrSges(3,1) * t688 + t693 * t663 - t694 * t662 - Ifges(4,6) * t675 - Ifges(4,5) * t676 + t667 * t645 - t668 * t644 + mrSges(3,3) * t671 - Ifges(5,6) * t641 - Ifges(5,5) * t642 - mrSges(4,1) * t630 + mrSges(4,2) * t631 - mrSges(5,1) * t602 + mrSges(5,2) * t603 + mrSges(6,2) * t595 + mrSges(7,1) * t593 - mrSges(6,1) * t594 - mrSges(7,3) * t592 - pkin(4) * t582 - pkin(3) * t574;
t556 = -mrSges(4,1) * t655 + mrSges(4,3) * t631 + Ifges(4,4) * t676 + Ifges(4,2) * t675 - Ifges(4,6) * t700 - pkin(3) * t723 + pkin(8) * t728 + t717 * t569 + t714 * t570 - t694 * t661 - t663 * t735;
t555 = mrSges(3,2) * t688 - mrSges(3,3) * t670 + Ifges(3,1) * t699 + Ifges(3,4) * t700 + Ifges(3,5) * qJDD(2) - qJ(3) * t568 - qJD(2) * t686 - t556 * t711 + t558 * t712 + t685 * t735;
t554 = Ifges(2,6) * qJDD(1) + t721 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t704 - Ifges(3,5) * t699 - Ifges(3,6) * t700 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t670 + mrSges(3,2) * t671 - t711 * t558 - t712 * t556 - pkin(2) * t722 - qJ(3) * t729 - pkin(1) * t562 + (-t686 * t715 + t687 * t718) * qJD(1);
t553 = -mrSges(2,2) * g(3) - mrSges(2,3) * t703 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t721 - pkin(7) * t562 + t555 * t718 - t557 * t715;
t1 = [-m(1) * g(1) + t731; -m(1) * g(2) + t741; (-m(1) - m(2)) * g(3) + t562; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t741 + t719 * t553 - t716 * t554; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t731 + t716 * t553 + t719 * t554; -mrSges(1,1) * g(2) + mrSges(2,1) * t703 + mrSges(1,2) * g(1) - mrSges(2,2) * t704 + Ifges(2,3) * qJDD(1) + pkin(1) * t724 + pkin(7) * t730 + t715 * t555 + t718 * t557;];
tauB  = t1;

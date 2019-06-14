% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:42:40
% EndTime: 2019-05-05 10:43:02
% DurationCPUTime: 21.16s
% Computational Cost: add. (375854->342), mult. (734844->441), div. (0->0), fcn. (540641->14), ass. (0->143)
t682 = sin(pkin(12));
t684 = cos(pkin(12));
t670 = g(1) * t682 - g(2) * t684;
t671 = -g(1) * t684 - g(2) * t682;
t681 = -g(3) + qJDD(1);
t683 = sin(pkin(6));
t685 = cos(pkin(6));
t690 = sin(qJ(2));
t695 = cos(qJ(2));
t639 = -t690 * t671 + (t670 * t685 + t681 * t683) * t695;
t696 = qJD(2) ^ 2;
t700 = -qJDD(2) * pkin(2) - t639;
t634 = -pkin(8) * t696 + t700;
t689 = sin(qJ(3));
t694 = cos(qJ(3));
t710 = qJD(2) * qJD(3);
t709 = t694 * t710;
t668 = qJDD(2) * t689 + t709;
t669 = qJDD(2) * t694 - t689 * t710;
t712 = qJD(2) * t689;
t672 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t712;
t711 = qJD(2) * t694;
t673 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t711;
t714 = t685 * t690;
t715 = t683 * t690;
t640 = t670 * t714 + t695 * t671 + t681 * t715;
t635 = -pkin(2) * t696 + qJDD(2) * pkin(8) + t640;
t651 = -t670 * t683 + t681 * t685;
t619 = -t635 * t689 + t694 * t651;
t607 = (-t668 + t709) * pkin(9) + (t689 * t694 * t696 + qJDD(3)) * pkin(3) + t619;
t620 = t694 * t635 + t689 * t651;
t675 = qJD(3) * pkin(3) - pkin(9) * t712;
t680 = t694 ^ 2;
t608 = -pkin(3) * t680 * t696 + pkin(9) * t669 - qJD(3) * t675 + t620;
t688 = sin(qJ(4));
t693 = cos(qJ(4));
t598 = t688 * t607 + t693 * t608;
t659 = -t688 * t712 + t693 * t711;
t660 = (t688 * t694 + t689 * t693) * qJD(2);
t643 = -pkin(4) * t659 - pkin(10) * t660;
t679 = qJD(3) + qJD(4);
t677 = t679 ^ 2;
t678 = qJDD(3) + qJDD(4);
t590 = -pkin(4) * t677 + pkin(10) * t678 + t643 * t659 + t598;
t614 = -pkin(3) * t669 + t675 * t712 + (-pkin(9) * t680 - pkin(8)) * t696 + t700;
t629 = -t660 * qJD(4) - t688 * t668 + t669 * t693;
t630 = qJD(4) * t659 + t668 * t693 + t669 * t688;
t596 = (-t659 * t679 - t630) * pkin(10) + (t660 * t679 - t629) * pkin(4) + t614;
t687 = sin(qJ(5));
t692 = cos(qJ(5));
t585 = -t590 * t687 + t692 * t596;
t645 = -t660 * t687 + t679 * t692;
t611 = qJD(5) * t645 + t630 * t692 + t678 * t687;
t627 = qJDD(5) - t629;
t646 = t660 * t692 + t679 * t687;
t654 = qJD(5) - t659;
t583 = (t645 * t654 - t611) * pkin(11) + (t645 * t646 + t627) * pkin(5) + t585;
t586 = t692 * t590 + t687 * t596;
t610 = -qJD(5) * t646 - t630 * t687 + t678 * t692;
t633 = pkin(5) * t654 - pkin(11) * t646;
t644 = t645 ^ 2;
t584 = -pkin(5) * t644 + pkin(11) * t610 - t633 * t654 + t586;
t686 = sin(qJ(6));
t691 = cos(qJ(6));
t581 = t583 * t691 - t584 * t686;
t621 = t645 * t691 - t646 * t686;
t593 = qJD(6) * t621 + t610 * t686 + t611 * t691;
t622 = t645 * t686 + t646 * t691;
t603 = -mrSges(7,1) * t621 + mrSges(7,2) * t622;
t652 = qJD(6) + t654;
t612 = -mrSges(7,2) * t652 + mrSges(7,3) * t621;
t624 = qJDD(6) + t627;
t579 = m(7) * t581 + mrSges(7,1) * t624 - mrSges(7,3) * t593 - t603 * t622 + t612 * t652;
t582 = t583 * t686 + t584 * t691;
t592 = -qJD(6) * t622 + t610 * t691 - t611 * t686;
t613 = mrSges(7,1) * t652 - mrSges(7,3) * t622;
t580 = m(7) * t582 - mrSges(7,2) * t624 + mrSges(7,3) * t592 + t603 * t621 - t613 * t652;
t571 = t691 * t579 + t686 * t580;
t623 = -mrSges(6,1) * t645 + mrSges(6,2) * t646;
t631 = -mrSges(6,2) * t654 + mrSges(6,3) * t645;
t569 = m(6) * t585 + mrSges(6,1) * t627 - mrSges(6,3) * t611 - t623 * t646 + t631 * t654 + t571;
t632 = mrSges(6,1) * t654 - mrSges(6,3) * t646;
t704 = -t579 * t686 + t691 * t580;
t570 = m(6) * t586 - mrSges(6,2) * t627 + mrSges(6,3) * t610 + t623 * t645 - t632 * t654 + t704;
t565 = t692 * t569 + t687 * t570;
t649 = -mrSges(5,2) * t679 + mrSges(5,3) * t659;
t650 = mrSges(5,1) * t679 - mrSges(5,3) * t660;
t699 = m(5) * t614 - t629 * mrSges(5,1) + mrSges(5,2) * t630 - t659 * t649 + t650 * t660 + t565;
t697 = -m(4) * t634 + t669 * mrSges(4,1) - mrSges(4,2) * t668 - t672 * t712 + t673 * t711 - t699;
t561 = m(3) * t639 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t696 + t697;
t716 = t561 * t695;
t642 = -mrSges(5,1) * t659 + mrSges(5,2) * t660;
t705 = -t569 * t687 + t692 * t570;
t564 = m(5) * t598 - mrSges(5,2) * t678 + mrSges(5,3) * t629 + t642 * t659 - t650 * t679 + t705;
t597 = t607 * t693 - t688 * t608;
t589 = -pkin(4) * t678 - pkin(10) * t677 + t660 * t643 - t597;
t587 = -pkin(5) * t610 - pkin(11) * t644 + t633 * t646 + t589;
t701 = m(7) * t587 - t592 * mrSges(7,1) + mrSges(7,2) * t593 - t621 * t612 + t613 * t622;
t698 = -m(6) * t589 + t610 * mrSges(6,1) - mrSges(6,2) * t611 + t645 * t631 - t632 * t646 - t701;
t575 = m(5) * t597 + mrSges(5,1) * t678 - mrSges(5,3) * t630 - t642 * t660 + t649 * t679 + t698;
t556 = t688 * t564 + t693 * t575;
t667 = (-mrSges(4,1) * t694 + mrSges(4,2) * t689) * qJD(2);
t554 = m(4) * t619 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t668 + qJD(3) * t673 - t667 * t712 + t556;
t706 = t693 * t564 - t575 * t688;
t555 = m(4) * t620 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t669 - qJD(3) * t672 + t667 * t711 + t706;
t707 = -t554 * t689 + t694 * t555;
t546 = m(3) * t640 - mrSges(3,1) * t696 - qJDD(2) * mrSges(3,2) + t707;
t549 = t694 * t554 + t689 * t555;
t548 = m(3) * t651 + t549;
t537 = t546 * t714 - t548 * t683 + t685 * t716;
t535 = m(2) * t670 + t537;
t541 = t695 * t546 - t561 * t690;
t540 = m(2) * t671 + t541;
t713 = t684 * t535 + t682 * t540;
t536 = t546 * t715 + t685 * t548 + t683 * t716;
t708 = -t535 * t682 + t684 * t540;
t599 = Ifges(7,5) * t622 + Ifges(7,6) * t621 + Ifges(7,3) * t652;
t601 = Ifges(7,1) * t622 + Ifges(7,4) * t621 + Ifges(7,5) * t652;
t572 = -mrSges(7,1) * t587 + mrSges(7,3) * t582 + Ifges(7,4) * t593 + Ifges(7,2) * t592 + Ifges(7,6) * t624 - t599 * t622 + t601 * t652;
t600 = Ifges(7,4) * t622 + Ifges(7,2) * t621 + Ifges(7,6) * t652;
t573 = mrSges(7,2) * t587 - mrSges(7,3) * t581 + Ifges(7,1) * t593 + Ifges(7,4) * t592 + Ifges(7,5) * t624 + t599 * t621 - t600 * t652;
t615 = Ifges(6,5) * t646 + Ifges(6,6) * t645 + Ifges(6,3) * t654;
t617 = Ifges(6,1) * t646 + Ifges(6,4) * t645 + Ifges(6,5) * t654;
t557 = -mrSges(6,1) * t589 + mrSges(6,3) * t586 + Ifges(6,4) * t611 + Ifges(6,2) * t610 + Ifges(6,6) * t627 - pkin(5) * t701 + pkin(11) * t704 + t691 * t572 + t686 * t573 - t646 * t615 + t654 * t617;
t616 = Ifges(6,4) * t646 + Ifges(6,2) * t645 + Ifges(6,6) * t654;
t558 = mrSges(6,2) * t589 - mrSges(6,3) * t585 + Ifges(6,1) * t611 + Ifges(6,4) * t610 + Ifges(6,5) * t627 - pkin(11) * t571 - t572 * t686 + t573 * t691 + t615 * t645 - t616 * t654;
t636 = Ifges(5,5) * t660 + Ifges(5,6) * t659 + Ifges(5,3) * t679;
t637 = Ifges(5,4) * t660 + Ifges(5,2) * t659 + Ifges(5,6) * t679;
t542 = mrSges(5,2) * t614 - mrSges(5,3) * t597 + Ifges(5,1) * t630 + Ifges(5,4) * t629 + Ifges(5,5) * t678 - pkin(10) * t565 - t557 * t687 + t558 * t692 + t636 * t659 - t637 * t679;
t638 = Ifges(5,1) * t660 + Ifges(5,4) * t659 + Ifges(5,5) * t679;
t550 = Ifges(5,4) * t630 + Ifges(5,2) * t629 + Ifges(5,6) * t678 - t660 * t636 + t679 * t638 - mrSges(5,1) * t614 + mrSges(5,3) * t598 - Ifges(6,5) * t611 - Ifges(6,6) * t610 - Ifges(6,3) * t627 - t646 * t616 + t645 * t617 - mrSges(6,1) * t585 + mrSges(6,2) * t586 - Ifges(7,5) * t593 - Ifges(7,6) * t592 - Ifges(7,3) * t624 - t622 * t600 + t621 * t601 - mrSges(7,1) * t581 + mrSges(7,2) * t582 - pkin(5) * t571 - pkin(4) * t565;
t656 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t689 + Ifges(4,6) * t694) * qJD(2);
t658 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t689 + Ifges(4,4) * t694) * qJD(2);
t531 = -mrSges(4,1) * t634 + mrSges(4,3) * t620 + Ifges(4,4) * t668 + Ifges(4,2) * t669 + Ifges(4,6) * qJDD(3) - pkin(3) * t699 + pkin(9) * t706 + qJD(3) * t658 + t688 * t542 + t693 * t550 - t656 * t712;
t657 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t689 + Ifges(4,2) * t694) * qJD(2);
t533 = mrSges(4,2) * t634 - mrSges(4,3) * t619 + Ifges(4,1) * t668 + Ifges(4,4) * t669 + Ifges(4,5) * qJDD(3) - pkin(9) * t556 - qJD(3) * t657 + t542 * t693 - t550 * t688 + t656 * t711;
t530 = mrSges(3,2) * t651 - mrSges(3,3) * t639 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t696 - pkin(8) * t549 - t531 * t689 + t533 * t694;
t532 = t696 * Ifges(3,5) - t660 * t637 + t659 * t638 - mrSges(3,1) * t651 + Ifges(3,6) * qJDD(2) - pkin(2) * t549 + mrSges(3,3) * t640 - pkin(3) * t556 - Ifges(4,5) * t668 - Ifges(4,6) * t669 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t619 + mrSges(4,2) * t620 - t687 * t558 - t692 * t557 - pkin(4) * t698 - pkin(10) * t705 - Ifges(5,5) * t630 - Ifges(5,6) * t629 - Ifges(5,3) * t678 - mrSges(5,1) * t597 + mrSges(5,2) * t598 + (-t657 * t689 + t658 * t694) * qJD(2);
t702 = pkin(7) * t541 + t530 * t690 + t532 * t695;
t529 = mrSges(3,1) * t639 - mrSges(3,2) * t640 + Ifges(3,3) * qJDD(2) + pkin(2) * t697 + pkin(8) * t707 + t694 * t531 + t689 * t533;
t528 = mrSges(2,2) * t681 - mrSges(2,3) * t670 + t530 * t695 - t532 * t690 + (-t536 * t683 - t537 * t685) * pkin(7);
t527 = -mrSges(2,1) * t681 + mrSges(2,3) * t671 - pkin(1) * t536 - t529 * t683 + t702 * t685;
t1 = [-m(1) * g(1) + t708; -m(1) * g(2) + t713; -m(1) * g(3) + m(2) * t681 + t536; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t713 - t682 * t527 + t684 * t528; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t708 + t684 * t527 + t682 * t528; -mrSges(1,1) * g(2) + mrSges(2,1) * t670 + mrSges(1,2) * g(1) - mrSges(2,2) * t671 + pkin(1) * t537 + t529 * t685 + t683 * t702;];
tauB  = t1;

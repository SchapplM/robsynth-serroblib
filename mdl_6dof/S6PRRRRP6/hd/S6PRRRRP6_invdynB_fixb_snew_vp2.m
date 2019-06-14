% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:12:10
% EndTime: 2019-05-05 10:12:32
% DurationCPUTime: 21.49s
% Computational Cost: add. (362363->329), mult. (741886->425), div. (0->0), fcn. (587815->14), ass. (0->146)
t727 = Ifges(6,1) + Ifges(7,1);
t720 = Ifges(6,4) - Ifges(7,5);
t726 = -Ifges(6,5) - Ifges(7,4);
t725 = Ifges(6,2) + Ifges(7,3);
t718 = Ifges(6,6) - Ifges(7,6);
t724 = -Ifges(6,3) - Ifges(7,2);
t678 = sin(pkin(7));
t685 = sin(qJ(3));
t688 = cos(qJ(3));
t704 = qJD(2) * qJD(3);
t662 = (-qJDD(2) * t688 + t685 * t704) * t678;
t677 = sin(pkin(12));
t680 = cos(pkin(12));
t669 = g(1) * t677 - g(2) * t680;
t670 = -g(1) * t680 - g(2) * t677;
t676 = -g(3) + qJDD(1);
t686 = sin(qJ(2));
t682 = cos(pkin(6));
t689 = cos(qJ(2));
t711 = t682 * t689;
t679 = sin(pkin(6));
t714 = t679 * t689;
t638 = t669 * t711 - t670 * t686 + t676 * t714;
t690 = qJD(2) ^ 2;
t722 = pkin(9) * t678;
t634 = qJDD(2) * pkin(2) + t690 * t722 + t638;
t712 = t682 * t686;
t715 = t679 * t686;
t639 = t669 * t712 + t689 * t670 + t676 * t715;
t635 = -pkin(2) * t690 + qJDD(2) * t722 + t639;
t655 = -t669 * t679 + t676 * t682;
t681 = cos(pkin(7));
t596 = -t685 * t635 + (t634 * t681 + t655 * t678) * t688;
t723 = cos(qJ(5));
t721 = -mrSges(6,3) - mrSges(7,2);
t675 = qJD(2) * t681 + qJD(3);
t705 = qJD(2) * t678;
t701 = t688 * t705;
t658 = -mrSges(4,2) * t675 + mrSges(4,3) * t701;
t659 = (-mrSges(4,1) * t688 + mrSges(4,2) * t685) * t705;
t661 = (qJDD(2) * t685 + t688 * t704) * t678;
t674 = qJDD(2) * t681 + qJDD(3);
t713 = t681 * t685;
t716 = t678 * t685;
t597 = t634 * t713 + t688 * t635 + t655 * t716;
t660 = (-pkin(3) * t688 - pkin(10) * t685) * t705;
t673 = t675 ^ 2;
t593 = -pkin(3) * t673 + pkin(10) * t674 + t660 * t701 + t597;
t649 = t681 * t655;
t595 = t662 * pkin(3) - t661 * pkin(10) + t649 + (-t634 + (pkin(3) * t685 - pkin(10) * t688) * t675 * qJD(2)) * t678;
t684 = sin(qJ(4));
t687 = cos(qJ(4));
t589 = t687 * t593 + t684 * t595;
t702 = t685 * t705;
t653 = t675 * t687 - t684 * t702;
t654 = t675 * t684 + t687 * t702;
t637 = -pkin(4) * t653 - pkin(11) * t654;
t656 = qJDD(4) + t662;
t668 = qJD(4) - t701;
t667 = t668 ^ 2;
t585 = -pkin(4) * t667 + pkin(11) * t656 + t637 * t653 + t589;
t592 = -t674 * pkin(3) - t673 * pkin(10) + t660 * t702 - t596;
t629 = -qJD(4) * t654 - t661 * t684 + t674 * t687;
t630 = qJD(4) * t653 + t661 * t687 + t674 * t684;
t587 = (-t653 * t668 - t630) * pkin(11) + (t654 * t668 - t629) * pkin(4) + t592;
t683 = sin(qJ(5));
t582 = t723 * t585 + t683 * t587;
t641 = t723 * t654 + t683 * t668;
t600 = qJD(5) * t641 + t630 * t683 - t723 * t656;
t651 = qJD(5) - t653;
t619 = mrSges(6,1) * t651 - mrSges(6,3) * t641;
t627 = qJDD(5) - t629;
t640 = t654 * t683 - t723 * t668;
t612 = pkin(5) * t640 - qJ(6) * t641;
t650 = t651 ^ 2;
t578 = -pkin(5) * t650 + qJ(6) * t627 + 0.2e1 * qJD(6) * t651 - t612 * t640 + t582;
t620 = -mrSges(7,1) * t651 + mrSges(7,2) * t641;
t703 = m(7) * t578 + t627 * mrSges(7,3) + t651 * t620;
t613 = mrSges(7,1) * t640 - mrSges(7,3) * t641;
t706 = -mrSges(6,1) * t640 - mrSges(6,2) * t641 - t613;
t574 = m(6) * t582 - t627 * mrSges(6,2) + t721 * t600 - t651 * t619 + t706 * t640 + t703;
t581 = -t683 * t585 + t723 * t587;
t601 = -t640 * qJD(5) + t723 * t630 + t683 * t656;
t618 = -mrSges(6,2) * t651 - mrSges(6,3) * t640;
t579 = -t627 * pkin(5) - t650 * qJ(6) + t641 * t612 + qJDD(6) - t581;
t617 = -mrSges(7,2) * t640 + mrSges(7,3) * t651;
t697 = -m(7) * t579 + t627 * mrSges(7,1) + t651 * t617;
t575 = m(6) * t581 + t627 * mrSges(6,1) + t721 * t601 + t651 * t618 + t706 * t641 + t697;
t570 = t683 * t574 + t723 * t575;
t642 = -mrSges(5,2) * t668 + mrSges(5,3) * t653;
t643 = mrSges(5,1) * t668 - mrSges(5,3) * t654;
t692 = -m(5) * t592 + t629 * mrSges(5,1) - t630 * mrSges(5,2) + t653 * t642 - t654 * t643 - t570;
t564 = m(4) * t596 + t674 * mrSges(4,1) - t661 * mrSges(4,3) + t675 * t658 - t659 * t702 + t692;
t717 = t564 * t688;
t657 = mrSges(4,1) * t675 - mrSges(4,3) * t702;
t636 = -mrSges(5,1) * t653 + mrSges(5,2) * t654;
t698 = t723 * t574 - t575 * t683;
t567 = m(5) * t589 - mrSges(5,2) * t656 + mrSges(5,3) * t629 + t636 * t653 - t643 * t668 + t698;
t588 = -t684 * t593 + t687 * t595;
t584 = -t656 * pkin(4) - t667 * pkin(11) + t654 * t637 - t588;
t580 = -0.2e1 * qJD(6) * t641 + (t640 * t651 - t601) * qJ(6) + (t641 * t651 + t600) * pkin(5) + t584;
t576 = m(7) * t580 + mrSges(7,1) * t600 - t601 * mrSges(7,3) + t617 * t640 - t641 * t620;
t691 = -m(6) * t584 - t600 * mrSges(6,1) - mrSges(6,2) * t601 - t640 * t618 - t619 * t641 - t576;
t572 = m(5) * t588 + mrSges(5,1) * t656 - mrSges(5,3) * t630 - t636 * t654 + t642 * t668 + t691;
t699 = t687 * t567 - t572 * t684;
t558 = m(4) * t597 - mrSges(4,2) * t674 - mrSges(4,3) * t662 - t657 * t675 + t659 * t701 + t699;
t561 = t684 * t567 + t687 * t572;
t615 = -t678 * t634 + t649;
t560 = m(4) * t615 + t662 * mrSges(4,1) + t661 * mrSges(4,2) + (t657 * t685 - t658 * t688) * t705 + t561;
t547 = t558 * t713 - t560 * t678 + t681 * t717;
t543 = m(3) * t638 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t690 + t547;
t546 = t558 * t716 + t681 * t560 + t678 * t717;
t545 = m(3) * t655 + t546;
t553 = t688 * t558 - t564 * t685;
t552 = m(3) * t639 - mrSges(3,1) * t690 - qJDD(2) * mrSges(3,2) + t553;
t533 = t543 * t711 - t545 * t679 + t552 * t712;
t531 = m(2) * t669 + t533;
t539 = -t543 * t686 + t689 * t552;
t538 = m(2) * t670 + t539;
t710 = t680 * t531 + t677 * t538;
t709 = t725 * t640 - t720 * t641 - t718 * t651;
t708 = t718 * t640 + t726 * t641 + t724 * t651;
t707 = -t720 * t640 + t727 * t641 - t726 * t651;
t532 = t543 * t714 + t682 * t545 + t552 * t715;
t700 = -t531 * t677 + t680 * t538;
t568 = -mrSges(6,1) * t584 - mrSges(7,1) * t580 + mrSges(7,2) * t578 + mrSges(6,3) * t582 - pkin(5) * t576 - t725 * t600 + t720 * t601 + t718 * t627 + t708 * t641 + t707 * t651;
t569 = mrSges(6,2) * t584 + mrSges(7,2) * t579 - mrSges(6,3) * t581 - mrSges(7,3) * t580 - qJ(6) * t576 - t720 * t600 + t727 * t601 - t627 * t726 + t708 * t640 + t709 * t651;
t623 = Ifges(5,5) * t654 + Ifges(5,6) * t653 + Ifges(5,3) * t668;
t624 = Ifges(5,4) * t654 + Ifges(5,2) * t653 + Ifges(5,6) * t668;
t548 = mrSges(5,2) * t592 - mrSges(5,3) * t588 + Ifges(5,1) * t630 + Ifges(5,4) * t629 + Ifges(5,5) * t656 - pkin(11) * t570 - t683 * t568 + t723 * t569 + t653 * t623 - t668 * t624;
t625 = Ifges(5,1) * t654 + Ifges(5,4) * t653 + Ifges(5,5) * t668;
t554 = Ifges(5,4) * t630 + Ifges(5,2) * t629 + Ifges(5,6) * t656 - t654 * t623 + t668 * t625 - mrSges(5,1) * t592 + mrSges(5,3) * t589 - mrSges(6,1) * t581 + mrSges(6,2) * t582 + mrSges(7,1) * t579 - mrSges(7,3) * t578 - pkin(5) * t697 - qJ(6) * t703 - pkin(4) * t570 + (pkin(5) * t613 + t709) * t641 + (qJ(6) * t613 - t707) * t640 + t724 * t627 + (mrSges(7,2) * pkin(5) + t726) * t601 + (mrSges(7,2) * qJ(6) + t718) * t600;
t646 = Ifges(4,6) * t675 + (Ifges(4,4) * t685 + Ifges(4,2) * t688) * t705;
t647 = Ifges(4,5) * t675 + (Ifges(4,1) * t685 + Ifges(4,4) * t688) * t705;
t534 = Ifges(4,5) * t661 - Ifges(4,6) * t662 + Ifges(4,3) * t674 + mrSges(4,1) * t596 - mrSges(4,2) * t597 + t684 * t548 + t687 * t554 + pkin(3) * t692 + pkin(10) * t699 + (t646 * t685 - t647 * t688) * t705;
t645 = Ifges(4,3) * t675 + (Ifges(4,5) * t685 + Ifges(4,6) * t688) * t705;
t535 = mrSges(4,2) * t615 - mrSges(4,3) * t596 + Ifges(4,1) * t661 - Ifges(4,4) * t662 + Ifges(4,5) * t674 - pkin(10) * t561 + t548 * t687 - t554 * t684 + t645 * t701 - t646 * t675;
t540 = Ifges(4,4) * t661 - Ifges(4,2) * t662 + Ifges(4,6) * t674 - t645 * t702 + t675 * t647 - mrSges(4,1) * t615 + mrSges(4,3) * t597 - Ifges(5,5) * t630 - Ifges(5,6) * t629 - Ifges(5,3) * t656 - t654 * t624 + t653 * t625 - mrSges(5,1) * t588 + mrSges(5,2) * t589 - t683 * t569 - t723 * t568 - pkin(4) * t691 - pkin(11) * t698 - pkin(3) * t561;
t693 = pkin(9) * t553 + t535 * t685 + t540 * t688;
t528 = -mrSges(3,1) * t655 + mrSges(3,3) * t639 + t690 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t546 - t678 * t534 + t693 * t681;
t529 = mrSges(3,2) * t655 - mrSges(3,3) * t638 + Ifges(3,5) * qJDD(2) - t690 * Ifges(3,6) + t688 * t535 - t685 * t540 + (-t546 * t678 - t547 * t681) * pkin(9);
t694 = pkin(8) * t539 + t528 * t689 + t529 * t686;
t527 = mrSges(3,1) * t638 - mrSges(3,2) * t639 + Ifges(3,3) * qJDD(2) + pkin(2) * t547 + t681 * t534 + t693 * t678;
t526 = mrSges(2,2) * t676 - mrSges(2,3) * t669 - t686 * t528 + t689 * t529 + (-t532 * t679 - t533 * t682) * pkin(8);
t525 = -mrSges(2,1) * t676 + mrSges(2,3) * t670 - pkin(1) * t532 - t679 * t527 + t694 * t682;
t1 = [-m(1) * g(1) + t700; -m(1) * g(2) + t710; -m(1) * g(3) + m(2) * t676 + t532; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t710 - t677 * t525 + t680 * t526; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t700 + t680 * t525 + t677 * t526; -mrSges(1,1) * g(2) + mrSges(2,1) * t669 + mrSges(1,2) * g(1) - mrSges(2,2) * t670 + pkin(1) * t533 + t682 * t527 + t679 * t694;];
tauB  = t1;

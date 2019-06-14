% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRR1
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
% Datum: 2019-05-05 10:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:24:37
% EndTime: 2019-05-05 10:25:01
% DurationCPUTime: 23.44s
% Computational Cost: add. (394428->342), mult. (820878->442), div. (0->0), fcn. (617130->14), ass. (0->143)
t678 = sin(pkin(12));
t680 = cos(pkin(12));
t665 = g(1) * t678 - g(2) * t680;
t666 = -g(1) * t680 - g(2) * t678;
t677 = -g(3) + qJDD(1);
t679 = sin(pkin(6));
t681 = cos(pkin(6));
t686 = sin(qJ(2));
t691 = cos(qJ(2));
t637 = -t686 * t666 + (t665 * t681 + t677 * t679) * t691;
t692 = qJD(2) ^ 2;
t695 = -qJDD(2) * pkin(2) - t637;
t631 = -t692 * pkin(8) + t695;
t685 = sin(qJ(3));
t690 = cos(qJ(3));
t706 = qJD(2) * qJD(3);
t705 = t690 * t706;
t663 = qJDD(2) * t685 + t705;
t664 = qJDD(2) * t690 - t685 * t706;
t708 = qJD(2) * t685;
t667 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t708;
t707 = qJD(2) * t690;
t668 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t707;
t670 = qJD(3) * pkin(3) - pkin(9) * t708;
t676 = t690 ^ 2;
t620 = -t664 * pkin(3) + t670 * t708 + (-pkin(9) * t676 - pkin(8)) * t692 + t695;
t684 = sin(qJ(4));
t689 = cos(qJ(4));
t656 = (t684 * t690 + t685 * t689) * qJD(2);
t627 = -qJD(4) * t656 - t663 * t684 + t664 * t689;
t655 = (-t684 * t685 + t689 * t690) * qJD(2);
t628 = qJD(4) * t655 + t663 * t689 + t664 * t684;
t675 = qJD(3) + qJD(4);
t646 = -mrSges(5,2) * t675 + mrSges(5,3) * t655;
t647 = mrSges(5,1) * t675 - mrSges(5,3) * t656;
t710 = t681 * t686;
t711 = t679 * t686;
t638 = t665 * t710 + t691 * t666 + t677 * t711;
t632 = -pkin(2) * t692 + qJDD(2) * pkin(8) + t638;
t649 = -t665 * t679 + t677 * t681;
t621 = -t685 * t632 + t690 * t649;
t608 = (-t663 + t705) * pkin(9) + (t685 * t690 * t692 + qJDD(3)) * pkin(3) + t621;
t622 = t690 * t632 + t685 * t649;
t610 = -pkin(3) * t676 * t692 + pkin(9) * t664 - qJD(3) * t670 + t622;
t589 = t689 * t608 - t684 * t610;
t674 = qJDD(3) + qJDD(4);
t586 = (t655 * t675 - t628) * pkin(10) + (t655 * t656 + t674) * pkin(4) + t589;
t590 = t684 * t608 + t689 * t610;
t648 = pkin(4) * t675 - pkin(10) * t656;
t651 = t655 ^ 2;
t588 = -pkin(4) * t651 + pkin(10) * t627 - t648 * t675 + t590;
t683 = sin(qJ(5));
t688 = cos(qJ(5));
t583 = t683 * t586 + t688 * t588;
t641 = t655 * t688 - t656 * t683;
t642 = t655 * t683 + t656 * t688;
t619 = -pkin(5) * t641 - pkin(11) * t642;
t673 = qJD(5) + t675;
t671 = t673 ^ 2;
t672 = qJDD(5) + t674;
t581 = -pkin(5) * t671 + pkin(11) * t672 + t619 * t641 + t583;
t595 = -t627 * pkin(4) - t651 * pkin(10) + t656 * t648 + t620;
t603 = -qJD(5) * t642 + t627 * t688 - t628 * t683;
t604 = qJD(5) * t641 + t627 * t683 + t628 * t688;
t584 = (-t641 * t673 - t604) * pkin(11) + (t642 * t673 - t603) * pkin(5) + t595;
t682 = sin(qJ(6));
t687 = cos(qJ(6));
t578 = -t581 * t682 + t584 * t687;
t623 = -t642 * t682 + t673 * t687;
t593 = qJD(6) * t623 + t604 * t687 + t672 * t682;
t601 = qJDD(6) - t603;
t624 = t642 * t687 + t673 * t682;
t609 = -mrSges(7,1) * t623 + mrSges(7,2) * t624;
t636 = qJD(6) - t641;
t611 = -mrSges(7,2) * t636 + mrSges(7,3) * t623;
t576 = m(7) * t578 + mrSges(7,1) * t601 - mrSges(7,3) * t593 - t609 * t624 + t611 * t636;
t579 = t581 * t687 + t584 * t682;
t592 = -qJD(6) * t624 - t604 * t682 + t672 * t687;
t612 = mrSges(7,1) * t636 - mrSges(7,3) * t624;
t577 = m(7) * t579 - mrSges(7,2) * t601 + mrSges(7,3) * t592 + t609 * t623 - t612 * t636;
t568 = t687 * t576 + t682 * t577;
t629 = -mrSges(6,2) * t673 + mrSges(6,3) * t641;
t630 = mrSges(6,1) * t673 - mrSges(6,3) * t642;
t698 = m(6) * t595 - t603 * mrSges(6,1) + t604 * mrSges(6,2) - t641 * t629 + t642 * t630 + t568;
t694 = m(5) * t620 - t627 * mrSges(5,1) + mrSges(5,2) * t628 - t655 * t646 + t647 * t656 + t698;
t693 = -m(4) * t631 + t664 * mrSges(4,1) - mrSges(4,2) * t663 - t667 * t708 + t668 * t707 - t694;
t564 = m(3) * t637 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t692 + t693;
t712 = t564 * t691;
t618 = -mrSges(6,1) * t641 + mrSges(6,2) * t642;
t700 = -t576 * t682 + t687 * t577;
t567 = m(6) * t583 - mrSges(6,2) * t672 + mrSges(6,3) * t603 + t618 * t641 - t630 * t673 + t700;
t582 = t586 * t688 - t588 * t683;
t580 = -pkin(5) * t672 - pkin(11) * t671 + t619 * t642 - t582;
t696 = -m(7) * t580 + t592 * mrSges(7,1) - mrSges(7,2) * t593 + t623 * t611 - t612 * t624;
t572 = m(6) * t582 + mrSges(6,1) * t672 - mrSges(6,3) * t604 - t618 * t642 + t629 * t673 + t696;
t561 = t683 * t567 + t688 * t572;
t643 = -mrSges(5,1) * t655 + mrSges(5,2) * t656;
t559 = m(5) * t589 + mrSges(5,1) * t674 - mrSges(5,3) * t628 - t643 * t656 + t646 * t675 + t561;
t701 = t688 * t567 - t572 * t683;
t560 = m(5) * t590 - mrSges(5,2) * t674 + mrSges(5,3) * t627 + t643 * t655 - t647 * t675 + t701;
t553 = t689 * t559 + t684 * t560;
t662 = (-mrSges(4,1) * t690 + mrSges(4,2) * t685) * qJD(2);
t551 = m(4) * t621 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t663 + qJD(3) * t668 - t662 * t708 + t553;
t702 = -t559 * t684 + t689 * t560;
t552 = m(4) * t622 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t664 - qJD(3) * t667 + t662 * t707 + t702;
t703 = -t551 * t685 + t690 * t552;
t543 = m(3) * t638 - mrSges(3,1) * t692 - qJDD(2) * mrSges(3,2) + t703;
t546 = t690 * t551 + t685 * t552;
t545 = m(3) * t649 + t546;
t534 = t543 * t710 - t545 * t679 + t681 * t712;
t532 = m(2) * t665 + t534;
t538 = t691 * t543 - t564 * t686;
t537 = m(2) * t666 + t538;
t709 = t680 * t532 + t678 * t537;
t533 = t543 * t711 + t681 * t545 + t679 * t712;
t704 = -t532 * t678 + t680 * t537;
t596 = Ifges(7,5) * t624 + Ifges(7,6) * t623 + Ifges(7,3) * t636;
t598 = Ifges(7,1) * t624 + Ifges(7,4) * t623 + Ifges(7,5) * t636;
t569 = -mrSges(7,1) * t580 + mrSges(7,3) * t579 + Ifges(7,4) * t593 + Ifges(7,2) * t592 + Ifges(7,6) * t601 - t596 * t624 + t598 * t636;
t597 = Ifges(7,4) * t624 + Ifges(7,2) * t623 + Ifges(7,6) * t636;
t570 = mrSges(7,2) * t580 - mrSges(7,3) * t578 + Ifges(7,1) * t593 + Ifges(7,4) * t592 + Ifges(7,5) * t601 + t596 * t623 - t597 * t636;
t613 = Ifges(6,5) * t642 + Ifges(6,6) * t641 + Ifges(6,3) * t673;
t614 = Ifges(6,4) * t642 + Ifges(6,2) * t641 + Ifges(6,6) * t673;
t554 = mrSges(6,2) * t595 - mrSges(6,3) * t582 + Ifges(6,1) * t604 + Ifges(6,4) * t603 + Ifges(6,5) * t672 - pkin(11) * t568 - t569 * t682 + t570 * t687 + t613 * t641 - t614 * t673;
t615 = Ifges(6,1) * t642 + Ifges(6,4) * t641 + Ifges(6,5) * t673;
t555 = -mrSges(6,1) * t595 - mrSges(7,1) * t578 + mrSges(7,2) * t579 + mrSges(6,3) * t583 + Ifges(6,4) * t604 - Ifges(7,5) * t593 + Ifges(6,2) * t603 + Ifges(6,6) * t672 - Ifges(7,6) * t592 - Ifges(7,3) * t601 - pkin(5) * t568 - t597 * t624 + t598 * t623 - t613 * t642 + t615 * t673;
t633 = Ifges(5,5) * t656 + Ifges(5,6) * t655 + Ifges(5,3) * t675;
t635 = Ifges(5,1) * t656 + Ifges(5,4) * t655 + Ifges(5,5) * t675;
t539 = -mrSges(5,1) * t620 + mrSges(5,3) * t590 + Ifges(5,4) * t628 + Ifges(5,2) * t627 + Ifges(5,6) * t674 - pkin(4) * t698 + pkin(10) * t701 + t683 * t554 + t688 * t555 - t656 * t633 + t675 * t635;
t634 = Ifges(5,4) * t656 + Ifges(5,2) * t655 + Ifges(5,6) * t675;
t547 = mrSges(5,2) * t620 - mrSges(5,3) * t589 + Ifges(5,1) * t628 + Ifges(5,4) * t627 + Ifges(5,5) * t674 - pkin(10) * t561 + t554 * t688 - t555 * t683 + t633 * t655 - t634 * t675;
t652 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t685 + Ifges(4,6) * t690) * qJD(2);
t654 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t685 + Ifges(4,4) * t690) * qJD(2);
t528 = -mrSges(4,1) * t631 + mrSges(4,3) * t622 + Ifges(4,4) * t663 + Ifges(4,2) * t664 + Ifges(4,6) * qJDD(3) - pkin(3) * t694 + pkin(9) * t702 + qJD(3) * t654 + t689 * t539 + t684 * t547 - t652 * t708;
t653 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t685 + Ifges(4,2) * t690) * qJD(2);
t529 = mrSges(4,2) * t631 - mrSges(4,3) * t621 + Ifges(4,1) * t663 + Ifges(4,4) * t664 + Ifges(4,5) * qJDD(3) - pkin(9) * t553 - qJD(3) * t653 - t539 * t684 + t547 * t689 + t652 * t707;
t527 = mrSges(3,2) * t649 - mrSges(3,3) * t637 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t692 - pkin(8) * t546 - t528 * t685 + t529 * t690;
t530 = (-t653 * t685 + t654 * t690) * qJD(2) + t692 * Ifges(3,5) - t687 * t569 - t682 * t570 - Ifges(6,3) * t672 - Ifges(5,3) * t674 - t656 * t634 - Ifges(4,5) * t663 - Ifges(4,6) * t664 - mrSges(3,1) * t649 + t655 * t635 + t641 * t615 - t642 * t614 - pkin(11) * t700 + mrSges(3,3) * t638 - Ifges(5,5) * t628 - Ifges(5,6) * t627 - mrSges(4,1) * t621 + mrSges(4,2) * t622 - Ifges(6,6) * t603 - Ifges(6,5) * t604 - mrSges(5,1) * t589 + mrSges(5,2) * t590 - mrSges(6,1) * t582 + mrSges(6,2) * t583 - pkin(4) * t561 - pkin(3) * t553 - pkin(5) * t696 - pkin(2) * t546 + Ifges(3,6) * qJDD(2) - Ifges(4,3) * qJDD(3);
t697 = pkin(7) * t538 + t527 * t686 + t530 * t691;
t526 = mrSges(3,1) * t637 - mrSges(3,2) * t638 + Ifges(3,3) * qJDD(2) + pkin(2) * t693 + pkin(8) * t703 + t690 * t528 + t685 * t529;
t525 = mrSges(2,2) * t677 - mrSges(2,3) * t665 + t691 * t527 - t686 * t530 + (-t533 * t679 - t534 * t681) * pkin(7);
t524 = -mrSges(2,1) * t677 + mrSges(2,3) * t666 - pkin(1) * t533 - t679 * t526 + t681 * t697;
t1 = [-m(1) * g(1) + t704; -m(1) * g(2) + t709; -m(1) * g(3) + m(2) * t677 + t533; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t709 - t678 * t524 + t680 * t525; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t704 + t680 * t524 + t678 * t525; -mrSges(1,1) * g(2) + mrSges(2,1) * t665 + mrSges(1,2) * g(1) - mrSges(2,2) * t666 + pkin(1) * t534 + t681 * t526 + t679 * t697;];
tauB  = t1;

% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR15_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:25:23
% EndTime: 2024-09-27 22:25:47
% DurationCPUTime: 11.50s
% Computational Cost: add. (574661->337), mult. (1616007->453), div. (0->0), fcn. (1309214->14), ass. (0->149)
t682 = sin(pkin(5));
t688 = sin(qJ(2));
t693 = cos(qJ(2));
t710 = qJD(1) * qJD(2);
t666 = (qJDD(1) * t688 + t693 * t710) * t682;
t684 = cos(pkin(5));
t674 = t684 * qJDD(1) + qJDD(2);
t675 = t684 * qJD(1) + qJD(2);
t689 = sin(qJ(1));
t694 = cos(qJ(1));
t671 = t689 * g(1) - t694 * g(2);
t695 = qJD(1) ^ 2;
t722 = pkin(8) * t682;
t662 = qJDD(1) * pkin(1) + t695 * t722 + t671;
t672 = -t694 * g(1) - t689 * g(2);
t663 = -t695 * pkin(1) + qJDD(1) * t722 + t672;
t713 = t684 * t693;
t703 = t662 * t713 - t688 * t663;
t717 = t682 ^ 2 * t695;
t622 = t674 * pkin(2) - t666 * pkin(9) + (pkin(2) * t688 * t717 + (pkin(9) * qJD(1) * t675 - g(3)) * t682) * t693 + t703;
t714 = t684 * t688;
t716 = t682 * t688;
t642 = -g(3) * t716 + t662 * t714 + t693 * t663;
t711 = qJD(1) * t682;
t708 = t688 * t711;
t665 = t675 * pkin(2) - pkin(9) * t708;
t667 = (qJDD(1) * t693 - t688 * t710) * t682;
t709 = t693 ^ 2 * t717;
t623 = -pkin(2) * t709 + t667 * pkin(9) - t675 * t665 + t642;
t687 = sin(qJ(3));
t692 = cos(qJ(3));
t605 = t692 * t622 - t687 * t623;
t658 = (-t687 * t688 + t692 * t693) * t711;
t633 = t658 * qJD(3) + t692 * t666 + t687 * t667;
t659 = (t687 * t693 + t688 * t692) * t711;
t670 = qJDD(3) + t674;
t673 = qJD(3) + t675;
t596 = (t658 * t673 - t633) * pkin(10) + (t658 * t659 + t670) * pkin(3) + t605;
t606 = t687 * t622 + t692 * t623;
t632 = -t659 * qJD(3) - t687 * t666 + t692 * t667;
t648 = t673 * pkin(3) - t659 * pkin(10);
t656 = t658 ^ 2;
t598 = -t656 * pkin(3) + t632 * pkin(10) - t673 * t648 + t606;
t686 = sin(qJ(4));
t691 = cos(qJ(4));
t591 = t686 * t596 + t691 * t598;
t643 = t691 * t658 - t686 * t659;
t644 = t686 * t658 + t691 * t659;
t681 = sin(pkin(6));
t721 = pkin(11) * t681;
t624 = -t643 * pkin(4) - t644 * t721;
t669 = qJD(4) + t673;
t683 = cos(pkin(6));
t720 = pkin(11) * t683;
t631 = t669 * pkin(4) - t644 * t720;
t611 = -t644 * qJD(4) + t691 * t632 - t686 * t633;
t668 = qJDD(4) + t670;
t701 = t611 * t683 + t668 * t681;
t588 = t701 * pkin(11) + t643 * t624 - t669 * t631 + t591;
t685 = sin(qJ(5));
t690 = cos(qJ(5));
t590 = t691 * t596 - t686 * t598;
t612 = t643 * qJD(4) + t686 * t632 + t691 * t633;
t700 = t643 * t683 + t669 * t681;
t628 = t700 * pkin(11);
t587 = t668 * pkin(4) - t612 * t720 - t644 * t624 + t669 * t628 + t590;
t652 = -t684 * g(3) - t682 * t662;
t627 = -t667 * pkin(2) - pkin(9) * t709 + t665 * t708 + t652;
t604 = -t632 * pkin(3) - t656 * pkin(10) + t659 * t648 + t627;
t589 = -t611 * pkin(4) - t612 * t721 - t643 * t628 + t644 * t631 + t604;
t702 = t587 * t683 + t589 * t681;
t584 = -t685 * t588 + t702 * t690;
t613 = -t685 * t644 + t700 * t690;
t593 = t613 * qJD(5) + t690 * t612 + t701 * t685;
t614 = t690 * t644 + t700 * t685;
t602 = -t613 * mrSges(6,1) + t614 * mrSges(6,2);
t607 = -t681 * t611 + t683 * t668 + qJDD(5);
t629 = -t681 * t643 + t683 * t669 + qJD(5);
t608 = -t629 * mrSges(6,2) + t613 * mrSges(6,3);
t580 = m(6) * t584 + t607 * mrSges(6,1) - t593 * mrSges(6,3) - t614 * t602 + t629 * t608;
t585 = t690 * t588 + t702 * t685;
t592 = -t614 * qJD(5) - t685 * t612 + t701 * t690;
t609 = t629 * mrSges(6,1) - t614 * mrSges(6,3);
t581 = m(6) * t585 - t607 * mrSges(6,2) + t592 * mrSges(6,3) + t613 * t602 - t629 * t609;
t723 = t580 * t690 + t581 * t685;
t715 = t682 * t693;
t586 = -t681 * t587 + t683 * t589;
t583 = m(6) * t586 - t592 * mrSges(6,1) + t593 * mrSges(6,2) - t613 * t608 + t614 * t609;
t567 = -t681 * t583 + t723 * t683;
t625 = -t643 * mrSges(5,1) + t644 * mrSges(5,2);
t634 = -t669 * mrSges(5,2) + t643 * mrSges(5,3);
t565 = m(5) * t590 + t668 * mrSges(5,1) - t612 * mrSges(5,3) - t644 * t625 + t669 * t634 + t567;
t571 = -t685 * t580 + t690 * t581;
t635 = t669 * mrSges(5,1) - t644 * mrSges(5,3);
t570 = m(5) * t591 - t668 * mrSges(5,2) + t611 * mrSges(5,3) + t643 * t625 - t669 * t635 + t571;
t561 = t691 * t565 + t686 * t570;
t645 = -t658 * mrSges(4,1) + t659 * mrSges(4,2);
t646 = -t673 * mrSges(4,2) + t658 * mrSges(4,3);
t559 = m(4) * t605 + t670 * mrSges(4,1) - t633 * mrSges(4,3) - t659 * t645 + t673 * t646 + t561;
t647 = t673 * mrSges(4,1) - t659 * mrSges(4,3);
t704 = -t686 * t565 + t691 * t570;
t560 = m(4) * t606 - t670 * mrSges(4,2) + t632 * mrSges(4,3) + t658 * t645 - t673 * t647 + t704;
t553 = t692 * t559 + t687 * t560;
t641 = -g(3) * t715 + t703;
t707 = t693 * t711;
t661 = -t675 * mrSges(3,2) + mrSges(3,3) * t707;
t664 = (-mrSges(3,1) * t693 + mrSges(3,2) * t688) * t711;
t551 = m(3) * t641 + t674 * mrSges(3,1) - t666 * mrSges(3,3) + t675 * t661 - t664 * t708 + t553;
t660 = t675 * mrSges(3,1) - mrSges(3,3) * t708;
t705 = -t687 * t559 + t692 * t560;
t552 = m(3) * t642 - t674 * mrSges(3,2) + t667 * mrSges(3,3) - t675 * t660 + t664 * t707 + t705;
t566 = t683 * t583 + t723 * t681;
t697 = m(5) * t604 - t611 * mrSges(5,1) + t612 * mrSges(5,2) - t643 * t634 + t644 * t635 + t566;
t696 = m(4) * t627 - t632 * mrSges(4,1) + t633 * mrSges(4,2) - t658 * t646 + t659 * t647 + t697;
t563 = t666 * mrSges(3,2) - t667 * mrSges(3,1) + m(3) * t652 + (t660 * t688 - t661 * t693) * t711 + t696;
t541 = t551 * t713 + t552 * t714 - t682 * t563;
t539 = m(2) * t671 + qJDD(1) * mrSges(2,1) - t695 * mrSges(2,2) + t541;
t545 = -t688 * t551 + t693 * t552;
t544 = m(2) * t672 - t695 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t545;
t712 = t694 * t539 + t689 * t544;
t540 = t551 * t715 + t552 * t716 + t684 * t563;
t706 = -t689 * t539 + t694 * t544;
t600 = Ifges(6,4) * t614 + Ifges(6,2) * t613 + Ifges(6,6) * t629;
t601 = Ifges(6,1) * t614 + Ifges(6,4) * t613 + Ifges(6,5) * t629;
t572 = mrSges(6,1) * t584 - mrSges(6,2) * t585 + Ifges(6,5) * t593 + Ifges(6,6) * t592 + Ifges(6,3) * t607 + t614 * t600 - t613 * t601;
t615 = Ifges(5,5) * t644 + Ifges(5,6) * t643 + Ifges(5,3) * t669;
t617 = Ifges(5,1) * t644 + Ifges(5,4) * t643 + Ifges(5,5) * t669;
t599 = Ifges(6,5) * t614 + Ifges(6,6) * t613 + Ifges(6,3) * t629;
t573 = -mrSges(6,1) * t586 + mrSges(6,3) * t585 + Ifges(6,4) * t593 + Ifges(6,2) * t592 + Ifges(6,6) * t607 - t614 * t599 + t629 * t601;
t574 = mrSges(6,2) * t586 - mrSges(6,3) * t584 + Ifges(6,1) * t593 + Ifges(6,4) * t592 + Ifges(6,5) * t607 + t613 * t599 - t629 * t600;
t698 = pkin(11) * t571 + t690 * t573 + t685 * t574;
t554 = -mrSges(5,1) * t604 + mrSges(5,3) * t591 + Ifges(5,4) * t612 + Ifges(5,2) * t611 + Ifges(5,6) * t668 - pkin(4) * t566 - t681 * t572 - t644 * t615 + t669 * t617 + t698 * t683;
t616 = Ifges(5,4) * t644 + Ifges(5,2) * t643 + Ifges(5,6) * t669;
t555 = mrSges(5,2) * t604 - mrSges(5,3) * t590 + Ifges(5,1) * t612 + Ifges(5,4) * t611 + Ifges(5,5) * t668 - t685 * t573 + t690 * t574 + t643 * t615 - t669 * t616 + (-t566 * t681 - t567 * t683) * pkin(11);
t636 = Ifges(4,5) * t659 + Ifges(4,6) * t658 + Ifges(4,3) * t673;
t638 = Ifges(4,1) * t659 + Ifges(4,4) * t658 + Ifges(4,5) * t673;
t535 = -mrSges(4,1) * t627 + mrSges(4,3) * t606 + Ifges(4,4) * t633 + Ifges(4,2) * t632 + Ifges(4,6) * t670 - pkin(3) * t697 + pkin(10) * t704 + t691 * t554 + t686 * t555 - t659 * t636 + t673 * t638;
t637 = Ifges(4,4) * t659 + Ifges(4,2) * t658 + Ifges(4,6) * t673;
t537 = mrSges(4,2) * t627 - mrSges(4,3) * t605 + Ifges(4,1) * t633 + Ifges(4,4) * t632 + Ifges(4,5) * t670 - pkin(10) * t561 - t686 * t554 + t691 * t555 + t658 * t636 - t673 * t637;
t649 = Ifges(3,3) * t675 + (Ifges(3,5) * t688 + Ifges(3,6) * t693) * t711;
t651 = Ifges(3,5) * t675 + (Ifges(3,1) * t688 + Ifges(3,4) * t693) * t711;
t533 = -mrSges(3,1) * t652 + mrSges(3,3) * t642 + Ifges(3,4) * t666 + Ifges(3,2) * t667 + Ifges(3,6) * t674 - pkin(2) * t696 + pkin(9) * t705 + t692 * t535 + t687 * t537 - t649 * t708 + t675 * t651;
t650 = Ifges(3,6) * t675 + (Ifges(3,4) * t688 + Ifges(3,2) * t693) * t711;
t534 = mrSges(3,2) * t652 - mrSges(3,3) * t641 + Ifges(3,1) * t666 + Ifges(3,4) * t667 + Ifges(3,5) * t674 - pkin(9) * t553 - t687 * t535 + t692 * t537 + t649 * t707 - t675 * t650;
t699 = pkin(8) * t545 + t533 * t693 + t534 * t688;
t536 = pkin(3) * t561 + Ifges(3,3) * t674 + t683 * t572 + Ifges(5,3) * t668 + Ifges(4,3) * t670 + t659 * t637 + Ifges(3,5) * t666 + Ifges(3,6) * t667 - t658 * t638 - t643 * t617 + t644 * t616 + Ifges(4,6) * t632 + Ifges(4,5) * t633 + mrSges(3,1) * t641 - mrSges(3,2) * t642 + Ifges(5,6) * t611 + Ifges(5,5) * t612 + mrSges(4,1) * t605 - mrSges(4,2) * t606 + mrSges(5,1) * t590 - mrSges(5,2) * t591 + pkin(2) * t553 + pkin(4) * t567 + t698 * t681 + (t650 * t688 - t651 * t693) * t711;
t532 = -mrSges(2,2) * g(3) - mrSges(2,3) * t671 + Ifges(2,5) * qJDD(1) - t695 * Ifges(2,6) - t688 * t533 + t693 * t534 + (-t540 * t682 - t541 * t684) * pkin(8);
t531 = mrSges(2,1) * g(3) + mrSges(2,3) * t672 + t695 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t540 - t682 * t536 + t699 * t684;
t1 = [-m(1) * g(1) + t706; -m(1) * g(2) + t712; (-m(1) - m(2)) * g(3) + t540; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t712 - t689 * t531 + t694 * t532; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t706 + t694 * t531 + t689 * t532; -mrSges(1,1) * g(2) + mrSges(2,1) * t671 + mrSges(1,2) * g(1) - mrSges(2,2) * t672 + Ifges(2,3) * qJDD(1) + pkin(1) * t541 + t684 * t536 + t699 * t682;];
tauB = t1;

% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:29:12
% EndTime: 2019-05-05 04:29:32
% DurationCPUTime: 19.41s
% Computational Cost: add. (318626->341), mult. (689537->442), div. (0->0), fcn. (502049->14), ass. (0->141)
t674 = sin(pkin(11));
t677 = cos(pkin(11));
t664 = g(1) * t674 - g(2) * t677;
t665 = -g(1) * t677 - g(2) * t674;
t672 = -g(3) + qJDD(1);
t686 = cos(qJ(2));
t678 = cos(pkin(6));
t682 = sin(qJ(2));
t708 = t678 * t682;
t675 = sin(pkin(6));
t709 = t675 * t682;
t628 = t664 * t708 + t686 * t665 + t672 * t709;
t688 = qJD(2) ^ 2;
t623 = -pkin(2) * t688 + qJDD(2) * pkin(8) + t628;
t644 = -t664 * t675 + t672 * t678;
t681 = sin(qJ(3));
t685 = cos(qJ(3));
t605 = -t681 * t623 + t685 * t644;
t703 = qJD(2) * qJD(3);
t702 = t685 * t703;
t662 = qJDD(2) * t681 + t702;
t600 = (-t662 + t702) * qJ(4) + (t681 * t685 * t688 + qJDD(3)) * pkin(3) + t605;
t606 = t685 * t623 + t681 * t644;
t663 = qJDD(2) * t685 - t681 * t703;
t706 = qJD(2) * t681;
t666 = qJD(3) * pkin(3) - qJ(4) * t706;
t671 = t685 ^ 2;
t601 = -pkin(3) * t671 * t688 + qJ(4) * t663 - qJD(3) * t666 + t606;
t673 = sin(pkin(12));
t676 = cos(pkin(12));
t650 = (t673 * t685 + t676 * t681) * qJD(2);
t584 = -0.2e1 * qJD(4) * t650 + t600 * t676 - t673 * t601;
t627 = -t682 * t665 + (t664 * t678 + t672 * t675) * t686;
t692 = -qJDD(2) * pkin(2) - t627;
t622 = -t688 * pkin(8) + t692;
t667 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t706;
t705 = qJD(2) * t685;
t668 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t705;
t649 = -t673 * t706 + t676 * t705;
t585 = 0.2e1 * qJD(4) * t649 + t673 * t600 + t676 * t601;
t631 = -pkin(4) * t649 - pkin(9) * t650;
t687 = qJD(3) ^ 2;
t583 = -pkin(4) * t687 + qJDD(3) * pkin(9) + t631 * t649 + t585;
t602 = -t663 * pkin(3) + qJDD(4) + t666 * t706 + (-qJ(4) * t671 - pkin(8)) * t688 + t692;
t635 = -t673 * t662 + t663 * t676;
t636 = t662 * t676 + t663 * t673;
t591 = (-qJD(3) * t649 - t636) * pkin(9) + (qJD(3) * t650 - t635) * pkin(4) + t602;
t680 = sin(qJ(5));
t684 = cos(qJ(5));
t578 = -t680 * t583 + t684 * t591;
t638 = qJD(3) * t684 - t650 * t680;
t613 = qJD(5) * t638 + qJDD(3) * t680 + t636 * t684;
t634 = qJDD(5) - t635;
t639 = qJD(3) * t680 + t650 * t684;
t648 = qJD(5) - t649;
t576 = (t638 * t648 - t613) * pkin(10) + (t638 * t639 + t634) * pkin(5) + t578;
t579 = t684 * t583 + t680 * t591;
t612 = -qJD(5) * t639 + qJDD(3) * t684 - t636 * t680;
t621 = pkin(5) * t648 - pkin(10) * t639;
t637 = t638 ^ 2;
t577 = -pkin(5) * t637 + pkin(10) * t612 - t621 * t648 + t579;
t679 = sin(qJ(6));
t683 = cos(qJ(6));
t574 = t576 * t683 - t577 * t679;
t614 = t638 * t683 - t639 * t679;
t588 = qJD(6) * t614 + t612 * t679 + t613 * t683;
t615 = t638 * t679 + t639 * t683;
t596 = -mrSges(7,1) * t614 + mrSges(7,2) * t615;
t645 = qJD(6) + t648;
t603 = -mrSges(7,2) * t645 + mrSges(7,3) * t614;
t632 = qJDD(6) + t634;
t572 = m(7) * t574 + mrSges(7,1) * t632 - mrSges(7,3) * t588 - t596 * t615 + t603 * t645;
t575 = t576 * t679 + t577 * t683;
t587 = -qJD(6) * t615 + t612 * t683 - t613 * t679;
t604 = mrSges(7,1) * t645 - mrSges(7,3) * t615;
t573 = m(7) * t575 - mrSges(7,2) * t632 + mrSges(7,3) * t587 + t596 * t614 - t604 * t645;
t564 = t683 * t572 + t679 * t573;
t616 = -mrSges(6,1) * t638 + mrSges(6,2) * t639;
t619 = -mrSges(6,2) * t648 + mrSges(6,3) * t638;
t562 = m(6) * t578 + mrSges(6,1) * t634 - mrSges(6,3) * t613 - t616 * t639 + t619 * t648 + t564;
t620 = mrSges(6,1) * t648 - mrSges(6,3) * t639;
t697 = -t572 * t679 + t683 * t573;
t563 = m(6) * t579 - mrSges(6,2) * t634 + mrSges(6,3) * t612 + t616 * t638 - t620 * t648 + t697;
t558 = t684 * t562 + t680 * t563;
t642 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t649;
t643 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t650;
t691 = m(5) * t602 - t635 * mrSges(5,1) + mrSges(5,2) * t636 - t649 * t642 + t643 * t650 + t558;
t689 = -m(4) * t622 + t663 * mrSges(4,1) - mrSges(4,2) * t662 - t667 * t706 + t668 * t705 - t691;
t554 = m(3) * t627 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t688 + t689;
t710 = t554 * t686;
t630 = -mrSges(5,1) * t649 + mrSges(5,2) * t650;
t698 = -t562 * t680 + t684 * t563;
t557 = m(5) * t585 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t635 - qJD(3) * t643 + t630 * t649 + t698;
t582 = -qJDD(3) * pkin(4) - pkin(9) * t687 + t650 * t631 - t584;
t580 = -pkin(5) * t612 - pkin(10) * t637 + t621 * t639 + t582;
t693 = m(7) * t580 - t587 * mrSges(7,1) + mrSges(7,2) * t588 - t614 * t603 + t604 * t615;
t690 = -m(6) * t582 + t612 * mrSges(6,1) - mrSges(6,2) * t613 + t638 * t619 - t620 * t639 - t693;
t568 = m(5) * t584 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t636 + qJD(3) * t642 - t630 * t650 + t690;
t549 = t673 * t557 + t676 * t568;
t661 = (-mrSges(4,1) * t685 + mrSges(4,2) * t681) * qJD(2);
t547 = m(4) * t605 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t662 + qJD(3) * t668 - t661 * t706 + t549;
t699 = t676 * t557 - t568 * t673;
t548 = m(4) * t606 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t663 - qJD(3) * t667 + t661 * t705 + t699;
t700 = -t547 * t681 + t685 * t548;
t539 = m(3) * t628 - mrSges(3,1) * t688 - qJDD(2) * mrSges(3,2) + t700;
t542 = t685 * t547 + t681 * t548;
t541 = m(3) * t644 + t542;
t530 = t539 * t708 - t541 * t675 + t678 * t710;
t528 = m(2) * t664 + t530;
t534 = t686 * t539 - t554 * t682;
t533 = m(2) * t665 + t534;
t707 = t677 * t528 + t674 * t533;
t529 = t539 * t709 + t678 * t541 + t675 * t710;
t701 = -t528 * t674 + t677 * t533;
t592 = Ifges(7,5) * t615 + Ifges(7,6) * t614 + Ifges(7,3) * t645;
t594 = Ifges(7,1) * t615 + Ifges(7,4) * t614 + Ifges(7,5) * t645;
t565 = -mrSges(7,1) * t580 + mrSges(7,3) * t575 + Ifges(7,4) * t588 + Ifges(7,2) * t587 + Ifges(7,6) * t632 - t592 * t615 + t594 * t645;
t593 = Ifges(7,4) * t615 + Ifges(7,2) * t614 + Ifges(7,6) * t645;
t566 = mrSges(7,2) * t580 - mrSges(7,3) * t574 + Ifges(7,1) * t588 + Ifges(7,4) * t587 + Ifges(7,5) * t632 + t592 * t614 - t593 * t645;
t607 = Ifges(6,5) * t639 + Ifges(6,6) * t638 + Ifges(6,3) * t648;
t609 = Ifges(6,1) * t639 + Ifges(6,4) * t638 + Ifges(6,5) * t648;
t550 = -mrSges(6,1) * t582 + mrSges(6,3) * t579 + Ifges(6,4) * t613 + Ifges(6,2) * t612 + Ifges(6,6) * t634 - pkin(5) * t693 + pkin(10) * t697 + t683 * t565 + t679 * t566 - t639 * t607 + t648 * t609;
t608 = Ifges(6,4) * t639 + Ifges(6,2) * t638 + Ifges(6,6) * t648;
t551 = mrSges(6,2) * t582 - mrSges(6,3) * t578 + Ifges(6,1) * t613 + Ifges(6,4) * t612 + Ifges(6,5) * t634 - pkin(10) * t564 - t565 * t679 + t566 * t683 + t607 * t638 - t608 * t648;
t624 = Ifges(5,5) * t650 + Ifges(5,6) * t649 + Ifges(5,3) * qJD(3);
t625 = Ifges(5,4) * t650 + Ifges(5,2) * t649 + Ifges(5,6) * qJD(3);
t535 = mrSges(5,2) * t602 - mrSges(5,3) * t584 + Ifges(5,1) * t636 + Ifges(5,4) * t635 + Ifges(5,5) * qJDD(3) - pkin(9) * t558 - qJD(3) * t625 - t550 * t680 + t551 * t684 + t624 * t649;
t626 = Ifges(5,1) * t650 + Ifges(5,4) * t649 + Ifges(5,5) * qJD(3);
t543 = Ifges(5,4) * t636 + Ifges(5,2) * t635 + Ifges(5,6) * qJDD(3) - t650 * t624 + qJD(3) * t626 - mrSges(5,1) * t602 + mrSges(5,3) * t585 - Ifges(6,5) * t613 - Ifges(6,6) * t612 - Ifges(6,3) * t634 - t639 * t608 + t638 * t609 - mrSges(6,1) * t578 + mrSges(6,2) * t579 - Ifges(7,5) * t588 - Ifges(7,6) * t587 - Ifges(7,3) * t632 - t615 * t593 + t614 * t594 - mrSges(7,1) * t574 + mrSges(7,2) * t575 - pkin(5) * t564 - pkin(4) * t558;
t652 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t681 + Ifges(4,6) * t685) * qJD(2);
t654 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t681 + Ifges(4,4) * t685) * qJD(2);
t524 = -mrSges(4,1) * t622 + mrSges(4,3) * t606 + Ifges(4,4) * t662 + Ifges(4,2) * t663 + Ifges(4,6) * qJDD(3) - pkin(3) * t691 + qJ(4) * t699 + qJD(3) * t654 + t673 * t535 + t676 * t543 - t652 * t706;
t653 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t681 + Ifges(4,2) * t685) * qJD(2);
t526 = mrSges(4,2) * t622 - mrSges(4,3) * t605 + Ifges(4,1) * t662 + Ifges(4,4) * t663 + Ifges(4,5) * qJDD(3) - qJ(4) * t549 - qJD(3) * t653 + t535 * t676 - t543 * t673 + t652 * t705;
t523 = mrSges(3,2) * t644 - mrSges(3,3) * t627 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t688 - pkin(8) * t542 - t524 * t681 + t526 * t685;
t525 = Ifges(3,6) * qJDD(2) - pkin(2) * t542 + mrSges(3,3) * t628 - mrSges(3,1) * t644 - pkin(3) * t549 - Ifges(4,5) * t662 - Ifges(4,6) * t663 - mrSges(4,1) * t605 + mrSges(4,2) * t606 - pkin(4) * t690 - pkin(9) * t698 - Ifges(5,6) * t635 - mrSges(5,1) * t584 + mrSges(5,2) * t585 - t680 * t551 - t684 * t550 - Ifges(5,5) * t636 - t650 * t625 + t649 * t626 + t688 * Ifges(3,5) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t653 * t681 + t654 * t685) * qJD(2);
t694 = pkin(7) * t534 + t523 * t682 + t525 * t686;
t522 = mrSges(3,1) * t627 - mrSges(3,2) * t628 + Ifges(3,3) * qJDD(2) + pkin(2) * t689 + pkin(8) * t700 + t685 * t524 + t681 * t526;
t521 = mrSges(2,2) * t672 - mrSges(2,3) * t664 + t686 * t523 - t682 * t525 + (-t529 * t675 - t530 * t678) * pkin(7);
t520 = -mrSges(2,1) * t672 + mrSges(2,3) * t665 - pkin(1) * t529 - t675 * t522 + t678 * t694;
t1 = [-m(1) * g(1) + t701; -m(1) * g(2) + t707; -m(1) * g(3) + m(2) * t672 + t529; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t707 - t674 * t520 + t677 * t521; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t701 + t677 * t520 + t674 * t521; -mrSges(1,1) * g(2) + mrSges(2,1) * t664 + mrSges(1,2) * g(1) - mrSges(2,2) * t665 + pkin(1) * t530 + t678 * t522 + t675 * t694;];
tauB  = t1;

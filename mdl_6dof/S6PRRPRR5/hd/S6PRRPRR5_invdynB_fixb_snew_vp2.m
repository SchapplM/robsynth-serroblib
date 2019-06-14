% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR5
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
% Datum: 2019-05-05 05:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:23:05
% EndTime: 2019-05-05 05:23:29
% DurationCPUTime: 23.33s
% Computational Cost: add. (395160->340), mult. (824921->438), div. (0->0), fcn. (603240->14), ass. (0->141)
t657 = sin(pkin(11));
t660 = cos(pkin(11));
t647 = g(1) * t657 - g(2) * t660;
t648 = -g(1) * t660 - g(2) * t657;
t655 = -g(3) + qJDD(1);
t658 = sin(pkin(6));
t661 = cos(pkin(6));
t665 = sin(qJ(2));
t669 = cos(qJ(2));
t608 = -t665 * t648 + (t647 * t661 + t655 * t658) * t669;
t671 = qJD(2) ^ 2;
t688 = t661 * t665;
t689 = t658 * t665;
t609 = t647 * t688 + t648 * t669 + t655 * t689;
t604 = -pkin(2) * t671 + qJDD(2) * pkin(8) + t609;
t626 = -t647 * t658 + t655 * t661;
t664 = sin(qJ(3));
t668 = cos(qJ(3));
t599 = t604 * t668 + t626 * t664;
t643 = (-pkin(3) * t668 - qJ(4) * t664) * qJD(2);
t670 = qJD(3) ^ 2;
t685 = qJD(2) * t668;
t583 = -pkin(3) * t670 + qJDD(3) * qJ(4) + t643 * t685 + t599;
t603 = -qJDD(2) * pkin(2) - t671 * pkin(8) - t608;
t684 = qJD(2) * qJD(3);
t683 = t668 * t684;
t645 = qJDD(2) * t664 + t683;
t654 = t664 * t684;
t646 = qJDD(2) * t668 - t654;
t589 = (-t645 - t683) * qJ(4) + (-t646 + t654) * pkin(3) + t603;
t656 = sin(pkin(12));
t659 = cos(pkin(12));
t686 = qJD(2) * t664;
t638 = qJD(3) * t656 + t659 * t686;
t572 = -0.2e1 * qJD(4) * t638 - t656 * t583 + t589 * t659;
t624 = qJDD(3) * t656 + t645 * t659;
t637 = qJD(3) * t659 - t656 * t686;
t565 = (-t637 * t685 - t624) * pkin(9) + (t637 * t638 - t646) * pkin(4) + t572;
t573 = 0.2e1 * qJD(4) * t637 + t583 * t659 + t589 * t656;
t623 = qJDD(3) * t659 - t645 * t656;
t625 = -pkin(4) * t685 - pkin(9) * t638;
t636 = t637 ^ 2;
t567 = -pkin(4) * t636 + pkin(9) * t623 + t625 * t685 + t573;
t663 = sin(qJ(5));
t667 = cos(qJ(5));
t559 = t565 * t667 - t663 * t567;
t616 = t637 * t667 - t638 * t663;
t591 = qJD(5) * t616 + t623 * t663 + t624 * t667;
t617 = t637 * t663 + t638 * t667;
t640 = qJDD(5) - t646;
t653 = qJD(5) - t685;
t557 = (t616 * t653 - t591) * pkin(10) + (t616 * t617 + t640) * pkin(5) + t559;
t560 = t565 * t663 + t567 * t667;
t590 = -qJD(5) * t617 + t623 * t667 - t624 * t663;
t607 = pkin(5) * t653 - pkin(10) * t617;
t615 = t616 ^ 2;
t558 = -pkin(5) * t615 + pkin(10) * t590 - t607 * t653 + t560;
t662 = sin(qJ(6));
t666 = cos(qJ(6));
t555 = t557 * t666 - t558 * t662;
t596 = t616 * t666 - t617 * t662;
t571 = qJD(6) * t596 + t590 * t662 + t591 * t666;
t597 = t616 * t662 + t617 * t666;
t580 = -mrSges(7,1) * t596 + mrSges(7,2) * t597;
t652 = qJD(6) + t653;
t586 = -mrSges(7,2) * t652 + mrSges(7,3) * t596;
t634 = qJDD(6) + t640;
t551 = m(7) * t555 + mrSges(7,1) * t634 - mrSges(7,3) * t571 - t580 * t597 + t586 * t652;
t556 = t557 * t662 + t558 * t666;
t570 = -qJD(6) * t597 + t590 * t666 - t591 * t662;
t587 = mrSges(7,1) * t652 - mrSges(7,3) * t597;
t552 = m(7) * t556 - mrSges(7,2) * t634 + mrSges(7,3) * t570 + t580 * t596 - t587 * t652;
t545 = t551 * t666 + t552 * t662;
t600 = -mrSges(6,1) * t616 + mrSges(6,2) * t617;
t605 = -mrSges(6,2) * t653 + mrSges(6,3) * t616;
t543 = m(6) * t559 + mrSges(6,1) * t640 - mrSges(6,3) * t591 - t600 * t617 + t605 * t653 + t545;
t606 = mrSges(6,1) * t653 - mrSges(6,3) * t617;
t678 = -t551 * t662 + t552 * t666;
t544 = m(6) * t560 - mrSges(6,2) * t640 + mrSges(6,3) * t590 + t600 * t616 - t606 * t653 + t678;
t539 = t543 * t667 + t544 * t663;
t618 = -mrSges(5,1) * t637 + mrSges(5,2) * t638;
t621 = mrSges(5,2) * t685 + mrSges(5,3) * t637;
t537 = m(5) * t572 - mrSges(5,1) * t646 - mrSges(5,3) * t624 - t618 * t638 - t621 * t685 + t539;
t622 = -mrSges(5,1) * t685 - mrSges(5,3) * t638;
t679 = -t543 * t663 + t544 * t667;
t538 = m(5) * t573 + mrSges(5,2) * t646 + mrSges(5,3) * t623 + t618 * t637 + t622 * t685 + t679;
t533 = t537 * t659 + t538 * t656;
t649 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t686;
t650 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t685;
t673 = -m(4) * t603 + mrSges(4,1) * t646 - mrSges(4,2) * t645 - t649 * t686 + t650 * t685 - t533;
t529 = m(3) * t608 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t671 + t673;
t690 = t529 * t669;
t644 = (-mrSges(4,1) * t668 + mrSges(4,2) * t664) * qJD(2);
t680 = -t537 * t656 + t538 * t659;
t532 = m(4) * t599 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t646 - qJD(3) * t649 + t644 * t685 + t680;
t598 = -t604 * t664 + t626 * t668;
t582 = -qJDD(3) * pkin(3) - qJ(4) * t670 + t643 * t686 + qJDD(4) - t598;
t574 = -pkin(4) * t623 - pkin(9) * t636 + t625 * t638 + t582;
t562 = -pkin(5) * t590 - pkin(10) * t615 + t607 * t617 + t574;
t677 = m(7) * t562 - mrSges(7,1) * t570 + mrSges(7,2) * t571 - t586 * t596 + t587 * t597;
t674 = m(6) * t574 - mrSges(6,1) * t590 + mrSges(6,2) * t591 - t605 * t616 + t606 * t617 + t677;
t672 = -m(5) * t582 + mrSges(5,1) * t623 - mrSges(5,2) * t624 + t621 * t637 - t622 * t638 - t674;
t554 = m(4) * t598 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t645 + qJD(3) * t650 - t644 * t686 + t672;
t681 = t532 * t668 - t554 * t664;
t523 = m(3) * t609 - mrSges(3,1) * t671 - qJDD(2) * mrSges(3,2) + t681;
t526 = t532 * t664 + t554 * t668;
t525 = m(3) * t626 + t526;
t512 = t523 * t688 - t525 * t658 + t661 * t690;
t510 = m(2) * t647 + t512;
t516 = t523 * t669 - t529 * t665;
t515 = m(2) * t648 + t516;
t687 = t510 * t660 + t515 * t657;
t511 = t523 * t689 + t525 * t661 + t658 * t690;
t682 = -t510 * t657 + t515 * t660;
t575 = Ifges(7,5) * t597 + Ifges(7,6) * t596 + Ifges(7,3) * t652;
t577 = Ifges(7,1) * t597 + Ifges(7,4) * t596 + Ifges(7,5) * t652;
t546 = -mrSges(7,1) * t562 + mrSges(7,3) * t556 + Ifges(7,4) * t571 + Ifges(7,2) * t570 + Ifges(7,6) * t634 - t575 * t597 + t577 * t652;
t576 = Ifges(7,4) * t597 + Ifges(7,2) * t596 + Ifges(7,6) * t652;
t547 = mrSges(7,2) * t562 - mrSges(7,3) * t555 + Ifges(7,1) * t571 + Ifges(7,4) * t570 + Ifges(7,5) * t634 + t575 * t596 - t576 * t652;
t592 = Ifges(6,5) * t617 + Ifges(6,6) * t616 + Ifges(6,3) * t653;
t594 = Ifges(6,1) * t617 + Ifges(6,4) * t616 + Ifges(6,5) * t653;
t534 = -mrSges(6,1) * t574 + mrSges(6,3) * t560 + Ifges(6,4) * t591 + Ifges(6,2) * t590 + Ifges(6,6) * t640 - pkin(5) * t677 + pkin(10) * t678 + t666 * t546 + t662 * t547 - t617 * t592 + t653 * t594;
t593 = Ifges(6,4) * t617 + Ifges(6,2) * t616 + Ifges(6,6) * t653;
t535 = mrSges(6,2) * t574 - mrSges(6,3) * t559 + Ifges(6,1) * t591 + Ifges(6,4) * t590 + Ifges(6,5) * t640 - pkin(10) * t545 - t546 * t662 + t547 * t666 + t592 * t616 - t593 * t653;
t610 = Ifges(5,5) * t638 + Ifges(5,6) * t637 - Ifges(5,3) * t685;
t612 = Ifges(5,1) * t638 + Ifges(5,4) * t637 - Ifges(5,5) * t685;
t518 = -mrSges(5,1) * t582 + mrSges(5,3) * t573 + Ifges(5,4) * t624 + Ifges(5,2) * t623 - Ifges(5,6) * t646 - pkin(4) * t674 + pkin(9) * t679 + t667 * t534 + t663 * t535 - t638 * t610 - t612 * t685;
t611 = Ifges(5,4) * t638 + Ifges(5,2) * t637 - Ifges(5,6) * t685;
t519 = mrSges(5,2) * t582 - mrSges(5,3) * t572 + Ifges(5,1) * t624 + Ifges(5,4) * t623 - Ifges(5,5) * t646 - pkin(9) * t539 - t534 * t663 + t535 * t667 + t610 * t637 + t611 * t685;
t631 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t664 + Ifges(4,6) * t668) * qJD(2);
t632 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t664 + Ifges(4,2) * t668) * qJD(2);
t508 = mrSges(4,2) * t603 - mrSges(4,3) * t598 + Ifges(4,1) * t645 + Ifges(4,4) * t646 + Ifges(4,5) * qJDD(3) - qJ(4) * t533 - qJD(3) * t632 - t518 * t656 + t519 * t659 + t631 * t685;
t633 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t664 + Ifges(4,4) * t668) * qJD(2);
t517 = Ifges(4,4) * t645 - Ifges(7,3) * t634 + t637 * t612 - t638 * t611 - Ifges(6,3) * t640 - Ifges(5,6) * t623 - Ifges(5,5) * t624 + qJD(3) * t633 + t616 * t594 - t617 * t593 - t597 * t576 + mrSges(4,3) * t599 - mrSges(4,1) * t603 - Ifges(6,6) * t590 - Ifges(6,5) * t591 + t596 * t577 - Ifges(7,6) * t570 - Ifges(7,5) * t571 - mrSges(5,1) * t572 + mrSges(5,2) * t573 + mrSges(6,2) * t560 - mrSges(6,1) * t559 - mrSges(7,1) * t555 + mrSges(7,2) * t556 + Ifges(4,6) * qJDD(3) - pkin(5) * t545 - pkin(4) * t539 - pkin(3) * t533 - t631 * t686 + (Ifges(4,2) + Ifges(5,3)) * t646;
t506 = mrSges(3,2) * t626 - mrSges(3,3) * t608 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t671 - pkin(8) * t526 + t508 * t668 - t517 * t664;
t507 = Ifges(3,6) * qJDD(2) + t671 * Ifges(3,5) - mrSges(3,1) * t626 + mrSges(3,3) * t609 - Ifges(4,5) * t645 - Ifges(4,6) * t646 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t598 + mrSges(4,2) * t599 - t656 * t519 - t659 * t518 - pkin(3) * t672 - qJ(4) * t680 - pkin(2) * t526 + (-t632 * t664 + t633 * t668) * qJD(2);
t675 = pkin(7) * t516 + t506 * t665 + t507 * t669;
t505 = mrSges(3,1) * t608 - mrSges(3,2) * t609 + Ifges(3,3) * qJDD(2) + pkin(2) * t673 + pkin(8) * t681 + t664 * t508 + t668 * t517;
t504 = mrSges(2,2) * t655 - mrSges(2,3) * t647 + t669 * t506 - t665 * t507 + (-t511 * t658 - t512 * t661) * pkin(7);
t503 = -mrSges(2,1) * t655 + mrSges(2,3) * t648 - pkin(1) * t511 - t658 * t505 + t661 * t675;
t1 = [-m(1) * g(1) + t682; -m(1) * g(2) + t687; -m(1) * g(3) + m(2) * t655 + t511; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t687 - t503 * t657 + t504 * t660; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t682 + t660 * t503 + t657 * t504; -mrSges(1,1) * g(2) + mrSges(2,1) * t647 + mrSges(1,2) * g(1) - mrSges(2,2) * t648 + pkin(1) * t512 + t661 * t505 + t658 * t675;];
tauB  = t1;

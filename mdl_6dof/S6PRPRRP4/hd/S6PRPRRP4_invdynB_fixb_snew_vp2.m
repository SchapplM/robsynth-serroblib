% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:47:48
% EndTime: 2019-05-04 23:47:57
% DurationCPUTime: 8.82s
% Computational Cost: add. (131405->296), mult. (288427->367), div. (0->0), fcn. (213558->12), ass. (0->134)
t711 = Ifges(6,1) + Ifges(7,1);
t704 = Ifges(6,4) - Ifges(7,5);
t710 = -Ifges(6,5) - Ifges(7,4);
t709 = Ifges(6,2) + Ifges(7,3);
t702 = Ifges(6,6) - Ifges(7,6);
t708 = -Ifges(6,3) - Ifges(7,2);
t667 = qJD(2) ^ 2;
t655 = sin(pkin(11));
t658 = cos(pkin(11));
t662 = sin(qJ(4));
t664 = cos(qJ(4));
t674 = t655 * t662 - t658 * t664;
t635 = t674 * qJD(2);
t656 = sin(pkin(10));
t659 = cos(pkin(10));
t642 = g(1) * t656 - g(2) * t659;
t643 = -g(1) * t659 - g(2) * t656;
t654 = -g(3) + qJDD(1);
t657 = sin(pkin(6));
t660 = cos(pkin(6));
t663 = sin(qJ(2));
t665 = cos(qJ(2));
t615 = -t663 * t643 + (t642 * t660 + t654 * t657) * t665;
t675 = t655 * t664 + t658 * t662;
t636 = t675 * qJD(2);
t688 = t636 * qJD(4);
t624 = -t674 * qJDD(2) - t688;
t707 = cos(qJ(5));
t706 = pkin(3) * t658;
t705 = -mrSges(6,3) - mrSges(7,2);
t701 = mrSges(4,2) * t655;
t671 = qJDD(3) - t615;
t606 = -qJDD(2) * pkin(2) - t667 * qJ(3) + t671;
t652 = t655 ^ 2;
t697 = t660 * t663;
t698 = t657 * t663;
t616 = t642 * t697 + t665 * t643 + t654 * t698;
t611 = -pkin(2) * t667 + qJDD(2) * qJ(3) + t616;
t632 = -t642 * t657 + t654 * t660;
t687 = qJD(2) * qJD(3);
t691 = t658 * t632 - 0.2e1 * t655 * t687;
t582 = (-pkin(8) * qJDD(2) + t667 * t706 - t611) * t655 + t691;
t585 = t655 * t632 + (t611 + 0.2e1 * t687) * t658;
t686 = qJDD(2) * t658;
t653 = t658 ^ 2;
t699 = t653 * t667;
t583 = -pkin(3) * t699 + pkin(8) * t686 + t585;
t576 = t662 * t582 + t664 * t583;
t623 = pkin(4) * t635 - pkin(9) * t636;
t666 = qJD(4) ^ 2;
t574 = -pkin(4) * t666 + qJDD(4) * pkin(9) - t623 * t635 + t576;
t598 = (-pkin(2) - t706) * qJDD(2) + (-qJ(3) + (-t652 - t653) * pkin(8)) * t667 + t671;
t689 = t635 * qJD(4);
t625 = t675 * qJDD(2) - t689;
t578 = (-t625 + t689) * pkin(9) + (-t624 + t688) * pkin(4) + t598;
t661 = sin(qJ(5));
t571 = t707 * t574 + t661 * t578;
t627 = t661 * qJD(4) + t707 * t636;
t596 = qJD(5) * t627 - t707 * qJDD(4) + t625 * t661;
t634 = qJD(5) + t635;
t609 = mrSges(6,1) * t634 - mrSges(6,3) * t627;
t622 = qJDD(5) - t624;
t626 = -t707 * qJD(4) + t636 * t661;
t601 = pkin(5) * t626 - qJ(6) * t627;
t633 = t634 ^ 2;
t567 = -pkin(5) * t633 + qJ(6) * t622 + 0.2e1 * qJD(6) * t634 - t601 * t626 + t571;
t610 = -mrSges(7,1) * t634 + mrSges(7,2) * t627;
t685 = m(7) * t567 + t622 * mrSges(7,3) + t634 * t610;
t602 = mrSges(7,1) * t626 - mrSges(7,3) * t627;
t692 = -mrSges(6,1) * t626 - mrSges(6,2) * t627 - t602;
t562 = m(6) * t571 - t622 * mrSges(6,2) + t705 * t596 - t634 * t609 + t692 * t626 + t685;
t570 = -t661 * t574 + t707 * t578;
t597 = -t626 * qJD(5) + t661 * qJDD(4) + t707 * t625;
t608 = -mrSges(6,2) * t634 - mrSges(6,3) * t626;
t568 = -t622 * pkin(5) - t633 * qJ(6) + t627 * t601 + qJDD(6) - t570;
t607 = -mrSges(7,2) * t626 + mrSges(7,3) * t634;
t680 = -m(7) * t568 + t622 * mrSges(7,1) + t634 * t607;
t564 = m(6) * t570 + t622 * mrSges(6,1) + t705 * t597 + t634 * t608 + t692 * t627 + t680;
t557 = t661 * t562 + t707 * t564;
t630 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t635;
t631 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t636;
t670 = m(5) * t598 - t624 * mrSges(5,1) + t625 * mrSges(5,2) + t635 * t630 + t636 * t631 + t557;
t669 = -m(4) * t606 + mrSges(4,1) * t686 - t670 + (t652 * t667 + t699) * mrSges(4,3);
t551 = t669 - t667 * mrSges(3,2) + m(3) * t615 + (mrSges(3,1) - t701) * qJDD(2);
t700 = t551 * t665;
t620 = mrSges(5,1) * t635 + mrSges(5,2) * t636;
t681 = t707 * t562 - t564 * t661;
t554 = m(5) * t576 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t624 - qJD(4) * t631 - t620 * t635 + t681;
t575 = t664 * t582 - t662 * t583;
t573 = -qJDD(4) * pkin(4) - t666 * pkin(9) + t636 * t623 - t575;
t569 = -0.2e1 * qJD(6) * t627 + (t626 * t634 - t597) * qJ(6) + (t627 * t634 + t596) * pkin(5) + t573;
t565 = m(7) * t569 + mrSges(7,1) * t596 - t597 * mrSges(7,3) + t607 * t626 - t627 * t610;
t668 = -m(6) * t573 - t596 * mrSges(6,1) - mrSges(6,2) * t597 - t626 * t608 - t609 * t627 - t565;
t559 = m(5) * t575 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t625 + qJD(4) * t630 - t620 * t636 + t668;
t548 = t662 * t554 + t664 * t559;
t584 = -t611 * t655 + t691;
t673 = mrSges(4,3) * qJDD(2) + t667 * (-mrSges(4,1) * t658 + t701);
t546 = m(4) * t584 - t673 * t655 + t548;
t682 = t664 * t554 - t662 * t559;
t547 = m(4) * t585 + t673 * t658 + t682;
t683 = -t546 * t655 + t658 * t547;
t537 = m(3) * t616 - mrSges(3,1) * t667 - qJDD(2) * mrSges(3,2) + t683;
t540 = t658 * t546 + t655 * t547;
t539 = m(3) * t632 + t540;
t528 = t537 * t697 - t539 * t657 + t660 * t700;
t526 = m(2) * t642 + t528;
t533 = t665 * t537 - t551 * t663;
t532 = m(2) * t643 + t533;
t696 = t659 * t526 + t656 * t532;
t695 = t626 * t709 - t627 * t704 - t634 * t702;
t694 = t626 * t702 + t627 * t710 + t634 * t708;
t693 = -t704 * t626 + t627 * t711 - t710 * t634;
t677 = Ifges(4,5) * t655 + Ifges(4,6) * t658;
t690 = t667 * t677;
t527 = t537 * t698 + t660 * t539 + t657 * t700;
t684 = -t526 * t656 + t659 * t532;
t679 = Ifges(4,1) * t655 + Ifges(4,4) * t658;
t678 = Ifges(4,4) * t655 + Ifges(4,2) * t658;
t555 = -mrSges(6,1) * t573 - mrSges(7,1) * t569 + mrSges(7,2) * t567 + mrSges(6,3) * t571 - pkin(5) * t565 - t596 * t709 + t704 * t597 + t702 * t622 + t694 * t627 + t693 * t634;
t556 = mrSges(6,2) * t573 + mrSges(7,2) * t568 - mrSges(6,3) * t570 - mrSges(7,3) * t569 - qJ(6) * t565 - t704 * t596 + t597 * t711 - t622 * t710 + t694 * t626 + t695 * t634;
t612 = Ifges(5,5) * t636 - Ifges(5,6) * t635 + Ifges(5,3) * qJD(4);
t613 = Ifges(5,4) * t636 - Ifges(5,2) * t635 + Ifges(5,6) * qJD(4);
t541 = mrSges(5,2) * t598 - mrSges(5,3) * t575 + Ifges(5,1) * t625 + Ifges(5,4) * t624 + Ifges(5,5) * qJDD(4) - pkin(9) * t557 - qJD(4) * t613 - t661 * t555 + t707 * t556 - t635 * t612;
t614 = Ifges(5,1) * t636 - Ifges(5,4) * t635 + Ifges(5,5) * qJD(4);
t542 = Ifges(5,4) * t625 + Ifges(5,2) * t624 + Ifges(5,6) * qJDD(4) - t636 * t612 + qJD(4) * t614 - mrSges(5,1) * t598 + mrSges(5,3) * t576 - mrSges(6,1) * t570 + mrSges(6,2) * t571 + mrSges(7,1) * t568 - mrSges(7,3) * t567 - pkin(5) * t680 - qJ(6) * t685 - pkin(4) * t557 + (pkin(5) * t602 + t695) * t627 + (qJ(6) * t602 - t693) * t626 + t708 * t622 + (mrSges(7,2) * pkin(5) + t710) * t597 + (mrSges(7,2) * qJ(6) + t702) * t596;
t524 = -mrSges(4,1) * t606 + mrSges(4,3) * t585 - pkin(3) * t670 + pkin(8) * t682 + t678 * qJDD(2) + t662 * t541 + t664 * t542 - t655 * t690;
t529 = mrSges(4,2) * t606 - mrSges(4,3) * t584 - pkin(8) * t548 + t679 * qJDD(2) + t664 * t541 - t662 * t542 + t658 * t690;
t522 = mrSges(3,2) * t632 - mrSges(3,3) * t615 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t667 - qJ(3) * t540 - t524 * t655 + t529 * t658;
t523 = -pkin(2) * t540 + mrSges(3,3) * t616 - mrSges(3,1) * t632 - Ifges(5,5) * t625 - Ifges(5,6) * t624 - Ifges(5,3) * qJDD(4) - t636 * t613 - t635 * t614 - mrSges(5,1) * t575 + mrSges(5,2) * t576 - t661 * t556 - t707 * t555 - pkin(4) * t668 - pkin(9) * t681 - mrSges(4,1) * t584 + mrSges(4,2) * t585 - pkin(3) * t548 + (Ifges(3,6) - t677) * qJDD(2) + (-t655 * t678 + t658 * t679 + Ifges(3,5)) * t667;
t672 = pkin(7) * t533 + t522 * t663 + t523 * t665;
t521 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t615 - mrSges(3,2) * t616 + t655 * t529 + t658 * t524 + pkin(2) * (-qJDD(2) * t701 + t669) + qJ(3) * t683;
t520 = mrSges(2,2) * t654 - mrSges(2,3) * t642 + t665 * t522 - t663 * t523 + (-t527 * t657 - t528 * t660) * pkin(7);
t519 = -mrSges(2,1) * t654 + mrSges(2,3) * t643 - pkin(1) * t527 - t657 * t521 + t672 * t660;
t1 = [-m(1) * g(1) + t684; -m(1) * g(2) + t696; -m(1) * g(3) + m(2) * t654 + t527; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t696 - t656 * t519 + t659 * t520; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t684 + t659 * t519 + t656 * t520; -mrSges(1,1) * g(2) + mrSges(2,1) * t642 + mrSges(1,2) * g(1) - mrSges(2,2) * t643 + pkin(1) * t528 + t660 * t521 + t672 * t657;];
tauB  = t1;

% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:42:36
% EndTime: 2019-05-05 07:43:00
% DurationCPUTime: 23.74s
% Computational Cost: add. (418023->340), mult. (841749->438), div. (0->0), fcn. (614154->14), ass. (0->141)
t659 = sin(pkin(11));
t662 = cos(pkin(11));
t649 = g(1) * t659 - g(2) * t662;
t650 = -g(1) * t662 - g(2) * t659;
t657 = -g(3) + qJDD(1);
t660 = sin(pkin(6));
t663 = cos(pkin(6));
t667 = sin(qJ(2));
t671 = cos(qJ(2));
t610 = -t667 * t650 + (t649 * t663 + t657 * t660) * t671;
t673 = qJD(2) ^ 2;
t690 = t663 * t667;
t691 = t660 * t667;
t611 = t649 * t690 + t671 * t650 + t657 * t691;
t606 = -pkin(2) * t673 + qJDD(2) * pkin(8) + t611;
t630 = -t649 * t660 + t657 * t663;
t666 = sin(qJ(3));
t670 = cos(qJ(3));
t601 = t670 * t606 + t666 * t630;
t646 = (-pkin(3) * t670 - pkin(9) * t666) * qJD(2);
t672 = qJD(3) ^ 2;
t687 = qJD(2) * t670;
t585 = -pkin(3) * t672 + qJDD(3) * pkin(9) + t646 * t687 + t601;
t605 = -qJDD(2) * pkin(2) - t673 * pkin(8) - t610;
t686 = qJD(2) * qJD(3);
t685 = t670 * t686;
t647 = qJDD(2) * t666 + t685;
t656 = t666 * t686;
t648 = qJDD(2) * t670 - t656;
t590 = (-t647 - t685) * pkin(9) + (-t648 + t656) * pkin(3) + t605;
t665 = sin(qJ(4));
t669 = cos(qJ(4));
t574 = -t665 * t585 + t669 * t590;
t688 = qJD(2) * t666;
t643 = qJD(3) * t669 - t665 * t688;
t621 = qJD(4) * t643 + qJDD(3) * t665 + t647 * t669;
t640 = qJDD(4) - t648;
t644 = qJD(3) * t665 + t669 * t688;
t655 = qJD(4) - t687;
t567 = (t643 * t655 - t621) * qJ(5) + (t643 * t644 + t640) * pkin(4) + t574;
t575 = t669 * t585 + t665 * t590;
t620 = -qJD(4) * t644 + qJDD(3) * t669 - t647 * t665;
t628 = pkin(4) * t655 - qJ(5) * t644;
t639 = t643 ^ 2;
t569 = -pkin(4) * t639 + qJ(5) * t620 - t628 * t655 + t575;
t658 = sin(pkin(12));
t661 = cos(pkin(12));
t624 = t643 * t658 + t644 * t661;
t561 = -0.2e1 * qJD(5) * t624 + t661 * t567 - t658 * t569;
t596 = t620 * t658 + t621 * t661;
t623 = t643 * t661 - t644 * t658;
t559 = (t623 * t655 - t596) * pkin(10) + (t623 * t624 + t640) * pkin(5) + t561;
t562 = 0.2e1 * qJD(5) * t623 + t658 * t567 + t661 * t569;
t595 = t620 * t661 - t621 * t658;
t609 = pkin(5) * t655 - pkin(10) * t624;
t622 = t623 ^ 2;
t560 = -pkin(5) * t622 + pkin(10) * t595 - t609 * t655 + t562;
t664 = sin(qJ(6));
t668 = cos(qJ(6));
t557 = t559 * t668 - t560 * t664;
t598 = t623 * t668 - t624 * t664;
t573 = qJD(6) * t598 + t595 * t664 + t596 * t668;
t599 = t623 * t664 + t624 * t668;
t582 = -mrSges(7,1) * t598 + mrSges(7,2) * t599;
t654 = qJD(6) + t655;
t588 = -mrSges(7,2) * t654 + mrSges(7,3) * t598;
t636 = qJDD(6) + t640;
t553 = m(7) * t557 + mrSges(7,1) * t636 - mrSges(7,3) * t573 - t582 * t599 + t588 * t654;
t558 = t559 * t664 + t560 * t668;
t572 = -qJD(6) * t599 + t595 * t668 - t596 * t664;
t589 = mrSges(7,1) * t654 - mrSges(7,3) * t599;
t554 = m(7) * t558 - mrSges(7,2) * t636 + mrSges(7,3) * t572 + t582 * t598 - t589 * t654;
t547 = t668 * t553 + t664 * t554;
t602 = -mrSges(6,1) * t623 + mrSges(6,2) * t624;
t607 = -mrSges(6,2) * t655 + mrSges(6,3) * t623;
t545 = m(6) * t561 + mrSges(6,1) * t640 - mrSges(6,3) * t596 - t602 * t624 + t607 * t655 + t547;
t608 = mrSges(6,1) * t655 - mrSges(6,3) * t624;
t680 = -t553 * t664 + t668 * t554;
t546 = m(6) * t562 - mrSges(6,2) * t640 + mrSges(6,3) * t595 + t602 * t623 - t608 * t655 + t680;
t541 = t661 * t545 + t658 * t546;
t625 = -mrSges(5,1) * t643 + mrSges(5,2) * t644;
t627 = -mrSges(5,2) * t655 + mrSges(5,3) * t643;
t539 = m(5) * t574 + mrSges(5,1) * t640 - mrSges(5,3) * t621 - t625 * t644 + t627 * t655 + t541;
t629 = mrSges(5,1) * t655 - mrSges(5,3) * t644;
t681 = -t545 * t658 + t661 * t546;
t540 = m(5) * t575 - mrSges(5,2) * t640 + mrSges(5,3) * t620 + t625 * t643 - t629 * t655 + t681;
t535 = t539 * t669 + t540 * t665;
t651 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t688;
t652 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t687;
t675 = -m(4) * t605 + t648 * mrSges(4,1) - mrSges(4,2) * t647 - t651 * t688 + t652 * t687 - t535;
t531 = m(3) * t610 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t673 + t675;
t692 = t531 * t671;
t645 = (-mrSges(4,1) * t670 + mrSges(4,2) * t666) * qJD(2);
t682 = -t539 * t665 + t669 * t540;
t534 = m(4) * t601 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t648 - qJD(3) * t651 + t645 * t687 + t682;
t600 = -t666 * t606 + t630 * t670;
t584 = -qJDD(3) * pkin(3) - pkin(9) * t672 + t646 * t688 - t600;
t576 = -pkin(4) * t620 - qJ(5) * t639 + t644 * t628 + qJDD(5) + t584;
t564 = -pkin(5) * t595 - pkin(10) * t622 + t609 * t624 + t576;
t679 = m(7) * t564 - t572 * mrSges(7,1) + t573 * mrSges(7,2) - t598 * t588 + t599 * t589;
t676 = m(6) * t576 - t595 * mrSges(6,1) + t596 * mrSges(6,2) - t623 * t607 + t624 * t608 + t679;
t674 = -m(5) * t584 + t620 * mrSges(5,1) - t621 * mrSges(5,2) + t643 * t627 - t644 * t629 - t676;
t556 = m(4) * t600 + qJDD(3) * mrSges(4,1) - t647 * mrSges(4,3) + qJD(3) * t652 - t645 * t688 + t674;
t683 = t670 * t534 - t556 * t666;
t525 = m(3) * t611 - mrSges(3,1) * t673 - qJDD(2) * mrSges(3,2) + t683;
t528 = t666 * t534 + t670 * t556;
t527 = m(3) * t630 + t528;
t514 = t525 * t690 - t527 * t660 + t663 * t692;
t512 = m(2) * t649 + t514;
t518 = t671 * t525 - t531 * t667;
t517 = m(2) * t650 + t518;
t689 = t662 * t512 + t659 * t517;
t513 = t525 * t691 + t663 * t527 + t660 * t692;
t684 = -t512 * t659 + t662 * t517;
t577 = Ifges(7,5) * t599 + Ifges(7,6) * t598 + Ifges(7,3) * t654;
t579 = Ifges(7,1) * t599 + Ifges(7,4) * t598 + Ifges(7,5) * t654;
t548 = -mrSges(7,1) * t564 + mrSges(7,3) * t558 + Ifges(7,4) * t573 + Ifges(7,2) * t572 + Ifges(7,6) * t636 - t577 * t599 + t579 * t654;
t578 = Ifges(7,4) * t599 + Ifges(7,2) * t598 + Ifges(7,6) * t654;
t549 = mrSges(7,2) * t564 - mrSges(7,3) * t557 + Ifges(7,1) * t573 + Ifges(7,4) * t572 + Ifges(7,5) * t636 + t577 * t598 - t578 * t654;
t592 = Ifges(6,5) * t624 + Ifges(6,6) * t623 + Ifges(6,3) * t655;
t594 = Ifges(6,1) * t624 + Ifges(6,4) * t623 + Ifges(6,5) * t655;
t536 = -mrSges(6,1) * t576 + mrSges(6,3) * t562 + Ifges(6,4) * t596 + Ifges(6,2) * t595 + Ifges(6,6) * t640 - pkin(5) * t679 + pkin(10) * t680 + t668 * t548 + t664 * t549 - t624 * t592 + t655 * t594;
t593 = Ifges(6,4) * t624 + Ifges(6,2) * t623 + Ifges(6,6) * t655;
t537 = mrSges(6,2) * t576 - mrSges(6,3) * t561 + Ifges(6,1) * t596 + Ifges(6,4) * t595 + Ifges(6,5) * t640 - pkin(10) * t547 - t548 * t664 + t549 * t668 + t592 * t623 - t593 * t655;
t612 = Ifges(5,5) * t644 + Ifges(5,6) * t643 + Ifges(5,3) * t655;
t614 = Ifges(5,1) * t644 + Ifges(5,4) * t643 + Ifges(5,5) * t655;
t520 = -mrSges(5,1) * t584 + mrSges(5,3) * t575 + Ifges(5,4) * t621 + Ifges(5,2) * t620 + Ifges(5,6) * t640 - pkin(4) * t676 + qJ(5) * t681 + t661 * t536 + t658 * t537 - t644 * t612 + t655 * t614;
t613 = Ifges(5,4) * t644 + Ifges(5,2) * t643 + Ifges(5,6) * t655;
t521 = mrSges(5,2) * t584 - mrSges(5,3) * t574 + Ifges(5,1) * t621 + Ifges(5,4) * t620 + Ifges(5,5) * t640 - qJ(5) * t541 - t536 * t658 + t537 * t661 + t612 * t643 - t613 * t655;
t633 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t666 + Ifges(4,6) * t670) * qJD(2);
t634 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t666 + Ifges(4,2) * t670) * qJD(2);
t510 = mrSges(4,2) * t605 - mrSges(4,3) * t600 + Ifges(4,1) * t647 + Ifges(4,4) * t648 + Ifges(4,5) * qJDD(3) - pkin(9) * t535 - qJD(3) * t634 - t520 * t665 + t521 * t669 + t633 * t687;
t635 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t666 + Ifges(4,4) * t670) * qJD(2);
t519 = Ifges(4,6) * qJDD(3) + (-Ifges(5,3) - Ifges(6,3)) * t640 + Ifges(4,4) * t647 + Ifges(4,2) * t648 + t643 * t614 - t644 * t613 - t624 * t593 + qJD(3) * t635 - Ifges(7,3) * t636 - t633 * t688 - Ifges(5,6) * t620 - Ifges(5,5) * t621 + t623 * t594 - t599 * t578 + mrSges(4,3) * t601 - mrSges(4,1) * t605 - Ifges(6,6) * t595 - Ifges(6,5) * t596 + t598 * t579 - Ifges(7,6) * t572 - Ifges(7,5) * t573 - mrSges(5,1) * t574 + mrSges(5,2) * t575 + mrSges(6,2) * t562 - mrSges(6,1) * t561 - mrSges(7,1) * t557 + mrSges(7,2) * t558 - pkin(5) * t547 - pkin(4) * t541 - pkin(3) * t535;
t508 = mrSges(3,2) * t630 - mrSges(3,3) * t610 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t673 - pkin(8) * t528 + t510 * t670 - t519 * t666;
t509 = Ifges(3,6) * qJDD(2) + t673 * Ifges(3,5) - mrSges(3,1) * t630 + mrSges(3,3) * t611 - Ifges(4,5) * t647 - Ifges(4,6) * t648 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t600 + mrSges(4,2) * t601 - t665 * t521 - t669 * t520 - pkin(3) * t674 - pkin(9) * t682 - pkin(2) * t528 + (-t634 * t666 + t635 * t670) * qJD(2);
t677 = pkin(7) * t518 + t508 * t667 + t509 * t671;
t507 = mrSges(3,1) * t610 - mrSges(3,2) * t611 + Ifges(3,3) * qJDD(2) + pkin(2) * t675 + pkin(8) * t683 + t666 * t510 + t670 * t519;
t506 = mrSges(2,2) * t657 - mrSges(2,3) * t649 + t671 * t508 - t667 * t509 + (-t513 * t660 - t514 * t663) * pkin(7);
t505 = -mrSges(2,1) * t657 + mrSges(2,3) * t650 - pkin(1) * t513 - t660 * t507 + t663 * t677;
t1 = [-m(1) * g(1) + t684; -m(1) * g(2) + t689; -m(1) * g(3) + m(2) * t657 + t513; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t689 - t659 * t505 + t662 * t506; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t684 + t662 * t505 + t659 * t506; -mrSges(1,1) * g(2) + mrSges(2,1) * t649 + mrSges(1,2) * g(1) - mrSges(2,2) * t650 + pkin(1) * t514 + t663 * t507 + t660 * t677;];
tauB  = t1;

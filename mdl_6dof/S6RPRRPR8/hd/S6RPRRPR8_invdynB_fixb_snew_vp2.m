% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:08:49
% EndTime: 2019-05-05 23:09:01
% DurationCPUTime: 8.36s
% Computational Cost: add. (135196->337), mult. (274399->413), div. (0->0), fcn. (182751->10), ass. (0->131)
t643 = sin(qJ(1));
t647 = cos(qJ(1));
t626 = -t647 * g(1) - t643 * g(2);
t672 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t626;
t671 = -pkin(1) - pkin(7);
t670 = mrSges(2,1) - mrSges(3,2);
t669 = -Ifges(3,4) + Ifges(2,5);
t668 = (Ifges(3,5) - Ifges(2,6));
t625 = t643 * g(1) - t647 * g(2);
t649 = qJD(1) ^ 2;
t655 = -t649 * qJ(2) + qJDD(2) - t625;
t604 = qJDD(1) * t671 + t655;
t642 = sin(qJ(3));
t646 = cos(qJ(3));
t597 = -t646 * g(3) + t642 * t604;
t619 = (mrSges(4,1) * t642 + mrSges(4,2) * t646) * qJD(1);
t665 = qJD(1) * qJD(3);
t629 = t646 * t665;
t621 = -t642 * qJDD(1) - t629;
t666 = qJD(1) * t646;
t624 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t666;
t631 = t642 * qJD(1);
t603 = t649 * t671 - t672;
t663 = t642 * t665;
t622 = t646 * qJDD(1) - t663;
t575 = (-t622 + t663) * pkin(8) + (-t621 + t629) * pkin(3) + t603;
t620 = (pkin(3) * t642 - pkin(8) * t646) * qJD(1);
t648 = qJD(3) ^ 2;
t578 = -t648 * pkin(3) + qJDD(3) * pkin(8) - t620 * t631 + t597;
t641 = sin(qJ(4));
t645 = cos(qJ(4));
t559 = t645 * t575 - t641 * t578;
t617 = t645 * qJD(3) - t641 * t666;
t591 = t617 * qJD(4) + t641 * qJDD(3) + t645 * t622;
t616 = qJDD(4) - t621;
t618 = t641 * qJD(3) + t645 * t666;
t628 = t631 + qJD(4);
t549 = (t617 * t628 - t591) * qJ(5) + (t617 * t618 + t616) * pkin(4) + t559;
t560 = t641 * t575 + t645 * t578;
t590 = -t618 * qJD(4) + t645 * qJDD(3) - t641 * t622;
t599 = t628 * pkin(4) - t618 * qJ(5);
t615 = t617 ^ 2;
t551 = -t615 * pkin(4) + t590 * qJ(5) - t628 * t599 + t560;
t638 = sin(pkin(10));
t639 = cos(pkin(10));
t594 = t638 * t617 + t639 * t618;
t539 = -0.2e1 * qJD(5) * t594 + t639 * t549 - t638 * t551;
t568 = t638 * t590 + t639 * t591;
t593 = t639 * t617 - t638 * t618;
t537 = (t593 * t628 - t568) * pkin(9) + (t593 * t594 + t616) * pkin(5) + t539;
t540 = 0.2e1 * qJD(5) * t593 + t638 * t549 + t639 * t551;
t567 = t639 * t590 - t638 * t591;
t581 = t628 * pkin(5) - t594 * pkin(9);
t592 = t593 ^ 2;
t538 = -t592 * pkin(5) + t567 * pkin(9) - t628 * t581 + t540;
t640 = sin(qJ(6));
t644 = cos(qJ(6));
t535 = t644 * t537 - t640 * t538;
t570 = t644 * t593 - t640 * t594;
t546 = t570 * qJD(6) + t640 * t567 + t644 * t568;
t571 = t640 * t593 + t644 * t594;
t557 = -t570 * mrSges(7,1) + t571 * mrSges(7,2);
t627 = qJD(6) + t628;
t561 = -t627 * mrSges(7,2) + t570 * mrSges(7,3);
t611 = qJDD(6) + t616;
t532 = m(7) * t535 + t611 * mrSges(7,1) - t546 * mrSges(7,3) - t571 * t557 + t627 * t561;
t536 = t640 * t537 + t644 * t538;
t545 = -t571 * qJD(6) + t644 * t567 - t640 * t568;
t562 = t627 * mrSges(7,1) - t571 * mrSges(7,3);
t533 = m(7) * t536 - t611 * mrSges(7,2) + t545 * mrSges(7,3) + t570 * t557 - t627 * t562;
t526 = t644 * t532 + t640 * t533;
t572 = -t593 * mrSges(6,1) + t594 * mrSges(6,2);
t579 = -t628 * mrSges(6,2) + t593 * mrSges(6,3);
t524 = m(6) * t539 + t616 * mrSges(6,1) - t568 * mrSges(6,3) - t594 * t572 + t628 * t579 + t526;
t580 = t628 * mrSges(6,1) - t594 * mrSges(6,3);
t658 = -t640 * t532 + t644 * t533;
t525 = m(6) * t540 - t616 * mrSges(6,2) + t567 * mrSges(6,3) + t593 * t572 - t628 * t580 + t658;
t520 = t639 * t524 + t638 * t525;
t595 = -t617 * mrSges(5,1) + t618 * mrSges(5,2);
t598 = -t628 * mrSges(5,2) + t617 * mrSges(5,3);
t518 = m(5) * t559 + t616 * mrSges(5,1) - t591 * mrSges(5,3) - t618 * t595 + t628 * t598 + t520;
t600 = t628 * mrSges(5,1) - t618 * mrSges(5,3);
t659 = -t638 * t524 + t639 * t525;
t519 = m(5) * t560 - t616 * mrSges(5,2) + t590 * mrSges(5,3) + t617 * t595 - t628 * t600 + t659;
t660 = -t641 * t518 + t645 * t519;
t511 = m(4) * t597 - qJDD(3) * mrSges(4,2) + t621 * mrSges(4,3) - qJD(3) * t624 - t619 * t631 + t660;
t596 = t642 * g(3) + t646 * t604;
t623 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t631;
t577 = -qJDD(3) * pkin(3) - t648 * pkin(8) + t620 * t666 - t596;
t558 = -t590 * pkin(4) - t615 * qJ(5) + t618 * t599 + qJDD(5) + t577;
t542 = -t567 * pkin(5) - t592 * pkin(9) + t594 * t581 + t558;
t657 = m(7) * t542 - t545 * mrSges(7,1) + t546 * mrSges(7,2) - t570 * t561 + t571 * t562;
t651 = m(6) * t558 - t567 * mrSges(6,1) + t568 * mrSges(6,2) - t593 * t579 + t594 * t580 + t657;
t650 = -m(5) * t577 + t590 * mrSges(5,1) - t591 * mrSges(5,2) + t617 * t598 - t618 * t600 - t651;
t534 = m(4) * t596 + qJDD(3) * mrSges(4,1) - t622 * mrSges(4,3) + qJD(3) * t623 - t619 * t666 + t650;
t506 = t642 * t511 + t646 * t534;
t606 = -qJDD(1) * pkin(1) + t655;
t654 = -m(3) * t606 + (t649 * mrSges(3,3)) - t506;
t504 = m(2) * t625 - (t649 * mrSges(2,2)) + qJDD(1) * t670 + t654;
t605 = t649 * pkin(1) + t672;
t512 = t645 * t518 + t641 * t519;
t653 = -m(4) * t603 + t621 * mrSges(4,1) - t622 * mrSges(4,2) - t623 * t631 - t624 * t666 - t512;
t652 = -m(3) * t605 + (t649 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t653;
t509 = m(2) * t626 - (t649 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t652;
t667 = t647 * t504 + t643 * t509;
t662 = -t643 * t504 + t647 * t509;
t661 = t646 * t511 - t642 * t534;
t610 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t646 - Ifges(4,4) * t642) * qJD(1);
t609 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t646 - Ifges(4,2) * t642) * qJD(1);
t608 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t646 - Ifges(4,6) * t642) * qJD(1);
t584 = Ifges(5,1) * t618 + Ifges(5,4) * t617 + Ifges(5,5) * t628;
t583 = Ifges(5,4) * t618 + Ifges(5,2) * t617 + Ifges(5,6) * t628;
t582 = Ifges(5,5) * t618 + Ifges(5,6) * t617 + Ifges(5,3) * t628;
t566 = Ifges(6,1) * t594 + Ifges(6,4) * t593 + Ifges(6,5) * t628;
t565 = Ifges(6,4) * t594 + Ifges(6,2) * t593 + Ifges(6,6) * t628;
t564 = Ifges(6,5) * t594 + Ifges(6,6) * t593 + Ifges(6,3) * t628;
t554 = Ifges(7,1) * t571 + Ifges(7,4) * t570 + Ifges(7,5) * t627;
t553 = Ifges(7,4) * t571 + Ifges(7,2) * t570 + Ifges(7,6) * t627;
t552 = Ifges(7,5) * t571 + Ifges(7,6) * t570 + Ifges(7,3) * t627;
t528 = mrSges(7,2) * t542 - mrSges(7,3) * t535 + Ifges(7,1) * t546 + Ifges(7,4) * t545 + Ifges(7,5) * t611 + t570 * t552 - t627 * t553;
t527 = -mrSges(7,1) * t542 + mrSges(7,3) * t536 + Ifges(7,4) * t546 + Ifges(7,2) * t545 + Ifges(7,6) * t611 - t571 * t552 + t627 * t554;
t514 = mrSges(6,2) * t558 - mrSges(6,3) * t539 + Ifges(6,1) * t568 + Ifges(6,4) * t567 + Ifges(6,5) * t616 - pkin(9) * t526 - t640 * t527 + t644 * t528 + t593 * t564 - t628 * t565;
t513 = -mrSges(6,1) * t558 + mrSges(6,3) * t540 + Ifges(6,4) * t568 + Ifges(6,2) * t567 + Ifges(6,6) * t616 - pkin(5) * t657 + pkin(9) * t658 + t644 * t527 + t640 * t528 - t594 * t564 + t628 * t566;
t505 = -m(3) * g(3) + t661;
t502 = mrSges(5,2) * t577 - mrSges(5,3) * t559 + Ifges(5,1) * t591 + Ifges(5,4) * t590 + Ifges(5,5) * t616 - qJ(5) * t520 - t638 * t513 + t639 * t514 + t617 * t582 - t628 * t583;
t501 = -mrSges(5,1) * t577 + mrSges(5,3) * t560 + Ifges(5,4) * t591 + Ifges(5,2) * t590 + Ifges(5,6) * t616 - pkin(4) * t651 + qJ(5) * t659 + t639 * t513 + t638 * t514 - t618 * t582 + t628 * t584;
t500 = Ifges(4,6) * qJDD(3) - t608 * t666 + (-Ifges(5,3) - Ifges(6,3)) * t616 + Ifges(4,4) * t622 + t617 * t584 - t618 * t583 + Ifges(4,2) * t621 - mrSges(4,1) * t603 + qJD(3) * t610 - Ifges(7,3) * t611 - Ifges(5,5) * t591 + t593 * t566 - t594 * t565 + mrSges(4,3) * t597 - Ifges(5,6) * t590 + t570 * t554 - t571 * t553 - Ifges(6,6) * t567 - Ifges(6,5) * t568 - mrSges(5,1) * t559 + mrSges(5,2) * t560 - Ifges(7,6) * t545 - Ifges(7,5) * t546 + mrSges(6,2) * t540 - mrSges(6,1) * t539 + mrSges(7,2) * t536 - mrSges(7,1) * t535 - pkin(5) * t526 - pkin(4) * t520 - pkin(3) * t512;
t499 = mrSges(4,2) * t603 - mrSges(4,3) * t596 + Ifges(4,1) * t622 + Ifges(4,4) * t621 + Ifges(4,5) * qJDD(3) - pkin(8) * t512 - qJD(3) * t609 - t641 * t501 + t645 * t502 - t608 * t631;
t498 = -qJ(2) * t505 - mrSges(2,3) * t625 + pkin(2) * t506 + mrSges(3,1) * t606 + t645 * t501 + pkin(3) * t650 + pkin(8) * t660 + t641 * t502 + Ifges(4,5) * t622 + Ifges(4,6) * t621 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t596 - mrSges(4,2) * t597 + (t668 * t649) + t669 * qJDD(1) + (t646 * t609 + t642 * t610) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t497 = -mrSges(3,1) * t605 + mrSges(2,3) * t626 - pkin(1) * t505 - pkin(2) * t653 - pkin(7) * t661 + g(3) * t670 - qJDD(1) * t668 - t642 * t499 - t646 * t500 + t649 * t669;
t1 = [-m(1) * g(1) + t662; -m(1) * g(2) + t667; (-m(1) - m(2) - m(3)) * g(3) + t661; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t667 - t643 * t497 + t647 * t498; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t662 + t647 * t497 + t643 * t498; pkin(1) * t654 + qJ(2) * t652 - pkin(7) * t506 + mrSges(2,1) * t625 - mrSges(2,2) * t626 + t646 * t499 - t642 * t500 + mrSges(3,2) * t606 - mrSges(3,3) * t605 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;

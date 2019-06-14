% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 17:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:57:06
% EndTime: 2019-05-05 17:57:13
% DurationCPUTime: 4.60s
% Computational Cost: add. (45525->316), mult. (98767->374), div. (0->0), fcn. (63246->8), ass. (0->122)
t688 = -2 * qJD(4);
t687 = Ifges(6,1) + Ifges(7,1);
t680 = Ifges(6,4) + Ifges(7,4);
t678 = Ifges(6,5) + Ifges(7,5);
t686 = Ifges(6,2) + Ifges(7,2);
t685 = Ifges(6,6) + Ifges(7,6);
t684 = Ifges(6,3) + Ifges(7,3);
t642 = sin(qJ(1));
t645 = cos(qJ(1));
t627 = g(1) * t642 - t645 * g(2);
t647 = qJD(1) ^ 2;
t653 = -qJ(2) * t647 + qJDD(2) - t627;
t683 = -pkin(1) - pkin(7);
t606 = qJDD(1) * t683 + t653;
t641 = sin(qJ(3));
t644 = cos(qJ(3));
t595 = t641 * g(3) + t644 * t606;
t667 = qJD(1) * qJD(3);
t663 = t641 * t667;
t623 = qJDD(1) * t644 - t663;
t572 = (-t623 - t663) * qJ(4) + (-t641 * t644 * t647 + qJDD(3)) * pkin(3) + t595;
t596 = -g(3) * t644 + t641 * t606;
t622 = -qJDD(1) * t641 - t644 * t667;
t669 = qJD(1) * t644;
t625 = qJD(3) * pkin(3) - qJ(4) * t669;
t635 = t641 ^ 2;
t573 = -pkin(3) * t635 * t647 + qJ(4) * t622 - qJD(3) * t625 + t596;
t638 = sin(pkin(9));
t639 = cos(pkin(9));
t670 = qJD(1) * t641;
t613 = -t638 * t670 + t639 * t669;
t552 = t572 * t639 - t638 * t573 + t613 * t688;
t628 = -g(1) * t645 - g(2) * t642;
t654 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t628;
t612 = (t638 * t644 + t639 * t641) * qJD(1);
t682 = mrSges(2,1) - mrSges(3,2);
t681 = -mrSges(6,2) - mrSges(7,2);
t679 = Ifges(2,5) - Ifges(3,4);
t677 = -Ifges(2,6) + Ifges(3,5);
t553 = t638 * t572 + t639 * t573 + t612 * t688;
t588 = mrSges(5,1) * t612 + mrSges(5,2) * t613;
t593 = t622 * t639 - t623 * t638;
t605 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t613;
t589 = pkin(4) * t612 - pkin(8) * t613;
t646 = qJD(3) ^ 2;
t548 = -pkin(4) * t646 + qJDD(3) * pkin(8) - t589 * t612 + t553;
t575 = -pkin(3) * t622 + qJDD(4) + t625 * t669 + (-qJ(4) * t635 + t683) * t647 + t654;
t594 = t622 * t638 + t623 * t639;
t551 = (qJD(3) * t612 - t594) * pkin(8) + (qJD(3) * t613 - t593) * pkin(4) + t575;
t640 = sin(qJ(5));
t643 = cos(qJ(5));
t544 = -t548 * t640 + t643 * t551;
t598 = qJD(3) * t643 - t613 * t640;
t567 = qJD(5) * t598 + qJDD(3) * t640 + t594 * t643;
t599 = qJD(3) * t640 + t613 * t643;
t576 = -mrSges(7,1) * t598 + mrSges(7,2) * t599;
t577 = -mrSges(6,1) * t598 + mrSges(6,2) * t599;
t610 = qJD(5) + t612;
t580 = -mrSges(6,2) * t610 + mrSges(6,3) * t598;
t592 = qJDD(5) - t593;
t540 = -0.2e1 * qJD(6) * t599 + (t598 * t610 - t567) * qJ(6) + (t598 * t599 + t592) * pkin(5) + t544;
t579 = -mrSges(7,2) * t610 + mrSges(7,3) * t598;
t665 = m(7) * t540 + t592 * mrSges(7,1) + t610 * t579;
t532 = m(6) * t544 + mrSges(6,1) * t592 + t580 * t610 + (-t576 - t577) * t599 + (-mrSges(6,3) - mrSges(7,3)) * t567 + t665;
t545 = t643 * t548 + t640 * t551;
t566 = -qJD(5) * t599 + qJDD(3) * t643 - t594 * t640;
t581 = pkin(5) * t610 - qJ(6) * t599;
t597 = t598 ^ 2;
t542 = -pkin(5) * t597 + qJ(6) * t566 + 0.2e1 * qJD(6) * t598 - t581 * t610 + t545;
t664 = m(7) * t542 + t566 * mrSges(7,3) + t598 * t576;
t582 = mrSges(7,1) * t610 - mrSges(7,3) * t599;
t671 = -mrSges(6,1) * t610 + mrSges(6,3) * t599 - t582;
t535 = m(6) * t545 + mrSges(6,3) * t566 + t577 * t598 + t592 * t681 + t610 * t671 + t664;
t659 = -t532 * t640 + t643 * t535;
t528 = m(5) * t553 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t593 - qJD(3) * t605 - t588 * t612 + t659;
t604 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t612;
t547 = -qJDD(3) * pkin(4) - pkin(8) * t646 + t613 * t589 - t552;
t543 = -pkin(5) * t566 - qJ(6) * t597 + t581 * t599 + qJDD(6) + t547;
t657 = m(7) * t543 - t566 * mrSges(7,1) - t598 * t579;
t649 = -m(6) * t547 + t566 * mrSges(6,1) + t567 * t681 + t598 * t580 + t599 * t671 - t657;
t537 = m(5) * t552 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t594 + qJD(3) * t604 - t588 * t613 + t649;
t521 = t638 * t528 + t639 * t537;
t621 = (mrSges(4,1) * t641 + mrSges(4,2) * t644) * qJD(1);
t624 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t670;
t519 = m(4) * t595 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t623 + qJD(3) * t624 - t621 * t669 + t521;
t626 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t669;
t660 = t639 * t528 - t537 * t638;
t520 = m(4) * t596 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t622 - qJD(3) * t626 - t621 * t670 + t660;
t516 = t519 * t644 + t520 * t641;
t611 = -qJDD(1) * pkin(1) + t653;
t652 = -m(3) * t611 + t647 * mrSges(3,3) - t516;
t514 = m(2) * t627 - mrSges(2,2) * t647 + qJDD(1) * t682 + t652;
t607 = pkin(1) * t647 - t654;
t603 = t647 * t683 + t654;
t530 = t643 * t532 + t640 * t535;
t651 = m(5) * t575 - mrSges(5,1) * t593 + t594 * mrSges(5,2) + t604 * t612 + t613 * t605 + t530;
t650 = -m(4) * t603 + mrSges(4,1) * t622 - t623 * mrSges(4,2) - t624 * t670 - t626 * t669 - t651;
t648 = -m(3) * t607 + t647 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t650;
t525 = m(2) * t628 - mrSges(2,1) * t647 - qJDD(1) * mrSges(2,2) + t648;
t675 = t645 * t514 + t642 * t525;
t674 = t598 * t685 + t599 * t678 + t610 * t684;
t673 = -t598 * t686 - t599 * t680 - t610 * t685;
t672 = t598 * t680 + t599 * t687 + t610 * t678;
t662 = -t514 * t642 + t645 * t525;
t661 = -t519 * t641 + t644 * t520;
t616 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t644 - Ifges(4,4) * t641) * qJD(1);
t615 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t644 - Ifges(4,2) * t641) * qJD(1);
t614 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t644 - Ifges(4,6) * t641) * qJD(1);
t586 = Ifges(5,1) * t613 - Ifges(5,4) * t612 + Ifges(5,5) * qJD(3);
t585 = Ifges(5,4) * t613 - Ifges(5,2) * t612 + Ifges(5,6) * qJD(3);
t584 = Ifges(5,5) * t613 - Ifges(5,6) * t612 + Ifges(5,3) * qJD(3);
t538 = -mrSges(7,3) * t567 - t576 * t599 + t665;
t529 = mrSges(6,2) * t547 + mrSges(7,2) * t543 - mrSges(6,3) * t544 - mrSges(7,3) * t540 - qJ(6) * t538 + t680 * t566 + t567 * t687 + t678 * t592 + t674 * t598 + t673 * t610;
t522 = -mrSges(6,1) * t547 + mrSges(6,3) * t545 - mrSges(7,1) * t543 + mrSges(7,3) * t542 - pkin(5) * t657 + qJ(6) * t664 + (-qJ(6) * t582 + t672) * t610 + (-pkin(5) * t582 - t674) * t599 + (-mrSges(7,2) * qJ(6) + t685) * t592 + (-mrSges(7,2) * pkin(5) + t680) * t567 + t686 * t566;
t517 = -mrSges(5,1) * t575 - mrSges(6,1) * t544 - mrSges(7,1) * t540 + mrSges(6,2) * t545 + mrSges(7,2) * t542 + mrSges(5,3) * t553 + Ifges(5,4) * t594 + Ifges(5,2) * t593 + Ifges(5,6) * qJDD(3) - pkin(4) * t530 - pkin(5) * t538 + qJD(3) * t586 - t613 * t584 + t673 * t599 + t672 * t598 - t684 * t592 - t678 * t567 - t685 * t566;
t515 = -m(3) * g(3) + t661;
t512 = mrSges(5,2) * t575 - mrSges(5,3) * t552 + Ifges(5,1) * t594 + Ifges(5,4) * t593 + Ifges(5,5) * qJDD(3) - pkin(8) * t530 - qJD(3) * t585 - t522 * t640 + t529 * t643 - t584 * t612;
t511 = mrSges(4,2) * t603 - mrSges(4,3) * t595 + Ifges(4,1) * t623 + Ifges(4,4) * t622 + Ifges(4,5) * qJDD(3) - qJ(4) * t521 - qJD(3) * t615 + t512 * t639 - t517 * t638 - t614 * t670;
t510 = -mrSges(4,1) * t603 + mrSges(4,3) * t596 + Ifges(4,4) * t623 + Ifges(4,2) * t622 + Ifges(4,6) * qJDD(3) - pkin(3) * t651 + qJ(4) * t660 + qJD(3) * t616 + t638 * t512 + t639 * t517 - t614 * t669;
t509 = pkin(8) * t659 + t640 * t529 + t643 * t522 - mrSges(2,3) * t627 + t613 * t585 + Ifges(4,6) * t622 + Ifges(4,5) * t623 + mrSges(3,1) * t611 + t612 * t586 + pkin(4) * t649 + Ifges(5,5) * t594 + mrSges(4,1) * t595 - mrSges(4,2) * t596 + Ifges(5,6) * t593 + mrSges(5,1) * t552 - mrSges(5,2) * t553 + pkin(3) * t521 + pkin(2) * t516 - qJ(2) * t515 + t677 * t647 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t679 * qJDD(1) + (t615 * t644 + t616 * t641) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t508 = -mrSges(3,1) * t607 + mrSges(2,3) * t628 - pkin(1) * t515 - pkin(2) * t650 - pkin(7) * t661 + g(3) * t682 - qJDD(1) * t677 - t644 * t510 - t641 * t511 + t647 * t679;
t1 = [-m(1) * g(1) + t662; -m(1) * g(2) + t675; (-m(1) - m(2) - m(3)) * g(3) + t661; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t675 - t642 * t508 + t645 * t509; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t662 + t645 * t508 + t642 * t509; pkin(1) * t652 + qJ(2) * t648 + t644 * t511 - t641 * t510 - pkin(7) * t516 + mrSges(2,1) * t627 - mrSges(2,2) * t628 + mrSges(3,2) * t611 - mrSges(3,3) * t607 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;

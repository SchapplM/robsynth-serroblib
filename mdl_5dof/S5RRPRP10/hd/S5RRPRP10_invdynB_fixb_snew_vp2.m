% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:37
% EndTime: 2019-12-31 20:09:40
% DurationCPUTime: 1.84s
% Computational Cost: add. (13084->279), mult. (26625->315), div. (0->0), fcn. (13911->6), ass. (0->109)
t643 = -2 * qJD(3);
t642 = Ifges(3,1) + Ifges(4,2);
t641 = Ifges(5,1) + Ifges(6,1);
t632 = Ifges(3,4) + Ifges(4,6);
t631 = Ifges(5,4) + Ifges(6,4);
t630 = Ifges(3,5) - Ifges(4,4);
t629 = Ifges(5,5) + Ifges(6,5);
t640 = Ifges(3,2) + Ifges(4,3);
t639 = Ifges(5,2) + Ifges(6,2);
t628 = Ifges(3,6) - Ifges(4,5);
t627 = Ifges(5,6) + Ifges(6,6);
t638 = Ifges(3,3) + Ifges(4,1);
t637 = Ifges(5,3) + Ifges(6,3);
t590 = sin(qJ(4));
t593 = cos(qJ(4));
t594 = cos(qJ(2));
t616 = qJD(1) * t594;
t565 = t593 * qJD(2) - t590 * t616;
t591 = sin(qJ(2));
t614 = qJD(1) * qJD(2);
t609 = t591 * t614;
t570 = t594 * qJDD(1) - t609;
t531 = -t565 * qJD(4) - t590 * qJDD(2) - t593 * t570;
t636 = (mrSges(5,1) + mrSges(6,1)) * t531;
t592 = sin(qJ(1));
t595 = cos(qJ(1));
t579 = -t595 * g(1) - t592 * g(2);
t597 = qJD(1) ^ 2;
t555 = -t597 * pkin(1) + qJDD(1) * pkin(6) + t579;
t542 = -t591 * g(3) + t594 * t555;
t566 = (-pkin(2) * t594 - qJ(3) * t591) * qJD(1);
t596 = qJD(2) ^ 2;
t517 = t596 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t643 - t566 * t616 - t542;
t564 = -t590 * qJD(2) - t593 * t616;
t615 = t591 * qJD(1);
t581 = qJD(4) + t615;
t536 = -t581 * mrSges(6,2) + t564 * mrSges(6,3);
t537 = -t581 * mrSges(5,2) + t564 * mrSges(5,3);
t635 = -t636 - (t536 + t537) * t564;
t634 = t597 * pkin(6);
t626 = t564 * t536;
t541 = -t594 * g(3) - t591 * t555;
t567 = (mrSges(4,2) * t594 - mrSges(4,3) * t591) * qJD(1);
t568 = (-mrSges(3,1) * t594 + mrSges(3,2) * t591) * qJD(1);
t608 = t594 * t614;
t569 = t591 * qJDD(1) + t608;
t574 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t616;
t575 = -mrSges(4,1) * t616 - qJD(2) * mrSges(4,3);
t577 = pkin(3) * t615 - qJD(2) * pkin(7);
t589 = t594 ^ 2;
t578 = t592 * g(1) - t595 * g(2);
t605 = -qJDD(1) * pkin(1) - t578;
t599 = pkin(2) * t609 + t615 * t643 + (-t569 - t608) * qJ(3) + t605;
t510 = -t577 * t615 + (-pkin(3) * t589 - pkin(6)) * t597 + (-pkin(2) - pkin(7)) * t570 + t599;
t518 = -qJDD(2) * pkin(2) - t596 * qJ(3) + t566 * t615 + qJDD(3) - t541;
t515 = (-t591 * t594 * t597 - qJDD(2)) * pkin(7) + (t569 - t608) * pkin(3) + t518;
t505 = -t590 * t510 + t593 * t515;
t532 = t564 * qJD(4) + t593 * qJDD(2) - t590 * t570;
t534 = -t564 * mrSges(6,1) + t565 * mrSges(6,2);
t535 = -t564 * mrSges(5,1) + t565 * mrSges(5,2);
t563 = qJDD(4) + t569;
t502 = -0.2e1 * qJD(5) * t565 + (t564 * t581 - t532) * qJ(5) + (t564 * t565 + t563) * pkin(4) + t505;
t612 = m(6) * t502 + t563 * mrSges(6,1) + t581 * t536;
t494 = m(5) * t505 + t563 * mrSges(5,1) + t581 * t537 + (-t534 - t535) * t565 + (-mrSges(5,3) - mrSges(6,3)) * t532 + t612;
t506 = t593 * t510 + t590 * t515;
t539 = t581 * mrSges(6,1) - t565 * mrSges(6,3);
t540 = t581 * mrSges(5,1) - t565 * mrSges(5,3);
t538 = t581 * pkin(4) - t565 * qJ(5);
t562 = t564 ^ 2;
t504 = -t562 * pkin(4) + t531 * qJ(5) + 0.2e1 * qJD(5) * t564 - t581 * t538 + t506;
t611 = m(6) * t504 + t531 * mrSges(6,3) + t564 * t534;
t499 = m(5) * t506 + t531 * mrSges(5,3) + t564 * t535 + (-t539 - t540) * t581 + (-mrSges(5,2) - mrSges(6,2)) * t563 + t611;
t492 = t593 * t494 + t590 * t499;
t601 = -m(4) * t518 - t569 * mrSges(4,1) - t492;
t489 = m(3) * t541 - t569 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t574 - t575) * qJD(2) + (-t567 - t568) * t615 + t601;
t573 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t615;
t576 = mrSges(4,1) * t615 + qJD(2) * mrSges(4,2);
t514 = -t589 * t597 * pkin(7) + t570 * pkin(3) + qJD(2) * t577 - t517;
t508 = -t531 * pkin(4) - t562 * qJ(5) + t565 * t538 + qJDD(5) + t514;
t610 = m(6) * t508 + t532 * mrSges(6,2) + t565 * t539;
t604 = -m(5) * t514 - t532 * mrSges(5,2) - t565 * t540 - t610;
t600 = -m(4) * t517 + qJDD(2) * mrSges(4,3) + qJD(2) * t576 + t567 * t616 - t604;
t498 = t568 * t616 + t600 + (mrSges(3,3) + mrSges(4,1)) * t570 - qJDD(2) * mrSges(3,2) - qJD(2) * t573 + m(3) * t542 + t635;
t606 = -t591 * t489 + t594 * t498;
t483 = m(2) * t579 - t597 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t606;
t554 = t605 - t634;
t516 = -t570 * pkin(2) + t599 - t634;
t624 = -t590 * t494 + t593 * t499;
t603 = -m(4) * t516 - t570 * mrSges(4,2) + t576 * t615 - t624;
t598 = -m(3) * t554 + t574 * t616 + t570 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t569 + (-t573 * t591 - t575 * t594) * qJD(1) + t603;
t487 = m(2) * t578 + qJDD(1) * mrSges(2,1) - t597 * mrSges(2,2) + t598;
t625 = t592 * t483 + t595 * t487;
t484 = t594 * t489 + t591 * t498;
t623 = -t627 * t564 - t629 * t565 - t637 * t581;
t622 = t639 * t564 + t631 * t565 + t627 * t581;
t621 = -t631 * t564 - t641 * t565 - t629 * t581;
t619 = t638 * qJD(2) + (t630 * t591 + t628 * t594) * qJD(1);
t618 = -t628 * qJD(2) + (-t632 * t591 - t640 * t594) * qJD(1);
t617 = t630 * qJD(2) + (t642 * t591 + t632 * t594) * qJD(1);
t607 = t595 * t483 - t592 * t487;
t500 = -t532 * mrSges(6,3) - t565 * t534 + t612;
t491 = mrSges(5,2) * t514 + mrSges(6,2) * t508 - mrSges(5,3) * t505 - mrSges(6,3) * t502 - qJ(5) * t500 + t631 * t531 + t641 * t532 + t629 * t563 - t623 * t564 - t622 * t581;
t490 = -t569 * mrSges(4,3) + t575 * t616 - t603;
t485 = -mrSges(5,1) * t514 + mrSges(5,3) * t506 - mrSges(6,1) * t508 + mrSges(6,3) * t504 - pkin(4) * (t610 - t626) + qJ(5) * t611 + (-qJ(5) * t539 - t621) * t581 + t623 * t565 + (-qJ(5) * mrSges(6,2) + t627) * t563 + t631 * t532 + (pkin(4) * mrSges(6,1) + t639) * t531;
t480 = mrSges(4,1) * t518 + mrSges(5,1) * t505 + mrSges(6,1) * t502 + mrSges(3,2) * t554 - mrSges(5,2) * t506 - mrSges(6,2) * t504 - mrSges(3,3) * t541 - mrSges(4,3) * t516 + pkin(3) * t492 + pkin(4) * t500 - qJ(3) * t490 + t632 * t570 + t642 * t569 + t622 * t565 + t621 * t564 + t637 * t563 + t629 * t532 + t627 * t531 + t630 * qJDD(2) + t618 * qJD(2) + t619 * t616;
t479 = -mrSges(3,1) * t554 + mrSges(3,3) * t542 - mrSges(4,1) * t517 + mrSges(4,2) * t516 - t590 * t491 - t593 * t485 - pkin(3) * (t604 - t635) - pkin(7) * t624 - pkin(2) * t490 + t640 * t570 + t632 * t569 + t628 * qJDD(2) + t617 * qJD(2) - t619 * t615;
t478 = -pkin(1) * t484 + mrSges(2,3) * t579 - pkin(2) * (-qJD(2) * t575 + t601) - qJ(3) * (-t564 * t537 + t600 - t626 - t636) + t590 * t485 + pkin(7) * t492 - mrSges(3,1) * t541 + mrSges(3,2) * t542 - mrSges(4,2) * t518 + mrSges(4,3) * t517 - t593 * t491 + t597 * Ifges(2,5) + mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) + (-qJ(3) * mrSges(4,1) - t628) * t570 - t630 * t569 + (pkin(2) * mrSges(4,2) - t638) * qJDD(2) + (t617 * t594 + (pkin(2) * t567 + t618) * t591) * qJD(1);
t477 = -mrSges(2,2) * g(3) - mrSges(2,3) * t578 + Ifges(2,5) * qJDD(1) - t597 * Ifges(2,6) - pkin(6) * t484 - t591 * t479 + t594 * t480;
t1 = [-m(1) * g(1) + t607; -m(1) * g(2) + t625; (-m(1) - m(2)) * g(3) + t484; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t625 + t595 * t477 - t592 * t478; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t607 + t592 * t477 + t595 * t478; -mrSges(1,1) * g(2) + mrSges(2,1) * t578 + mrSges(1,2) * g(1) - mrSges(2,2) * t579 + Ifges(2,3) * qJDD(1) + pkin(1) * t598 + pkin(6) * t606 + t594 * t479 + t591 * t480;];
tauB = t1;

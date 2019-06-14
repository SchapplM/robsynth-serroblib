% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-05 00:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:59:35
% EndTime: 2019-05-04 23:59:40
% DurationCPUTime: 3.93s
% Computational Cost: add. (44100->268), mult. (80350->324), div. (0->0), fcn. (49937->10), ass. (0->116)
t635 = Ifges(6,1) + Ifges(7,1);
t625 = Ifges(6,4) - Ifges(7,5);
t634 = -Ifges(6,5) - Ifges(7,4);
t633 = Ifges(6,2) + Ifges(7,3);
t622 = Ifges(6,6) - Ifges(7,6);
t632 = -Ifges(6,3) - Ifges(7,2);
t584 = sin(pkin(10));
t586 = cos(pkin(10));
t569 = g(1) * t584 - g(2) * t586;
t570 = -g(1) * t586 - g(2) * t584;
t581 = -g(3) + qJDD(1);
t592 = cos(qJ(2));
t587 = cos(pkin(6));
t590 = sin(qJ(2));
t619 = t587 * t590;
t585 = sin(pkin(6));
t620 = t585 * t590;
t524 = t569 * t619 + t592 * t570 + t581 * t620;
t631 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t524;
t523 = -t590 * t570 + (t569 * t587 + t581 * t585) * t592;
t630 = -pkin(2) - pkin(8);
t629 = cos(qJ(5));
t628 = mrSges(3,1) - mrSges(4,2);
t627 = -mrSges(6,3) - mrSges(7,2);
t626 = (-Ifges(4,4) + Ifges(3,5));
t623 = Ifges(4,5) - Ifges(3,6);
t594 = qJD(2) ^ 2;
t597 = -t594 * qJ(3) + qJDD(3) - t523;
t520 = qJDD(2) * t630 + t597;
t549 = -t569 * t585 + t581 * t587;
t589 = sin(qJ(4));
t591 = cos(qJ(4));
t516 = t589 * t520 + t591 * t549;
t565 = (mrSges(5,1) * t589 + mrSges(5,2) * t591) * qJD(2);
t611 = qJD(2) * qJD(4);
t607 = t591 * t611;
t567 = -qJDD(2) * t589 - t607;
t612 = qJD(2) * t591;
t572 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t612;
t566 = (pkin(4) * t589 - pkin(9) * t591) * qJD(2);
t593 = qJD(4) ^ 2;
t613 = qJD(2) * t589;
t512 = -pkin(4) * t593 + qJDD(4) * pkin(9) - t566 * t613 + t516;
t519 = t594 * t630 - t631;
t608 = t589 * t611;
t568 = qJDD(2) * t591 - t608;
t514 = (-t568 + t608) * pkin(9) + (-t567 + t607) * pkin(4) + t519;
t588 = sin(qJ(5));
t509 = t629 * t512 + t588 * t514;
t564 = t588 * qJD(4) + t612 * t629;
t535 = t564 * qJD(5) - qJDD(4) * t629 + t588 * t568;
t575 = qJD(5) + t613;
t545 = mrSges(6,1) * t575 - mrSges(6,3) * t564;
t560 = qJDD(5) - t567;
t563 = -qJD(4) * t629 + t588 * t612;
t539 = pkin(5) * t563 - qJ(6) * t564;
t574 = t575 ^ 2;
t505 = -pkin(5) * t574 + qJ(6) * t560 + 0.2e1 * qJD(6) * t575 - t539 * t563 + t509;
t546 = -mrSges(7,1) * t575 + mrSges(7,2) * t564;
t609 = m(7) * t505 + t560 * mrSges(7,3) + t575 * t546;
t540 = mrSges(7,1) * t563 - mrSges(7,3) * t564;
t614 = -mrSges(6,1) * t563 - mrSges(6,2) * t564 - t540;
t500 = m(6) * t509 - t560 * mrSges(6,2) + t535 * t627 - t575 * t545 + t563 * t614 + t609;
t508 = -t588 * t512 + t514 * t629;
t536 = -t563 * qJD(5) + t588 * qJDD(4) + t568 * t629;
t544 = -mrSges(6,2) * t575 - mrSges(6,3) * t563;
t506 = -t560 * pkin(5) - t574 * qJ(6) + t564 * t539 + qJDD(6) - t508;
t547 = -mrSges(7,2) * t563 + mrSges(7,3) * t575;
t602 = -m(7) * t506 + t560 * mrSges(7,1) + t575 * t547;
t502 = m(6) * t508 + t560 * mrSges(6,1) + t536 * t627 + t575 * t544 + t564 * t614 + t602;
t604 = t629 * t500 - t502 * t588;
t493 = m(5) * t516 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t567 - qJD(4) * t572 - t565 * t613 + t604;
t515 = t591 * t520 - t589 * t549;
t571 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t613;
t511 = -qJDD(4) * pkin(4) - t593 * pkin(9) + t566 * t612 - t515;
t507 = -0.2e1 * qJD(6) * t564 + (t563 * t575 - t536) * qJ(6) + (t564 * t575 + t535) * pkin(5) + t511;
t503 = m(7) * t507 + mrSges(7,1) * t535 - t536 * mrSges(7,3) - t564 * t546 + t547 * t563;
t595 = -m(6) * t511 - t535 * mrSges(6,1) - mrSges(6,2) * t536 - t563 * t544 - t545 * t564 - t503;
t497 = m(5) * t515 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t568 + qJD(4) * t571 - t565 * t612 + t595;
t487 = t589 * t493 + t591 * t497;
t522 = -qJDD(2) * pkin(2) + t597;
t599 = -m(4) * t522 + (t594 * mrSges(4,3)) - t487;
t483 = m(3) * t523 - (t594 * mrSges(3,2)) + qJDD(2) * t628 + t599;
t621 = t483 * t592;
t605 = t591 * t493 - t497 * t589;
t486 = m(4) * t549 + t605;
t485 = m(3) * t549 + t486;
t521 = t594 * pkin(2) + t631;
t496 = t588 * t500 + t629 * t502;
t598 = -m(5) * t519 + mrSges(5,1) * t567 - t568 * mrSges(5,2) - t571 * t613 - t572 * t612 - t496;
t596 = -m(4) * t521 + (t594 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t598;
t491 = m(3) * t524 - (mrSges(3,1) * t594) - qJDD(2) * mrSges(3,2) + t596;
t474 = -t485 * t585 + t491 * t619 + t587 * t621;
t472 = m(2) * t569 + t474;
t479 = -t483 * t590 + t592 * t491;
t478 = m(2) * t570 + t479;
t618 = t586 * t472 + t584 * t478;
t617 = t563 * t633 - t564 * t625 - t575 * t622;
t616 = t563 * t622 + t564 * t634 + t575 * t632;
t615 = -t625 * t563 + t564 * t635 - t634 * t575;
t473 = t587 * t485 + t491 * t620 + t585 * t621;
t606 = -t472 * t584 + t586 * t478;
t494 = -mrSges(6,1) * t511 - mrSges(7,1) * t507 + mrSges(7,2) * t505 + mrSges(6,3) * t509 - pkin(5) * t503 - t535 * t633 + t625 * t536 + t622 * t560 + t616 * t564 + t615 * t575;
t495 = mrSges(6,2) * t511 + mrSges(7,2) * t506 - mrSges(6,3) * t508 - mrSges(7,3) * t507 - qJ(6) * t503 - t625 * t535 + t536 * t635 - t560 * t634 + t616 * t563 + t617 * t575;
t552 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t591 - Ifges(5,6) * t589) * qJD(2);
t553 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t591 - Ifges(5,2) * t589) * qJD(2);
t475 = mrSges(5,2) * t519 - mrSges(5,3) * t515 + Ifges(5,1) * t568 + Ifges(5,4) * t567 + Ifges(5,5) * qJDD(4) - pkin(9) * t496 - qJD(4) * t553 - t588 * t494 + t495 * t629 - t552 * t613;
t554 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t591 - Ifges(5,4) * t589) * qJD(2);
t480 = Ifges(5,4) * t568 + Ifges(5,2) * t567 + Ifges(5,6) * qJDD(4) - t552 * t612 + qJD(4) * t554 - mrSges(5,1) * t519 + mrSges(5,3) * t516 - mrSges(6,1) * t508 + mrSges(6,2) * t509 + mrSges(7,1) * t506 - mrSges(7,3) * t505 - pkin(5) * t602 - qJ(6) * t609 - pkin(4) * t496 + (pkin(5) * t540 + t617) * t564 + (qJ(6) * t540 - t615) * t563 + t632 * t560 + (mrSges(7,2) * pkin(5) + t634) * t536 + (mrSges(7,2) * qJ(6) + t622) * t535;
t469 = -mrSges(4,1) * t521 + mrSges(3,3) * t524 - pkin(2) * t486 - pkin(3) * t598 - pkin(8) * t605 - qJDD(2) * t623 - t589 * t475 - t591 * t480 - t549 * t628 + (t594 * t626);
t470 = pkin(3) * t487 + mrSges(4,1) * t522 + t629 * t494 + pkin(4) * t595 + pkin(9) * t604 + mrSges(5,1) * t515 - mrSges(5,2) * t516 + t588 * t495 + Ifges(5,5) * t568 + Ifges(5,6) * t567 + Ifges(5,3) * qJDD(4) - qJ(3) * t486 - mrSges(3,3) * t523 + t623 * t594 + (mrSges(3,2) - mrSges(4,3)) * t549 + t626 * qJDD(2) + (t553 * t591 + t554 * t589) * qJD(2);
t600 = pkin(7) * t479 + t469 * t592 + t470 * t590;
t468 = mrSges(3,1) * t523 - mrSges(3,2) * t524 + mrSges(4,2) * t522 - mrSges(4,3) * t521 + t591 * t475 - t589 * t480 - pkin(8) * t487 + pkin(2) * t599 + qJ(3) * t596 + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t467 = mrSges(2,2) * t581 - mrSges(2,3) * t569 - t590 * t469 + t592 * t470 + (-t473 * t585 - t474 * t587) * pkin(7);
t466 = -mrSges(2,1) * t581 + mrSges(2,3) * t570 - pkin(1) * t473 - t585 * t468 + t587 * t600;
t1 = [-m(1) * g(1) + t606; -m(1) * g(2) + t618; -m(1) * g(3) + m(2) * t581 + t473; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t618 - t584 * t466 + t586 * t467; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t606 + t586 * t466 + t584 * t467; -mrSges(1,1) * g(2) + mrSges(2,1) * t569 + mrSges(1,2) * g(1) - mrSges(2,2) * t570 + pkin(1) * t474 + t587 * t468 + t585 * t600;];
tauB  = t1;

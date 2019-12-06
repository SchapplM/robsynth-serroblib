% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:45:13
% EndTime: 2019-12-05 18:45:19
% DurationCPUTime: 5.22s
% Computational Cost: add. (55397->292), mult. (121964->357), div. (0->0), fcn. (84274->8), ass. (0->113)
t628 = Ifges(5,1) + Ifges(6,1);
t624 = Ifges(5,4) + Ifges(6,4);
t623 = Ifges(5,5) + Ifges(6,5);
t627 = Ifges(5,2) + Ifges(6,2);
t622 = -Ifges(5,6) - Ifges(6,6);
t626 = -Ifges(5,3) - Ifges(6,3);
t602 = qJD(1) ^ 2;
t625 = pkin(2) * t602;
t597 = sin(qJ(1));
t601 = cos(qJ(1));
t586 = -t601 * g(1) - t597 * g(2);
t575 = -t602 * pkin(1) + qJDD(1) * pkin(6) + t586;
t596 = sin(qJ(2));
t621 = t596 * t575;
t600 = cos(qJ(2));
t614 = qJD(1) * qJD(2);
t580 = t596 * qJDD(1) + t600 * t614;
t541 = qJDD(2) * pkin(2) - t580 * pkin(7) - t621 + (pkin(7) * t614 + t596 * t625 - g(3)) * t600;
t563 = -t596 * g(3) + t600 * t575;
t581 = t600 * qJDD(1) - t596 * t614;
t616 = qJD(1) * t596;
t584 = qJD(2) * pkin(2) - pkin(7) * t616;
t593 = t600 ^ 2;
t542 = t581 * pkin(7) - qJD(2) * t584 - t593 * t625 + t563;
t595 = sin(qJ(3));
t599 = cos(qJ(3));
t524 = t599 * t541 - t595 * t542;
t572 = (-t595 * t596 + t599 * t600) * qJD(1);
t546 = t572 * qJD(3) + t599 * t580 + t595 * t581;
t573 = (t595 * t600 + t596 * t599) * qJD(1);
t591 = qJDD(2) + qJDD(3);
t592 = qJD(2) + qJD(3);
t511 = (t572 * t592 - t546) * pkin(8) + (t572 * t573 + t591) * pkin(3) + t524;
t525 = t595 * t541 + t599 * t542;
t545 = -t573 * qJD(3) - t595 * t580 + t599 * t581;
t566 = t592 * pkin(3) - t573 * pkin(8);
t568 = t572 ^ 2;
t513 = -t568 * pkin(3) + t545 * pkin(8) - t592 * t566 + t525;
t594 = sin(qJ(4));
t598 = cos(qJ(4));
t505 = t598 * t511 - t594 * t513;
t559 = t598 * t572 - t594 * t573;
t522 = t559 * qJD(4) + t594 * t545 + t598 * t546;
t560 = t594 * t572 + t598 * t573;
t536 = -t559 * mrSges(6,1) + t560 * mrSges(6,2);
t537 = -t559 * mrSges(5,1) + t560 * mrSges(5,2);
t589 = qJD(4) + t592;
t549 = -t589 * mrSges(5,2) + t559 * mrSges(5,3);
t588 = qJDD(4) + t591;
t502 = -0.2e1 * qJD(5) * t560 + (t559 * t589 - t522) * qJ(5) + (t559 * t560 + t588) * pkin(4) + t505;
t548 = -t589 * mrSges(6,2) + t559 * mrSges(6,3);
t613 = m(6) * t502 + t588 * mrSges(6,1) + t589 * t548;
t494 = m(5) * t505 + t588 * mrSges(5,1) + t589 * t549 + (-t536 - t537) * t560 + (-mrSges(5,3) - mrSges(6,3)) * t522 + t613;
t506 = t594 * t511 + t598 * t513;
t521 = -t560 * qJD(4) + t598 * t545 - t594 * t546;
t551 = t589 * mrSges(6,1) - t560 * mrSges(6,3);
t552 = t589 * mrSges(5,1) - t560 * mrSges(5,3);
t550 = t589 * pkin(4) - t560 * qJ(5);
t558 = t559 ^ 2;
t504 = -t558 * pkin(4) + t521 * qJ(5) + 0.2e1 * qJD(5) * t559 - t589 * t550 + t506;
t612 = m(6) * t504 + t521 * mrSges(6,3) + t559 * t536;
t499 = m(5) * t506 + t521 * mrSges(5,3) + t559 * t537 + (-t551 - t552) * t589 + (-mrSges(5,2) - mrSges(6,2)) * t588 + t612;
t492 = t598 * t494 + t594 * t499;
t561 = -t572 * mrSges(4,1) + t573 * mrSges(4,2);
t564 = -t592 * mrSges(4,2) + t572 * mrSges(4,3);
t489 = m(4) * t524 + t591 * mrSges(4,1) - t546 * mrSges(4,3) - t573 * t561 + t592 * t564 + t492;
t565 = t592 * mrSges(4,1) - t573 * mrSges(4,3);
t608 = -t594 * t494 + t598 * t499;
t490 = m(4) * t525 - t591 * mrSges(4,2) + t545 * mrSges(4,3) + t572 * t561 - t592 * t565 + t608;
t484 = t599 * t489 + t595 * t490;
t562 = -t600 * g(3) - t621;
t579 = (-mrSges(3,1) * t600 + mrSges(3,2) * t596) * qJD(1);
t615 = qJD(1) * t600;
t583 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t615;
t482 = m(3) * t562 + qJDD(2) * mrSges(3,1) - t580 * mrSges(3,3) + qJD(2) * t583 - t579 * t616 + t484;
t582 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t616;
t609 = -t595 * t489 + t599 * t490;
t483 = m(3) * t563 - qJDD(2) * mrSges(3,2) + t581 * mrSges(3,3) - qJD(2) * t582 + t579 * t615 + t609;
t610 = -t596 * t482 + t600 * t483;
t475 = m(2) * t586 - t602 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t610;
t585 = t597 * g(1) - t601 * g(2);
t606 = -qJDD(1) * pkin(1) - t585;
t574 = -t602 * pkin(6) + t606;
t547 = -t581 * pkin(2) + t584 * t616 + (-pkin(7) * t593 - pkin(6)) * t602 + t606;
t515 = -t545 * pkin(3) - t568 * pkin(8) + t573 * t566 + t547;
t508 = -t521 * pkin(4) - t558 * qJ(5) + t560 * t550 + qJDD(5) + t515;
t607 = m(6) * t508 - t521 * mrSges(6,1) + t522 * mrSges(6,2) - t559 * t548 + t560 * t551;
t605 = m(5) * t515 - t521 * mrSges(5,1) + t522 * mrSges(5,2) - t559 * t549 + t560 * t552 + t607;
t604 = m(4) * t547 - t545 * mrSges(4,1) + t546 * mrSges(4,2) - t572 * t564 + t573 * t565 + t605;
t603 = -m(3) * t574 + t581 * mrSges(3,1) - t580 * mrSges(3,2) - t582 * t616 + t583 * t615 - t604;
t496 = m(2) * t585 + qJDD(1) * mrSges(2,1) - t602 * mrSges(2,2) + t603;
t620 = t597 * t475 + t601 * t496;
t476 = t600 * t482 + t596 * t483;
t619 = t622 * t559 - t623 * t560 + t626 * t589;
t618 = -t627 * t559 - t624 * t560 + t622 * t589;
t617 = t624 * t559 + t628 * t560 + t623 * t589;
t611 = t601 * t475 - t597 * t496;
t571 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t596 + Ifges(3,4) * t600) * qJD(1);
t570 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t596 + Ifges(3,2) * t600) * qJD(1);
t569 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t596 + Ifges(3,6) * t600) * qJD(1);
t555 = Ifges(4,1) * t573 + Ifges(4,4) * t572 + Ifges(4,5) * t592;
t554 = Ifges(4,4) * t573 + Ifges(4,2) * t572 + Ifges(4,6) * t592;
t553 = Ifges(4,5) * t573 + Ifges(4,6) * t572 + Ifges(4,3) * t592;
t500 = -t522 * mrSges(6,3) - t560 * t536 + t613;
t491 = mrSges(5,2) * t515 + mrSges(6,2) * t508 - mrSges(5,3) * t505 - mrSges(6,3) * t502 - qJ(5) * t500 + t624 * t521 + t628 * t522 - t619 * t559 + t623 * t588 + t618 * t589;
t485 = -mrSges(5,1) * t515 + mrSges(5,3) * t506 - mrSges(6,1) * t508 + mrSges(6,3) * t504 - pkin(4) * t607 + qJ(5) * t612 + (-qJ(5) * t551 + t617) * t589 + (-qJ(5) * mrSges(6,2) - t622) * t588 + t619 * t560 + t624 * t522 + t627 * t521;
t478 = mrSges(4,2) * t547 - mrSges(4,3) * t524 + Ifges(4,1) * t546 + Ifges(4,4) * t545 + Ifges(4,5) * t591 - pkin(8) * t492 - t594 * t485 + t598 * t491 + t572 * t553 - t592 * t554;
t477 = -mrSges(4,1) * t547 + mrSges(4,3) * t525 + Ifges(4,4) * t546 + Ifges(4,2) * t545 + Ifges(4,6) * t591 - pkin(3) * t605 + pkin(8) * t608 + t598 * t485 + t594 * t491 - t573 * t553 + t592 * t555;
t472 = t626 * t588 + Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) + (-t596 * t570 + t600 * t571) * qJD(1) + mrSges(2,1) * g(3) + t602 * Ifges(2,5) - Ifges(4,3) * t591 + mrSges(2,3) * t586 - Ifges(3,5) * t580 - Ifges(3,6) * t581 + mrSges(3,2) * t563 + t572 * t555 - t573 * t554 - mrSges(3,1) * t562 - Ifges(4,6) * t545 - Ifges(4,5) * t546 - mrSges(4,1) * t524 + mrSges(4,2) * t525 - mrSges(5,1) * t505 + mrSges(5,2) * t506 - mrSges(6,1) * t502 + mrSges(6,2) * t504 - pkin(4) * t500 - pkin(3) * t492 - pkin(2) * t484 - pkin(1) * t476 + t617 * t559 + t618 * t560 + t622 * t521 - t623 * t522;
t471 = mrSges(3,2) * t574 - mrSges(3,3) * t562 + Ifges(3,1) * t580 + Ifges(3,4) * t581 + Ifges(3,5) * qJDD(2) - pkin(7) * t484 - qJD(2) * t570 - t595 * t477 + t599 * t478 + t569 * t615;
t470 = -mrSges(3,1) * t574 + mrSges(3,3) * t563 + Ifges(3,4) * t580 + Ifges(3,2) * t581 + Ifges(3,6) * qJDD(2) - pkin(2) * t604 + pkin(7) * t609 + qJD(2) * t571 + t599 * t477 + t595 * t478 - t569 * t616;
t469 = -mrSges(2,2) * g(3) - mrSges(2,3) * t585 + Ifges(2,5) * qJDD(1) - t602 * Ifges(2,6) - pkin(6) * t476 - t596 * t470 + t600 * t471;
t1 = [-m(1) * g(1) + t611; -m(1) * g(2) + t620; (-m(1) - m(2)) * g(3) + t476; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t620 + t601 * t469 - t597 * t472; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t611 + t597 * t469 + t601 * t472; -mrSges(1,1) * g(2) + mrSges(2,1) * t585 + mrSges(1,2) * g(1) - mrSges(2,2) * t586 + Ifges(2,3) * qJDD(1) + pkin(1) * t603 + pkin(6) * t610 + t600 * t470 + t596 * t471;];
tauB = t1;

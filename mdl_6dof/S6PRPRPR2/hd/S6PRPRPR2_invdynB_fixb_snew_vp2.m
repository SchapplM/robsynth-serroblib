% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:20:46
% EndTime: 2019-05-04 22:20:58
% DurationCPUTime: 11.29s
% Computational Cost: add. (197983->297), mult. (385233->383), div. (0->0), fcn. (271955->14), ass. (0->128)
t585 = sin(pkin(6));
t592 = sin(qJ(2));
t616 = t585 * t592;
t595 = cos(qJ(2));
t615 = t585 * t595;
t589 = cos(pkin(6));
t614 = t589 * t592;
t613 = t589 * t595;
t584 = sin(pkin(10));
t588 = cos(pkin(10));
t572 = g(1) * t584 - g(2) * t588;
t573 = -g(1) * t588 - g(2) * t584;
t581 = -g(3) + qJDD(1);
t534 = t572 * t613 - t573 * t592 + t581 * t615;
t532 = qJDD(2) * pkin(2) + t534;
t535 = t572 * t614 + t595 * t573 + t581 * t616;
t597 = qJD(2) ^ 2;
t533 = -pkin(2) * t597 + t535;
t583 = sin(pkin(11));
t587 = cos(pkin(11));
t518 = t583 * t532 + t587 * t533;
t516 = -pkin(3) * t597 + qJDD(2) * pkin(8) + t518;
t552 = -t572 * t585 + t589 * t581;
t546 = qJDD(3) + t552;
t591 = sin(qJ(4));
t594 = cos(qJ(4));
t512 = t594 * t516 + t591 * t546;
t569 = (-mrSges(5,1) * t594 + mrSges(5,2) * t591) * qJD(2);
t609 = qJD(2) * qJD(4);
t580 = t591 * t609;
t571 = qJDD(2) * t594 - t580;
t611 = qJD(2) * t591;
t574 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t611;
t568 = (-pkin(4) * t594 - qJ(5) * t591) * qJD(2);
t596 = qJD(4) ^ 2;
t610 = qJD(2) * t594;
t507 = -pkin(4) * t596 + qJDD(4) * qJ(5) + t568 * t610 + t512;
t517 = t587 * t532 - t583 * t533;
t515 = -qJDD(2) * pkin(3) - t597 * pkin(8) - t517;
t607 = t594 * t609;
t570 = qJDD(2) * t591 + t607;
t510 = (-t570 - t607) * qJ(5) + (-t571 + t580) * pkin(4) + t515;
t582 = sin(pkin(12));
t586 = cos(pkin(12));
t564 = qJD(4) * t582 + t586 * t611;
t502 = -0.2e1 * qJD(5) * t564 - t582 * t507 + t586 * t510;
t550 = qJDD(4) * t582 + t570 * t586;
t563 = qJD(4) * t586 - t582 * t611;
t500 = (-t563 * t610 - t550) * pkin(9) + (t563 * t564 - t571) * pkin(5) + t502;
t503 = 0.2e1 * qJD(5) * t563 + t586 * t507 + t582 * t510;
t549 = qJDD(4) * t586 - t570 * t582;
t551 = -pkin(5) * t610 - pkin(9) * t564;
t562 = t563 ^ 2;
t501 = -pkin(5) * t562 + pkin(9) * t549 + t551 * t610 + t503;
t590 = sin(qJ(6));
t593 = cos(qJ(6));
t498 = t500 * t593 - t501 * t590;
t540 = t563 * t593 - t564 * t590;
t521 = qJD(6) * t540 + t549 * t590 + t550 * t593;
t541 = t563 * t590 + t564 * t593;
t526 = -mrSges(7,1) * t540 + mrSges(7,2) * t541;
t578 = qJD(6) - t610;
t530 = -mrSges(7,2) * t578 + mrSges(7,3) * t540;
t566 = qJDD(6) - t571;
t494 = m(7) * t498 + mrSges(7,1) * t566 - mrSges(7,3) * t521 - t526 * t541 + t530 * t578;
t499 = t500 * t590 + t501 * t593;
t520 = -qJD(6) * t541 + t549 * t593 - t550 * t590;
t531 = mrSges(7,1) * t578 - mrSges(7,3) * t541;
t495 = m(7) * t499 - mrSges(7,2) * t566 + mrSges(7,3) * t520 + t526 * t540 - t531 * t578;
t488 = t593 * t494 + t590 * t495;
t542 = -mrSges(6,1) * t563 + mrSges(6,2) * t564;
t547 = mrSges(6,2) * t610 + mrSges(6,3) * t563;
t486 = m(6) * t502 - mrSges(6,1) * t571 - mrSges(6,3) * t550 - t542 * t564 - t547 * t610 + t488;
t548 = -mrSges(6,1) * t610 - mrSges(6,3) * t564;
t602 = -t494 * t590 + t593 * t495;
t487 = m(6) * t503 + mrSges(6,2) * t571 + mrSges(6,3) * t549 + t542 * t563 + t548 * t610 + t602;
t603 = -t486 * t582 + t586 * t487;
t483 = m(5) * t512 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t571 - qJD(4) * t574 + t569 * t610 + t603;
t511 = -t591 * t516 + t546 * t594;
t575 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t610;
t506 = -qJDD(4) * pkin(4) - qJ(5) * t596 + t568 * t611 + qJDD(5) - t511;
t504 = -pkin(5) * t549 - pkin(9) * t562 + t551 * t564 + t506;
t600 = m(7) * t504 - t520 * mrSges(7,1) + mrSges(7,2) * t521 - t540 * t530 + t531 * t541;
t598 = -m(6) * t506 + t549 * mrSges(6,1) - mrSges(6,2) * t550 + t563 * t547 - t548 * t564 - t600;
t497 = m(5) * t511 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t570 + qJD(4) * t575 - t569 * t611 + t598;
t604 = t594 * t483 - t497 * t591;
t473 = m(4) * t518 - mrSges(4,1) * t597 - qJDD(2) * mrSges(4,2) + t604;
t484 = t486 * t586 + t487 * t582;
t599 = -m(5) * t515 + t571 * mrSges(5,1) - mrSges(5,2) * t570 - t574 * t611 + t575 * t610 - t484;
t480 = m(4) * t517 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t597 + t599;
t469 = t583 * t473 + t587 * t480;
t467 = m(3) * t534 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t597 + t469;
t605 = t587 * t473 - t480 * t583;
t468 = m(3) * t535 - mrSges(3,1) * t597 - qJDD(2) * mrSges(3,2) + t605;
t476 = t591 * t483 + t594 * t497;
t608 = m(4) * t546 + t476;
t475 = m(3) * t552 + t608;
t455 = t467 * t613 + t468 * t614 - t475 * t585;
t453 = m(2) * t572 + t455;
t460 = -t467 * t592 + t595 * t468;
t459 = m(2) * t573 + t460;
t612 = t588 * t453 + t584 * t459;
t454 = t467 * t615 + t468 * t616 + t589 * t475;
t606 = -t453 * t584 + t588 * t459;
t522 = Ifges(7,5) * t541 + Ifges(7,6) * t540 + Ifges(7,3) * t578;
t524 = Ifges(7,1) * t541 + Ifges(7,4) * t540 + Ifges(7,5) * t578;
t489 = -mrSges(7,1) * t504 + mrSges(7,3) * t499 + Ifges(7,4) * t521 + Ifges(7,2) * t520 + Ifges(7,6) * t566 - t522 * t541 + t524 * t578;
t523 = Ifges(7,4) * t541 + Ifges(7,2) * t540 + Ifges(7,6) * t578;
t490 = mrSges(7,2) * t504 - mrSges(7,3) * t498 + Ifges(7,1) * t521 + Ifges(7,4) * t520 + Ifges(7,5) * t566 + t522 * t540 - t523 * t578;
t536 = Ifges(6,5) * t564 + Ifges(6,6) * t563 - Ifges(6,3) * t610;
t538 = Ifges(6,1) * t564 + Ifges(6,4) * t563 - Ifges(6,5) * t610;
t477 = -mrSges(6,1) * t506 + mrSges(6,3) * t503 + Ifges(6,4) * t550 + Ifges(6,2) * t549 - Ifges(6,6) * t571 - pkin(5) * t600 + pkin(9) * t602 + t593 * t489 + t590 * t490 - t564 * t536 - t538 * t610;
t537 = Ifges(6,4) * t564 + Ifges(6,2) * t563 - Ifges(6,6) * t610;
t478 = mrSges(6,2) * t506 - mrSges(6,3) * t502 + Ifges(6,1) * t550 + Ifges(6,4) * t549 - Ifges(6,5) * t571 - pkin(9) * t488 - t489 * t590 + t490 * t593 + t536 * t563 + t537 * t610;
t558 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t591 + Ifges(5,6) * t594) * qJD(2);
t559 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t591 + Ifges(5,2) * t594) * qJD(2);
t461 = mrSges(5,2) * t515 - mrSges(5,3) * t511 + Ifges(5,1) * t570 + Ifges(5,4) * t571 + Ifges(5,5) * qJDD(4) - qJ(5) * t484 - qJD(4) * t559 - t477 * t582 + t478 * t586 + t558 * t610;
t560 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t591 + Ifges(5,4) * t594) * qJD(2);
t470 = Ifges(5,4) * t570 + Ifges(5,6) * qJDD(4) - t558 * t611 + qJD(4) * t560 - mrSges(5,1) * t515 + mrSges(5,3) * t512 - Ifges(6,5) * t550 - Ifges(6,6) * t549 - t564 * t537 + t563 * t538 - mrSges(6,1) * t502 + mrSges(6,2) * t503 - Ifges(7,5) * t521 - Ifges(7,6) * t520 - Ifges(7,3) * t566 - t541 * t523 + t540 * t524 - mrSges(7,1) * t498 + mrSges(7,2) * t499 - pkin(5) * t488 - pkin(4) * t484 + (Ifges(5,2) + Ifges(6,3)) * t571;
t451 = mrSges(4,2) * t546 - mrSges(4,3) * t517 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t597 - pkin(8) * t476 + t461 * t594 - t470 * t591;
t456 = Ifges(4,6) * qJDD(2) + t597 * Ifges(4,5) - mrSges(4,1) * t546 + mrSges(4,3) * t518 - Ifges(5,5) * t570 - Ifges(5,6) * t571 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t511 + mrSges(5,2) * t512 - t582 * t478 - t586 * t477 - pkin(4) * t598 - qJ(5) * t603 - pkin(3) * t476 + (-t559 * t591 + t560 * t594) * qJD(2);
t448 = -mrSges(3,1) * t552 + mrSges(3,3) * t535 + t597 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t608 + qJ(3) * t605 + t583 * t451 + t587 * t456;
t449 = mrSges(3,2) * t552 - mrSges(3,3) * t534 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t597 - qJ(3) * t469 + t451 * t587 - t456 * t583;
t601 = pkin(7) * t460 + t448 * t595 + t449 * t592;
t450 = mrSges(3,1) * t534 - mrSges(3,2) * t535 + mrSges(4,1) * t517 - mrSges(4,2) * t518 + t591 * t461 + t594 * t470 + pkin(3) * t599 + pkin(8) * t604 + pkin(2) * t469 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t447 = mrSges(2,2) * t581 - mrSges(2,3) * t572 - t592 * t448 + t595 * t449 + (-t454 * t585 - t455 * t589) * pkin(7);
t446 = -mrSges(2,1) * t581 + mrSges(2,3) * t573 - pkin(1) * t454 - t585 * t450 + t589 * t601;
t1 = [-m(1) * g(1) + t606; -m(1) * g(2) + t612; -m(1) * g(3) + m(2) * t581 + t454; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t612 - t584 * t446 + t588 * t447; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t606 + t588 * t446 + t584 * t447; -mrSges(1,1) * g(2) + mrSges(2,1) * t572 + mrSges(1,2) * g(1) - mrSges(2,2) * t573 + pkin(1) * t455 + t589 * t450 + t585 * t601;];
tauB  = t1;

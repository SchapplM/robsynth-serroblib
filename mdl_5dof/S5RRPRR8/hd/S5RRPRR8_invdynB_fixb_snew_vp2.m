% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:17:22
% EndTime: 2019-12-31 20:17:28
% DurationCPUTime: 6.26s
% Computational Cost: add. (88201->313), mult. (204486->400), div. (0->0), fcn. (145132->10), ass. (0->121)
t575 = qJD(1) ^ 2;
t591 = pkin(2) * t575;
t570 = sin(qJ(1));
t574 = cos(qJ(1));
t559 = -t574 * g(1) - t570 * g(2);
t548 = -t575 * pkin(1) + qJDD(1) * pkin(6) + t559;
t569 = sin(qJ(2));
t590 = t569 * t548;
t573 = cos(qJ(2));
t586 = qJD(1) * qJD(2);
t553 = t569 * qJDD(1) + t573 * t586;
t514 = qJDD(2) * pkin(2) - t553 * qJ(3) - t590 + (qJ(3) * t586 + t569 * t591 - g(3)) * t573;
t534 = -t569 * g(3) + t573 * t548;
t554 = t573 * qJDD(1) - t569 * t586;
t588 = qJD(1) * t569;
t555 = qJD(2) * pkin(2) - qJ(3) * t588;
t564 = t573 ^ 2;
t515 = t554 * qJ(3) - qJD(2) * t555 - t564 * t591 + t534;
t565 = sin(pkin(9));
t566 = cos(pkin(9));
t543 = (t565 * t573 + t566 * t569) * qJD(1);
t494 = -0.2e1 * qJD(3) * t543 + t566 * t514 - t565 * t515;
t532 = t566 * t553 + t565 * t554;
t542 = (-t565 * t569 + t566 * t573) * qJD(1);
t482 = (qJD(2) * t542 - t532) * pkin(7) + (t542 * t543 + qJDD(2)) * pkin(3) + t494;
t495 = 0.2e1 * qJD(3) * t542 + t565 * t514 + t566 * t515;
t531 = -t565 * t553 + t566 * t554;
t537 = qJD(2) * pkin(3) - t543 * pkin(7);
t541 = t542 ^ 2;
t484 = -t541 * pkin(3) + t531 * pkin(7) - qJD(2) * t537 + t495;
t568 = sin(qJ(4));
t572 = cos(qJ(4));
t480 = t568 * t482 + t572 * t484;
t526 = t568 * t542 + t572 * t543;
t499 = -t526 * qJD(4) + t572 * t531 - t568 * t532;
t525 = t572 * t542 - t568 * t543;
t509 = -t525 * mrSges(5,1) + t526 * mrSges(5,2);
t563 = qJD(2) + qJD(4);
t520 = t563 * mrSges(5,1) - t526 * mrSges(5,3);
t562 = qJDD(2) + qJDD(4);
t510 = -t525 * pkin(4) - t526 * pkin(8);
t561 = t563 ^ 2;
t477 = -t561 * pkin(4) + t562 * pkin(8) + t525 * t510 + t480;
t558 = t570 * g(1) - t574 * g(2);
t580 = -qJDD(1) * pkin(1) - t558;
t516 = -t554 * pkin(2) + qJDD(3) + t555 * t588 + (-qJ(3) * t564 - pkin(6)) * t575 + t580;
t493 = -t531 * pkin(3) - t541 * pkin(7) + t543 * t537 + t516;
t500 = t525 * qJD(4) + t568 * t531 + t572 * t532;
t478 = (-t525 * t563 - t500) * pkin(8) + (t526 * t563 - t499) * pkin(4) + t493;
t567 = sin(qJ(5));
t571 = cos(qJ(5));
t474 = -t567 * t477 + t571 * t478;
t517 = -t567 * t526 + t571 * t563;
t487 = t517 * qJD(5) + t571 * t500 + t567 * t562;
t498 = qJDD(5) - t499;
t518 = t571 * t526 + t567 * t563;
t501 = -t517 * mrSges(6,1) + t518 * mrSges(6,2);
t521 = qJD(5) - t525;
t502 = -t521 * mrSges(6,2) + t517 * mrSges(6,3);
t472 = m(6) * t474 + t498 * mrSges(6,1) - t487 * mrSges(6,3) - t518 * t501 + t521 * t502;
t475 = t571 * t477 + t567 * t478;
t486 = -t518 * qJD(5) - t567 * t500 + t571 * t562;
t503 = t521 * mrSges(6,1) - t518 * mrSges(6,3);
t473 = m(6) * t475 - t498 * mrSges(6,2) + t486 * mrSges(6,3) + t517 * t501 - t521 * t503;
t581 = -t567 * t472 + t571 * t473;
t463 = m(5) * t480 - t562 * mrSges(5,2) + t499 * mrSges(5,3) + t525 * t509 - t563 * t520 + t581;
t479 = t572 * t482 - t568 * t484;
t519 = -t563 * mrSges(5,2) + t525 * mrSges(5,3);
t476 = -t562 * pkin(4) - t561 * pkin(8) + t526 * t510 - t479;
t578 = -m(6) * t476 + t486 * mrSges(6,1) - t487 * mrSges(6,2) + t517 * t502 - t518 * t503;
t468 = m(5) * t479 + t562 * mrSges(5,1) - t500 * mrSges(5,3) - t526 * t509 + t563 * t519 + t578;
t458 = t568 * t463 + t572 * t468;
t529 = -t542 * mrSges(4,1) + t543 * mrSges(4,2);
t535 = -qJD(2) * mrSges(4,2) + t542 * mrSges(4,3);
t456 = m(4) * t494 + qJDD(2) * mrSges(4,1) - t532 * mrSges(4,3) + qJD(2) * t535 - t543 * t529 + t458;
t536 = qJD(2) * mrSges(4,1) - t543 * mrSges(4,3);
t582 = t572 * t463 - t568 * t468;
t457 = m(4) * t495 - qJDD(2) * mrSges(4,2) + t531 * mrSges(4,3) - qJD(2) * t536 + t542 * t529 + t582;
t450 = t566 * t456 + t565 * t457;
t533 = -t573 * g(3) - t590;
t552 = (-mrSges(3,1) * t573 + mrSges(3,2) * t569) * qJD(1);
t587 = qJD(1) * t573;
t557 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t587;
t448 = m(3) * t533 + qJDD(2) * mrSges(3,1) - t553 * mrSges(3,3) + qJD(2) * t557 - t552 * t588 + t450;
t556 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t588;
t583 = -t565 * t456 + t566 * t457;
t449 = m(3) * t534 - qJDD(2) * mrSges(3,2) + t554 * mrSges(3,3) - qJD(2) * t556 + t552 * t587 + t583;
t584 = -t569 * t448 + t573 * t449;
t442 = m(2) * t559 - t575 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t584;
t547 = -t575 * pkin(6) + t580;
t464 = t571 * t472 + t567 * t473;
t579 = m(5) * t493 - t499 * mrSges(5,1) + t500 * mrSges(5,2) - t525 * t519 + t526 * t520 + t464;
t577 = m(4) * t516 - t531 * mrSges(4,1) + t532 * mrSges(4,2) - t542 * t535 + t543 * t536 + t579;
t576 = -m(3) * t547 + t554 * mrSges(3,1) - t553 * mrSges(3,2) - t556 * t588 + t557 * t587 - t577;
t460 = m(2) * t558 + qJDD(1) * mrSges(2,1) - t575 * mrSges(2,2) + t576;
t589 = t570 * t442 + t574 * t460;
t443 = t573 * t448 + t569 * t449;
t585 = t574 * t442 - t570 * t460;
t546 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t569 + Ifges(3,4) * t573) * qJD(1);
t545 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t569 + Ifges(3,2) * t573) * qJD(1);
t544 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t569 + Ifges(3,6) * t573) * qJD(1);
t524 = Ifges(4,1) * t543 + Ifges(4,4) * t542 + Ifges(4,5) * qJD(2);
t523 = Ifges(4,4) * t543 + Ifges(4,2) * t542 + Ifges(4,6) * qJD(2);
t522 = Ifges(4,5) * t543 + Ifges(4,6) * t542 + Ifges(4,3) * qJD(2);
t506 = Ifges(5,1) * t526 + Ifges(5,4) * t525 + Ifges(5,5) * t563;
t505 = Ifges(5,4) * t526 + Ifges(5,2) * t525 + Ifges(5,6) * t563;
t504 = Ifges(5,5) * t526 + Ifges(5,6) * t525 + Ifges(5,3) * t563;
t490 = Ifges(6,1) * t518 + Ifges(6,4) * t517 + Ifges(6,5) * t521;
t489 = Ifges(6,4) * t518 + Ifges(6,2) * t517 + Ifges(6,6) * t521;
t488 = Ifges(6,5) * t518 + Ifges(6,6) * t517 + Ifges(6,3) * t521;
t466 = mrSges(6,2) * t476 - mrSges(6,3) * t474 + Ifges(6,1) * t487 + Ifges(6,4) * t486 + Ifges(6,5) * t498 + t517 * t488 - t521 * t489;
t465 = -mrSges(6,1) * t476 + mrSges(6,3) * t475 + Ifges(6,4) * t487 + Ifges(6,2) * t486 + Ifges(6,6) * t498 - t518 * t488 + t521 * t490;
t452 = -mrSges(5,1) * t493 - mrSges(6,1) * t474 + mrSges(6,2) * t475 + mrSges(5,3) * t480 + Ifges(5,4) * t500 - Ifges(6,5) * t487 + Ifges(5,2) * t499 + Ifges(5,6) * t562 - Ifges(6,6) * t486 - Ifges(6,3) * t498 - pkin(4) * t464 - t518 * t489 + t517 * t490 - t526 * t504 + t563 * t506;
t451 = mrSges(5,2) * t493 - mrSges(5,3) * t479 + Ifges(5,1) * t500 + Ifges(5,4) * t499 + Ifges(5,5) * t562 - pkin(8) * t464 - t567 * t465 + t571 * t466 + t525 * t504 - t563 * t505;
t444 = mrSges(4,2) * t516 - mrSges(4,3) * t494 + Ifges(4,1) * t532 + Ifges(4,4) * t531 + Ifges(4,5) * qJDD(2) - pkin(7) * t458 - qJD(2) * t523 + t572 * t451 - t568 * t452 + t542 * t522;
t439 = -mrSges(4,1) * t516 + mrSges(4,3) * t495 + Ifges(4,4) * t532 + Ifges(4,2) * t531 + Ifges(4,6) * qJDD(2) - pkin(3) * t579 + pkin(7) * t582 + qJD(2) * t524 + t568 * t451 + t572 * t452 - t543 * t522;
t438 = -pkin(4) * t578 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t569 * t545 + t573 * t546) * qJD(1) - t571 * t465 + t575 * Ifges(2,5) - t567 * t466 + mrSges(2,3) * t559 - Ifges(5,3) * t562 - t543 * t523 - Ifges(3,5) * t553 - Ifges(3,6) * t554 + mrSges(3,2) * t534 + t542 * t524 - t526 * t505 - Ifges(4,6) * t531 - Ifges(4,5) * t532 - mrSges(3,1) * t533 + t525 * t506 - Ifges(5,6) * t499 - Ifges(5,5) * t500 - mrSges(4,1) * t494 + mrSges(4,2) * t495 - mrSges(5,1) * t479 + mrSges(5,2) * t480 - pkin(3) * t458 - pkin(2) * t450 - pkin(1) * t443 - pkin(8) * t581;
t437 = mrSges(3,2) * t547 - mrSges(3,3) * t533 + Ifges(3,1) * t553 + Ifges(3,4) * t554 + Ifges(3,5) * qJDD(2) - qJ(3) * t450 - qJD(2) * t545 - t565 * t439 + t566 * t444 + t544 * t587;
t436 = -mrSges(3,1) * t547 + mrSges(3,3) * t534 + Ifges(3,4) * t553 + Ifges(3,2) * t554 + Ifges(3,6) * qJDD(2) - pkin(2) * t577 + qJ(3) * t583 + qJD(2) * t546 + t566 * t439 + t565 * t444 - t544 * t588;
t435 = -mrSges(2,2) * g(3) - mrSges(2,3) * t558 + Ifges(2,5) * qJDD(1) - t575 * Ifges(2,6) - pkin(6) * t443 - t569 * t436 + t573 * t437;
t1 = [-m(1) * g(1) + t585; -m(1) * g(2) + t589; (-m(1) - m(2)) * g(3) + t443; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t589 + t574 * t435 - t570 * t438; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t585 + t570 * t435 + t574 * t438; -mrSges(1,1) * g(2) + mrSges(2,1) * t558 + mrSges(1,2) * g(1) - mrSges(2,2) * t559 + Ifges(2,3) * qJDD(1) + pkin(1) * t576 + pkin(6) * t584 + t573 * t436 + t569 * t437;];
tauB = t1;

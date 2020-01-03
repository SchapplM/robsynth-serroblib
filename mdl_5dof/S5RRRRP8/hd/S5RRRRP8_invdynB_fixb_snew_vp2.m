% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP8
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:00:06
% EndTime: 2019-12-31 22:00:14
% DurationCPUTime: 4.44s
% Computational Cost: add. (46540->290), mult. (93123->351), div. (0->0), fcn. (62358->8), ass. (0->112)
t602 = Ifges(5,1) + Ifges(6,1);
t599 = Ifges(5,4) + Ifges(6,4);
t598 = Ifges(5,5) + Ifges(6,5);
t601 = Ifges(5,2) + Ifges(6,2);
t597 = -Ifges(5,6) - Ifges(6,6);
t600 = -Ifges(5,3) - Ifges(6,3);
t572 = sin(qJ(1));
t576 = cos(qJ(1));
t563 = -t576 * g(1) - t572 * g(2);
t578 = qJD(1) ^ 2;
t548 = -t578 * pkin(1) + qJDD(1) * pkin(6) + t563;
t571 = sin(qJ(2));
t575 = cos(qJ(2));
t538 = -t571 * g(3) + t575 * t548;
t556 = (-mrSges(3,1) * t575 + mrSges(3,2) * t571) * qJD(1);
t590 = qJD(1) * qJD(2);
t566 = t571 * t590;
t559 = t575 * qJDD(1) - t566;
t592 = qJD(1) * t571;
t560 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t592;
t562 = t572 * g(1) - t576 * g(2);
t547 = -qJDD(1) * pkin(1) - t578 * pkin(6) - t562;
t587 = t575 * t590;
t558 = t571 * qJDD(1) + t587;
t513 = (-t558 - t587) * pkin(7) + (-t559 + t566) * pkin(2) + t547;
t557 = (-pkin(2) * t575 - pkin(7) * t571) * qJD(1);
t577 = qJD(2) ^ 2;
t591 = t575 * qJD(1);
t517 = -t577 * pkin(2) + qJDD(2) * pkin(7) + t557 * t591 + t538;
t570 = sin(qJ(3));
t574 = cos(qJ(3));
t497 = t574 * t513 - t570 * t517;
t554 = t574 * qJD(2) - t570 * t592;
t530 = t554 * qJD(3) + t570 * qJDD(2) + t574 * t558;
t553 = qJDD(3) - t559;
t555 = t570 * qJD(2) + t574 * t592;
t565 = qJD(3) - t591;
t487 = (t554 * t565 - t530) * pkin(8) + (t554 * t555 + t553) * pkin(3) + t497;
t498 = t570 * t513 + t574 * t517;
t529 = -t555 * qJD(3) + t574 * qJDD(2) - t570 * t558;
t539 = t565 * pkin(3) - t555 * pkin(8);
t552 = t554 ^ 2;
t489 = -t552 * pkin(3) + t529 * pkin(8) - t565 * t539 + t498;
t569 = sin(qJ(4));
t573 = cos(qJ(4));
t481 = t573 * t487 - t569 * t489;
t532 = t573 * t554 - t569 * t555;
t496 = t532 * qJD(4) + t569 * t529 + t573 * t530;
t533 = t569 * t554 + t573 * t555;
t509 = -t532 * mrSges(6,1) + t533 * mrSges(6,2);
t510 = -t532 * mrSges(5,1) + t533 * mrSges(5,2);
t564 = qJD(4) + t565;
t519 = -t564 * mrSges(5,2) + t532 * mrSges(5,3);
t549 = qJDD(4) + t553;
t478 = -0.2e1 * qJD(5) * t533 + (t532 * t564 - t496) * qJ(5) + (t532 * t533 + t549) * pkin(4) + t481;
t518 = -t564 * mrSges(6,2) + t532 * mrSges(6,3);
t589 = m(6) * t478 + t549 * mrSges(6,1) + t564 * t518;
t470 = m(5) * t481 + t549 * mrSges(5,1) + t564 * t519 + (-t509 - t510) * t533 + (-mrSges(5,3) - mrSges(6,3)) * t496 + t589;
t482 = t569 * t487 + t573 * t489;
t495 = -t533 * qJD(4) + t573 * t529 - t569 * t530;
t521 = t564 * mrSges(6,1) - t533 * mrSges(6,3);
t522 = t564 * mrSges(5,1) - t533 * mrSges(5,3);
t520 = t564 * pkin(4) - t533 * qJ(5);
t531 = t532 ^ 2;
t480 = -t531 * pkin(4) + t495 * qJ(5) + 0.2e1 * qJD(5) * t532 - t564 * t520 + t482;
t588 = m(6) * t480 + t495 * mrSges(6,3) + t532 * t509;
t473 = m(5) * t482 + t495 * mrSges(5,3) + t532 * t510 + (-t521 - t522) * t564 + (-mrSges(5,2) - mrSges(6,2)) * t549 + t588;
t468 = t573 * t470 + t569 * t473;
t534 = -t554 * mrSges(4,1) + t555 * mrSges(4,2);
t535 = -t565 * mrSges(4,2) + t554 * mrSges(4,3);
t465 = m(4) * t497 + t553 * mrSges(4,1) - t530 * mrSges(4,3) - t555 * t534 + t565 * t535 + t468;
t536 = t565 * mrSges(4,1) - t555 * mrSges(4,3);
t583 = -t569 * t470 + t573 * t473;
t466 = m(4) * t498 - t553 * mrSges(4,2) + t529 * mrSges(4,3) + t554 * t534 - t565 * t536 + t583;
t584 = -t570 * t465 + t574 * t466;
t461 = m(3) * t538 - qJDD(2) * mrSges(3,2) + t559 * mrSges(3,3) - qJD(2) * t560 + t556 * t591 + t584;
t537 = -t575 * g(3) - t571 * t548;
t561 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t591;
t516 = -qJDD(2) * pkin(2) - t577 * pkin(7) + t557 * t592 - t537;
t490 = -t529 * pkin(3) - t552 * pkin(8) + t555 * t539 + t516;
t484 = -t495 * pkin(4) - t531 * qJ(5) + t533 * t520 + qJDD(5) + t490;
t582 = m(6) * t484 - t495 * mrSges(6,1) + t496 * mrSges(6,2) - t532 * t518 + t533 * t521;
t581 = m(5) * t490 - t495 * mrSges(5,1) + t496 * mrSges(5,2) - t532 * t519 + t533 * t522 + t582;
t579 = -m(4) * t516 + t529 * mrSges(4,1) - t530 * mrSges(4,2) + t554 * t535 - t555 * t536 - t581;
t475 = m(3) * t537 + qJDD(2) * mrSges(3,1) - t558 * mrSges(3,3) + qJD(2) * t561 - t556 * t592 + t579;
t585 = t575 * t461 - t571 * t475;
t455 = m(2) * t563 - t578 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t585;
t462 = t574 * t465 + t570 * t466;
t580 = -m(3) * t547 + t559 * mrSges(3,1) - t558 * mrSges(3,2) - t560 * t592 + t561 * t591 - t462;
t458 = m(2) * t562 + qJDD(1) * mrSges(2,1) - t578 * mrSges(2,2) + t580;
t596 = t572 * t455 + t576 * t458;
t456 = t571 * t461 + t575 * t475;
t595 = t597 * t532 - t598 * t533 + t600 * t564;
t594 = -t601 * t532 - t599 * t533 + t597 * t564;
t593 = t599 * t532 + t602 * t533 + t598 * t564;
t586 = t576 * t455 - t572 * t458;
t546 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t571 + Ifges(3,4) * t575) * qJD(1);
t545 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t571 + Ifges(3,2) * t575) * qJD(1);
t544 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t571 + Ifges(3,6) * t575) * qJD(1);
t525 = Ifges(4,1) * t555 + Ifges(4,4) * t554 + Ifges(4,5) * t565;
t524 = Ifges(4,4) * t555 + Ifges(4,2) * t554 + Ifges(4,6) * t565;
t523 = Ifges(4,5) * t555 + Ifges(4,6) * t554 + Ifges(4,3) * t565;
t476 = -t496 * mrSges(6,3) - t533 * t509 + t589;
t467 = mrSges(5,2) * t490 + mrSges(6,2) * t484 - mrSges(5,3) * t481 - mrSges(6,3) * t478 - qJ(5) * t476 + t599 * t495 + t602 * t496 - t595 * t532 + t598 * t549 + t594 * t564;
t463 = -mrSges(5,1) * t490 + mrSges(5,3) * t482 - mrSges(6,1) * t484 + mrSges(6,3) * t480 - pkin(4) * t582 + qJ(5) * t588 + (-qJ(5) * t521 + t593) * t564 + (-qJ(5) * mrSges(6,2) - t597) * t549 + t595 * t533 + t599 * t496 + t601 * t495;
t452 = mrSges(4,2) * t516 - mrSges(4,3) * t497 + Ifges(4,1) * t530 + Ifges(4,4) * t529 + Ifges(4,5) * t553 - pkin(8) * t468 - t569 * t463 + t573 * t467 + t554 * t523 - t565 * t524;
t451 = -mrSges(4,1) * t516 + mrSges(4,3) * t498 + Ifges(4,4) * t530 + Ifges(4,2) * t529 + Ifges(4,6) * t553 - pkin(3) * t581 + pkin(8) * t583 + t573 * t463 + t569 * t467 - t555 * t523 + t565 * t525;
t450 = -Ifges(4,3) * t553 + t554 * t525 - t555 * t524 + Ifges(3,4) * t558 + Ifges(3,2) * t559 - Ifges(4,6) * t529 - Ifges(4,5) * t530 + mrSges(6,2) * t480 - mrSges(5,1) * t481 + mrSges(5,2) * t482 - pkin(4) * t476 - mrSges(6,1) * t478 - pkin(3) * t468 + t600 * t549 + Ifges(3,6) * qJDD(2) + mrSges(3,3) * t538 + qJD(2) * t546 - mrSges(3,1) * t547 + t597 * t495 - t598 * t496 - t544 * t592 + t593 * t532 + t594 * t533 - pkin(2) * t462 - mrSges(4,1) * t497 + mrSges(4,2) * t498;
t449 = mrSges(3,2) * t547 - mrSges(3,3) * t537 + Ifges(3,1) * t558 + Ifges(3,4) * t559 + Ifges(3,5) * qJDD(2) - pkin(7) * t462 - qJD(2) * t545 - t570 * t451 + t574 * t452 + t544 * t591;
t448 = Ifges(2,6) * qJDD(1) + t578 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t563 - Ifges(3,5) * t558 - Ifges(3,6) * t559 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t537 + mrSges(3,2) * t538 - t570 * t452 - t574 * t451 - pkin(2) * t579 - pkin(7) * t584 - pkin(1) * t456 + (-t571 * t545 + t575 * t546) * qJD(1);
t447 = -mrSges(2,2) * g(3) - mrSges(2,3) * t562 + Ifges(2,5) * qJDD(1) - t578 * Ifges(2,6) - pkin(6) * t456 + t575 * t449 - t571 * t450;
t1 = [-m(1) * g(1) + t586; -m(1) * g(2) + t596; (-m(1) - m(2)) * g(3) + t456; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t596 + t576 * t447 - t572 * t448; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t586 + t572 * t447 + t576 * t448; -mrSges(1,1) * g(2) + mrSges(2,1) * t562 + mrSges(1,2) * g(1) - mrSges(2,2) * t563 + Ifges(2,3) * qJDD(1) + pkin(1) * t580 + pkin(6) * t585 + t571 * t449 + t575 * t450;];
tauB = t1;

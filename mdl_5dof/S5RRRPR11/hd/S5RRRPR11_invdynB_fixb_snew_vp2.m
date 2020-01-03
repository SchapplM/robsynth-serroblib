% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:44
% EndTime: 2019-12-31 21:32:49
% DurationCPUTime: 3.24s
% Computational Cost: add. (30984->291), mult. (60766->351), div. (0->0), fcn. (38285->8), ass. (0->115)
t617 = Ifges(4,1) + Ifges(5,1);
t610 = Ifges(4,4) - Ifges(5,5);
t609 = Ifges(4,5) + Ifges(5,4);
t616 = Ifges(4,2) + Ifges(5,3);
t608 = Ifges(4,6) - Ifges(5,6);
t615 = -Ifges(4,3) - Ifges(5,2);
t578 = sin(qJ(3));
t579 = sin(qJ(2));
t601 = qJD(1) * t579;
t612 = cos(qJ(3));
t558 = -t612 * qJD(2) + t578 * t601;
t582 = cos(qJ(2));
t599 = qJD(1) * qJD(2);
t597 = t582 * t599;
t562 = t579 * qJDD(1) + t597;
t527 = -t558 * qJD(3) + t578 * qJDD(2) + t612 * t562;
t580 = sin(qJ(1));
t583 = cos(qJ(1));
t568 = -t583 * g(1) - t580 * g(2);
t585 = qJD(1) ^ 2;
t549 = -t585 * pkin(1) + qJDD(1) * pkin(6) + t568;
t539 = -t582 * g(3) - t579 * t549;
t561 = (-pkin(2) * t582 - pkin(7) * t579) * qJD(1);
t584 = qJD(2) ^ 2;
t589 = qJDD(2) * pkin(2) + t584 * pkin(7) - t561 * t601 + t539;
t600 = t582 * qJD(1);
t571 = qJD(3) - t600;
t607 = t558 * t571;
t614 = (-t527 + t607) * qJ(4) - t589;
t613 = 2 * qJD(4);
t611 = -mrSges(4,3) - mrSges(5,2);
t540 = -t579 * g(3) + t582 * t549;
t560 = (-mrSges(3,1) * t582 + mrSges(3,2) * t579) * qJD(1);
t598 = t579 * t599;
t563 = t582 * qJDD(1) - t598;
t564 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t601;
t567 = t580 * g(1) - t583 * g(2);
t548 = -qJDD(1) * pkin(1) - t585 * pkin(6) - t567;
t507 = (-t562 - t597) * pkin(7) + (-t563 + t598) * pkin(2) + t548;
t511 = -t584 * pkin(2) + qJDD(2) * pkin(7) + t561 * t600 + t540;
t499 = t578 * t507 + t612 * t511;
t559 = t578 * qJD(2) + t612 * t601;
t526 = t559 * qJD(3) - t612 * qJDD(2) + t578 * t562;
t536 = t571 * mrSges(4,1) - t559 * mrSges(4,3);
t557 = qJDD(3) - t563;
t532 = t558 * pkin(3) - t559 * qJ(4);
t570 = t571 ^ 2;
t491 = -t570 * pkin(3) + t557 * qJ(4) - t558 * t532 + t571 * t613 + t499;
t537 = -t571 * mrSges(5,1) + t559 * mrSges(5,2);
t498 = t612 * t507 - t578 * t511;
t493 = -t557 * pkin(3) - t570 * qJ(4) + t559 * t532 + qJDD(4) - t498;
t486 = (-t527 - t607) * pkin(8) + (t558 * t559 - t557) * pkin(4) + t493;
t541 = -t571 * pkin(4) - t559 * pkin(8);
t556 = t558 ^ 2;
t488 = -t556 * pkin(4) + t526 * pkin(8) + t571 * t541 + t491;
t577 = sin(qJ(5));
t581 = cos(qJ(5));
t484 = t581 * t486 - t577 * t488;
t528 = t581 * t558 - t577 * t559;
t497 = t528 * qJD(5) + t577 * t526 + t581 * t527;
t529 = t577 * t558 + t581 * t559;
t505 = -t528 * mrSges(6,1) + t529 * mrSges(6,2);
t569 = qJD(5) - t571;
t512 = -t569 * mrSges(6,2) + t528 * mrSges(6,3);
t553 = qJDD(5) - t557;
t482 = m(6) * t484 + t553 * mrSges(6,1) - t497 * mrSges(6,3) - t529 * t505 + t569 * t512;
t485 = t577 * t486 + t581 * t488;
t496 = -t529 * qJD(5) + t581 * t526 - t577 * t527;
t513 = t569 * mrSges(6,1) - t529 * mrSges(6,3);
t483 = m(6) * t485 - t553 * mrSges(6,2) + t496 * mrSges(6,3) + t528 * t505 - t569 * t513;
t593 = -t577 * t482 + t581 * t483;
t590 = m(5) * t491 + t557 * mrSges(5,3) + t571 * t537 + t593;
t533 = t558 * mrSges(5,1) - t559 * mrSges(5,3);
t602 = -t558 * mrSges(4,1) - t559 * mrSges(4,2) - t533;
t473 = m(4) * t499 - t557 * mrSges(4,2) + t611 * t526 - t571 * t536 + t602 * t558 + t590;
t535 = -t571 * mrSges(4,2) - t558 * mrSges(4,3);
t475 = t581 * t482 + t577 * t483;
t538 = -t558 * mrSges(5,2) + t571 * mrSges(5,3);
t588 = -m(5) * t493 + t557 * mrSges(5,1) + t571 * t538 - t475;
t474 = m(4) * t498 + t557 * mrSges(4,1) + t611 * t527 + t571 * t535 + t602 * t559 + t588;
t594 = t612 * t473 - t578 * t474;
t470 = m(3) * t540 - qJDD(2) * mrSges(3,2) + t563 * mrSges(3,3) - qJD(2) * t564 + t560 * t600 + t594;
t565 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t600;
t492 = -0.2e1 * qJD(4) * t559 + (t559 * t571 + t526) * pkin(3) + t614;
t489 = -t556 * pkin(8) + (-pkin(3) - pkin(4)) * t526 + (-pkin(3) * t571 + t541 + t613) * t559 - t614;
t591 = -m(6) * t489 + t496 * mrSges(6,1) - t497 * mrSges(6,2) + t528 * t512 - t529 * t513;
t480 = m(5) * t492 + t526 * mrSges(5,1) - t527 * mrSges(5,3) - t559 * t537 + t558 * t538 + t591;
t586 = m(4) * t589 - t526 * mrSges(4,1) - t527 * mrSges(4,2) - t558 * t535 - t559 * t536 - t480;
t479 = m(3) * t539 + qJDD(2) * mrSges(3,1) - t562 * mrSges(3,3) + qJD(2) * t565 - t560 * t601 + t586;
t595 = t582 * t470 - t579 * t479;
t464 = m(2) * t568 - t585 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t595;
t471 = t578 * t473 + t612 * t474;
t587 = -m(3) * t548 + t563 * mrSges(3,1) - t562 * mrSges(3,2) - t564 * t601 + t565 * t600 - t471;
t467 = m(2) * t567 + qJDD(1) * mrSges(2,1) - t585 * mrSges(2,2) + t587;
t606 = t580 * t464 + t583 * t467;
t465 = t579 * t470 + t582 * t479;
t605 = t616 * t558 - t610 * t559 - t608 * t571;
t604 = t608 * t558 - t609 * t559 + t615 * t571;
t603 = -t610 * t558 + t617 * t559 + t609 * t571;
t596 = t583 * t464 - t580 * t467;
t547 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t579 + Ifges(3,4) * t582) * qJD(1);
t546 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t579 + Ifges(3,2) * t582) * qJD(1);
t545 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t579 + Ifges(3,6) * t582) * qJD(1);
t502 = Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * t569;
t501 = Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * t569;
t500 = Ifges(6,5) * t529 + Ifges(6,6) * t528 + Ifges(6,3) * t569;
t477 = mrSges(6,2) * t489 - mrSges(6,3) * t484 + Ifges(6,1) * t497 + Ifges(6,4) * t496 + Ifges(6,5) * t553 + t528 * t500 - t569 * t501;
t476 = -mrSges(6,1) * t489 + mrSges(6,3) * t485 + Ifges(6,4) * t497 + Ifges(6,2) * t496 + Ifges(6,6) * t553 - t529 * t500 + t569 * t502;
t461 = -mrSges(4,2) * t589 + mrSges(5,2) * t493 - mrSges(4,3) * t498 - mrSges(5,3) * t492 - pkin(8) * t475 - qJ(4) * t480 - t577 * t476 + t581 * t477 - t610 * t526 + t617 * t527 + t609 * t557 + t604 * t558 + t605 * t571;
t460 = mrSges(4,1) * t589 - mrSges(5,1) * t492 + mrSges(5,2) * t491 + mrSges(4,3) * t499 - pkin(3) * t480 - pkin(4) * t591 - pkin(8) * t593 - t581 * t476 - t577 * t477 - t616 * t526 + t610 * t527 + t608 * t557 + t604 * t559 + t603 * t571;
t459 = -pkin(3) * t588 - qJ(4) * t590 + (qJ(4) * mrSges(5,2) + t608) * t526 + (pkin(3) * mrSges(5,2) - t609) * t527 - t545 * t601 + (qJ(4) * t533 - t603) * t558 + (pkin(3) * t533 + t605) * t559 + t615 * t557 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t562 + Ifges(3,2) * t563 + mrSges(3,3) * t540 + qJD(2) * t547 - mrSges(3,1) * t548 + Ifges(6,3) * t553 - t528 * t502 + t529 * t501 - mrSges(4,1) * t498 + mrSges(4,2) * t499 - mrSges(5,3) * t491 + mrSges(5,1) * t493 + Ifges(6,6) * t496 + Ifges(6,5) * t497 + mrSges(6,1) * t484 - mrSges(6,2) * t485 + pkin(4) * t475 - pkin(2) * t471;
t458 = mrSges(3,2) * t548 - mrSges(3,3) * t539 + Ifges(3,1) * t562 + Ifges(3,4) * t563 + Ifges(3,5) * qJDD(2) - pkin(7) * t471 - qJD(2) * t546 - t578 * t460 + t612 * t461 + t545 * t600;
t457 = Ifges(2,6) * qJDD(1) + t585 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t568 - Ifges(3,5) * t562 - Ifges(3,6) * t563 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t539 + mrSges(3,2) * t540 - t578 * t461 - t612 * t460 - pkin(2) * t586 - pkin(7) * t594 - pkin(1) * t465 + (-t579 * t546 + t582 * t547) * qJD(1);
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t567 + Ifges(2,5) * qJDD(1) - t585 * Ifges(2,6) - pkin(6) * t465 + t582 * t458 - t579 * t459;
t1 = [-m(1) * g(1) + t596; -m(1) * g(2) + t606; (-m(1) - m(2)) * g(3) + t465; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t606 + t583 * t456 - t580 * t457; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t596 + t580 * t456 + t583 * t457; -mrSges(1,1) * g(2) + mrSges(2,1) * t567 + mrSges(1,2) * g(1) - mrSges(2,2) * t568 + Ifges(2,3) * qJDD(1) + pkin(1) * t587 + pkin(6) * t595 + t579 * t458 + t582 * t459;];
tauB = t1;

% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:23
% EndTime: 2019-12-05 18:38:32
% DurationCPUTime: 9.18s
% Computational Cost: add. (131671->313), mult. (302368->400), div. (0->0), fcn. (215354->10), ass. (0->121)
t596 = qJD(1) ^ 2;
t612 = pkin(2) * t596;
t591 = sin(qJ(1));
t595 = cos(qJ(1));
t579 = -t595 * g(1) - t591 * g(2);
t568 = -t596 * pkin(1) + qJDD(1) * pkin(6) + t579;
t590 = sin(qJ(2));
t611 = t590 * t568;
t594 = cos(qJ(2));
t607 = qJD(1) * qJD(2);
t573 = t590 * qJDD(1) + t594 * t607;
t535 = qJDD(2) * pkin(2) - t573 * pkin(7) - t611 + (pkin(7) * t607 + t590 * t612 - g(3)) * t594;
t556 = -t590 * g(3) + t594 * t568;
t574 = t594 * qJDD(1) - t590 * t607;
t609 = qJD(1) * t590;
t577 = qJD(2) * pkin(2) - pkin(7) * t609;
t585 = t594 ^ 2;
t536 = t574 * pkin(7) - qJD(2) * t577 - t585 * t612 + t556;
t589 = sin(qJ(3));
t593 = cos(qJ(3));
t516 = t593 * t535 - t589 * t536;
t565 = (-t589 * t590 + t593 * t594) * qJD(1);
t539 = t565 * qJD(3) + t593 * t573 + t589 * t574;
t566 = (t589 * t594 + t590 * t593) * qJD(1);
t583 = qJDD(2) + qJDD(3);
t584 = qJD(2) + qJD(3);
t505 = (t565 * t584 - t539) * qJ(4) + (t565 * t566 + t583) * pkin(3) + t516;
t517 = t589 * t535 + t593 * t536;
t538 = -t566 * qJD(3) - t589 * t573 + t593 * t574;
t558 = t584 * pkin(3) - t566 * qJ(4);
t561 = t565 ^ 2;
t507 = -t561 * pkin(3) + t538 * qJ(4) - t584 * t558 + t517;
t586 = sin(pkin(9));
t587 = cos(pkin(9));
t553 = t586 * t565 + t587 * t566;
t495 = -0.2e1 * qJD(4) * t553 + t587 * t505 - t586 * t507;
t521 = t586 * t538 + t587 * t539;
t552 = t587 * t565 - t586 * t566;
t493 = (t552 * t584 - t521) * pkin(8) + (t552 * t553 + t583) * pkin(4) + t495;
t496 = 0.2e1 * qJD(4) * t552 + t586 * t505 + t587 * t507;
t520 = t587 * t538 - t586 * t539;
t543 = t584 * pkin(4) - t553 * pkin(8);
t549 = t552 ^ 2;
t494 = -t549 * pkin(4) + t520 * pkin(8) - t584 * t543 + t496;
t588 = sin(qJ(5));
t592 = cos(qJ(5));
t491 = t592 * t493 - t588 * t494;
t529 = t592 * t552 - t588 * t553;
t502 = t529 * qJD(5) + t588 * t520 + t592 * t521;
t530 = t588 * t552 + t592 * t553;
t513 = -t529 * mrSges(6,1) + t530 * mrSges(6,2);
t581 = qJD(5) + t584;
t522 = -t581 * mrSges(6,2) + t529 * mrSges(6,3);
t580 = qJDD(5) + t583;
t489 = m(6) * t491 + t580 * mrSges(6,1) - t502 * mrSges(6,3) - t530 * t513 + t581 * t522;
t492 = t588 * t493 + t592 * t494;
t501 = -t530 * qJD(5) + t592 * t520 - t588 * t521;
t523 = t581 * mrSges(6,1) - t530 * mrSges(6,3);
t490 = m(6) * t492 - t580 * mrSges(6,2) + t501 * mrSges(6,3) + t529 * t513 - t581 * t523;
t481 = t592 * t489 + t588 * t490;
t531 = -t552 * mrSges(5,1) + t553 * mrSges(5,2);
t541 = -t584 * mrSges(5,2) + t552 * mrSges(5,3);
t479 = m(5) * t495 + t583 * mrSges(5,1) - t521 * mrSges(5,3) - t553 * t531 + t584 * t541 + t481;
t542 = t584 * mrSges(5,1) - t553 * mrSges(5,3);
t602 = -t588 * t489 + t592 * t490;
t480 = m(5) * t496 - t583 * mrSges(5,2) + t520 * mrSges(5,3) + t552 * t531 - t584 * t542 + t602;
t475 = t587 * t479 + t586 * t480;
t554 = -t565 * mrSges(4,1) + t566 * mrSges(4,2);
t557 = -t584 * mrSges(4,2) + t565 * mrSges(4,3);
t473 = m(4) * t516 + t583 * mrSges(4,1) - t539 * mrSges(4,3) - t566 * t554 + t584 * t557 + t475;
t559 = t584 * mrSges(4,1) - t566 * mrSges(4,3);
t603 = -t586 * t479 + t587 * t480;
t474 = m(4) * t517 - t583 * mrSges(4,2) + t538 * mrSges(4,3) + t565 * t554 - t584 * t559 + t603;
t467 = t593 * t473 + t589 * t474;
t555 = -t594 * g(3) - t611;
t572 = (-mrSges(3,1) * t594 + mrSges(3,2) * t590) * qJD(1);
t608 = qJD(1) * t594;
t576 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t608;
t465 = m(3) * t555 + qJDD(2) * mrSges(3,1) - t573 * mrSges(3,3) + qJD(2) * t576 - t572 * t609 + t467;
t575 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t609;
t604 = -t589 * t473 + t593 * t474;
t466 = m(3) * t556 - qJDD(2) * mrSges(3,2) + t574 * mrSges(3,3) - qJD(2) * t575 + t572 * t608 + t604;
t605 = -t590 * t465 + t594 * t466;
t458 = m(2) * t579 - t596 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t605;
t578 = t591 * g(1) - t595 * g(2);
t600 = -qJDD(1) * pkin(1) - t578;
t567 = -t596 * pkin(6) + t600;
t540 = -t574 * pkin(2) + t577 * t609 + (-pkin(7) * t585 - pkin(6)) * t596 + t600;
t515 = -t538 * pkin(3) - t561 * qJ(4) + t566 * t558 + qJDD(4) + t540;
t498 = -t520 * pkin(4) - t549 * pkin(8) + t553 * t543 + t515;
t601 = m(6) * t498 - t501 * mrSges(6,1) + t502 * mrSges(6,2) - t529 * t522 + t530 * t523;
t599 = m(5) * t515 - t520 * mrSges(5,1) + t521 * mrSges(5,2) - t552 * t541 + t553 * t542 + t601;
t598 = m(4) * t540 - t538 * mrSges(4,1) + t539 * mrSges(4,2) - t565 * t557 + t566 * t559 + t599;
t597 = -m(3) * t567 + t574 * mrSges(3,1) - t573 * mrSges(3,2) - t575 * t609 + t576 * t608 - t598;
t485 = m(2) * t578 + qJDD(1) * mrSges(2,1) - t596 * mrSges(2,2) + t597;
t610 = t591 * t458 + t595 * t485;
t459 = t594 * t465 + t590 * t466;
t606 = t595 * t458 - t591 * t485;
t564 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t590 + Ifges(3,4) * t594) * qJD(1);
t563 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t590 + Ifges(3,2) * t594) * qJD(1);
t562 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t590 + Ifges(3,6) * t594) * qJD(1);
t548 = Ifges(4,1) * t566 + Ifges(4,4) * t565 + Ifges(4,5) * t584;
t547 = Ifges(4,4) * t566 + Ifges(4,2) * t565 + Ifges(4,6) * t584;
t546 = Ifges(4,5) * t566 + Ifges(4,6) * t565 + Ifges(4,3) * t584;
t526 = Ifges(5,1) * t553 + Ifges(5,4) * t552 + Ifges(5,5) * t584;
t525 = Ifges(5,4) * t553 + Ifges(5,2) * t552 + Ifges(5,6) * t584;
t524 = Ifges(5,5) * t553 + Ifges(5,6) * t552 + Ifges(5,3) * t584;
t510 = Ifges(6,1) * t530 + Ifges(6,4) * t529 + Ifges(6,5) * t581;
t509 = Ifges(6,4) * t530 + Ifges(6,2) * t529 + Ifges(6,6) * t581;
t508 = Ifges(6,5) * t530 + Ifges(6,6) * t529 + Ifges(6,3) * t581;
t483 = mrSges(6,2) * t498 - mrSges(6,3) * t491 + Ifges(6,1) * t502 + Ifges(6,4) * t501 + Ifges(6,5) * t580 + t529 * t508 - t581 * t509;
t482 = -mrSges(6,1) * t498 + mrSges(6,3) * t492 + Ifges(6,4) * t502 + Ifges(6,2) * t501 + Ifges(6,6) * t580 - t530 * t508 + t581 * t510;
t469 = mrSges(5,2) * t515 - mrSges(5,3) * t495 + Ifges(5,1) * t521 + Ifges(5,4) * t520 + Ifges(5,5) * t583 - pkin(8) * t481 - t588 * t482 + t592 * t483 + t552 * t524 - t584 * t525;
t468 = -mrSges(5,1) * t515 + mrSges(5,3) * t496 + Ifges(5,4) * t521 + Ifges(5,2) * t520 + Ifges(5,6) * t583 - pkin(4) * t601 + pkin(8) * t602 + t592 * t482 + t588 * t483 - t553 * t524 + t584 * t526;
t461 = mrSges(4,2) * t540 - mrSges(4,3) * t516 + Ifges(4,1) * t539 + Ifges(4,4) * t538 + Ifges(4,5) * t583 - qJ(4) * t475 - t586 * t468 + t587 * t469 + t565 * t546 - t584 * t547;
t460 = -mrSges(4,1) * t540 + mrSges(4,3) * t517 + Ifges(4,4) * t539 + Ifges(4,2) * t538 + Ifges(4,6) * t583 - pkin(3) * t599 + qJ(4) * t603 + t587 * t468 + t586 * t469 - t566 * t546 + t584 * t548;
t455 = -pkin(1) * t459 + Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) + (-Ifges(4,3) - Ifges(5,3)) * t583 + (-t590 * t563 + t594 * t564) * qJD(1) + t596 * Ifges(2,5) - Ifges(6,3) * t580 - Ifges(3,5) * t573 - Ifges(3,6) * t574 + mrSges(2,3) * t579 + t565 * t548 - t566 * t547 + t552 * t526 - t553 * t525 - mrSges(3,1) * t555 + mrSges(3,2) * t556 - Ifges(4,6) * t538 - Ifges(4,5) * t539 + t529 * t510 - t530 * t509 - Ifges(5,6) * t520 - Ifges(5,5) * t521 - mrSges(4,1) * t516 + mrSges(4,2) * t517 - Ifges(6,6) * t501 - Ifges(6,5) * t502 - mrSges(5,1) * t495 + mrSges(5,2) * t496 - mrSges(6,1) * t491 + mrSges(6,2) * t492 - pkin(4) * t481 - pkin(3) * t475 - pkin(2) * t467;
t454 = mrSges(3,2) * t567 - mrSges(3,3) * t555 + Ifges(3,1) * t573 + Ifges(3,4) * t574 + Ifges(3,5) * qJDD(2) - pkin(7) * t467 - qJD(2) * t563 - t589 * t460 + t593 * t461 + t562 * t608;
t453 = -mrSges(3,1) * t567 + mrSges(3,3) * t556 + Ifges(3,4) * t573 + Ifges(3,2) * t574 + Ifges(3,6) * qJDD(2) - pkin(2) * t598 + pkin(7) * t604 + qJD(2) * t564 + t593 * t460 + t589 * t461 - t562 * t609;
t452 = -mrSges(2,2) * g(3) - mrSges(2,3) * t578 + Ifges(2,5) * qJDD(1) - t596 * Ifges(2,6) - pkin(6) * t459 - t590 * t453 + t594 * t454;
t1 = [-m(1) * g(1) + t606; -m(1) * g(2) + t610; (-m(1) - m(2)) * g(3) + t459; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t610 + t595 * t452 - t591 * t455; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t606 + t591 * t452 + t595 * t455; -mrSges(1,1) * g(2) + mrSges(2,1) * t578 + mrSges(1,2) * g(1) - mrSges(2,2) * t579 + Ifges(2,3) * qJDD(1) + pkin(1) * t597 + pkin(6) * t605 + t594 * t453 + t590 * t454;];
tauB = t1;

% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:32:04
% EndTime: 2019-12-31 19:32:11
% DurationCPUTime: 7.22s
% Computational Cost: add. (76025->312), mult. (178968->399), div. (0->0), fcn. (122106->10), ass. (0->120)
t596 = -2 * qJD(3);
t570 = sin(qJ(2));
t573 = cos(qJ(2));
t588 = qJD(1) * qJD(2);
t556 = t570 * qJDD(1) + t573 * t588;
t571 = sin(qJ(1));
t574 = cos(qJ(1));
t562 = -t574 * g(1) - t571 * g(2);
t576 = qJD(1) ^ 2;
t551 = -t576 * pkin(1) + qJDD(1) * pkin(6) + t562;
t593 = t570 * t551;
t595 = pkin(2) * t576;
t510 = qJDD(2) * pkin(2) - t556 * qJ(3) - t593 + (qJ(3) * t588 + t570 * t595 - g(3)) * t573;
t537 = -t570 * g(3) + t573 * t551;
t557 = t573 * qJDD(1) - t570 * t588;
t591 = qJD(1) * t570;
t558 = qJD(2) * pkin(2) - qJ(3) * t591;
t565 = t573 ^ 2;
t512 = t557 * qJ(3) - qJD(2) * t558 - t565 * t595 + t537;
t567 = sin(pkin(8));
t594 = cos(pkin(8));
t545 = (t567 * t573 + t594 * t570) * qJD(1);
t496 = t594 * t510 - t567 * t512 + t545 * t596;
t590 = qJD(1) * t573;
t544 = t567 * t591 - t594 * t590;
t497 = t567 * t510 + t594 * t512 + t544 * t596;
t526 = t544 * mrSges(4,1) + t545 * mrSges(4,2);
t529 = t567 * t556 - t594 * t557;
t539 = qJD(2) * mrSges(4,1) - t545 * mrSges(4,3);
t525 = t544 * pkin(3) - t545 * qJ(4);
t575 = qJD(2) ^ 2;
t485 = -t575 * pkin(3) + qJDD(2) * qJ(4) - t544 * t525 + t497;
t561 = t571 * g(1) - t574 * g(2);
t581 = -qJDD(1) * pkin(1) - t561;
t514 = -t557 * pkin(2) + qJDD(3) + t558 * t591 + (-qJ(3) * t565 - pkin(6)) * t576 + t581;
t530 = t594 * t556 + t567 * t557;
t488 = (qJD(2) * t544 - t530) * qJ(4) + (qJD(2) * t545 + t529) * pkin(3) + t514;
t566 = sin(pkin(9));
t568 = cos(pkin(9));
t535 = t566 * qJD(2) + t568 * t545;
t480 = -0.2e1 * qJD(4) * t535 - t566 * t485 + t568 * t488;
t520 = t566 * qJDD(2) + t568 * t530;
t534 = t568 * qJD(2) - t566 * t545;
t478 = (t534 * t544 - t520) * pkin(7) + (t534 * t535 + t529) * pkin(4) + t480;
t481 = 0.2e1 * qJD(4) * t534 + t568 * t485 + t566 * t488;
t517 = t544 * pkin(4) - t535 * pkin(7);
t519 = t568 * qJDD(2) - t566 * t530;
t533 = t534 ^ 2;
t479 = -t533 * pkin(4) + t519 * pkin(7) - t544 * t517 + t481;
t569 = sin(qJ(5));
t572 = cos(qJ(5));
t476 = t572 * t478 - t569 * t479;
t506 = t572 * t534 - t569 * t535;
t491 = t506 * qJD(5) + t569 * t519 + t572 * t520;
t507 = t569 * t534 + t572 * t535;
t498 = -t506 * mrSges(6,1) + t507 * mrSges(6,2);
t543 = qJD(5) + t544;
t499 = -t543 * mrSges(6,2) + t506 * mrSges(6,3);
t528 = qJDD(5) + t529;
t474 = m(6) * t476 + t528 * mrSges(6,1) - t491 * mrSges(6,3) - t507 * t498 + t543 * t499;
t477 = t569 * t478 + t572 * t479;
t490 = -t507 * qJD(5) + t572 * t519 - t569 * t520;
t500 = t543 * mrSges(6,1) - t507 * mrSges(6,3);
t475 = m(6) * t477 - t528 * mrSges(6,2) + t490 * mrSges(6,3) + t506 * t498 - t543 * t500;
t466 = t572 * t474 + t569 * t475;
t511 = -t534 * mrSges(5,1) + t535 * mrSges(5,2);
t515 = -t544 * mrSges(5,2) + t534 * mrSges(5,3);
t464 = m(5) * t480 + t529 * mrSges(5,1) - t520 * mrSges(5,3) - t535 * t511 + t544 * t515 + t466;
t516 = t544 * mrSges(5,1) - t535 * mrSges(5,3);
t583 = -t569 * t474 + t572 * t475;
t465 = m(5) * t481 - t529 * mrSges(5,2) + t519 * mrSges(5,3) + t534 * t511 - t544 * t516 + t583;
t584 = -t566 * t464 + t568 * t465;
t459 = m(4) * t497 - qJDD(2) * mrSges(4,2) - t529 * mrSges(4,3) - qJD(2) * t539 - t544 * t526 + t584;
t538 = -qJD(2) * mrSges(4,2) - t544 * mrSges(4,3);
t484 = -qJDD(2) * pkin(3) - t575 * qJ(4) + t545 * t525 + qJDD(4) - t496;
t482 = -t519 * pkin(4) - t533 * pkin(7) + t535 * t517 + t484;
t580 = m(6) * t482 - t490 * mrSges(6,1) + t491 * mrSges(6,2) - t506 * t499 + t507 * t500;
t578 = -m(5) * t484 + t519 * mrSges(5,1) - t520 * mrSges(5,2) + t534 * t515 - t535 * t516 - t580;
t470 = m(4) * t496 + qJDD(2) * mrSges(4,1) - t530 * mrSges(4,3) + qJD(2) * t538 - t545 * t526 + t578;
t452 = t567 * t459 + t594 * t470;
t536 = -t573 * g(3) - t593;
t555 = (-mrSges(3,1) * t573 + mrSges(3,2) * t570) * qJD(1);
t560 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t590;
t450 = m(3) * t536 + qJDD(2) * mrSges(3,1) - t556 * mrSges(3,3) + qJD(2) * t560 - t555 * t591 + t452;
t559 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t591;
t585 = t594 * t459 - t567 * t470;
t451 = m(3) * t537 - qJDD(2) * mrSges(3,2) + t557 * mrSges(3,3) - qJD(2) * t559 + t555 * t590 + t585;
t586 = -t570 * t450 + t573 * t451;
t444 = m(2) * t562 - t576 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t586;
t550 = -t576 * pkin(6) + t581;
t460 = t568 * t464 + t566 * t465;
t579 = m(4) * t514 + t529 * mrSges(4,1) + t530 * mrSges(4,2) + t544 * t538 + t545 * t539 + t460;
t577 = -m(3) * t550 + t557 * mrSges(3,1) - t556 * mrSges(3,2) - t559 * t591 + t560 * t590 - t579;
t456 = m(2) * t561 + qJDD(1) * mrSges(2,1) - t576 * mrSges(2,2) + t577;
t592 = t571 * t444 + t574 * t456;
t445 = t573 * t450 + t570 * t451;
t587 = t574 * t444 - t571 * t456;
t548 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t570 + Ifges(3,4) * t573) * qJD(1);
t547 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t570 + Ifges(3,2) * t573) * qJD(1);
t546 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t570 + Ifges(3,6) * t573) * qJD(1);
t523 = Ifges(4,1) * t545 - Ifges(4,4) * t544 + Ifges(4,5) * qJD(2);
t522 = Ifges(4,4) * t545 - Ifges(4,2) * t544 + Ifges(4,6) * qJD(2);
t521 = Ifges(4,5) * t545 - Ifges(4,6) * t544 + Ifges(4,3) * qJD(2);
t503 = Ifges(5,1) * t535 + Ifges(5,4) * t534 + Ifges(5,5) * t544;
t502 = Ifges(5,4) * t535 + Ifges(5,2) * t534 + Ifges(5,6) * t544;
t501 = Ifges(5,5) * t535 + Ifges(5,6) * t534 + Ifges(5,3) * t544;
t494 = Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * t543;
t493 = Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * t543;
t492 = Ifges(6,5) * t507 + Ifges(6,6) * t506 + Ifges(6,3) * t543;
t468 = mrSges(6,2) * t482 - mrSges(6,3) * t476 + Ifges(6,1) * t491 + Ifges(6,4) * t490 + Ifges(6,5) * t528 + t506 * t492 - t543 * t493;
t467 = -mrSges(6,1) * t482 + mrSges(6,3) * t477 + Ifges(6,4) * t491 + Ifges(6,2) * t490 + Ifges(6,6) * t528 - t507 * t492 + t543 * t494;
t454 = mrSges(5,2) * t484 - mrSges(5,3) * t480 + Ifges(5,1) * t520 + Ifges(5,4) * t519 + Ifges(5,5) * t529 - pkin(7) * t466 - t569 * t467 + t572 * t468 + t534 * t501 - t544 * t502;
t453 = -mrSges(5,1) * t484 + mrSges(5,3) * t481 + Ifges(5,4) * t520 + Ifges(5,2) * t519 + Ifges(5,6) * t529 - pkin(4) * t580 + pkin(7) * t583 + t572 * t467 + t569 * t468 - t535 * t501 + t544 * t503;
t446 = Ifges(4,4) * t530 + Ifges(4,6) * qJDD(2) - t545 * t521 + qJD(2) * t523 - mrSges(4,1) * t514 + mrSges(4,3) * t497 - Ifges(5,5) * t520 - Ifges(5,6) * t519 - t535 * t502 + t534 * t503 - mrSges(5,1) * t480 + mrSges(5,2) * t481 - Ifges(6,5) * t491 - Ifges(6,6) * t490 - Ifges(6,3) * t528 - t507 * t493 + t506 * t494 - mrSges(6,1) * t476 + mrSges(6,2) * t477 - pkin(4) * t466 - pkin(3) * t460 + (-Ifges(4,2) - Ifges(5,3)) * t529;
t441 = mrSges(4,2) * t514 - mrSges(4,3) * t496 + Ifges(4,1) * t530 - Ifges(4,4) * t529 + Ifges(4,5) * qJDD(2) - qJ(4) * t460 - qJD(2) * t522 - t566 * t453 + t568 * t454 - t544 * t521;
t440 = mrSges(3,2) * t550 - mrSges(3,3) * t536 + Ifges(3,1) * t556 + Ifges(3,4) * t557 + Ifges(3,5) * qJDD(2) - qJ(3) * t452 - qJD(2) * t547 + t594 * t441 - t567 * t446 + t546 * t590;
t439 = -pkin(1) * t445 + mrSges(2,3) * t562 - pkin(2) * t452 - Ifges(3,5) * t556 - Ifges(3,6) * t557 - mrSges(3,1) * t536 + mrSges(3,2) * t537 - t568 * t453 - pkin(3) * t578 - qJ(4) * t584 - Ifges(4,5) * t530 + Ifges(4,6) * t529 - mrSges(4,1) * t496 + mrSges(4,2) * t497 - t566 * t454 + mrSges(2,1) * g(3) - t545 * t522 - t544 * t523 + t576 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t570 * t547 + t573 * t548) * qJD(1);
t438 = -mrSges(3,1) * t550 + mrSges(3,3) * t537 + Ifges(3,4) * t556 + Ifges(3,2) * t557 + Ifges(3,6) * qJDD(2) - pkin(2) * t579 + qJ(3) * t585 + qJD(2) * t548 + t567 * t441 + t594 * t446 - t546 * t591;
t437 = -mrSges(2,2) * g(3) - mrSges(2,3) * t561 + Ifges(2,5) * qJDD(1) - t576 * Ifges(2,6) - pkin(6) * t445 - t570 * t438 + t573 * t440;
t1 = [-m(1) * g(1) + t587; -m(1) * g(2) + t592; (-m(1) - m(2)) * g(3) + t445; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t592 + t574 * t437 - t571 * t439; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t587 + t571 * t437 + t574 * t439; -mrSges(1,1) * g(2) + mrSges(2,1) * t561 + mrSges(1,2) * g(1) - mrSges(2,2) * t562 + Ifges(2,3) * qJDD(1) + pkin(1) * t577 + pkin(6) * t586 + t573 * t438 + t570 * t440;];
tauB = t1;

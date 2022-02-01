% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:25:20
% EndTime: 2022-01-23 09:25:27
% DurationCPUTime: 6.93s
% Computational Cost: add. (65060->288), mult. (167253->380), div. (0->0), fcn. (114048->10), ass. (0->121)
t575 = sin(qJ(1));
t578 = cos(qJ(1));
t556 = -t578 * g(1) - t575 * g(2);
t579 = qJD(1) ^ 2;
t610 = -t579 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t556;
t570 = sin(pkin(8));
t572 = cos(pkin(8));
t527 = -t572 * g(3) - t610 * t570;
t600 = t572 * qJD(1);
t559 = qJD(3) - t600;
t574 = sin(qJ(3));
t601 = t570 * qJD(1);
t595 = t574 * t601;
t542 = -t559 * mrSges(4,2) - mrSges(4,3) * t595;
t577 = cos(qJ(3));
t594 = t577 * t601;
t544 = t559 * mrSges(4,1) - mrSges(4,3) * t594;
t609 = -t542 * t574 - t544 * t577;
t608 = mrSges(3,2) * t570;
t605 = t570 ^ 2 * t579;
t528 = -t570 * g(3) + t610 * t572;
t551 = (-mrSges(3,1) * t572 + t608) * qJD(1);
t587 = -pkin(2) * t572 - pkin(6) * t570;
t553 = t587 * qJD(1);
t517 = t553 * t600 + t528;
t555 = t575 * g(1) - t578 * g(2);
t583 = -t579 * qJ(2) + qJDD(2) - t555;
t529 = (-pkin(1) + t587) * qJDD(1) + t583;
t526 = t577 * t529;
t598 = qJD(1) * qJD(3);
t547 = (qJDD(1) * t577 - t574 * t598) * t570;
t597 = t572 * qJDD(1);
t558 = qJDD(3) - t597;
t496 = t558 * pkin(3) - t547 * qJ(4) + t526 + (-pkin(3) * t577 * t605 - qJ(4) * t559 * t601 - t517) * t574;
t506 = t577 * t517 + t574 * t529;
t543 = t559 * pkin(3) - qJ(4) * t594;
t546 = (-qJDD(1) * t574 - t577 * t598) * t570;
t596 = t574 ^ 2 * t605;
t497 = -pkin(3) * t596 + t546 * qJ(4) - t559 * t543 + t506;
t569 = sin(pkin(9));
t571 = cos(pkin(9));
t538 = (-t569 * t574 + t571 * t577) * t601;
t485 = -0.2e1 * qJD(4) * t538 + t571 * t496 - t569 * t497;
t521 = t569 * t546 + t571 * t547;
t537 = (-t569 * t577 - t571 * t574) * t601;
t483 = (t537 * t559 - t521) * pkin(7) + (t537 * t538 + t558) * pkin(4) + t485;
t486 = 0.2e1 * qJD(4) * t537 + t569 * t496 + t571 * t497;
t520 = t571 * t546 - t569 * t547;
t524 = t559 * pkin(4) - t538 * pkin(7);
t536 = t537 ^ 2;
t484 = -t536 * pkin(4) + t520 * pkin(7) - t559 * t524 + t486;
t573 = sin(qJ(5));
t576 = cos(qJ(5));
t481 = t576 * t483 - t573 * t484;
t514 = t576 * t537 - t573 * t538;
t492 = t514 * qJD(5) + t573 * t520 + t576 * t521;
t515 = t573 * t537 + t576 * t538;
t503 = -t514 * mrSges(6,1) + t515 * mrSges(6,2);
t557 = qJD(5) + t559;
t507 = -t557 * mrSges(6,2) + t514 * mrSges(6,3);
t554 = qJDD(5) + t558;
t479 = m(6) * t481 + t554 * mrSges(6,1) - t492 * mrSges(6,3) - t515 * t503 + t557 * t507;
t482 = t573 * t483 + t576 * t484;
t491 = -t515 * qJD(5) + t576 * t520 - t573 * t521;
t508 = t557 * mrSges(6,1) - t515 * mrSges(6,3);
t480 = m(6) * t482 - t554 * mrSges(6,2) + t491 * mrSges(6,3) + t514 * t503 - t557 * t508;
t471 = t576 * t479 + t573 * t480;
t518 = -t537 * mrSges(5,1) + t538 * mrSges(5,2);
t522 = -t559 * mrSges(5,2) + t537 * mrSges(5,3);
t469 = m(5) * t485 + t558 * mrSges(5,1) - t521 * mrSges(5,3) - t538 * t518 + t559 * t522 + t471;
t523 = t559 * mrSges(5,1) - t538 * mrSges(5,3);
t588 = -t573 * t479 + t576 * t480;
t470 = m(5) * t486 - t558 * mrSges(5,2) + t520 * mrSges(5,3) + t537 * t518 - t559 * t523 + t588;
t465 = t571 * t469 + t569 * t470;
t505 = -t574 * t517 + t526;
t545 = (mrSges(4,1) * t574 + mrSges(4,2) * t577) * t601;
t463 = m(4) * t505 + t558 * mrSges(4,1) - t547 * mrSges(4,3) + t559 * t542 - t545 * t594 + t465;
t589 = -t569 * t469 + t571 * t470;
t464 = m(4) * t506 - t558 * mrSges(4,2) + t546 * mrSges(4,3) - t559 * t544 - t545 * t595 + t589;
t590 = -t574 * t463 + t577 * t464;
t602 = qJDD(1) * mrSges(3,3);
t458 = m(3) * t528 + (qJD(1) * t551 + t602) * t572 + t590;
t516 = t553 * t601 - t527;
t504 = -t546 * pkin(3) - qJ(4) * t596 + t543 * t594 + qJDD(4) + t516;
t488 = -t520 * pkin(4) - t536 * pkin(7) + t538 * t524 + t504;
t584 = m(6) * t488 - t491 * mrSges(6,1) + t492 * mrSges(6,2) - t514 * t507 + t515 * t508;
t581 = m(5) * t504 - t520 * mrSges(5,1) + t521 * mrSges(5,2) - t537 * t522 + t538 * t523 + t584;
t580 = -m(4) * t516 + t546 * mrSges(4,1) - t547 * mrSges(4,2) - t581;
t478 = t580 + (-t602 + (-t551 + t609) * qJD(1)) * t570 + m(3) * t527;
t591 = t572 * t458 - t570 * t478;
t452 = m(2) * t556 - t579 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t591;
t459 = t577 * t463 + t574 * t464;
t549 = -qJDD(1) * pkin(1) + t583;
t582 = -m(3) * t549 + mrSges(3,1) * t597 - t459 + (t572 ^ 2 * t579 + t605) * mrSges(3,3);
t455 = m(2) * t555 - t579 * mrSges(2,2) + (mrSges(2,1) - t608) * qJDD(1) + t582;
t604 = t575 * t452 + t578 * t455;
t453 = t570 * t458 + t572 * t478;
t592 = t578 * t452 - t575 * t455;
t586 = Ifges(3,1) * t570 + Ifges(3,4) * t572;
t585 = Ifges(3,5) * t570 + Ifges(3,6) * t572;
t552 = t585 * qJD(1);
t532 = Ifges(4,5) * t559 + (Ifges(4,1) * t577 - Ifges(4,4) * t574) * t601;
t531 = Ifges(4,6) * t559 + (Ifges(4,4) * t577 - Ifges(4,2) * t574) * t601;
t530 = Ifges(4,3) * t559 + (Ifges(4,5) * t577 - Ifges(4,6) * t574) * t601;
t511 = Ifges(5,1) * t538 + Ifges(5,4) * t537 + Ifges(5,5) * t559;
t510 = Ifges(5,4) * t538 + Ifges(5,2) * t537 + Ifges(5,6) * t559;
t509 = Ifges(5,5) * t538 + Ifges(5,6) * t537 + Ifges(5,3) * t559;
t500 = Ifges(6,1) * t515 + Ifges(6,4) * t514 + Ifges(6,5) * t557;
t499 = Ifges(6,4) * t515 + Ifges(6,2) * t514 + Ifges(6,6) * t557;
t498 = Ifges(6,5) * t515 + Ifges(6,6) * t514 + Ifges(6,3) * t557;
t473 = mrSges(6,2) * t488 - mrSges(6,3) * t481 + Ifges(6,1) * t492 + Ifges(6,4) * t491 + Ifges(6,5) * t554 + t514 * t498 - t557 * t499;
t472 = -mrSges(6,1) * t488 + mrSges(6,3) * t482 + Ifges(6,4) * t492 + Ifges(6,2) * t491 + Ifges(6,6) * t554 - t515 * t498 + t557 * t500;
t461 = mrSges(5,2) * t504 - mrSges(5,3) * t485 + Ifges(5,1) * t521 + Ifges(5,4) * t520 + Ifges(5,5) * t558 - pkin(7) * t471 - t573 * t472 + t576 * t473 + t537 * t509 - t559 * t510;
t460 = -mrSges(5,1) * t504 + mrSges(5,3) * t486 + Ifges(5,4) * t521 + Ifges(5,2) * t520 + Ifges(5,6) * t558 - pkin(4) * t584 + pkin(7) * t588 + t576 * t472 + t573 * t473 - t538 * t509 + t559 * t511;
t449 = mrSges(4,2) * t516 - mrSges(4,3) * t505 + Ifges(4,1) * t547 + Ifges(4,4) * t546 + Ifges(4,5) * t558 - qJ(4) * t465 - t569 * t460 + t571 * t461 - t530 * t595 - t559 * t531;
t448 = -mrSges(4,1) * t516 + mrSges(4,3) * t506 + Ifges(4,4) * t547 + Ifges(4,2) * t546 + Ifges(4,6) * t558 - pkin(3) * t581 + qJ(4) * t589 + t571 * t460 + t569 * t461 - t530 * t594 + t559 * t532;
t447 = -pkin(2) * t459 + Ifges(3,2) * t597 + (Ifges(3,4) * qJDD(1) + (-t531 * t577 - t532 * t574 - t552) * qJD(1)) * t570 + (-Ifges(4,3) - Ifges(5,3)) * t558 - Ifges(4,5) * t547 - mrSges(3,1) * t549 - Ifges(6,3) * t554 - Ifges(4,6) * t546 + t537 * t511 - t538 * t510 - Ifges(5,5) * t521 + mrSges(3,3) * t528 + t514 * t500 - t515 * t499 - Ifges(5,6) * t520 - mrSges(4,1) * t505 + mrSges(4,2) * t506 - Ifges(6,6) * t491 - Ifges(6,5) * t492 - mrSges(5,1) * t485 + mrSges(5,2) * t486 + mrSges(6,2) * t482 - mrSges(6,1) * t481 - pkin(3) * t465 - pkin(4) * t471;
t446 = mrSges(3,2) * t549 - mrSges(3,3) * t527 - pkin(6) * t459 + t586 * qJDD(1) - t574 * t448 + t577 * t449 + t552 * t600;
t445 = t579 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t556 - mrSges(3,1) * t527 + mrSges(3,2) * t528 - t574 * t449 - t577 * t448 - pkin(2) * t580 - pkin(6) * t590 - pkin(1) * t453 + (Ifges(2,6) - t585) * qJDD(1) + (-pkin(2) * t609 * t570 + (-t570 * (Ifges(3,4) * t570 + Ifges(3,2) * t572) + t572 * t586) * qJD(1)) * qJD(1);
t444 = -mrSges(2,2) * g(3) - mrSges(2,3) * t555 + Ifges(2,5) * qJDD(1) - t579 * Ifges(2,6) - qJ(2) * t453 + t572 * t446 - t570 * t447;
t1 = [-m(1) * g(1) + t592; -m(1) * g(2) + t604; (-m(1) - m(2)) * g(3) + t453; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t604 + t578 * t444 - t575 * t445; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t592 + t575 * t444 + t578 * t445; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t555 - mrSges(2,2) * t556 + t570 * t446 + t572 * t447 + pkin(1) * (-qJDD(1) * t608 + t582) + qJ(2) * t591;];
tauB = t1;

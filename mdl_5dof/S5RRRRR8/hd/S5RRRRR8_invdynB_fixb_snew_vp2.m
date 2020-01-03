% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:56
% EndTime: 2019-12-31 22:25:03
% DurationCPUTime: 6.30s
% Computational Cost: add. (98813->314), mult. (198688->398), div. (0->0), fcn. (139534->10), ass. (0->123)
t588 = qJD(1) ^ 2;
t604 = pkin(2) * t588;
t582 = sin(qJ(1));
t587 = cos(qJ(1));
t571 = -g(1) * t587 - g(2) * t582;
t559 = -pkin(1) * t588 + qJDD(1) * pkin(6) + t571;
t581 = sin(qJ(2));
t603 = t581 * t559;
t586 = cos(qJ(2));
t599 = qJD(1) * qJD(2);
t565 = qJDD(1) * t581 + t586 * t599;
t524 = qJDD(2) * pkin(2) - t565 * pkin(7) - t603 + (pkin(7) * t599 + t581 * t604 - g(3)) * t586;
t546 = -g(3) * t581 + t586 * t559;
t566 = qJDD(1) * t586 - t581 * t599;
t601 = qJD(1) * t581;
t569 = qJD(2) * pkin(2) - pkin(7) * t601;
t577 = t586 ^ 2;
t525 = pkin(7) * t566 - qJD(2) * t569 - t577 * t604 + t546;
t580 = sin(qJ(3));
t585 = cos(qJ(3));
t508 = t580 * t524 + t585 * t525;
t557 = (t580 * t586 + t581 * t585) * qJD(1);
t530 = -t557 * qJD(3) - t580 * t565 + t566 * t585;
t600 = qJD(1) * t586;
t556 = -t580 * t601 + t585 * t600;
t540 = -mrSges(4,1) * t556 + mrSges(4,2) * t557;
t576 = qJD(2) + qJD(3);
t548 = mrSges(4,1) * t576 - mrSges(4,3) * t557;
t575 = qJDD(2) + qJDD(3);
t531 = qJD(3) * t556 + t565 * t585 + t566 * t580;
t570 = t582 * g(1) - t587 * g(2);
t593 = -qJDD(1) * pkin(1) - t570;
t532 = -t566 * pkin(2) + t569 * t601 + (-pkin(7) * t577 - pkin(6)) * t588 + t593;
t498 = (-t556 * t576 - t531) * pkin(8) + (t557 * t576 - t530) * pkin(3) + t532;
t541 = -pkin(3) * t556 - pkin(8) * t557;
t574 = t576 ^ 2;
t501 = -pkin(3) * t574 + pkin(8) * t575 + t541 * t556 + t508;
t579 = sin(qJ(4));
t584 = cos(qJ(4));
t490 = t584 * t498 - t579 * t501;
t543 = -t557 * t579 + t576 * t584;
t511 = qJD(4) * t543 + t531 * t584 + t575 * t579;
t529 = qJDD(4) - t530;
t544 = t557 * t584 + t576 * t579;
t552 = qJD(4) - t556;
t488 = (t543 * t552 - t511) * pkin(9) + (t543 * t544 + t529) * pkin(4) + t490;
t491 = t579 * t498 + t584 * t501;
t510 = -qJD(4) * t544 - t531 * t579 + t575 * t584;
t535 = pkin(4) * t552 - pkin(9) * t544;
t542 = t543 ^ 2;
t489 = -pkin(4) * t542 + pkin(9) * t510 - t535 * t552 + t491;
t578 = sin(qJ(5));
t583 = cos(qJ(5));
t486 = t488 * t583 - t489 * t578;
t518 = t543 * t583 - t544 * t578;
t495 = qJD(5) * t518 + t510 * t578 + t511 * t583;
t519 = t543 * t578 + t544 * t583;
t506 = -mrSges(6,1) * t518 + mrSges(6,2) * t519;
t550 = qJD(5) + t552;
t512 = -mrSges(6,2) * t550 + mrSges(6,3) * t518;
t526 = qJDD(5) + t529;
t484 = m(6) * t486 + mrSges(6,1) * t526 - mrSges(6,3) * t495 - t506 * t519 + t512 * t550;
t487 = t488 * t578 + t489 * t583;
t494 = -qJD(5) * t519 + t510 * t583 - t511 * t578;
t513 = mrSges(6,1) * t550 - mrSges(6,3) * t519;
t485 = m(6) * t487 - mrSges(6,2) * t526 + mrSges(6,3) * t494 + t506 * t518 - t513 * t550;
t476 = t583 * t484 + t578 * t485;
t523 = -mrSges(5,1) * t543 + mrSges(5,2) * t544;
t533 = -mrSges(5,2) * t552 + mrSges(5,3) * t543;
t474 = m(5) * t490 + mrSges(5,1) * t529 - mrSges(5,3) * t511 - t523 * t544 + t533 * t552 + t476;
t534 = mrSges(5,1) * t552 - mrSges(5,3) * t544;
t594 = -t484 * t578 + t583 * t485;
t475 = m(5) * t491 - mrSges(5,2) * t529 + mrSges(5,3) * t510 + t523 * t543 - t534 * t552 + t594;
t595 = -t474 * t579 + t475 * t584;
t469 = m(4) * t508 - mrSges(4,2) * t575 + mrSges(4,3) * t530 + t540 * t556 - t548 * t576 + t595;
t507 = t524 * t585 - t580 * t525;
t547 = -mrSges(4,2) * t576 + mrSges(4,3) * t556;
t500 = -pkin(3) * t575 - pkin(8) * t574 + t557 * t541 - t507;
t492 = -pkin(4) * t510 - pkin(9) * t542 + t535 * t544 + t500;
t592 = m(6) * t492 - t494 * mrSges(6,1) + mrSges(6,2) * t495 - t518 * t512 + t513 * t519;
t590 = -m(5) * t500 + t510 * mrSges(5,1) - mrSges(5,2) * t511 + t543 * t533 - t534 * t544 - t592;
t480 = m(4) * t507 + mrSges(4,1) * t575 - mrSges(4,3) * t531 - t540 * t557 + t547 * t576 + t590;
t463 = t469 * t580 + t585 * t480;
t545 = -t586 * g(3) - t603;
t564 = (-mrSges(3,1) * t586 + mrSges(3,2) * t581) * qJD(1);
t568 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t600;
t460 = m(3) * t545 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t565 + qJD(2) * t568 - t564 * t601 + t463;
t567 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t601;
t596 = t469 * t585 - t480 * t580;
t461 = m(3) * t546 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t566 - qJD(2) * t567 + t564 * t600 + t596;
t597 = -t460 * t581 + t461 * t586;
t454 = m(2) * t571 - mrSges(2,1) * t588 - qJDD(1) * mrSges(2,2) + t597;
t558 = -t588 * pkin(6) + t593;
t470 = t474 * t584 + t475 * t579;
t591 = m(4) * t532 - t530 * mrSges(4,1) + mrSges(4,2) * t531 - t556 * t547 + t548 * t557 + t470;
t589 = -m(3) * t558 + t566 * mrSges(3,1) - mrSges(3,2) * t565 - t567 * t601 + t568 * t600 - t591;
t466 = m(2) * t570 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t588 + t589;
t602 = t454 * t582 + t466 * t587;
t455 = t460 * t586 + t461 * t581;
t598 = t454 * t587 - t466 * t582;
t555 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t581 + Ifges(3,4) * t586) * qJD(1);
t554 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t581 + Ifges(3,2) * t586) * qJD(1);
t553 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t581 + Ifges(3,6) * t586) * qJD(1);
t538 = Ifges(4,1) * t557 + Ifges(4,4) * t556 + Ifges(4,5) * t576;
t537 = Ifges(4,4) * t557 + Ifges(4,2) * t556 + Ifges(4,6) * t576;
t536 = Ifges(4,5) * t557 + Ifges(4,6) * t556 + Ifges(4,3) * t576;
t516 = Ifges(5,1) * t544 + Ifges(5,4) * t543 + Ifges(5,5) * t552;
t515 = Ifges(5,4) * t544 + Ifges(5,2) * t543 + Ifges(5,6) * t552;
t514 = Ifges(5,5) * t544 + Ifges(5,6) * t543 + Ifges(5,3) * t552;
t504 = Ifges(6,1) * t519 + Ifges(6,4) * t518 + Ifges(6,5) * t550;
t503 = Ifges(6,4) * t519 + Ifges(6,2) * t518 + Ifges(6,6) * t550;
t502 = Ifges(6,5) * t519 + Ifges(6,6) * t518 + Ifges(6,3) * t550;
t478 = mrSges(6,2) * t492 - mrSges(6,3) * t486 + Ifges(6,1) * t495 + Ifges(6,4) * t494 + Ifges(6,5) * t526 + t502 * t518 - t503 * t550;
t477 = -mrSges(6,1) * t492 + mrSges(6,3) * t487 + Ifges(6,4) * t495 + Ifges(6,2) * t494 + Ifges(6,6) * t526 - t502 * t519 + t504 * t550;
t464 = mrSges(5,2) * t500 - mrSges(5,3) * t490 + Ifges(5,1) * t511 + Ifges(5,4) * t510 + Ifges(5,5) * t529 - pkin(9) * t476 - t477 * t578 + t478 * t583 + t514 * t543 - t515 * t552;
t462 = -mrSges(5,1) * t500 + mrSges(5,3) * t491 + Ifges(5,4) * t511 + Ifges(5,2) * t510 + Ifges(5,6) * t529 - pkin(4) * t592 + pkin(9) * t594 + t583 * t477 + t578 * t478 - t544 * t514 + t552 * t516;
t456 = Ifges(4,4) * t531 + Ifges(4,2) * t530 + Ifges(4,6) * t575 - t557 * t536 + t576 * t538 - mrSges(4,1) * t532 + mrSges(4,3) * t508 - Ifges(5,5) * t511 - Ifges(5,6) * t510 - Ifges(5,3) * t529 - t544 * t515 + t543 * t516 - mrSges(5,1) * t490 + mrSges(5,2) * t491 - Ifges(6,5) * t495 - Ifges(6,6) * t494 - Ifges(6,3) * t526 - t519 * t503 + t518 * t504 - mrSges(6,1) * t486 + mrSges(6,2) * t487 - pkin(4) * t476 - pkin(3) * t470;
t451 = mrSges(4,2) * t532 - mrSges(4,3) * t507 + Ifges(4,1) * t531 + Ifges(4,4) * t530 + Ifges(4,5) * t575 - pkin(8) * t470 - t462 * t579 + t464 * t584 + t536 * t556 - t537 * t576;
t450 = mrSges(3,2) * t558 - mrSges(3,3) * t545 + Ifges(3,1) * t565 + Ifges(3,4) * t566 + Ifges(3,5) * qJDD(2) - pkin(7) * t463 - qJD(2) * t554 + t451 * t585 - t456 * t580 + t553 * t600;
t449 = -pkin(1) * t455 + mrSges(2,3) * t571 - pkin(2) * t463 - Ifges(3,5) * t565 - Ifges(3,6) * t566 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t545 + mrSges(3,2) * t546 - t579 * t464 - t584 * t462 - pkin(3) * t590 - pkin(8) * t595 - Ifges(4,5) * t531 - Ifges(4,6) * t530 - Ifges(4,3) * t575 - mrSges(4,1) * t507 + mrSges(4,2) * t508 + t588 * Ifges(2,5) + mrSges(2,1) * g(3) - t557 * t537 + t556 * t538 + Ifges(2,6) * qJDD(1) + (-t554 * t581 + t555 * t586) * qJD(1);
t448 = -mrSges(3,1) * t558 + mrSges(3,3) * t546 + Ifges(3,4) * t565 + Ifges(3,2) * t566 + Ifges(3,6) * qJDD(2) - pkin(2) * t591 + pkin(7) * t596 + qJD(2) * t555 + t580 * t451 + t585 * t456 - t553 * t601;
t447 = -mrSges(2,2) * g(3) - mrSges(2,3) * t570 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t588 - pkin(6) * t455 - t448 * t581 + t450 * t586;
t1 = [-m(1) * g(1) + t598; -m(1) * g(2) + t602; (-m(1) - m(2)) * g(3) + t455; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t602 + t447 * t587 - t449 * t582; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t598 + t582 * t447 + t587 * t449; -mrSges(1,1) * g(2) + mrSges(2,1) * t570 + mrSges(1,2) * g(1) - mrSges(2,2) * t571 + Ifges(2,3) * qJDD(1) + pkin(1) * t589 + pkin(6) * t597 + t586 * t448 + t581 * t450;];
tauB = t1;

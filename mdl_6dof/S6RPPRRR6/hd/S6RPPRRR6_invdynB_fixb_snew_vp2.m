% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:57:19
% EndTime: 2019-05-05 15:57:25
% DurationCPUTime: 2.88s
% Computational Cost: add. (32676->289), mult. (61555->342), div. (0->0), fcn. (35363->8), ass. (0->108)
t570 = sin(qJ(1));
t574 = cos(qJ(1));
t546 = t570 * g(1) - t574 * g(2);
t576 = qJD(1) ^ 2;
t526 = -qJDD(1) * pkin(1) - t576 * qJ(2) + qJDD(2) - t546;
t520 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t526;
t547 = -t574 * g(1) - t570 * g(2);
t599 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t547;
t598 = -m(3) - m(4);
t597 = mrSges(2,1) - mrSges(3,2);
t525 = t576 * pkin(1) - t599;
t521 = qJDD(3) + (-pkin(1) - qJ(3)) * t576 + t599;
t518 = -qJDD(1) * pkin(7) + t521;
t569 = sin(qJ(4));
t573 = cos(qJ(4));
t511 = -t573 * g(3) + t569 * t518;
t540 = (mrSges(5,1) * t569 + mrSges(5,2) * t573) * qJD(1);
t594 = qJD(1) * qJD(4);
t550 = t573 * t594;
t542 = -t569 * qJDD(1) - t550;
t595 = qJD(1) * t573;
t545 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t595;
t552 = t569 * qJD(1);
t517 = -t576 * pkin(7) - t520;
t589 = t569 * t594;
t543 = t573 * qJDD(1) - t589;
t497 = (-t543 + t589) * pkin(8) + (-t542 + t550) * pkin(4) + t517;
t541 = (pkin(4) * t569 - pkin(8) * t573) * qJD(1);
t575 = qJD(4) ^ 2;
t500 = -t575 * pkin(4) + qJDD(4) * pkin(8) - t541 * t552 + t511;
t568 = sin(qJ(5));
t572 = cos(qJ(5));
t485 = t572 * t497 - t568 * t500;
t538 = t572 * qJD(4) - t568 * t595;
t509 = t538 * qJD(5) + t568 * qJDD(4) + t572 * t543;
t537 = qJDD(5) - t542;
t539 = t568 * qJD(4) + t572 * t595;
t549 = t552 + qJD(5);
t482 = (t538 * t549 - t509) * pkin(9) + (t538 * t539 + t537) * pkin(5) + t485;
t486 = t568 * t497 + t572 * t500;
t508 = -t539 * qJD(5) + t572 * qJDD(4) - t568 * t543;
t524 = t549 * pkin(5) - t539 * pkin(9);
t536 = t538 ^ 2;
t483 = -t536 * pkin(5) + t508 * pkin(9) - t549 * t524 + t486;
t567 = sin(qJ(6));
t571 = cos(qJ(6));
t480 = t571 * t482 - t567 * t483;
t512 = t571 * t538 - t567 * t539;
t489 = t512 * qJD(6) + t567 * t508 + t571 * t509;
t513 = t567 * t538 + t571 * t539;
t496 = -t512 * mrSges(7,1) + t513 * mrSges(7,2);
t548 = qJD(6) + t549;
t501 = -t548 * mrSges(7,2) + t512 * mrSges(7,3);
t531 = qJDD(6) + t537;
t478 = m(7) * t480 + t531 * mrSges(7,1) - t489 * mrSges(7,3) - t513 * t496 + t548 * t501;
t481 = t567 * t482 + t571 * t483;
t488 = -t513 * qJD(6) + t571 * t508 - t567 * t509;
t502 = t548 * mrSges(7,1) - t513 * mrSges(7,3);
t479 = m(7) * t481 - t531 * mrSges(7,2) + t488 * mrSges(7,3) + t512 * t496 - t548 * t502;
t470 = t571 * t478 + t567 * t479;
t515 = -t538 * mrSges(6,1) + t539 * mrSges(6,2);
t522 = -t549 * mrSges(6,2) + t538 * mrSges(6,3);
t468 = m(6) * t485 + t537 * mrSges(6,1) - t509 * mrSges(6,3) - t539 * t515 + t549 * t522 + t470;
t523 = t549 * mrSges(6,1) - t539 * mrSges(6,3);
t585 = -t567 * t478 + t571 * t479;
t469 = m(6) * t486 - t537 * mrSges(6,2) + t508 * mrSges(6,3) + t538 * t515 - t549 * t523 + t585;
t586 = -t568 * t468 + t572 * t469;
t463 = m(5) * t511 - qJDD(4) * mrSges(5,2) + t542 * mrSges(5,3) - qJD(4) * t545 - t540 * t552 + t586;
t510 = t569 * g(3) + t573 * t518;
t544 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t552;
t499 = -qJDD(4) * pkin(4) - t575 * pkin(8) + t541 * t595 - t510;
t484 = -t508 * pkin(5) - t536 * pkin(9) + t539 * t524 + t499;
t579 = m(7) * t484 - t488 * mrSges(7,1) + t489 * mrSges(7,2) - t512 * t501 + t513 * t502;
t577 = -m(6) * t499 + t508 * mrSges(6,1) - t509 * mrSges(6,2) + t538 * t522 - t539 * t523 - t579;
t474 = m(5) * t510 + qJDD(4) * mrSges(5,1) - t543 * mrSges(5,3) + qJD(4) * t544 - t540 * t595 + t577;
t456 = t569 * t463 + t573 * t474;
t584 = -m(4) * t521 - qJDD(1) * mrSges(4,2) - t456;
t580 = -m(3) * t525 + t576 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t584;
t454 = m(2) * t547 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t576 + t580;
t464 = t572 * t468 + t568 * t469;
t581 = -m(5) * t517 + t542 * mrSges(5,1) - t543 * mrSges(5,2) - t544 * t552 - t545 * t595 - t464;
t460 = m(4) * t520 - t576 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t581;
t578 = -m(3) * t526 + t576 * mrSges(3,3) - t460;
t459 = m(2) * t546 - t576 * mrSges(2,2) + t597 * qJDD(1) + t578;
t596 = t570 * t454 + t574 * t459;
t591 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t590 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t588 = t574 * t454 - t570 * t459;
t587 = t573 * t463 - t569 * t474;
t530 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t573 - Ifges(5,4) * t569) * qJD(1);
t529 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t573 - Ifges(5,2) * t569) * qJD(1);
t528 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t573 - Ifges(5,6) * t569) * qJD(1);
t505 = Ifges(6,1) * t539 + Ifges(6,4) * t538 + Ifges(6,5) * t549;
t504 = Ifges(6,4) * t539 + Ifges(6,2) * t538 + Ifges(6,6) * t549;
t503 = Ifges(6,5) * t539 + Ifges(6,6) * t538 + Ifges(6,3) * t549;
t492 = Ifges(7,1) * t513 + Ifges(7,4) * t512 + Ifges(7,5) * t548;
t491 = Ifges(7,4) * t513 + Ifges(7,2) * t512 + Ifges(7,6) * t548;
t490 = Ifges(7,5) * t513 + Ifges(7,6) * t512 + Ifges(7,3) * t548;
t472 = mrSges(7,2) * t484 - mrSges(7,3) * t480 + Ifges(7,1) * t489 + Ifges(7,4) * t488 + Ifges(7,5) * t531 + t512 * t490 - t548 * t491;
t471 = -mrSges(7,1) * t484 + mrSges(7,3) * t481 + Ifges(7,4) * t489 + Ifges(7,2) * t488 + Ifges(7,6) * t531 - t513 * t490 + t548 * t492;
t457 = mrSges(6,2) * t499 - mrSges(6,3) * t485 + Ifges(6,1) * t509 + Ifges(6,4) * t508 + Ifges(6,5) * t537 - pkin(9) * t470 - t567 * t471 + t571 * t472 + t538 * t503 - t549 * t504;
t455 = t598 * g(3) + t587;
t451 = -mrSges(6,1) * t499 + mrSges(6,3) * t486 + Ifges(6,4) * t509 + Ifges(6,2) * t508 + Ifges(6,6) * t537 - pkin(5) * t579 + pkin(9) * t585 + t571 * t471 + t567 * t472 - t539 * t503 + t549 * t505;
t450 = Ifges(5,4) * t543 + Ifges(5,2) * t542 + Ifges(5,6) * qJDD(4) - t528 * t595 + qJD(4) * t530 - mrSges(5,1) * t517 + mrSges(5,3) * t511 - Ifges(6,5) * t509 - Ifges(6,6) * t508 - Ifges(6,3) * t537 - t539 * t504 + t538 * t505 - mrSges(6,1) * t485 + mrSges(6,2) * t486 - Ifges(7,5) * t489 - Ifges(7,6) * t488 - Ifges(7,3) * t531 - t513 * t491 + t512 * t492 - mrSges(7,1) * t480 + mrSges(7,2) * t481 - pkin(5) * t470 - pkin(4) * t464;
t449 = mrSges(5,2) * t517 - mrSges(5,3) * t510 + Ifges(5,1) * t543 + Ifges(5,4) * t542 + Ifges(5,5) * qJDD(4) - pkin(8) * t464 - qJD(4) * t529 - t568 * t451 + t572 * t457 - t528 * t552;
t448 = -pkin(2) * t584 + pkin(8) * t586 + t568 * t457 - qJ(3) * t587 + t572 * t451 + mrSges(2,3) * t547 + pkin(4) * t577 + Ifges(5,6) * t542 + Ifges(5,5) * t543 - mrSges(3,1) * t525 + mrSges(5,1) * t510 - mrSges(5,2) * t511 + mrSges(4,1) * t521 + Ifges(5,3) * qJDD(4) + pkin(3) * t456 - pkin(1) * t455 + (t573 * t529 + t569 * t530) * qJD(1) + (-pkin(2) * mrSges(4,3) + t591) * t576 + t590 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t597) * g(3);
t447 = -qJ(2) * t455 - mrSges(2,3) * t546 + pkin(2) * t460 + mrSges(3,1) * t526 + pkin(3) * t581 + pkin(7) * t587 + t569 * t449 + t573 * t450 + mrSges(4,1) * t520 - t590 * t576 + t591 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t588; -m(1) * g(2) + t596; (-m(1) - m(2) + t598) * g(3) + t587; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t596 + t574 * t447 - t570 * t448; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t588 + t570 * t447 + t574 * t448; qJ(2) * (-t576 * mrSges(4,3) + t580) + pkin(1) * t578 + mrSges(2,1) * t546 - mrSges(2,2) * t547 - qJ(3) * t460 + mrSges(3,2) * t526 - mrSges(3,3) * t525 + t573 * t449 - t569 * t450 - pkin(7) * t456 + mrSges(4,2) * t521 - mrSges(4,3) * t520 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;

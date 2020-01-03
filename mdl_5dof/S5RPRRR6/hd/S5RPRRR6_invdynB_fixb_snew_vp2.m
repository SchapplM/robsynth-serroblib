% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:13
% EndTime: 2019-12-31 19:01:17
% DurationCPUTime: 3.86s
% Computational Cost: add. (46073->270), mult. (90089->343), div. (0->0), fcn. (57131->10), ass. (0->108)
t501 = sin(qJ(4));
t502 = sin(qJ(3));
t505 = cos(qJ(4));
t506 = cos(qJ(3));
t474 = (t501 * t502 - t505 * t506) * qJD(1);
t503 = sin(qJ(1));
t507 = cos(qJ(1));
t487 = t503 * g(1) - g(2) * t507;
t479 = qJDD(1) * pkin(1) + t487;
t488 = -g(1) * t507 - g(2) * t503;
t508 = qJD(1) ^ 2;
t481 = -pkin(1) * t508 + t488;
t498 = sin(pkin(9));
t499 = cos(pkin(9));
t462 = t498 * t479 + t499 * t481;
t458 = -pkin(2) * t508 + qJDD(1) * pkin(6) + t462;
t497 = -g(3) + qJDD(2);
t446 = -t502 * t458 + t506 * t497;
t521 = qJD(1) * qJD(3);
t519 = t506 * t521;
t482 = qJDD(1) * t502 + t519;
t436 = (-t482 + t519) * pkin(7) + (t502 * t506 * t508 + qJDD(3)) * pkin(3) + t446;
t447 = t506 * t458 + t502 * t497;
t483 = qJDD(1) * t506 - t502 * t521;
t523 = qJD(1) * t502;
t486 = qJD(3) * pkin(3) - pkin(7) * t523;
t496 = t506 ^ 2;
t437 = -pkin(3) * t496 * t508 + pkin(7) * t483 - qJD(3) * t486 + t447;
t430 = t501 * t436 + t505 * t437;
t475 = (t501 * t506 + t502 * t505) * qJD(1);
t448 = -qJD(4) * t475 - t482 * t501 + t483 * t505;
t459 = mrSges(5,1) * t474 + mrSges(5,2) * t475;
t495 = qJD(3) + qJD(4);
t466 = mrSges(5,1) * t495 - mrSges(5,3) * t475;
t494 = qJDD(3) + qJDD(4);
t460 = pkin(4) * t474 - pkin(8) * t475;
t493 = t495 ^ 2;
t427 = -pkin(4) * t493 + pkin(8) * t494 - t460 * t474 + t430;
t461 = t479 * t499 - t498 * t481;
t512 = -qJDD(1) * pkin(2) - t461;
t442 = -pkin(3) * t483 + t486 * t523 + (-pkin(7) * t496 - pkin(6)) * t508 + t512;
t449 = -qJD(4) * t474 + t482 * t505 + t483 * t501;
t428 = (t474 * t495 - t449) * pkin(8) + (t475 * t495 - t448) * pkin(4) + t442;
t500 = sin(qJ(5));
t504 = cos(qJ(5));
t424 = -t427 * t500 + t428 * t504;
t463 = -t475 * t500 + t495 * t504;
t433 = qJD(5) * t463 + t449 * t504 + t494 * t500;
t464 = t475 * t504 + t495 * t500;
t443 = -mrSges(6,1) * t463 + mrSges(6,2) * t464;
t445 = qJDD(5) - t448;
t467 = qJD(5) + t474;
t450 = -mrSges(6,2) * t467 + mrSges(6,3) * t463;
t422 = m(6) * t424 + mrSges(6,1) * t445 - mrSges(6,3) * t433 - t443 * t464 + t450 * t467;
t425 = t427 * t504 + t428 * t500;
t432 = -qJD(5) * t464 - t449 * t500 + t494 * t504;
t451 = mrSges(6,1) * t467 - mrSges(6,3) * t464;
t423 = m(6) * t425 - mrSges(6,2) * t445 + mrSges(6,3) * t432 + t443 * t463 - t451 * t467;
t514 = -t422 * t500 + t504 * t423;
t413 = m(5) * t430 - mrSges(5,2) * t494 + mrSges(5,3) * t448 - t459 * t474 - t466 * t495 + t514;
t429 = t436 * t505 - t437 * t501;
t465 = -mrSges(5,2) * t495 - mrSges(5,3) * t474;
t426 = -pkin(4) * t494 - pkin(8) * t493 + t460 * t475 - t429;
t511 = -m(6) * t426 + t432 * mrSges(6,1) - mrSges(6,2) * t433 + t463 * t450 - t451 * t464;
t418 = m(5) * t429 + mrSges(5,1) * t494 - mrSges(5,3) * t449 - t459 * t475 + t465 * t495 + t511;
t408 = t501 * t413 + t505 * t418;
t480 = (-mrSges(4,1) * t506 + mrSges(4,2) * t502) * qJD(1);
t522 = qJD(1) * t506;
t485 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t522;
t406 = m(4) * t446 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t482 + qJD(3) * t485 - t480 * t523 + t408;
t484 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t523;
t515 = t505 * t413 - t418 * t501;
t407 = m(4) * t447 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t483 - qJD(3) * t484 + t480 * t522 + t515;
t516 = -t406 * t502 + t506 * t407;
t399 = m(3) * t462 - mrSges(3,1) * t508 - qJDD(1) * mrSges(3,2) + t516;
t457 = -pkin(6) * t508 + t512;
t414 = t504 * t422 + t500 * t423;
t510 = m(5) * t442 - t448 * mrSges(5,1) + mrSges(5,2) * t449 + t474 * t465 + t466 * t475 + t414;
t509 = -m(4) * t457 + t483 * mrSges(4,1) - mrSges(4,2) * t482 - t484 * t523 + t485 * t522 - t510;
t410 = m(3) * t461 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t508 + t509;
t396 = t498 * t399 + t499 * t410;
t394 = m(2) * t487 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t508 + t396;
t517 = t499 * t399 - t498 * t410;
t395 = m(2) * t488 - mrSges(2,1) * t508 - qJDD(1) * mrSges(2,2) + t517;
t524 = t507 * t394 + t503 * t395;
t400 = t506 * t406 + t502 * t407;
t520 = m(3) * t497 + t400;
t518 = -t394 * t503 + t507 * t395;
t473 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t502 + Ifges(4,4) * t506) * qJD(1);
t472 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t502 + Ifges(4,2) * t506) * qJD(1);
t471 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t502 + Ifges(4,6) * t506) * qJD(1);
t455 = Ifges(5,1) * t475 - Ifges(5,4) * t474 + Ifges(5,5) * t495;
t454 = Ifges(5,4) * t475 - Ifges(5,2) * t474 + Ifges(5,6) * t495;
t453 = Ifges(5,5) * t475 - Ifges(5,6) * t474 + Ifges(5,3) * t495;
t440 = Ifges(6,1) * t464 + Ifges(6,4) * t463 + Ifges(6,5) * t467;
t439 = Ifges(6,4) * t464 + Ifges(6,2) * t463 + Ifges(6,6) * t467;
t438 = Ifges(6,5) * t464 + Ifges(6,6) * t463 + Ifges(6,3) * t467;
t416 = mrSges(6,2) * t426 - mrSges(6,3) * t424 + Ifges(6,1) * t433 + Ifges(6,4) * t432 + Ifges(6,5) * t445 + t438 * t463 - t439 * t467;
t415 = -mrSges(6,1) * t426 + mrSges(6,3) * t425 + Ifges(6,4) * t433 + Ifges(6,2) * t432 + Ifges(6,6) * t445 - t438 * t464 + t440 * t467;
t402 = -mrSges(5,1) * t442 - mrSges(6,1) * t424 + mrSges(6,2) * t425 + mrSges(5,3) * t430 + Ifges(5,4) * t449 - Ifges(6,5) * t433 + Ifges(5,2) * t448 + Ifges(5,6) * t494 - Ifges(6,6) * t432 - Ifges(6,3) * t445 - pkin(4) * t414 - t439 * t464 + t440 * t463 - t453 * t475 + t455 * t495;
t401 = mrSges(5,2) * t442 - mrSges(5,3) * t429 + Ifges(5,1) * t449 + Ifges(5,4) * t448 + Ifges(5,5) * t494 - pkin(8) * t414 - t415 * t500 + t416 * t504 - t453 * t474 - t454 * t495;
t390 = mrSges(4,2) * t457 - mrSges(4,3) * t446 + Ifges(4,1) * t482 + Ifges(4,4) * t483 + Ifges(4,5) * qJDD(3) - pkin(7) * t408 - qJD(3) * t472 + t401 * t505 - t402 * t501 + t471 * t522;
t389 = -mrSges(4,1) * t457 + mrSges(4,3) * t447 + Ifges(4,4) * t482 + Ifges(4,2) * t483 + Ifges(4,6) * qJDD(3) - pkin(3) * t510 + pkin(7) * t515 + qJD(3) * t473 + t501 * t401 + t505 * t402 - t471 * t523;
t388 = -pkin(2) * t400 - mrSges(3,1) * t497 + mrSges(3,3) * t462 - pkin(3) * t408 - Ifges(4,5) * t482 - Ifges(4,6) * t483 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t446 + mrSges(4,2) * t447 - Ifges(5,5) * t449 - Ifges(5,6) * t448 - Ifges(5,3) * t494 - mrSges(5,1) * t429 + mrSges(5,2) * t430 - t500 * t416 - t504 * t415 - pkin(4) * t511 - pkin(8) * t514 - t475 * t454 - t474 * t455 + t508 * Ifges(3,5) + Ifges(3,6) * qJDD(1) + (-t472 * t502 + t473 * t506) * qJD(1);
t387 = mrSges(3,2) * t497 - mrSges(3,3) * t461 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t508 - pkin(6) * t400 - t389 * t502 + t390 * t506;
t386 = -mrSges(2,2) * g(3) - mrSges(2,3) * t487 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t508 - qJ(2) * t396 + t387 * t499 - t388 * t498;
t385 = mrSges(2,1) * g(3) + mrSges(2,3) * t488 + t508 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t520 + qJ(2) * t517 + t498 * t387 + t499 * t388;
t1 = [-m(1) * g(1) + t518; -m(1) * g(2) + t524; (-m(1) - m(2)) * g(3) + t520; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t524 - t503 * t385 + t507 * t386; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t518 + t507 * t385 + t503 * t386; pkin(1) * t396 + mrSges(2,1) * t487 - mrSges(2,2) * t488 + t502 * t390 + t506 * t389 + pkin(2) * t509 + pkin(6) * t516 + mrSges(3,1) * t461 - mrSges(3,2) * t462 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;

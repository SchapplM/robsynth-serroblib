% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:28:19
% EndTime: 2019-05-05 09:28:23
% DurationCPUTime: 2.26s
% Computational Cost: add. (20216->269), mult. (39514->333), div. (0->0), fcn. (28211->12), ass. (0->115)
t523 = Ifges(6,1) + Ifges(7,1);
t517 = Ifges(6,4) + Ifges(7,4);
t516 = Ifges(6,5) + Ifges(7,5);
t522 = Ifges(6,2) + Ifges(7,2);
t515 = Ifges(6,6) + Ifges(7,6);
t521 = Ifges(6,3) + Ifges(7,3);
t477 = sin(pkin(11));
t479 = cos(pkin(11));
t465 = g(1) * t477 - g(2) * t479;
t476 = -g(3) + qJDD(1);
t478 = sin(pkin(6));
t480 = cos(pkin(6));
t520 = t465 * t480 + t476 * t478;
t482 = sin(qJ(4));
t483 = sin(qJ(3));
t486 = cos(qJ(4));
t487 = cos(qJ(3));
t455 = (t482 * t483 - t486 * t487) * qJD(2);
t466 = -g(1) * t479 - g(2) * t477;
t484 = sin(qJ(2));
t488 = cos(qJ(2));
t437 = -t484 * t466 + t488 * t520;
t505 = qJD(2) * qJD(3);
t502 = t487 * t505;
t463 = qJDD(2) * t483 + t502;
t464 = qJDD(2) * t487 - t483 * t505;
t426 = -qJD(4) * t455 + t463 * t486 + t464 * t482;
t456 = (t482 * t487 + t483 * t486) * qJD(2);
t474 = qJD(3) + qJD(4);
t481 = sin(qJ(5));
t485 = cos(qJ(5));
t443 = -t456 * t481 + t474 * t485;
t473 = qJDD(3) + qJDD(4);
t401 = qJD(5) * t443 + t426 * t485 + t473 * t481;
t444 = t456 * t485 + t474 * t481;
t415 = -mrSges(7,1) * t443 + mrSges(7,2) * t444;
t438 = t488 * t466 + t484 * t520;
t489 = qJD(2) ^ 2;
t433 = -pkin(2) * t489 + qJDD(2) * pkin(8) + t438;
t449 = -t465 * t478 + t476 * t480;
t412 = -t483 * t433 + t487 * t449;
t395 = (-t463 + t502) * pkin(9) + (t483 * t487 * t489 + qJDD(3)) * pkin(3) + t412;
t413 = t487 * t433 + t483 * t449;
t507 = qJD(2) * t483;
t470 = qJD(3) * pkin(3) - pkin(9) * t507;
t475 = t487 ^ 2;
t396 = -pkin(3) * t475 * t489 + pkin(9) * t464 - qJD(3) * t470 + t413;
t391 = t482 * t395 + t486 * t396;
t441 = pkin(4) * t455 - pkin(10) * t456;
t472 = t474 ^ 2;
t385 = -pkin(4) * t472 + pkin(10) * t473 - t441 * t455 + t391;
t495 = -qJDD(2) * pkin(2) - t437;
t403 = -t464 * pkin(3) + t470 * t507 + (-pkin(9) * t475 - pkin(8)) * t489 + t495;
t425 = -qJD(4) * t456 - t463 * t482 + t464 * t486;
t388 = (t455 * t474 - t426) * pkin(10) + (t456 * t474 - t425) * pkin(4) + t403;
t380 = -t481 * t385 + t485 * t388;
t423 = qJDD(5) - t425;
t450 = qJD(5) + t455;
t377 = -0.2e1 * qJD(6) * t444 + (t443 * t450 - t401) * qJ(6) + (t443 * t444 + t423) * pkin(5) + t380;
t427 = -mrSges(7,2) * t450 + mrSges(7,3) * t443;
t504 = m(7) * t377 + t423 * mrSges(7,1) + t450 * t427;
t374 = -t401 * mrSges(7,3) - t444 * t415 + t504;
t381 = t485 * t385 + t481 * t388;
t400 = -qJD(5) * t444 - t426 * t481 + t473 * t485;
t429 = pkin(5) * t450 - qJ(6) * t444;
t442 = t443 ^ 2;
t379 = -pkin(5) * t442 + qJ(6) * t400 + 0.2e1 * qJD(6) * t443 - t429 * t450 + t381;
t509 = t517 * t443 + t444 * t523 + t516 * t450;
t510 = -t443 * t522 - t444 * t517 - t450 * t515;
t519 = mrSges(6,1) * t380 + mrSges(7,1) * t377 - mrSges(6,2) * t381 - mrSges(7,2) * t379 + pkin(5) * t374 + t400 * t515 + t401 * t516 + t423 * t521 - t443 * t509 - t444 * t510;
t518 = -mrSges(6,2) - mrSges(7,2);
t440 = mrSges(5,1) * t455 + mrSges(5,2) * t456;
t448 = mrSges(5,1) * t474 - mrSges(5,3) * t456;
t416 = -mrSges(6,1) * t443 + mrSges(6,2) * t444;
t428 = -mrSges(6,2) * t450 + mrSges(6,3) * t443;
t367 = m(6) * t380 + t423 * mrSges(6,1) + t450 * t428 + (-t415 - t416) * t444 + (-mrSges(6,3) - mrSges(7,3)) * t401 + t504;
t503 = m(7) * t379 + t400 * mrSges(7,3) + t443 * t415;
t430 = mrSges(7,1) * t450 - mrSges(7,3) * t444;
t508 = -mrSges(6,1) * t450 + mrSges(6,3) * t444 - t430;
t370 = m(6) * t381 + t400 * mrSges(6,3) + t443 * t416 + t423 * t518 + t508 * t450 + t503;
t499 = -t367 * t481 + t485 * t370;
t361 = m(5) * t391 - mrSges(5,2) * t473 + mrSges(5,3) * t425 - t440 * t455 - t448 * t474 + t499;
t390 = t395 * t486 - t482 * t396;
t447 = -mrSges(5,2) * t474 - mrSges(5,3) * t455;
t384 = -pkin(4) * t473 - pkin(10) * t472 + t456 * t441 - t390;
t382 = -pkin(5) * t400 - qJ(6) * t442 + t429 * t444 + qJDD(6) + t384;
t498 = -m(7) * t382 + t400 * mrSges(7,1) + t443 * t427;
t491 = -m(6) * t384 + t400 * mrSges(6,1) + t401 * t518 + t443 * t428 + t508 * t444 + t498;
t372 = m(5) * t390 + t473 * mrSges(5,1) - t426 * mrSges(5,3) - t456 * t440 + t474 * t447 + t491;
t356 = t482 * t361 + t486 * t372;
t365 = t485 * t367 + t481 * t370;
t511 = -t443 * t515 - t444 * t516 - t450 * t521;
t506 = qJD(2) * t487;
t462 = (-mrSges(4,1) * t487 + mrSges(4,2) * t483) * qJD(2);
t468 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t506;
t354 = m(4) * t412 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t463 + qJD(3) * t468 - t462 * t507 + t356;
t467 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t507;
t500 = t486 * t361 - t372 * t482;
t355 = m(4) * t413 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t464 - qJD(3) * t467 + t462 * t506 + t500;
t501 = -t354 * t483 + t487 * t355;
t494 = m(5) * t403 - t425 * mrSges(5,1) + mrSges(5,2) * t426 + t455 * t447 + t448 * t456 + t365;
t375 = t401 * mrSges(7,2) + t444 * t430 - t498;
t358 = -mrSges(6,1) * t384 + mrSges(6,3) * t381 - mrSges(7,1) * t382 + mrSges(7,3) * t379 - pkin(5) * t375 + qJ(6) * t503 + (-qJ(6) * t430 + t509) * t450 + t511 * t444 + (-mrSges(7,2) * qJ(6) + t515) * t423 + t517 * t401 + t522 * t400;
t363 = mrSges(6,2) * t384 + mrSges(7,2) * t382 - mrSges(6,3) * t380 - mrSges(7,3) * t377 - qJ(6) * t374 + t517 * t400 + t401 * t523 + t516 * t423 - t511 * t443 + t510 * t450;
t435 = Ifges(5,4) * t456 - Ifges(5,2) * t455 + Ifges(5,6) * t474;
t436 = Ifges(5,1) * t456 - Ifges(5,4) * t455 + Ifges(5,5) * t474;
t492 = mrSges(5,1) * t390 - mrSges(5,2) * t391 + Ifges(5,5) * t426 + Ifges(5,6) * t425 + Ifges(5,3) * t473 + pkin(4) * t491 + pkin(10) * t499 + t485 * t358 + t481 * t363 + t456 * t435 + t455 * t436;
t432 = -t489 * pkin(8) + t495;
t490 = -m(4) * t432 + t464 * mrSges(4,1) - mrSges(4,2) * t463 - t467 * t507 + t468 * t506 - t494;
t454 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t483 + Ifges(4,4) * t487) * qJD(2);
t453 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t483 + Ifges(4,2) * t487) * qJD(2);
t434 = Ifges(5,5) * t456 - Ifges(5,6) * t455 + Ifges(5,3) * t474;
t352 = -mrSges(5,1) * t403 + mrSges(5,3) * t391 + Ifges(5,4) * t426 + Ifges(5,2) * t425 + Ifges(5,6) * t473 - pkin(4) * t365 - t456 * t434 + t474 * t436 - t519;
t351 = mrSges(5,2) * t403 - mrSges(5,3) * t390 + Ifges(5,1) * t426 + Ifges(5,4) * t425 + Ifges(5,5) * t473 - pkin(10) * t365 - t358 * t481 + t363 * t485 - t434 * t455 - t435 * t474;
t1 = [m(2) * t476 + t480 * (m(3) * t449 + t354 * t487 + t355 * t483) + (t484 * (m(3) * t438 - mrSges(3,1) * t489 - qJDD(2) * mrSges(3,2) + t501) + t488 * (m(3) * t437 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t489 + t490)) * t478; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t437 - mrSges(3,2) * t438 + t483 * (mrSges(4,2) * t432 - mrSges(4,3) * t412 + Ifges(4,1) * t463 + Ifges(4,4) * t464 + Ifges(4,5) * qJDD(3) - pkin(9) * t356 - qJD(3) * t453 + t351 * t486 - t352 * t482) + t487 * (-mrSges(4,1) * t432 + mrSges(4,3) * t413 + Ifges(4,4) * t463 + Ifges(4,2) * t464 + Ifges(4,6) * qJDD(3) - pkin(3) * t494 + pkin(9) * t500 + qJD(3) * t454 + t482 * t351 + t486 * t352) + pkin(2) * t490 + pkin(8) * t501; Ifges(4,3) * qJDD(3) + (t453 * t483 - t454 * t487) * qJD(2) + t492 + Ifges(4,5) * t463 + Ifges(4,6) * t464 + mrSges(4,1) * t412 - mrSges(4,2) * t413 + pkin(3) * t356; t492; t519; t375;];
tauJ  = t1;

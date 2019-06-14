% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:43:59
% EndTime: 2019-05-05 03:44:01
% DurationCPUTime: 2.07s
% Computational Cost: add. (15912->266), mult. (33985->333), div. (0->0), fcn. (23876->12), ass. (0->112)
t520 = -2 * qJD(4);
t519 = Ifges(6,1) + Ifges(7,1);
t512 = Ifges(6,4) - Ifges(7,5);
t511 = -Ifges(6,5) - Ifges(7,4);
t518 = Ifges(6,2) + Ifges(7,3);
t510 = Ifges(6,6) - Ifges(7,6);
t517 = -Ifges(6,3) - Ifges(7,2);
t474 = sin(pkin(10));
t477 = cos(pkin(10));
t463 = t474 * g(1) - t477 * g(2);
t472 = -g(3) + qJDD(1);
t475 = sin(pkin(6));
t478 = cos(pkin(6));
t516 = t463 * t478 + t472 * t475;
t464 = -t477 * g(1) - t474 * g(2);
t481 = sin(qJ(2));
t483 = cos(qJ(2));
t428 = t483 * t464 + t516 * t481;
t485 = qJD(2) ^ 2;
t423 = -t485 * pkin(2) + qJDD(2) * pkin(8) + t428;
t444 = -t475 * t463 + t478 * t472;
t480 = sin(qJ(3));
t482 = cos(qJ(3));
t397 = -t480 * t423 + t482 * t444;
t499 = qJD(2) * qJD(3);
t497 = t482 * t499;
t461 = t480 * qJDD(2) + t497;
t394 = (-t461 + t497) * qJ(4) + (t480 * t482 * t485 + qJDD(3)) * pkin(3) + t397;
t398 = t482 * t423 + t480 * t444;
t462 = t482 * qJDD(2) - t480 * t499;
t502 = qJD(2) * t480;
t465 = qJD(3) * pkin(3) - qJ(4) * t502;
t471 = t482 ^ 2;
t395 = -t471 * t485 * pkin(3) + t462 * qJ(4) - qJD(3) * t465 + t398;
t473 = sin(pkin(11));
t476 = cos(pkin(11));
t450 = (t482 * t473 + t480 * t476) * qJD(2);
t387 = t476 * t394 - t473 * t395 + t450 * t520;
t449 = (t480 * t473 - t482 * t476) * qJD(2);
t427 = -t481 * t464 + t483 * t516;
t437 = t476 * t461 + t473 * t462;
t479 = sin(qJ(5));
t514 = cos(qJ(5));
t438 = -qJD(3) * t514 + t479 * t450;
t410 = -t438 * qJD(5) + t479 * qJDD(3) + t437 * t514;
t439 = t479 * qJD(3) + t450 * t514;
t414 = t438 * mrSges(7,1) - t439 * mrSges(7,3);
t388 = t473 * t394 + t476 * t395 + t449 * t520;
t431 = t449 * pkin(4) - t450 * pkin(9);
t484 = qJD(3) ^ 2;
t386 = -t484 * pkin(4) + qJDD(3) * pkin(9) - t449 * t431 + t388;
t489 = -qJDD(2) * pkin(2) - t427;
t396 = -t462 * pkin(3) + qJDD(4) + t465 * t502 + (-qJ(4) * t471 - pkin(8)) * t485 + t489;
t436 = -t473 * t461 + t476 * t462;
t390 = (qJD(3) * t449 - t437) * pkin(9) + (qJD(3) * t450 - t436) * pkin(4) + t396;
t382 = -t479 * t386 + t390 * t514;
t413 = t438 * pkin(5) - t439 * qJ(6);
t435 = qJDD(5) - t436;
t448 = qJD(5) + t449;
t447 = t448 ^ 2;
t380 = -t435 * pkin(5) - t447 * qJ(6) + t439 * t413 + qJDD(6) - t382;
t418 = -t438 * mrSges(7,2) + t448 * mrSges(7,3);
t492 = -m(7) * t380 + t435 * mrSges(7,1) + t448 * t418;
t376 = t410 * mrSges(7,2) + t439 * t414 - t492;
t383 = t386 * t514 + t479 * t390;
t379 = -t447 * pkin(5) + t435 * qJ(6) + 0.2e1 * qJD(6) * t448 - t438 * t413 + t383;
t409 = t439 * qJD(5) - qJDD(3) * t514 + t479 * t437;
t421 = -t448 * mrSges(7,1) + t439 * mrSges(7,2);
t498 = m(7) * t379 + t435 * mrSges(7,3) + t448 * t421;
t504 = t512 * t438 - t519 * t439 + t511 * t448;
t505 = t518 * t438 - t512 * t439 - t510 * t448;
t515 = -t409 * t510 - t410 * t511 - t517 * t435 - t438 * t504 - t439 * t505 + mrSges(6,1) * t382 - mrSges(7,1) * t380 - mrSges(6,2) * t383 + mrSges(7,3) * t379 - pkin(5) * t376 + qJ(6) * (-t409 * mrSges(7,2) - t438 * t414 + t498);
t513 = -mrSges(6,3) - mrSges(7,2);
t430 = t449 * mrSges(5,1) + t450 * mrSges(5,2);
t443 = qJD(3) * mrSges(5,1) - t450 * mrSges(5,3);
t420 = t448 * mrSges(6,1) - t439 * mrSges(6,3);
t503 = -t438 * mrSges(6,1) - t439 * mrSges(6,2) - t414;
t372 = m(6) * t383 - t435 * mrSges(6,2) + t409 * t513 - t448 * t420 + t438 * t503 + t498;
t419 = -t448 * mrSges(6,2) - t438 * mrSges(6,3);
t374 = m(6) * t382 + t435 * mrSges(6,1) + t410 * t513 + t448 * t419 + t439 * t503 + t492;
t494 = t372 * t514 - t479 * t374;
t363 = m(5) * t388 - qJDD(3) * mrSges(5,2) + t436 * mrSges(5,3) - qJD(3) * t443 - t449 * t430 + t494;
t442 = -qJD(3) * mrSges(5,2) - t449 * mrSges(5,3);
t385 = -qJDD(3) * pkin(4) - t484 * pkin(9) + t450 * t431 - t387;
t381 = -0.2e1 * qJD(6) * t439 + (t438 * t448 - t410) * qJ(6) + (t439 * t448 + t409) * pkin(5) + t385;
t377 = m(7) * t381 + t409 * mrSges(7,1) - t410 * mrSges(7,3) + t438 * t418 - t439 * t421;
t487 = -m(6) * t385 - t409 * mrSges(6,1) - t410 * mrSges(6,2) - t438 * t419 - t439 * t420 - t377;
t369 = m(5) * t387 + qJDD(3) * mrSges(5,1) - t437 * mrSges(5,3) + qJD(3) * t442 - t450 * t430 + t487;
t360 = t473 * t363 + t476 * t369;
t367 = t479 * t372 + t374 * t514;
t506 = t510 * t438 + t511 * t439 + t517 * t448;
t501 = qJD(2) * t482;
t460 = (-t482 * mrSges(4,1) + t480 * mrSges(4,2)) * qJD(2);
t467 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t501;
t358 = m(4) * t397 + qJDD(3) * mrSges(4,1) - t461 * mrSges(4,3) + qJD(3) * t467 - t460 * t502 + t360;
t466 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t502;
t495 = t476 * t363 - t473 * t369;
t359 = m(4) * t398 - qJDD(3) * mrSges(4,2) + t462 * mrSges(4,3) - qJD(3) * t466 + t460 * t501 + t495;
t496 = -t480 * t358 + t482 * t359;
t366 = m(5) * t396 - t436 * mrSges(5,1) + t437 * mrSges(5,2) + t449 * t442 + t450 * t443 + t367;
t422 = -t485 * pkin(8) + t489;
t486 = -m(4) * t422 + t462 * mrSges(4,1) - t461 * mrSges(4,2) - t466 * t502 + t467 * t501 - t366;
t454 = Ifges(4,5) * qJD(3) + (t480 * Ifges(4,1) + t482 * Ifges(4,4)) * qJD(2);
t453 = Ifges(4,6) * qJD(3) + (t480 * Ifges(4,4) + t482 * Ifges(4,2)) * qJD(2);
t426 = Ifges(5,1) * t450 - Ifges(5,4) * t449 + Ifges(5,5) * qJD(3);
t425 = Ifges(5,4) * t450 - Ifges(5,2) * t449 + Ifges(5,6) * qJD(3);
t424 = Ifges(5,5) * t450 - Ifges(5,6) * t449 + Ifges(5,3) * qJD(3);
t365 = mrSges(6,2) * t385 + mrSges(7,2) * t380 - mrSges(6,3) * t382 - mrSges(7,3) * t381 - qJ(6) * t377 - t512 * t409 + t519 * t410 - t511 * t435 + t506 * t438 + t505 * t448;
t364 = -mrSges(6,1) * t385 - mrSges(7,1) * t381 + mrSges(7,2) * t379 + mrSges(6,3) * t383 - pkin(5) * t377 - t518 * t409 + t512 * t410 + t510 * t435 + t506 * t439 - t504 * t448;
t356 = -mrSges(5,1) * t396 + mrSges(5,3) * t388 + Ifges(5,4) * t437 + Ifges(5,2) * t436 + Ifges(5,6) * qJDD(3) - pkin(4) * t367 + qJD(3) * t426 - t450 * t424 - t515;
t355 = mrSges(5,2) * t396 - mrSges(5,3) * t387 + Ifges(5,1) * t437 + Ifges(5,4) * t436 + Ifges(5,5) * qJDD(3) - pkin(9) * t367 - qJD(3) * t425 - t479 * t364 + t365 * t514 - t449 * t424;
t1 = [m(2) * t472 + t478 * (m(3) * t444 + t482 * t358 + t480 * t359) + (t481 * (m(3) * t428 - t485 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t496) + t483 * (m(3) * t427 + qJDD(2) * mrSges(3,1) - t485 * mrSges(3,2) + t486)) * t475; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t427 - mrSges(3,2) * t428 + t480 * (mrSges(4,2) * t422 - mrSges(4,3) * t397 + Ifges(4,1) * t461 + Ifges(4,4) * t462 + Ifges(4,5) * qJDD(3) - qJ(4) * t360 - qJD(3) * t453 + t476 * t355 - t473 * t356) + t482 * (-mrSges(4,1) * t422 + mrSges(4,3) * t398 + Ifges(4,4) * t461 + Ifges(4,2) * t462 + Ifges(4,6) * qJDD(3) - pkin(3) * t366 + qJ(4) * t495 + qJD(3) * t454 + t473 * t355 + t476 * t356) + pkin(2) * t486 + pkin(8) * t496; Ifges(4,5) * t461 + Ifges(4,6) * t462 + mrSges(4,1) * t397 - mrSges(4,2) * t398 + Ifges(5,5) * t437 + Ifges(5,6) * t436 + t450 * t425 + t449 * t426 + mrSges(5,1) * t387 - mrSges(5,2) * t388 + t479 * t365 + t514 * t364 + pkin(4) * t487 + pkin(9) * t494 + pkin(3) * t360 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t480 * t453 - t482 * t454) * qJD(2); t366; t515; t376;];
tauJ  = t1;

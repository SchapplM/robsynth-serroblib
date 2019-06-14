% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:42:05
% EndTime: 2019-05-04 23:42:06
% DurationCPUTime: 1.83s
% Computational Cost: add. (14181->241), mult. (31318->295), div. (0->0), fcn. (23235->12), ass. (0->114)
t520 = Ifges(6,1) + Ifges(7,1);
t512 = Ifges(6,4) + Ifges(7,4);
t511 = Ifges(6,5) + Ifges(7,5);
t519 = Ifges(6,2) + Ifges(7,2);
t510 = Ifges(6,6) + Ifges(7,6);
t518 = Ifges(6,3) + Ifges(7,3);
t467 = sin(pkin(10));
t470 = cos(pkin(10));
t455 = g(1) * t467 - g(2) * t470;
t465 = -g(3) + qJDD(1);
t468 = sin(pkin(6));
t471 = cos(pkin(6));
t517 = t455 * t471 + t465 * t468;
t479 = qJD(2) ^ 2;
t466 = sin(pkin(11));
t469 = cos(pkin(11));
t473 = sin(qJ(4));
t476 = cos(qJ(4));
t486 = t466 * t473 - t469 * t476;
t448 = t486 * qJD(2);
t456 = -g(1) * t470 - g(2) * t467;
t474 = sin(qJ(2));
t477 = cos(qJ(2));
t429 = -t474 * t456 + t477 * t517;
t487 = t466 * t476 + t469 * t473;
t449 = t487 * qJD(2);
t497 = t449 * qJD(4);
t437 = -t486 * qJDD(2) - t497;
t498 = t448 * qJD(4);
t438 = t487 * qJDD(2) - t498;
t472 = sin(qJ(5));
t475 = cos(qJ(5));
t440 = qJD(4) * t475 - t449 * t472;
t412 = qJD(5) * t440 + qJDD(4) * t472 + t438 * t475;
t441 = qJD(4) * t472 + t449 * t475;
t415 = -mrSges(7,1) * t440 + mrSges(7,2) * t441;
t430 = t477 * t456 + t474 * t517;
t425 = -pkin(2) * t479 + qJDD(2) * qJ(3) + t430;
t446 = -t455 * t468 + t465 * t471;
t496 = qJD(2) * qJD(3);
t499 = t469 * t446 - 0.2e1 * t466 * t496;
t514 = pkin(3) * t469;
t395 = (-pkin(8) * qJDD(2) + t479 * t514 - t425) * t466 + t499;
t515 = 0.2e1 * t469;
t398 = t469 * t425 + t466 * t446 + t496 * t515;
t495 = qJDD(2) * t469;
t464 = t469 ^ 2;
t506 = t464 * t479;
t396 = -pkin(3) * t506 + pkin(8) * t495 + t398;
t388 = t473 * t395 + t476 * t396;
t436 = pkin(4) * t448 - pkin(9) * t449;
t478 = qJD(4) ^ 2;
t386 = -pkin(4) * t478 + qJDD(4) * pkin(9) - t436 * t448 + t388;
t463 = t466 ^ 2;
t484 = qJDD(3) - t429;
t413 = (-pkin(2) - t514) * qJDD(2) + (-qJ(3) + (-t463 - t464) * pkin(8)) * t479 + t484;
t391 = (-t438 + t498) * pkin(9) + (-t437 + t497) * pkin(4) + t413;
t381 = -t386 * t472 + t475 * t391;
t435 = qJDD(5) - t437;
t447 = qJD(5) + t448;
t378 = -0.2e1 * qJD(6) * t441 + (t440 * t447 - t412) * qJ(6) + (t440 * t441 + t435) * pkin(5) + t381;
t420 = -mrSges(7,2) * t447 + mrSges(7,3) * t440;
t494 = m(7) * t378 + t435 * mrSges(7,1) + t447 * t420;
t375 = -mrSges(7,3) * t412 - t415 * t441 + t494;
t382 = t475 * t386 + t472 * t391;
t411 = -qJD(5) * t441 + qJDD(4) * t475 - t438 * t472;
t422 = pkin(5) * t447 - qJ(6) * t441;
t439 = t440 ^ 2;
t380 = -pkin(5) * t439 + qJ(6) * t411 + 0.2e1 * qJD(6) * t440 - t422 * t447 + t382;
t501 = t512 * t440 + t441 * t520 + t511 * t447;
t502 = -t440 * t519 - t441 * t512 - t447 * t510;
t516 = mrSges(6,1) * t381 + mrSges(7,1) * t378 - mrSges(6,2) * t382 - mrSges(7,2) * t380 + pkin(5) * t375 + t510 * t411 + t511 * t412 + t435 * t518 - t501 * t440 - t502 * t441;
t513 = -mrSges(6,2) - mrSges(7,2);
t508 = mrSges(4,2) * t466;
t433 = mrSges(5,1) * t448 + mrSges(5,2) * t449;
t445 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t449;
t416 = -mrSges(6,1) * t440 + mrSges(6,2) * t441;
t421 = -mrSges(6,2) * t447 + mrSges(6,3) * t440;
t369 = m(6) * t381 + mrSges(6,1) * t435 + t421 * t447 + (-t415 - t416) * t441 + (-mrSges(6,3) - mrSges(7,3)) * t412 + t494;
t493 = m(7) * t380 + t411 * mrSges(7,3) + t440 * t415;
t423 = mrSges(7,1) * t447 - mrSges(7,3) * t441;
t500 = -mrSges(6,1) * t447 + mrSges(6,3) * t441 - t423;
t373 = m(6) * t382 + mrSges(6,3) * t411 + t416 * t440 + t513 * t435 + t500 * t447 + t493;
t490 = -t369 * t472 + t475 * t373;
t364 = m(5) * t388 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t437 - qJD(4) * t445 - t433 * t448 + t490;
t387 = t395 * t476 - t473 * t396;
t444 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t448;
t385 = -qJDD(4) * pkin(4) - pkin(9) * t478 + t449 * t436 - t387;
t383 = -pkin(5) * t411 - qJ(6) * t439 + t422 * t441 + qJDD(6) + t385;
t489 = -m(7) * t383 + t411 * mrSges(7,1) + t440 * t420;
t481 = -m(6) * t385 + t411 * mrSges(6,1) + t513 * t412 + t440 * t421 + t500 * t441 + t489;
t374 = m(5) * t387 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t438 + qJD(4) * t444 - t433 * t449 + t481;
t504 = t473 * t364 + t476 * t374;
t367 = t475 * t369 + t472 * t373;
t503 = -t440 * t510 - t441 * t511 - t447 * t518;
t397 = -t425 * t466 + t499;
t485 = mrSges(4,3) * qJDD(2) + t479 * (-mrSges(4,1) * t469 + t508);
t359 = m(4) * t397 - t485 * t466 + t504;
t491 = t476 * t364 - t374 * t473;
t360 = m(4) * t398 + t485 * t469 + t491;
t492 = -t359 * t466 + t469 * t360;
t483 = m(5) * t413 - t437 * mrSges(5,1) + mrSges(5,2) * t438 + t448 * t444 + t445 * t449 + t367;
t419 = -qJDD(2) * pkin(2) - qJ(3) * t479 + t484;
t480 = -m(4) * t419 + mrSges(4,1) * t495 - t483 + (t463 * t479 + t506) * mrSges(4,3);
t428 = Ifges(5,1) * t449 - Ifges(5,4) * t448 + Ifges(5,5) * qJD(4);
t427 = Ifges(5,4) * t449 - Ifges(5,2) * t448 + Ifges(5,6) * qJD(4);
t426 = Ifges(5,5) * t449 - Ifges(5,6) * t448 + Ifges(5,3) * qJD(4);
t376 = mrSges(7,2) * t412 + t423 * t441 - t489;
t366 = mrSges(6,2) * t385 + mrSges(7,2) * t383 - mrSges(6,3) * t381 - mrSges(7,3) * t378 - qJ(6) * t375 + t512 * t411 + t412 * t520 + t511 * t435 - t503 * t440 + t502 * t447;
t365 = qJDD(2) * t508 - t480;
t361 = -mrSges(6,1) * t385 + mrSges(6,3) * t382 - mrSges(7,1) * t383 + mrSges(7,3) * t380 - pkin(5) * t376 + qJ(6) * t493 + (-qJ(6) * t423 + t501) * t447 + t503 * t441 + (-mrSges(7,2) * qJ(6) + t510) * t435 + t512 * t412 + t519 * t411;
t357 = -mrSges(5,1) * t413 + mrSges(5,3) * t388 + Ifges(5,4) * t438 + Ifges(5,2) * t437 + Ifges(5,6) * qJDD(4) - pkin(4) * t367 + qJD(4) * t428 - t449 * t426 - t516;
t356 = mrSges(5,2) * t413 - mrSges(5,3) * t387 + Ifges(5,1) * t438 + Ifges(5,4) * t437 + Ifges(5,5) * qJDD(4) - pkin(9) * t367 - qJD(4) * t427 - t361 * t472 + t366 * t475 - t426 * t448;
t1 = [m(2) * t465 + t471 * (m(3) * t446 + t359 * t469 + t360 * t466) + (t474 * (m(3) * t430 - mrSges(3,1) * t479 - qJDD(2) * mrSges(3,2) + t492) + t477 * (t480 + (mrSges(3,1) - t508) * qJDD(2) + m(3) * t429 - mrSges(3,2) * t479)) * t468; mrSges(3,1) * t429 - mrSges(3,2) * t430 + t466 * (mrSges(4,2) * t419 - mrSges(4,3) * t397 - pkin(8) * t504 + t476 * t356 - t473 * t357) + t469 * (-mrSges(4,1) * t419 + mrSges(4,3) * t398 - pkin(3) * t483 + pkin(8) * t491 + t473 * t356 + t476 * t357) - pkin(2) * t365 + qJ(3) * t492 + (Ifges(4,2) * t464 + Ifges(3,3) + (Ifges(4,1) * t466 + Ifges(4,4) * t515) * t466) * qJDD(2); t365; mrSges(5,1) * t387 - mrSges(5,2) * t388 + Ifges(5,5) * t438 + Ifges(5,6) * t437 + Ifges(5,3) * qJDD(4) + pkin(4) * t481 + pkin(9) * t490 + t475 * t361 + t472 * t366 + t449 * t427 + t448 * t428; t516; t376;];
tauJ  = t1;

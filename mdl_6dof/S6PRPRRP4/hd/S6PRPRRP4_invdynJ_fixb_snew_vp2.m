% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRRP4
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
% Datum: 2019-05-04 23:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:47:46
% EndTime: 2019-05-04 23:47:48
% DurationCPUTime: 1.85s
% Computational Cost: add. (13995->239), mult. (30758->295), div. (0->0), fcn. (22765->12), ass. (0->113)
t518 = Ifges(6,1) + Ifges(7,1);
t509 = Ifges(6,4) - Ifges(7,5);
t508 = -Ifges(6,5) - Ifges(7,4);
t517 = Ifges(6,2) + Ifges(7,3);
t507 = Ifges(6,6) - Ifges(7,6);
t516 = -Ifges(6,3) - Ifges(7,2);
t466 = sin(pkin(10));
t469 = cos(pkin(10));
t452 = g(1) * t466 - g(2) * t469;
t464 = -g(3) + qJDD(1);
t467 = sin(pkin(6));
t470 = cos(pkin(6));
t515 = t452 * t470 + t464 * t467;
t477 = qJD(2) ^ 2;
t465 = sin(pkin(11));
t468 = cos(pkin(11));
t472 = sin(qJ(4));
t474 = cos(qJ(4));
t484 = t465 * t472 - t468 * t474;
t445 = t484 * qJD(2);
t453 = -g(1) * t469 - g(2) * t466;
t473 = sin(qJ(2));
t475 = cos(qJ(2));
t425 = -t473 * t453 + t475 * t515;
t485 = t465 * t474 + t468 * t472;
t446 = t485 * qJD(2);
t494 = t446 * qJD(4);
t434 = -qJDD(2) * t484 - t494;
t495 = t445 * qJD(4);
t435 = qJDD(2) * t485 - t495;
t471 = sin(qJ(5));
t512 = cos(qJ(5));
t436 = -qJD(4) * t512 + t446 * t471;
t407 = -t436 * qJD(5) + t471 * qJDD(4) + t435 * t512;
t437 = t471 * qJD(4) + t446 * t512;
t412 = mrSges(7,1) * t436 - mrSges(7,3) * t437;
t426 = t475 * t453 + t515 * t473;
t421 = -pkin(2) * t477 + qJDD(2) * qJ(3) + t426;
t442 = -t452 * t467 + t464 * t470;
t493 = qJD(2) * qJD(3);
t496 = t468 * t442 - 0.2e1 * t465 * t493;
t511 = pkin(3) * t468;
t392 = (-pkin(8) * qJDD(2) + t477 * t511 - t421) * t465 + t496;
t513 = 0.2e1 * t468;
t395 = t468 * t421 + t465 * t442 + t493 * t513;
t492 = qJDD(2) * t468;
t463 = t468 ^ 2;
t503 = t463 * t477;
t393 = -pkin(3) * t503 + pkin(8) * t492 + t395;
t386 = t472 * t392 + t474 * t393;
t433 = pkin(4) * t445 - pkin(9) * t446;
t476 = qJD(4) ^ 2;
t384 = -pkin(4) * t476 + qJDD(4) * pkin(9) - t433 * t445 + t386;
t462 = t465 ^ 2;
t482 = qJDD(3) - t425;
t408 = (-pkin(2) - t511) * qJDD(2) + (-qJ(3) + (-t462 - t463) * pkin(8)) * t477 + t482;
t388 = (-t435 + t495) * pkin(9) + (-t434 + t494) * pkin(4) + t408;
t380 = -t471 * t384 + t388 * t512;
t411 = pkin(5) * t436 - qJ(6) * t437;
t432 = qJDD(5) - t434;
t444 = qJD(5) + t445;
t443 = t444 ^ 2;
t378 = -t432 * pkin(5) - t443 * qJ(6) + t437 * t411 + qJDD(6) - t380;
t417 = -mrSges(7,2) * t436 + mrSges(7,3) * t444;
t487 = -m(7) * t378 + t432 * mrSges(7,1) + t444 * t417;
t374 = t407 * mrSges(7,2) + t437 * t412 - t487;
t381 = t512 * t384 + t471 * t388;
t377 = -pkin(5) * t443 + qJ(6) * t432 + 0.2e1 * qJD(6) * t444 - t411 * t436 + t381;
t406 = qJD(5) * t437 - qJDD(4) * t512 + t435 * t471;
t420 = -mrSges(7,1) * t444 + mrSges(7,2) * t437;
t491 = m(7) * t377 + t432 * mrSges(7,3) + t444 * t420;
t498 = t509 * t436 - t518 * t437 + t508 * t444;
t499 = t517 * t436 - t509 * t437 - t507 * t444;
t514 = -t406 * t507 - t407 * t508 - t516 * t432 - t436 * t498 - t437 * t499 + mrSges(6,1) * t380 - mrSges(7,1) * t378 - mrSges(6,2) * t381 + mrSges(7,3) * t377 - pkin(5) * t374 + qJ(6) * (-t406 * mrSges(7,2) - t436 * t412 + t491);
t510 = -mrSges(6,3) - mrSges(7,2);
t505 = mrSges(4,2) * t465;
t430 = mrSges(5,1) * t445 + mrSges(5,2) * t446;
t441 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t446;
t419 = mrSges(6,1) * t444 - mrSges(6,3) * t437;
t497 = -mrSges(6,1) * t436 - mrSges(6,2) * t437 - t412;
t370 = m(6) * t381 - t432 * mrSges(6,2) + t406 * t510 - t444 * t419 + t436 * t497 + t491;
t418 = -mrSges(6,2) * t444 - mrSges(6,3) * t436;
t372 = m(6) * t380 + t432 * mrSges(6,1) + t407 * t510 + t444 * t418 + t437 * t497 + t487;
t488 = t512 * t370 - t372 * t471;
t361 = m(5) * t386 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t434 - qJD(4) * t441 - t430 * t445 + t488;
t385 = t474 * t392 - t472 * t393;
t440 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t445;
t383 = -qJDD(4) * pkin(4) - t476 * pkin(9) + t446 * t433 - t385;
t379 = -0.2e1 * qJD(6) * t437 + (t436 * t444 - t407) * qJ(6) + (t437 * t444 + t406) * pkin(5) + t383;
t375 = m(7) * t379 + mrSges(7,1) * t406 - t407 * mrSges(7,3) + t417 * t436 - t437 * t420;
t478 = -m(6) * t383 - t406 * mrSges(6,1) - mrSges(6,2) * t407 - t436 * t418 - t419 * t437 - t375;
t367 = m(5) * t385 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t435 + qJD(4) * t440 - t430 * t446 + t478;
t501 = t472 * t361 + t474 * t367;
t365 = t471 * t370 + t512 * t372;
t500 = t507 * t436 + t508 * t437 + t516 * t444;
t394 = -t421 * t465 + t496;
t483 = mrSges(4,3) * qJDD(2) + t477 * (-mrSges(4,1) * t468 + t505);
t357 = m(4) * t394 - t465 * t483 + t501;
t489 = t474 * t361 - t472 * t367;
t358 = m(4) * t395 + t468 * t483 + t489;
t490 = -t357 * t465 + t468 * t358;
t481 = m(5) * t408 - t434 * mrSges(5,1) + t435 * mrSges(5,2) + t445 * t440 + t446 * t441 + t365;
t416 = -qJDD(2) * pkin(2) - t477 * qJ(3) + t482;
t479 = -m(4) * t416 + mrSges(4,1) * t492 - t481 + (t462 * t477 + t503) * mrSges(4,3);
t424 = Ifges(5,1) * t446 - Ifges(5,4) * t445 + Ifges(5,5) * qJD(4);
t423 = Ifges(5,4) * t446 - Ifges(5,2) * t445 + Ifges(5,6) * qJD(4);
t422 = Ifges(5,5) * t446 - Ifges(5,6) * t445 + Ifges(5,3) * qJD(4);
t364 = mrSges(6,2) * t383 + mrSges(7,2) * t378 - mrSges(6,3) * t380 - mrSges(7,3) * t379 - qJ(6) * t375 - t509 * t406 + t518 * t407 - t508 * t432 + t500 * t436 + t499 * t444;
t363 = -mrSges(6,1) * t383 - mrSges(7,1) * t379 + mrSges(7,2) * t377 + mrSges(6,3) * t381 - pkin(5) * t375 - t517 * t406 + t509 * t407 + t507 * t432 + t500 * t437 - t498 * t444;
t362 = qJDD(2) * t505 - t479;
t355 = -mrSges(5,1) * t408 + mrSges(5,3) * t386 + Ifges(5,4) * t435 + Ifges(5,2) * t434 + Ifges(5,6) * qJDD(4) - pkin(4) * t365 + qJD(4) * t424 - t446 * t422 - t514;
t354 = mrSges(5,2) * t408 - mrSges(5,3) * t385 + Ifges(5,1) * t435 + Ifges(5,4) * t434 + Ifges(5,5) * qJDD(4) - pkin(9) * t365 - qJD(4) * t423 - t471 * t363 + t364 * t512 - t445 * t422;
t1 = [m(2) * t464 + t470 * (m(3) * t442 + t357 * t468 + t358 * t465) + (t473 * (m(3) * t426 - mrSges(3,1) * t477 - qJDD(2) * mrSges(3,2) + t490) + t475 * (-t477 * mrSges(3,2) + m(3) * t425 + t479 + (mrSges(3,1) - t505) * qJDD(2))) * t467; mrSges(3,1) * t425 - mrSges(3,2) * t426 + t465 * (mrSges(4,2) * t416 - mrSges(4,3) * t394 - pkin(8) * t501 + t474 * t354 - t472 * t355) + t468 * (-mrSges(4,1) * t416 + mrSges(4,3) * t395 - pkin(3) * t481 + pkin(8) * t489 + t472 * t354 + t474 * t355) - pkin(2) * t362 + qJ(3) * t490 + (Ifges(4,2) * t463 + Ifges(3,3) + (Ifges(4,1) * t465 + Ifges(4,4) * t513) * t465) * qJDD(2); t362; mrSges(5,1) * t385 - mrSges(5,2) * t386 + Ifges(5,5) * t435 + Ifges(5,6) * t434 + Ifges(5,3) * qJDD(4) + pkin(4) * t478 + pkin(9) * t488 + t363 * t512 + t471 * t364 + t446 * t423 + t445 * t424; t514; t374;];
tauJ  = t1;

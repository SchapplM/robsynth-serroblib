% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-05-05 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:45:00
% EndTime: 2019-05-05 20:45:03
% DurationCPUTime: 2.04s
% Computational Cost: add. (12675->271), mult. (25335->325), div. (0->0), fcn. (14103->8), ass. (0->116)
t482 = sin(qJ(3));
t486 = cos(qJ(3));
t519 = -Ifges(5,6) - Ifges(4,4);
t532 = t486 * (Ifges(4,1) + Ifges(5,2)) + t482 * t519;
t531 = t486 * t519 - t482 * (-Ifges(4,2) - Ifges(5,3));
t530 = -2 * qJD(4);
t521 = Ifges(4,5) - Ifges(5,4);
t520 = Ifges(4,6) - Ifges(5,5);
t527 = (-t531 * qJD(1) + t520 * qJD(3)) * t486;
t489 = qJD(1) ^ 2;
t483 = sin(qJ(1));
t487 = cos(qJ(1));
t508 = g(1) * t483 - t487 * g(2);
t499 = -qJ(2) * t489 + qJDD(2) - t508;
t526 = -pkin(1) - pkin(7);
t435 = qJDD(1) * t526 + t499;
t429 = -g(3) * t486 + t482 * t435;
t454 = (pkin(3) * t482 - qJ(4) * t486) * qJD(1);
t488 = qJD(3) ^ 2;
t514 = qJD(1) * t482;
t411 = pkin(3) * t488 - qJDD(3) * qJ(4) + qJD(3) * t530 + t454 * t514 - t429;
t504 = -g(1) * t487 - g(2) * t483;
t500 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t504;
t525 = pkin(3) + pkin(8);
t524 = pkin(8) * t489;
t523 = g(3) * t482;
t522 = mrSges(4,1) - mrSges(5,2);
t518 = t435 * t486;
t513 = qJD(1) * qJD(3);
t509 = t486 * t513;
t457 = qJDD(1) * t482 + t509;
t472 = t486 * qJD(1);
t465 = pkin(4) * t472 - qJD(3) * pkin(8);
t479 = t482 ^ 2;
t471 = t482 * t513;
t458 = qJDD(1) * t486 - t471;
t493 = pkin(3) * t509 + t472 * t530 + t500 + (-t458 + t471) * qJ(4);
t398 = -t465 * t472 + t525 * t457 + (-pkin(4) * t479 + t526) * t489 + t493;
t498 = -qJ(4) * t488 + t454 * t472 + qJDD(4) - t518;
t403 = pkin(4) * t458 - t525 * qJDD(3) + (pkin(4) * t513 + t486 * t524 - g(3)) * t482 + t498;
t481 = sin(qJ(5));
t485 = cos(qJ(5));
t387 = -t398 * t481 + t485 * t403;
t452 = -qJD(3) * t481 + t485 * t514;
t424 = qJD(5) * t452 + qJDD(3) * t485 + t457 * t481;
t451 = qJDD(5) + t458;
t453 = qJD(3) * t485 + t481 * t514;
t469 = t472 + qJD(5);
t384 = (t452 * t469 - t424) * pkin(9) + (t452 * t453 + t451) * pkin(5) + t387;
t388 = t485 * t398 + t481 * t403;
t423 = -qJD(5) * t453 - qJDD(3) * t481 + t457 * t485;
t433 = pkin(5) * t469 - pkin(9) * t453;
t450 = t452 ^ 2;
t385 = -pkin(5) * t450 + pkin(9) * t423 - t433 * t469 + t388;
t480 = sin(qJ(6));
t484 = cos(qJ(6));
t382 = t384 * t484 - t385 * t480;
t425 = t452 * t484 - t453 * t480;
t396 = qJD(6) * t425 + t423 * t480 + t424 * t484;
t426 = t452 * t480 + t453 * t484;
t409 = -mrSges(7,1) * t425 + mrSges(7,2) * t426;
t466 = qJD(6) + t469;
t414 = -mrSges(7,2) * t466 + mrSges(7,3) * t425;
t447 = qJDD(6) + t451;
t379 = m(7) * t382 + mrSges(7,1) * t447 - mrSges(7,3) * t396 - t409 * t426 + t414 * t466;
t383 = t384 * t480 + t385 * t484;
t395 = -qJD(6) * t426 + t423 * t484 - t424 * t480;
t415 = mrSges(7,1) * t466 - mrSges(7,3) * t426;
t380 = m(7) * t383 - mrSges(7,2) * t447 + mrSges(7,3) * t395 + t409 * t425 - t415 * t466;
t372 = t484 * t379 + t480 * t380;
t516 = t532 * qJD(1) + t521 * qJD(3);
t463 = mrSges(5,1) * t514 - qJD(3) * mrSges(5,3);
t515 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t514 - t463;
t510 = t526 * t489;
t427 = -mrSges(6,1) * t452 + mrSges(6,2) * t453;
t430 = -mrSges(6,2) * t469 + mrSges(6,3) * t452;
t369 = m(6) * t387 + mrSges(6,1) * t451 - mrSges(6,3) * t424 - t427 * t453 + t430 * t469 + t372;
t431 = mrSges(6,1) * t469 - mrSges(6,3) * t453;
t506 = -t379 * t480 + t484 * t380;
t370 = m(6) * t388 - mrSges(6,2) * t451 + mrSges(6,3) * t423 + t427 * t452 - t431 * t469 + t506;
t507 = -t369 * t481 + t485 * t370;
t455 = (-mrSges(5,2) * t482 - mrSges(5,3) * t486) * qJD(1);
t505 = qJD(1) * (-t455 - (mrSges(4,1) * t482 + mrSges(4,2) * t486) * qJD(1));
t428 = t518 + t523;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t472;
t464 = mrSges(5,1) * t472 + qJD(3) * mrSges(5,2);
t402 = -pkin(4) * t457 + qJD(3) * t465 - t479 * t524 - t411;
t390 = -pkin(5) * t423 - pkin(9) * t450 + t433 * t453 + t402;
t497 = m(7) * t390 - mrSges(7,1) * t395 + t396 * mrSges(7,2) - t414 * t425 + t426 * t415;
t492 = -m(6) * t402 + mrSges(6,1) * t423 - t424 * mrSges(6,2) + t430 * t452 - t453 * t431 - t497;
t491 = -m(5) * t411 + qJDD(3) * mrSges(5,3) + qJD(3) * t464 - t492;
t367 = t369 * t485 + t370 * t481;
t412 = -qJDD(3) * pkin(3) + t498 - t523;
t496 = -m(5) * t412 - t458 * mrSges(5,1) - t367;
t503 = (m(4) * t428 - mrSges(4,3) * t458 + qJD(3) * t515 + qJDD(3) * t522 + t486 * t505 + t496) * t486 + (-qJDD(3) * mrSges(4,2) + t491 + t482 * t505 + m(4) * t429 - qJD(3) * t462 + (-mrSges(4,3) - mrSges(5,1)) * t457) * t482;
t410 = pkin(3) * t457 + t493 + t510;
t501 = m(5) * t410 - t458 * mrSges(5,3) + t507;
t405 = Ifges(7,4) * t426 + Ifges(7,2) * t425 + Ifges(7,6) * t466;
t406 = Ifges(7,1) * t426 + Ifges(7,4) * t425 + Ifges(7,5) * t466;
t495 = mrSges(7,1) * t382 - mrSges(7,2) * t383 + Ifges(7,5) * t396 + Ifges(7,6) * t395 + Ifges(7,3) * t447 + t426 * t405 - t425 * t406;
t417 = Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t469;
t418 = Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t469;
t490 = mrSges(6,1) * t387 - mrSges(6,2) * t388 + Ifges(6,5) * t424 + Ifges(6,6) * t423 + Ifges(6,3) * t451 + pkin(5) * t372 + t453 * t417 - t452 * t418 + t495;
t437 = -qJDD(1) * pkin(1) + t499;
t436 = pkin(1) * t489 - t500;
t434 = t510 + t500;
t416 = Ifges(6,5) * t453 + Ifges(6,6) * t452 + Ifges(6,3) * t469;
t404 = Ifges(7,5) * t426 + Ifges(7,6) * t425 + Ifges(7,3) * t466;
t374 = mrSges(7,2) * t390 - mrSges(7,3) * t382 + Ifges(7,1) * t396 + Ifges(7,4) * t395 + Ifges(7,5) * t447 + t404 * t425 - t405 * t466;
t373 = -mrSges(7,1) * t390 + mrSges(7,3) * t383 + Ifges(7,4) * t396 + Ifges(7,2) * t395 + Ifges(7,6) * t447 - t404 * t426 + t406 * t466;
t366 = qJDD(3) * mrSges(5,2) + qJD(3) * t463 + t455 * t472 - t496;
t365 = -mrSges(5,2) * t457 + (-t463 * t482 - t464 * t486) * qJD(1) + t501;
t363 = mrSges(6,2) * t402 - mrSges(6,3) * t387 + Ifges(6,1) * t424 + Ifges(6,4) * t423 + Ifges(6,5) * t451 - pkin(9) * t372 - t373 * t480 + t374 * t484 + t416 * t452 - t417 * t469;
t362 = -mrSges(6,1) * t402 + mrSges(6,3) * t388 + Ifges(6,4) * t424 + Ifges(6,2) * t423 + Ifges(6,6) * t451 - pkin(5) * t497 + pkin(9) * t506 + t484 * t373 + t480 * t374 - t453 * t416 + t469 * t418;
t361 = m(3) * t437 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t489 + t503;
t1 = [mrSges(2,1) * t508 - mrSges(2,2) * t504 + mrSges(3,2) * t437 - mrSges(3,3) * t436 + t486 * (mrSges(5,1) * t412 + mrSges(4,2) * t434 - mrSges(4,3) * t428 - mrSges(5,3) * t410 + pkin(4) * t367 - qJ(4) * t365 + t490) - t482 * (-mrSges(4,1) * t434 - mrSges(5,1) * t411 + mrSges(5,2) * t410 + mrSges(4,3) * t429 - pkin(3) * t365 - pkin(4) * t492 - pkin(8) * t507 - t485 * t362 - t481 * t363) - pkin(7) * t503 - pkin(1) * t361 + t532 * t458 + (-t482 * t520 + t486 * t521) * qJDD(3) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t482 * t516 - t527) * qJD(3) + t531 * t457 + (-m(3) * t436 + m(4) * t434 + mrSges(3,2) * t489 + t501 + mrSges(4,2) * t458 + qJDD(1) * mrSges(3,3) + t522 * t457 + ((t462 - t464) * t486 + t515 * t482) * qJD(1)) * qJ(2); t361; mrSges(4,1) * t428 - mrSges(4,2) * t429 + mrSges(5,2) * t412 - mrSges(5,3) * t411 + t485 * t363 - t481 * t362 - pkin(8) * t367 - pkin(3) * t366 + qJ(4) * t491 + t521 * t458 + (-mrSges(5,1) * qJ(4) - t520) * t457 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + (t527 + (-qJ(4) * t455 + t516) * t482) * qJD(1); t366; t490; t495;];
tauJ  = t1;

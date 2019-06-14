% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:24:32
% EndTime: 2019-05-05 10:24:37
% DurationCPUTime: 4.19s
% Computational Cost: add. (45299->291), mult. (94366->375), div. (0->0), fcn. (70989->14), ass. (0->124)
t489 = sin(pkin(12));
t491 = cos(pkin(12));
t474 = g(1) * t489 - g(2) * t491;
t488 = -g(3) + qJDD(1);
t490 = sin(pkin(6));
t492 = cos(pkin(6));
t523 = t474 * t492 + t488 * t490;
t475 = -g(1) * t491 - g(2) * t489;
t497 = sin(qJ(2));
t502 = cos(qJ(2));
t446 = -t497 * t475 + t502 * t523;
t447 = t502 * t475 + t523 * t497;
t503 = qJD(2) ^ 2;
t441 = -pkin(2) * t503 + qJDD(2) * pkin(8) + t447;
t458 = -t474 * t490 + t488 * t492;
t496 = sin(qJ(3));
t501 = cos(qJ(3));
t427 = -t441 * t496 + t501 * t458;
t518 = qJD(2) * qJD(3);
t517 = t501 * t518;
t472 = qJDD(2) * t496 + t517;
t414 = (-t472 + t517) * pkin(9) + (t496 * t501 * t503 + qJDD(3)) * pkin(3) + t427;
t428 = t501 * t441 + t496 * t458;
t473 = qJDD(2) * t501 - t496 * t518;
t519 = t496 * qJD(2);
t479 = qJD(3) * pkin(3) - pkin(9) * t519;
t487 = t501 ^ 2;
t416 = -pkin(3) * t487 * t503 + pkin(9) * t473 - qJD(3) * t479 + t428;
t495 = sin(qJ(4));
t500 = cos(qJ(4));
t392 = t500 * t414 - t416 * t495;
t464 = (-t496 * t495 + t501 * t500) * qJD(2);
t437 = qJD(4) * t464 + t472 * t500 + t473 * t495;
t465 = (t501 * t495 + t496 * t500) * qJD(2);
t485 = qJDD(3) + qJDD(4);
t486 = qJD(3) + qJD(4);
t388 = (t464 * t486 - t437) * pkin(10) + (t464 * t465 + t485) * pkin(4) + t392;
t393 = t495 * t414 + t500 * t416;
t436 = -qJD(4) * t465 - t472 * t495 + t473 * t500;
t457 = pkin(4) * t486 - pkin(10) * t465;
t460 = t464 ^ 2;
t390 = -pkin(4) * t460 + pkin(10) * t436 - t457 * t486 + t393;
t494 = sin(qJ(5));
t499 = cos(qJ(5));
t385 = t494 * t388 + t499 * t390;
t451 = t464 * t494 + t465 * t499;
t408 = -qJD(5) * t451 + t436 * t499 - t437 * t494;
t450 = t464 * t499 - t465 * t494;
t424 = -mrSges(6,1) * t450 + mrSges(6,2) * t451;
t484 = qJD(5) + t486;
t439 = mrSges(6,1) * t484 - mrSges(6,3) * t451;
t483 = qJDD(5) + t485;
t425 = -pkin(5) * t450 - pkin(11) * t451;
t482 = t484 ^ 2;
t382 = -pkin(5) * t482 + pkin(11) * t483 + t425 * t450 + t385;
t509 = -qJDD(2) * pkin(2) - t446;
t426 = -pkin(3) * t473 + t479 * t519 + (-pkin(9) * t487 - pkin(8)) * t503 + t509;
t398 = -pkin(4) * t436 - pkin(10) * t460 + t465 * t457 + t426;
t409 = qJD(5) * t450 + t436 * t494 + t437 * t499;
t386 = (-t450 * t484 - t409) * pkin(11) + (t451 * t484 - t408) * pkin(5) + t398;
t493 = sin(qJ(6));
t498 = cos(qJ(6));
t379 = -t382 * t493 + t386 * t498;
t430 = -t451 * t493 + t484 * t498;
t396 = qJD(6) * t430 + t409 * t498 + t483 * t493;
t406 = qJDD(6) - t408;
t431 = t451 * t498 + t484 * t493;
t415 = -mrSges(7,1) * t430 + mrSges(7,2) * t431;
t445 = qJD(6) - t450;
t417 = -mrSges(7,2) * t445 + mrSges(7,3) * t430;
t376 = m(7) * t379 + mrSges(7,1) * t406 - mrSges(7,3) * t396 - t415 * t431 + t417 * t445;
t380 = t382 * t498 + t386 * t493;
t395 = -qJD(6) * t431 - t409 * t493 + t483 * t498;
t418 = mrSges(7,1) * t445 - mrSges(7,3) * t431;
t377 = m(7) * t380 - mrSges(7,2) * t406 + mrSges(7,3) * t395 + t415 * t430 - t418 * t445;
t513 = -t376 * t493 + t498 * t377;
t364 = m(6) * t385 - mrSges(6,2) * t483 + mrSges(6,3) * t408 + t424 * t450 - t439 * t484 + t513;
t384 = t388 * t499 - t390 * t494;
t438 = -mrSges(6,2) * t484 + mrSges(6,3) * t450;
t381 = -pkin(5) * t483 - pkin(11) * t482 + t425 * t451 - t384;
t510 = -m(7) * t381 + t395 * mrSges(7,1) - mrSges(7,2) * t396 + t430 * t417 - t418 * t431;
t372 = m(6) * t384 + mrSges(6,1) * t483 - mrSges(6,3) * t409 - t424 * t451 + t438 * t484 + t510;
t361 = t494 * t364 + t499 * t372;
t452 = -mrSges(5,1) * t464 + mrSges(5,2) * t465;
t455 = -mrSges(5,2) * t486 + mrSges(5,3) * t464;
t358 = m(5) * t392 + mrSges(5,1) * t485 - mrSges(5,3) * t437 - t452 * t465 + t455 * t486 + t361;
t456 = mrSges(5,1) * t486 - mrSges(5,3) * t465;
t514 = t499 * t364 - t372 * t494;
t359 = m(5) * t393 - mrSges(5,2) * t485 + mrSges(5,3) * t436 + t452 * t464 - t456 * t486 + t514;
t352 = t500 * t358 + t495 * t359;
t366 = t498 * t376 + t493 * t377;
t520 = qJD(2) * t501;
t471 = (-t501 * mrSges(4,1) + t496 * mrSges(4,2)) * qJD(2);
t477 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t520;
t350 = m(4) * t427 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t472 + qJD(3) * t477 - t471 * t519 + t352;
t476 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t519;
t515 = -t358 * t495 + t500 * t359;
t351 = m(4) * t428 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t473 - qJD(3) * t476 + t471 * t520 + t515;
t516 = -t496 * t350 + t501 * t351;
t511 = m(6) * t398 - t408 * mrSges(6,1) + t409 * mrSges(6,2) - t450 * t438 + t451 * t439 + t366;
t399 = Ifges(7,5) * t431 + Ifges(7,6) * t430 + Ifges(7,3) * t445;
t401 = Ifges(7,1) * t431 + Ifges(7,4) * t430 + Ifges(7,5) * t445;
t369 = -mrSges(7,1) * t381 + mrSges(7,3) * t380 + Ifges(7,4) * t396 + Ifges(7,2) * t395 + Ifges(7,6) * t406 - t399 * t431 + t401 * t445;
t400 = Ifges(7,4) * t431 + Ifges(7,2) * t430 + Ifges(7,6) * t445;
t370 = mrSges(7,2) * t381 - mrSges(7,3) * t379 + Ifges(7,1) * t396 + Ifges(7,4) * t395 + Ifges(7,5) * t406 + t399 * t430 - t400 * t445;
t420 = Ifges(6,4) * t451 + Ifges(6,2) * t450 + Ifges(6,6) * t484;
t421 = Ifges(6,1) * t451 + Ifges(6,4) * t450 + Ifges(6,5) * t484;
t508 = mrSges(6,1) * t384 - mrSges(6,2) * t385 + Ifges(6,5) * t409 + Ifges(6,6) * t408 + Ifges(6,3) * t483 + pkin(5) * t510 + pkin(11) * t513 + t498 * t369 + t493 * t370 + t451 * t420 - t450 * t421;
t507 = m(5) * t426 - t436 * mrSges(5,1) + t437 * mrSges(5,2) - t464 * t455 + t465 * t456 + t511;
t506 = mrSges(7,1) * t379 - mrSges(7,2) * t380 + Ifges(7,5) * t396 + Ifges(7,6) * t395 + Ifges(7,3) * t406 + t400 * t431 - t401 * t430;
t443 = Ifges(5,4) * t465 + Ifges(5,2) * t464 + Ifges(5,6) * t486;
t444 = Ifges(5,1) * t465 + Ifges(5,4) * t464 + Ifges(5,5) * t486;
t505 = mrSges(5,1) * t392 - mrSges(5,2) * t393 + Ifges(5,5) * t437 + Ifges(5,6) * t436 + Ifges(5,3) * t485 + pkin(4) * t361 + t465 * t443 - t464 * t444 + t508;
t440 = -pkin(8) * t503 + t509;
t504 = -m(4) * t440 + t473 * mrSges(4,1) - t472 * mrSges(4,2) - t476 * t519 + t477 * t520 - t507;
t463 = Ifges(4,5) * qJD(3) + (t496 * Ifges(4,1) + t501 * Ifges(4,4)) * qJD(2);
t462 = Ifges(4,6) * qJD(3) + (t496 * Ifges(4,4) + t501 * Ifges(4,2)) * qJD(2);
t442 = Ifges(5,5) * t465 + Ifges(5,6) * t464 + Ifges(5,3) * t486;
t419 = Ifges(6,5) * t451 + Ifges(6,6) * t450 + Ifges(6,3) * t484;
t354 = -mrSges(6,1) * t398 + mrSges(6,3) * t385 + Ifges(6,4) * t409 + Ifges(6,2) * t408 + Ifges(6,6) * t483 - pkin(5) * t366 - t419 * t451 + t421 * t484 - t506;
t353 = mrSges(6,2) * t398 - mrSges(6,3) * t384 + Ifges(6,1) * t409 + Ifges(6,4) * t408 + Ifges(6,5) * t483 - pkin(11) * t366 - t369 * t493 + t370 * t498 + t419 * t450 - t420 * t484;
t348 = mrSges(5,2) * t426 - mrSges(5,3) * t392 + Ifges(5,1) * t437 + Ifges(5,4) * t436 + Ifges(5,5) * t485 - pkin(10) * t361 + t353 * t499 - t354 * t494 + t442 * t464 - t443 * t486;
t347 = -mrSges(5,1) * t426 + mrSges(5,3) * t393 + Ifges(5,4) * t437 + Ifges(5,2) * t436 + Ifges(5,6) * t485 - pkin(4) * t511 + pkin(10) * t514 + t494 * t353 + t499 * t354 - t465 * t442 + t486 * t444;
t1 = [m(2) * t488 + t492 * (m(3) * t458 + t501 * t350 + t496 * t351) + (t497 * (m(3) * t447 - t503 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t516) + t502 * (m(3) * t446 + qJDD(2) * mrSges(3,1) - t503 * mrSges(3,2) + t504)) * t490; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t446 - mrSges(3,2) * t447 + t496 * (mrSges(4,2) * t440 - mrSges(4,3) * t427 + Ifges(4,1) * t472 + Ifges(4,4) * t473 + Ifges(4,5) * qJDD(3) - pkin(9) * t352 - qJD(3) * t462 - t347 * t495 + t348 * t500) + t501 * (-mrSges(4,1) * t440 + mrSges(4,3) * t428 + Ifges(4,4) * t472 + Ifges(4,2) * t473 + Ifges(4,6) * qJDD(3) - pkin(3) * t507 + pkin(9) * t515 + qJD(3) * t463 + t500 * t347 + t495 * t348) + pkin(2) * t504 + pkin(8) * t516; t505 + pkin(3) * t352 + Ifges(4,5) * t472 + Ifges(4,6) * t473 + mrSges(4,1) * t427 - mrSges(4,2) * t428 + Ifges(4,3) * qJDD(3) + (t496 * t462 - t501 * t463) * qJD(2); t505; t508; t506;];
tauJ  = t1;

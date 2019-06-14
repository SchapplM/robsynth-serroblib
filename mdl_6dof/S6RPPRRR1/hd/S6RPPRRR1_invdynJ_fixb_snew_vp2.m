% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:10:18
% EndTime: 2019-05-05 15:10:21
% DurationCPUTime: 3.20s
% Computational Cost: add. (32229->266), mult. (73323->338), div. (0->0), fcn. (54212->12), ass. (0->119)
t479 = qJD(1) ^ 2;
t469 = cos(pkin(11));
t504 = pkin(3) * t469;
t467 = sin(pkin(11));
t503 = t467 * mrSges(4,2);
t464 = t469 ^ 2;
t502 = t464 * t479;
t474 = sin(qJ(1));
t478 = cos(qJ(1));
t495 = g(1) * t474 - g(2) * t478;
t450 = qJDD(1) * pkin(1) + t495;
t490 = -g(1) * t478 - g(2) * t474;
t451 = -pkin(1) * t479 + t490;
t468 = sin(pkin(10));
t470 = cos(pkin(10));
t438 = t450 * t468 + t451 * t470;
t430 = -pkin(2) * t479 + qJDD(1) * qJ(3) + t438;
t466 = -g(3) + qJDD(2);
t497 = qJD(1) * qJD(3);
t500 = t466 * t469 - 0.2e1 * t467 * t497;
t413 = (-pkin(7) * qJDD(1) + t479 * t504 - t430) * t467 + t500;
t417 = t467 * t466 + (t430 + 0.2e1 * t497) * t469;
t496 = qJDD(1) * t469;
t414 = -pkin(3) * t502 + pkin(7) * t496 + t417;
t473 = sin(qJ(4));
t477 = cos(qJ(4));
t389 = t413 * t477 - t414 * t473;
t488 = t467 * t477 + t469 * t473;
t487 = -t467 * t473 + t469 * t477;
t443 = t487 * qJD(1);
t498 = qJD(4) * t443;
t436 = qJDD(1) * t488 + t498;
t444 = t488 * qJD(1);
t378 = (-t436 + t498) * pkin(8) + (t443 * t444 + qJDD(4)) * pkin(4) + t389;
t390 = t413 * t473 + t414 * t477;
t435 = -qJD(4) * t444 + qJDD(1) * t487;
t441 = qJD(4) * pkin(4) - pkin(8) * t444;
t442 = t443 ^ 2;
t380 = -pkin(4) * t442 + pkin(8) * t435 - qJD(4) * t441 + t390;
t472 = sin(qJ(5));
t476 = cos(qJ(5));
t376 = t378 * t472 + t380 * t476;
t429 = t443 * t472 + t444 * t476;
t398 = -qJD(5) * t429 + t435 * t476 - t436 * t472;
t428 = t443 * t476 - t444 * t472;
t408 = -mrSges(6,1) * t428 + mrSges(6,2) * t429;
t465 = qJD(4) + qJD(5);
t421 = mrSges(6,1) * t465 - mrSges(6,3) * t429;
t462 = qJDD(4) + qJDD(5);
t409 = -pkin(5) * t428 - pkin(9) * t429;
t461 = t465 ^ 2;
t372 = -pkin(5) * t461 + pkin(9) * t462 + t409 * t428 + t376;
t463 = t467 ^ 2;
t437 = t450 * t470 - t451 * t468;
t489 = qJDD(3) - t437;
t415 = (-pkin(2) - t504) * qJDD(1) + (-qJ(3) + (-t463 - t464) * pkin(7)) * t479 + t489;
t385 = -pkin(4) * t435 - pkin(8) * t442 + t441 * t444 + t415;
t399 = qJD(5) * t428 + t435 * t472 + t436 * t476;
t373 = (-t428 * t465 - t399) * pkin(9) + (t429 * t465 - t398) * pkin(5) + t385;
t471 = sin(qJ(6));
t475 = cos(qJ(6));
t369 = -t372 * t471 + t373 * t475;
t418 = -t429 * t471 + t465 * t475;
t383 = qJD(6) * t418 + t399 * t475 + t462 * t471;
t397 = qJDD(6) - t398;
t419 = t429 * t475 + t465 * t471;
t400 = -mrSges(7,1) * t418 + mrSges(7,2) * t419;
t423 = qJD(6) - t428;
t401 = -mrSges(7,2) * t423 + mrSges(7,3) * t418;
t366 = m(7) * t369 + mrSges(7,1) * t397 - mrSges(7,3) * t383 - t400 * t419 + t401 * t423;
t370 = t372 * t475 + t373 * t471;
t382 = -qJD(6) * t419 - t399 * t471 + t462 * t475;
t402 = mrSges(7,1) * t423 - mrSges(7,3) * t419;
t367 = m(7) * t370 - mrSges(7,2) * t397 + mrSges(7,3) * t382 + t400 * t418 - t402 * t423;
t491 = -t366 * t471 + t367 * t475;
t354 = m(6) * t376 - mrSges(6,2) * t462 + mrSges(6,3) * t398 + t408 * t428 - t421 * t465 + t491;
t375 = t378 * t476 - t380 * t472;
t420 = -mrSges(6,2) * t465 + mrSges(6,3) * t428;
t371 = -pkin(5) * t462 - pkin(9) * t461 + t409 * t429 - t375;
t484 = -m(7) * t371 + t382 * mrSges(7,1) - mrSges(7,2) * t383 + t418 * t401 - t402 * t419;
t362 = m(6) * t375 + mrSges(6,1) * t462 - mrSges(6,3) * t399 - t408 * t429 + t420 * t465 + t484;
t350 = t354 * t472 + t362 * t476;
t433 = -mrSges(5,1) * t443 + mrSges(5,2) * t444;
t439 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t443;
t348 = m(5) * t389 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t436 + qJD(4) * t439 - t433 * t444 + t350;
t440 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t444;
t492 = t354 * t476 - t362 * t472;
t349 = m(5) * t390 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t435 - qJD(4) * t440 + t433 * t443 + t492;
t501 = t348 * t477 + t349 * t473;
t356 = t366 * t475 + t367 * t471;
t416 = -t430 * t467 + t500;
t486 = mrSges(4,3) * qJDD(1) + t479 * (-mrSges(4,1) * t469 + t503);
t341 = m(4) * t416 - t467 * t486 + t501;
t493 = -t348 * t473 + t349 * t477;
t342 = m(4) * t417 + t469 * t486 + t493;
t494 = -t341 * t467 + t342 * t469;
t485 = m(6) * t385 - mrSges(6,1) * t398 + mrSges(6,2) * t399 - t420 * t428 + t421 * t429 + t356;
t386 = Ifges(7,5) * t419 + Ifges(7,6) * t418 + Ifges(7,3) * t423;
t388 = Ifges(7,1) * t419 + Ifges(7,4) * t418 + Ifges(7,5) * t423;
t359 = -mrSges(7,1) * t371 + mrSges(7,3) * t370 + Ifges(7,4) * t383 + Ifges(7,2) * t382 + Ifges(7,6) * t397 - t386 * t419 + t388 * t423;
t387 = Ifges(7,4) * t419 + Ifges(7,2) * t418 + Ifges(7,6) * t423;
t360 = mrSges(7,2) * t371 - mrSges(7,3) * t369 + Ifges(7,1) * t383 + Ifges(7,4) * t382 + Ifges(7,5) * t397 + t386 * t418 - t387 * t423;
t404 = Ifges(6,4) * t429 + Ifges(6,2) * t428 + Ifges(6,6) * t465;
t405 = Ifges(6,1) * t429 + Ifges(6,4) * t428 + Ifges(6,5) * t465;
t483 = mrSges(6,1) * t375 - mrSges(6,2) * t376 + Ifges(6,5) * t399 + Ifges(6,6) * t398 + Ifges(6,3) * t462 + pkin(5) * t484 + pkin(9) * t491 + t359 * t475 + t360 * t471 + t404 * t429 - t405 * t428;
t482 = m(5) * t415 - mrSges(5,1) * t435 + mrSges(5,2) * t436 - t439 * t443 + t440 * t444 + t485;
t481 = mrSges(7,1) * t369 - mrSges(7,2) * t370 + Ifges(7,5) * t383 + Ifges(7,6) * t382 + Ifges(7,3) * t397 + t387 * t419 - t388 * t418;
t424 = -qJDD(1) * pkin(2) - qJ(3) * t479 + t489;
t480 = -m(4) * t424 + mrSges(4,1) * t496 - t482 + (t463 * t479 + t502) * mrSges(4,3);
t427 = Ifges(5,1) * t444 + Ifges(5,4) * t443 + Ifges(5,5) * qJD(4);
t426 = Ifges(5,4) * t444 + Ifges(5,2) * t443 + Ifges(5,6) * qJD(4);
t425 = Ifges(5,5) * t444 + Ifges(5,6) * t443 + Ifges(5,3) * qJD(4);
t403 = Ifges(6,5) * t429 + Ifges(6,6) * t428 + Ifges(6,3) * t465;
t351 = qJDD(1) * t503 - t480;
t344 = -mrSges(6,1) * t385 + mrSges(6,3) * t376 + Ifges(6,4) * t399 + Ifges(6,2) * t398 + Ifges(6,6) * t462 - pkin(5) * t356 - t403 * t429 + t405 * t465 - t481;
t343 = mrSges(6,2) * t385 - mrSges(6,3) * t375 + Ifges(6,1) * t399 + Ifges(6,4) * t398 + Ifges(6,5) * t462 - pkin(9) * t356 - t359 * t471 + t360 * t475 + t403 * t428 - t404 * t465;
t339 = mrSges(5,2) * t415 - mrSges(5,3) * t389 + Ifges(5,1) * t436 + Ifges(5,4) * t435 + Ifges(5,5) * qJDD(4) - pkin(8) * t350 - qJD(4) * t426 + t343 * t476 - t344 * t472 + t425 * t443;
t338 = -mrSges(5,1) * t415 + mrSges(5,3) * t390 + Ifges(5,4) * t436 + Ifges(5,2) * t435 + Ifges(5,6) * qJDD(4) - pkin(4) * t485 + pkin(8) * t492 + qJD(4) * t427 + t472 * t343 + t476 * t344 - t444 * t425;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t495 - mrSges(2,2) * t490 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t437 - mrSges(3,2) * t438 + t467 * (mrSges(4,2) * t424 - mrSges(4,3) * t416 + t477 * t339 - t473 * t338 - pkin(7) * t501 + (Ifges(4,1) * t467 + Ifges(4,4) * t469) * qJDD(1)) + t469 * (-mrSges(4,1) * t424 + mrSges(4,3) * t417 + t473 * t339 + t477 * t338 - pkin(3) * t482 + pkin(7) * t493 + (Ifges(4,4) * t467 + Ifges(4,2) * t469) * qJDD(1)) - pkin(2) * t351 + qJ(3) * t494 + pkin(1) * (t468 * (m(3) * t438 - mrSges(3,1) * t479 - qJDD(1) * mrSges(3,2) + t494) + t470 * (t480 + m(3) * t437 - mrSges(3,2) * t479 + (mrSges(3,1) - t503) * qJDD(1))); m(3) * t466 + t341 * t469 + t342 * t467; t351; mrSges(5,1) * t389 - mrSges(5,2) * t390 + Ifges(5,5) * t436 + Ifges(5,6) * t435 + Ifges(5,3) * qJDD(4) + pkin(4) * t350 + t426 * t444 - t427 * t443 + t483; t483; t481;];
tauJ  = t1;

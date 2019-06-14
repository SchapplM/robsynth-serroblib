% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:43:01
% EndTime: 2019-05-05 18:43:05
% DurationCPUTime: 2.38s
% Computational Cost: add. (14412->275), mult. (28634->333), div. (0->0), fcn. (16530->10), ass. (0->120)
t485 = sin(qJ(3));
t489 = cos(qJ(3));
t521 = Ifges(4,4) + Ifges(5,6);
t533 = t485 * t521 + t489 * (Ifges(4,2) + Ifges(5,3));
t532 = t485 * (Ifges(4,1) + Ifges(5,2)) + t489 * t521;
t531 = -2 * qJD(4);
t520 = (Ifges(4,5) - Ifges(5,4));
t519 = (Ifges(4,6) - Ifges(5,5));
t482 = cos(pkin(10));
t528 = pkin(1) * t482 + pkin(2);
t486 = sin(qJ(1));
t490 = cos(qJ(1));
t509 = t486 * g(1) - g(2) * t490;
t454 = qJDD(1) * pkin(1) + t509;
t492 = qJD(1) ^ 2;
t506 = -g(1) * t490 - g(2) * t486;
t458 = -pkin(1) * t492 + t506;
t481 = sin(pkin(10));
t425 = t481 * t454 + t482 * t458;
t412 = -pkin(2) * t492 + qJDD(1) * pkin(7) + t425;
t480 = -g(3) + qJDD(2);
t406 = t489 * t412 + t485 * t480;
t455 = (-pkin(3) * t489 - qJ(4) * t485) * qJD(1);
t491 = qJD(3) ^ 2;
t514 = qJD(1) * t489;
t401 = pkin(3) * t491 - qJDD(3) * qJ(4) + (qJD(3) * t531) - t455 * t514 - t406;
t527 = t485 * (t533 * qJD(1) + (t519 * qJD(3))) - t489 * (t532 * qJD(1) + (t520 * qJD(3)));
t526 = -pkin(3) - pkin(8);
t524 = pkin(7) * t492;
t523 = pkin(8) * t492;
t518 = t480 * t489;
t513 = qJD(1) * qJD(3);
t510 = t485 * t513;
t460 = qJDD(1) * t489 - t510;
t474 = t485 * qJD(1);
t467 = pkin(4) * t474 - qJD(3) * pkin(8);
t479 = t489 ^ 2;
t511 = t489 * t513;
t459 = qJDD(1) * t485 + t511;
t424 = t454 * t482 - t481 * t458;
t503 = -qJDD(1) * pkin(2) - t424;
t496 = pkin(3) * t510 + t474 * t531 + (-t459 - t511) * qJ(4) + t503;
t384 = -t467 * t474 + (-pkin(4) * t479 - pkin(7)) * t492 + t526 * t460 + t496;
t409 = t485 * t412;
t504 = -qJ(4) * t491 + t455 * t474 + qJDD(4) + t409;
t395 = pkin(4) * t459 + t526 * qJDD(3) + (-pkin(4) * t513 - t485 * t523 - t480) * t489 + t504;
t484 = sin(qJ(5));
t488 = cos(qJ(5));
t379 = -t384 * t484 + t488 * t395;
t452 = -qJD(3) * t484 - t488 * t514;
t421 = qJD(5) * t452 + qJDD(3) * t488 - t460 * t484;
t451 = qJDD(5) + t459;
t453 = qJD(3) * t488 - t484 * t514;
t470 = t474 + qJD(5);
t376 = (t452 * t470 - t421) * pkin(9) + (t452 * t453 + t451) * pkin(5) + t379;
t380 = t488 * t384 + t484 * t395;
t420 = -qJD(5) * t453 - qJDD(3) * t484 - t460 * t488;
t429 = pkin(5) * t470 - pkin(9) * t453;
t450 = t452 ^ 2;
t377 = -pkin(5) * t450 + pkin(9) * t420 - t429 * t470 + t380;
t483 = sin(qJ(6));
t487 = cos(qJ(6));
t374 = t376 * t487 - t377 * t483;
t422 = t452 * t487 - t453 * t483;
t390 = qJD(6) * t422 + t420 * t483 + t421 * t487;
t423 = t452 * t483 + t453 * t487;
t403 = -mrSges(7,1) * t422 + mrSges(7,2) * t423;
t468 = qJD(6) + t470;
t407 = -mrSges(7,2) * t468 + mrSges(7,3) * t422;
t444 = qJDD(6) + t451;
t371 = m(7) * t374 + mrSges(7,1) * t444 - mrSges(7,3) * t390 - t403 * t423 + t407 * t468;
t375 = t376 * t483 + t377 * t487;
t389 = -qJD(6) * t423 + t420 * t487 - t421 * t483;
t408 = mrSges(7,1) * t468 - mrSges(7,3) * t423;
t372 = m(7) * t375 - mrSges(7,2) * t444 + mrSges(7,3) * t389 + t403 * t422 - t408 * t468;
t363 = t487 * t371 + t483 * t372;
t426 = -mrSges(6,1) * t452 + mrSges(6,2) * t453;
t427 = -mrSges(6,2) * t470 + mrSges(6,3) * t452;
t360 = m(6) * t379 + mrSges(6,1) * t451 - mrSges(6,3) * t421 - t426 * t453 + t427 * t470 + t363;
t428 = mrSges(6,1) * t470 - mrSges(6,3) * t453;
t507 = -t371 * t483 + t487 * t372;
t361 = m(6) * t380 - mrSges(6,2) * t451 + mrSges(6,3) * t420 + t426 * t452 - t428 * t470 + t507;
t517 = -t484 * t360 + t488 * t361;
t405 = -t409 + t518;
t456 = (mrSges(5,2) * t489 - mrSges(5,3) * t485) * qJD(1);
t457 = (-mrSges(4,1) * t489 + mrSges(4,2) * t485) * qJD(1);
t464 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t514;
t465 = -mrSges(5,1) * t514 - qJD(3) * mrSges(5,3);
t357 = t360 * t488 + t361 * t484;
t402 = -qJDD(3) * pkin(3) + t504 - t518;
t499 = -m(5) * t402 - t459 * mrSges(5,1) - t357;
t354 = m(4) * t405 - mrSges(4,3) * t459 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t464 - t465) * qJD(3) + (-t456 - t457) * t474 + t499;
t463 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t474;
t466 = mrSges(5,1) * t474 + qJD(3) * mrSges(5,2);
t394 = pkin(4) * t460 + qJD(3) * t467 - t479 * t523 - t401;
t382 = -pkin(5) * t420 - pkin(9) * t450 + t429 * t453 + t394;
t500 = m(7) * t382 - mrSges(7,1) * t389 + t390 * mrSges(7,2) - t422 * t407 + t423 * t408;
t495 = -m(6) * t394 + mrSges(6,1) * t420 - t421 * mrSges(6,2) + t427 * t452 - t453 * t428 - t500;
t493 = -m(5) * t401 + qJDD(3) * mrSges(5,3) + qJD(3) * t466 + t456 * t514 - t495;
t367 = -qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t460 + m(4) * t406 - qJD(3) * t463 + t457 * t514 + t493;
t508 = -t354 * t485 + t489 * t367;
t396 = -pkin(3) * t460 + t496 - t524;
t502 = m(5) * t396 + t460 * mrSges(5,2) - t466 * t474 + t517;
t398 = Ifges(7,4) * t423 + Ifges(7,2) * t422 + Ifges(7,6) * t468;
t399 = Ifges(7,1) * t423 + Ifges(7,4) * t422 + Ifges(7,5) * t468;
t498 = mrSges(7,1) * t374 - mrSges(7,2) * t375 + Ifges(7,5) * t390 + Ifges(7,6) * t389 + Ifges(7,3) * t444 + t423 * t398 - t422 * t399;
t411 = t503 - t524;
t497 = -m(4) * t411 + t460 * mrSges(4,1) + t464 * t514 - t502;
t414 = Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t470;
t415 = Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t470;
t494 = mrSges(6,1) * t379 - mrSges(6,2) * t380 + Ifges(6,5) * t421 + Ifges(6,6) * t420 + Ifges(6,3) * t451 + pkin(5) * t363 + t453 * t414 - t452 * t415 + t498;
t413 = Ifges(6,5) * t453 + Ifges(6,6) * t452 + Ifges(6,3) * t470;
t397 = Ifges(7,5) * t423 + Ifges(7,6) * t422 + Ifges(7,3) * t468;
t365 = mrSges(7,2) * t382 - mrSges(7,3) * t374 + Ifges(7,1) * t390 + Ifges(7,4) * t389 + Ifges(7,5) * t444 + t397 * t422 - t398 * t468;
t364 = -mrSges(7,1) * t382 + mrSges(7,3) * t375 + Ifges(7,4) * t390 + Ifges(7,2) * t389 + Ifges(7,6) * t444 - t397 * t423 + t399 * t468;
t356 = qJDD(3) * mrSges(5,2) + qJD(3) * t465 + t456 * t474 - t499;
t355 = -mrSges(5,3) * t459 + t465 * t514 + t502;
t353 = mrSges(6,2) * t394 - mrSges(6,3) * t379 + Ifges(6,1) * t421 + Ifges(6,4) * t420 + Ifges(6,5) * t451 - pkin(9) * t363 - t364 * t483 + t365 * t487 + t413 * t452 - t414 * t470;
t352 = -mrSges(6,1) * t394 + mrSges(6,3) * t380 + Ifges(6,4) * t421 + Ifges(6,2) * t420 + Ifges(6,6) * t451 - pkin(5) * t500 + pkin(9) * t507 + t487 * t364 + t483 * t365 - t453 * t413 + t470 * t415;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t509 - mrSges(2,2) * t506 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t424 - mrSges(3,2) * t425 + t485 * (mrSges(5,1) * t402 + mrSges(4,2) * t411 - mrSges(4,3) * t405 - mrSges(5,3) * t396 + pkin(4) * t357 - qJ(4) * t355 + t494) + t489 * (-mrSges(4,1) * t411 - mrSges(5,1) * t401 + mrSges(5,2) * t396 + mrSges(4,3) * t406 - pkin(3) * t355 - pkin(4) * t495 - pkin(8) * t517 - t488 * t352 - t484 * t353) + pkin(2) * t497 + pkin(7) * t508 + pkin(1) * (t481 * (m(3) * t425 - mrSges(3,1) * t492 - qJDD(1) * mrSges(3,2) + t508) + t482 * (m(3) * t424 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t492 + t497)) + t533 * t460 + (t485 * t520 + t489 * t519) * qJDD(3) - t527 * qJD(3) + (t528 * (-mrSges(4,2) + mrSges(5,3)) + t532) * t459 + t528 * qJD(1) * (-t463 * t485 - t465 * t489); m(3) * t480 + t354 * t489 + t367 * t485; mrSges(4,1) * t405 - mrSges(4,2) * t406 + mrSges(5,2) * t402 - mrSges(5,3) * t401 + t488 * t353 - t484 * t352 - pkin(8) * t357 - pkin(3) * t356 + qJ(4) * t493 + (mrSges(5,1) * qJ(4) + t519) * t460 + t520 * t459 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + t527 * qJD(1); t356; t494; t498;];
tauJ  = t1;

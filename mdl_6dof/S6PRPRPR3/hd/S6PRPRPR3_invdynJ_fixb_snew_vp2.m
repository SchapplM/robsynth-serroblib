% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 22:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:33:16
% EndTime: 2019-05-04 22:33:18
% DurationCPUTime: 1.30s
% Computational Cost: add. (5393->218), mult. (10081->265), div. (0->0), fcn. (6338->12), ass. (0->106)
t437 = sin(qJ(4));
t440 = cos(qJ(4));
t473 = Ifges(5,4) + Ifges(6,6);
t484 = t437 * t473 + t440 * (Ifges(5,2) + Ifges(6,3));
t483 = t437 * (Ifges(5,1) + Ifges(6,2)) + t440 * t473;
t472 = Ifges(5,5) - Ifges(6,4);
t471 = Ifges(5,6) - Ifges(6,5);
t431 = sin(pkin(10));
t434 = cos(pkin(10));
t414 = -g(1) * t434 - g(2) * t431;
t429 = -g(3) + qJDD(1);
t438 = sin(qJ(2));
t441 = cos(qJ(2));
t432 = sin(pkin(6));
t467 = t432 * t441;
t413 = g(1) * t431 - g(2) * t434;
t435 = cos(pkin(6));
t469 = t413 * t435;
t370 = -t414 * t438 + t429 * t467 + t441 * t469;
t368 = qJDD(2) * pkin(2) + t370;
t468 = t432 * t438;
t371 = t441 * t414 + t429 * t468 + t438 * t469;
t443 = qJD(2) ^ 2;
t369 = -pkin(2) * t443 + t371;
t430 = sin(pkin(11));
t433 = cos(pkin(11));
t362 = t368 * t433 - t430 * t369;
t453 = -qJDD(2) * pkin(3) - t362;
t476 = pkin(8) * t443;
t360 = t453 - t476;
t460 = qJD(2) * qJD(4);
t458 = t440 * t460;
t410 = qJDD(2) * t437 + t458;
t457 = t437 * t460;
t411 = qJDD(2) * t440 - t457;
t462 = qJD(2) * t437;
t415 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t462;
t461 = qJD(2) * t440;
t416 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t461;
t417 = -mrSges(6,1) * t461 - qJD(4) * mrSges(6,3);
t478 = -2 * qJD(5);
t444 = pkin(4) * t457 + t462 * t478 + (-t410 - t458) * qJ(5) + t453;
t355 = -pkin(4) * t411 + t444 - t476;
t418 = mrSges(6,1) * t462 + qJD(4) * mrSges(6,2);
t455 = -t413 * t432 + t435 * t429;
t386 = qJDD(3) + t455;
t363 = t430 * t368 + t433 * t369;
t361 = -pkin(3) * t443 + qJDD(2) * pkin(8) + t363;
t358 = t437 * t361;
t407 = (-t440 * pkin(4) - qJ(5) * t437) * qJD(2);
t442 = qJD(4) ^ 2;
t454 = -qJ(5) * t442 + t407 * t462 + qJDD(5) + t358;
t475 = pkin(9) * t443;
t477 = -pkin(4) - pkin(9);
t351 = pkin(5) * t410 + t477 * qJDD(4) + (-pkin(5) * t460 - t437 * t475 - t386) * t440 + t454;
t421 = pkin(5) * t462 - qJD(4) * pkin(9);
t428 = t440 ^ 2;
t352 = -t421 * t462 + (-pkin(5) * t428 - pkin(8)) * t443 + t477 * t411 + t444;
t436 = sin(qJ(6));
t439 = cos(qJ(6));
t347 = t351 * t439 - t352 * t436;
t405 = -qJD(4) * t436 - t439 * t461;
t380 = qJD(6) * t405 + qJDD(4) * t439 - t411 * t436;
t406 = qJD(4) * t439 - t436 * t461;
t381 = -mrSges(7,1) * t405 + mrSges(7,2) * t406;
t423 = qJD(6) + t462;
t384 = -mrSges(7,2) * t423 + mrSges(7,3) * t405;
t403 = qJDD(6) + t410;
t344 = m(7) * t347 + mrSges(7,1) * t403 - mrSges(7,3) * t380 - t381 * t406 + t384 * t423;
t348 = t351 * t436 + t352 * t439;
t379 = -qJD(6) * t406 - qJDD(4) * t436 - t411 * t439;
t385 = mrSges(7,1) * t423 - mrSges(7,3) * t406;
t345 = m(7) * t348 - mrSges(7,2) * t403 + mrSges(7,3) * t379 + t381 * t405 - t385 * t423;
t465 = -t436 * t344 + t439 * t345;
t452 = m(6) * t355 + t411 * mrSges(6,2) - t418 * t462 + t465;
t480 = -m(5) * t360 + t411 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t410 + t416 * t461 + (-t415 * t437 - t417 * t440) * qJD(2) - t452;
t479 = (t484 * qJD(2) + t471 * qJD(4)) * t437 - t440 * (t483 * qJD(2) + t472 * qJD(4));
t470 = t386 * t440;
t356 = -t358 + t470;
t408 = (mrSges(6,2) * t440 - mrSges(6,3) * t437) * qJD(2);
t409 = (-mrSges(5,1) * t440 + mrSges(5,2) * t437) * qJD(2);
t336 = t344 * t439 + t345 * t436;
t354 = -qJDD(4) * pkin(4) + t454 - t470;
t449 = -m(6) * t354 - t410 * mrSges(6,1) - t336;
t333 = m(5) * t356 - mrSges(5,3) * t410 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t416 - t417) * qJD(4) + (-t408 - t409) * t462 + t449;
t357 = t440 * t361 + t437 * t386;
t446 = -pkin(4) * t442 + qJDD(4) * qJ(5) + t407 * t461 + t357;
t353 = qJD(4) * t478 - t446;
t350 = -t428 * t475 + pkin(5) * t411 + ((2 * qJD(5)) + t421) * qJD(4) + t446;
t450 = -m(7) * t350 + t379 * mrSges(7,1) - t380 * mrSges(7,2) + t405 * t384 - t406 * t385;
t445 = -m(6) * t353 + qJDD(4) * mrSges(6,3) + qJD(4) * t418 + t408 * t461 - t450;
t341 = t409 * t461 + m(5) * t357 - qJDD(4) * mrSges(5,2) - qJD(4) * t415 + (mrSges(5,3) + mrSges(6,1)) * t411 + t445;
t456 = -t333 * t437 + t440 * t341;
t329 = m(4) * t363 - mrSges(4,1) * t443 - qJDD(2) * mrSges(4,2) + t456;
t331 = m(4) * t362 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t443 + t480;
t466 = t430 * t329 + t433 * t331;
t459 = m(4) * t386 + t440 * t333 + t437 * t341;
t373 = Ifges(7,4) * t406 + Ifges(7,2) * t405 + Ifges(7,6) * t423;
t374 = Ifges(7,1) * t406 + Ifges(7,4) * t405 + Ifges(7,5) * t423;
t448 = mrSges(7,1) * t347 - mrSges(7,2) * t348 + Ifges(7,5) * t380 + Ifges(7,6) * t379 + Ifges(7,3) * t403 + t406 * t373 - t405 * t374;
t372 = Ifges(7,5) * t406 + Ifges(7,6) * t405 + Ifges(7,3) * t423;
t338 = mrSges(7,2) * t350 - mrSges(7,3) * t347 + Ifges(7,1) * t380 + Ifges(7,4) * t379 + Ifges(7,5) * t403 + t372 * t405 - t373 * t423;
t337 = -mrSges(7,1) * t350 + mrSges(7,3) * t348 + Ifges(7,4) * t380 + Ifges(7,2) * t379 + Ifges(7,6) * t403 - t372 * t406 + t374 * t423;
t335 = qJDD(4) * mrSges(6,2) + qJD(4) * t417 + t408 * t462 - t449;
t334 = -mrSges(6,3) * t410 + t417 * t461 + t452;
t1 = [m(2) * t429 + (m(3) * t371 - mrSges(3,1) * t443 - qJDD(2) * mrSges(3,2) + t329 * t433 - t331 * t430) * t468 + (m(3) * t370 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t443 + t466) * t467 + t435 * (m(3) * t455 + t459); mrSges(3,1) * t370 - mrSges(3,2) * t371 + mrSges(4,1) * t362 - mrSges(4,2) * t363 + t437 * (mrSges(6,1) * t354 + mrSges(5,2) * t360 - mrSges(5,3) * t356 - mrSges(6,3) * t355 + pkin(5) * t336 - qJ(5) * t334 + t448) + t440 * (-mrSges(5,1) * t360 - mrSges(6,1) * t353 + mrSges(6,2) * t355 + mrSges(5,3) * t357 - pkin(4) * t334 - pkin(5) * t450 - pkin(9) * t465 - t439 * t337 - t436 * t338) + pkin(8) * t456 + pkin(2) * t466 + t484 * t411 + (t437 * t472 + t440 * t471) * qJDD(4) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) - t479 * qJD(4) + t483 * t410 + t480 * pkin(3); t459; mrSges(5,1) * t356 - mrSges(5,2) * t357 + mrSges(6,2) * t354 - mrSges(6,3) * t353 + t439 * t338 - t436 * t337 - pkin(9) * t336 - pkin(4) * t335 + qJ(5) * t445 + (qJ(5) * mrSges(6,1) + t471) * t411 + t472 * t410 + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + t479 * qJD(2); t335; t448;];
tauJ  = t1;

% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:50:10
% EndTime: 2019-05-05 15:50:12
% DurationCPUTime: 1.30s
% Computational Cost: add. (8578->232), mult. (16624->289), div. (0->0), fcn. (9779->8), ass. (0->91)
t464 = 2 * qJD(1);
t439 = qJD(1) ^ 2;
t434 = sin(qJ(1));
t438 = cos(qJ(1));
t460 = t434 * g(1) - t438 * g(2);
t402 = -qJDD(1) * pkin(1) - (t439 * qJ(2)) + qJDD(2) - t460;
t445 = qJDD(1) * qJ(3) + (qJD(3) * t464) - t402;
t392 = -(t439 * pkin(7)) + t445;
t433 = sin(qJ(4));
t437 = cos(qJ(4));
t457 = qJD(1) * qJD(4);
t411 = -t433 * qJDD(1) - t437 * t457;
t454 = t433 * t457;
t412 = t437 * qJDD(1) - t454;
t459 = qJD(1) * t433;
t413 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t459;
t458 = t437 * qJD(1);
t414 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t458;
t451 = -t438 * g(1) - t434 * g(2);
t462 = qJDD(1) * qJ(2) + (qJD(2) * t464) + t451;
t398 = qJDD(3) + (-pkin(1) - qJ(3)) * t439 + t462;
t393 = -qJDD(1) * pkin(7) + t398;
t386 = t433 * g(3) + t437 * t393;
t368 = (-t412 - t454) * pkin(8) + (-t433 * t437 * t439 + qJDD(4)) * pkin(4) + t386;
t387 = -t437 * g(3) + t433 * t393;
t415 = qJD(4) * pkin(4) - pkin(8) * t458;
t429 = t433 ^ 2;
t369 = -t429 * t439 * pkin(4) + t411 * pkin(8) - qJD(4) * t415 + t387;
t432 = sin(qJ(5));
t436 = cos(qJ(5));
t358 = t432 * t368 + t436 * t369;
t407 = (t437 * t432 + t433 * t436) * qJD(1);
t408 = (-t433 * t432 + t437 * t436) * qJD(1);
t389 = t407 * pkin(5) - t408 * pkin(9);
t424 = qJD(4) + qJD(5);
t422 = t424 ^ 2;
t423 = qJDD(4) + qJDD(5);
t354 = -t422 * pkin(5) + t423 * pkin(9) - t407 * t389 + t358;
t371 = -t411 * pkin(4) + t415 * t458 + (-pkin(8) * t429 - pkin(7)) * t439 + t445;
t378 = -t408 * qJD(5) + t436 * t411 - t432 * t412;
t379 = -t407 * qJD(5) + t432 * t411 + t436 * t412;
t355 = (t407 * t424 - t379) * pkin(9) + (t408 * t424 - t378) * pkin(5) + t371;
t431 = sin(qJ(6));
t435 = cos(qJ(6));
t351 = -t431 * t354 + t435 * t355;
t394 = -t431 * t408 + t435 * t424;
t361 = t394 * qJD(6) + t435 * t379 + t431 * t423;
t395 = t435 * t408 + t431 * t424;
t372 = -t394 * mrSges(7,1) + t395 * mrSges(7,2);
t377 = qJDD(6) - t378;
t403 = qJD(6) + t407;
t380 = -t403 * mrSges(7,2) + t394 * mrSges(7,3);
t348 = m(7) * t351 + t377 * mrSges(7,1) - t361 * mrSges(7,3) - t395 * t372 + t403 * t380;
t352 = t435 * t354 + t431 * t355;
t360 = -t395 * qJD(6) - t431 * t379 + t435 * t423;
t381 = t403 * mrSges(7,1) - t395 * mrSges(7,3);
t349 = m(7) * t352 - t377 * mrSges(7,2) + t360 * mrSges(7,3) + t394 * t372 - t403 * t381;
t338 = t435 * t348 + t431 * t349;
t399 = -(t424 * mrSges(6,2)) - t407 * mrSges(6,3);
t400 = (t424 * mrSges(6,1)) - t408 * mrSges(6,3);
t443 = -m(6) * t371 + t378 * mrSges(6,1) - t379 * mrSges(6,2) - t407 * t399 - t408 * t400 - t338;
t463 = -m(4) * t445 - m(5) * t392 + t411 * mrSges(5,1) - t412 * mrSges(5,2) + t443 + (-t413 * t433 - t414 * t437) * qJD(1);
t388 = t407 * mrSges(6,1) + t408 * mrSges(6,2);
t452 = -t431 * t348 + t435 * t349;
t336 = m(6) * t358 - t423 * mrSges(6,2) + t378 * mrSges(6,3) - t407 * t388 - t424 * t400 + t452;
t357 = t436 * t368 - t432 * t369;
t353 = -t423 * pkin(5) - t422 * pkin(9) + t408 * t389 - t357;
t444 = -m(7) * t353 + t360 * mrSges(7,1) - t361 * mrSges(7,2) + t394 * t380 - t395 * t381;
t344 = m(6) * t357 + t423 * mrSges(6,1) - t379 * mrSges(6,3) - t408 * t388 + t424 * t399 + t444;
t332 = t432 * t336 + t436 * t344;
t410 = (t433 * mrSges(5,1) + t437 * mrSges(5,2)) * qJD(1);
t453 = t436 * t336 - t432 * t344;
t461 = t433 * (m(5) * t387 - qJDD(4) * mrSges(5,2) + t411 * mrSges(5,3) - qJD(4) * t414 - t410 * t459 + t453) + t437 * (m(5) * t386 + qJDD(4) * mrSges(5,1) - t412 * mrSges(5,3) + qJD(4) * t413 - t410 * t458 + t332);
t447 = m(4) * t398 + qJDD(1) * mrSges(4,2) - (t439 * mrSges(4,3)) + t461;
t362 = Ifges(7,5) * t395 + Ifges(7,6) * t394 + Ifges(7,3) * t403;
t364 = Ifges(7,1) * t395 + Ifges(7,4) * t394 + Ifges(7,5) * t403;
t341 = -mrSges(7,1) * t353 + mrSges(7,3) * t352 + Ifges(7,4) * t361 + Ifges(7,2) * t360 + Ifges(7,6) * t377 - t395 * t362 + t403 * t364;
t363 = Ifges(7,4) * t395 + Ifges(7,2) * t394 + Ifges(7,6) * t403;
t342 = mrSges(7,2) * t353 - mrSges(7,3) * t351 + Ifges(7,1) * t361 + Ifges(7,4) * t360 + Ifges(7,5) * t377 + t394 * t362 - t403 * t363;
t383 = Ifges(6,4) * t408 - Ifges(6,2) * t407 + (Ifges(6,6) * t424);
t384 = Ifges(6,1) * t408 - Ifges(6,4) * t407 + (Ifges(6,5) * t424);
t442 = mrSges(6,1) * t357 - mrSges(6,2) * t358 + Ifges(6,5) * t379 + Ifges(6,6) * t378 + Ifges(6,3) * t423 + pkin(5) * t444 + pkin(9) * t452 + t435 * t341 + t431 * t342 + t408 * t383 + t407 * t384;
t441 = mrSges(7,1) * t351 - mrSges(7,2) * t352 + Ifges(7,5) * t361 + Ifges(7,6) * t360 + Ifges(7,3) * t377 + t395 * t363 - t394 * t364;
t406 = (Ifges(5,5) * qJD(4)) + (t437 * Ifges(5,1) - t433 * Ifges(5,4)) * qJD(1);
t405 = (Ifges(5,6) * qJD(4)) + (t437 * Ifges(5,4) - t433 * Ifges(5,2)) * qJD(1);
t401 = t439 * pkin(1) - t462;
t382 = Ifges(6,5) * t408 - Ifges(6,6) * t407 + Ifges(6,3) * t424;
t333 = (-mrSges(4,2) - mrSges(3,3)) * t439 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + m(3) * t402 + t463;
t329 = -mrSges(6,1) * t371 + mrSges(6,3) * t358 + Ifges(6,4) * t379 + Ifges(6,2) * t378 + Ifges(6,6) * t423 - pkin(5) * t338 - t408 * t382 + t424 * t384 - t441;
t328 = mrSges(6,2) * t371 - mrSges(6,3) * t357 + Ifges(6,1) * t379 + Ifges(6,4) * t378 + Ifges(6,5) * t423 - pkin(9) * t338 - t431 * t341 + t435 * t342 - t407 * t382 - t424 * t383;
t1 = [-pkin(1) * t333 + mrSges(2,1) * t460 - mrSges(2,2) * t451 + qJ(2) * (-m(3) * t401 + (t439 * mrSges(3,2)) + t447) + mrSges(3,2) * t402 - mrSges(3,3) * t401 + t437 * (mrSges(5,2) * t392 - mrSges(5,3) * t386 + Ifges(5,1) * t412 + Ifges(5,4) * t411 + Ifges(5,5) * qJDD(4) - pkin(8) * t332 - qJD(4) * t405 + t436 * t328 - t432 * t329) - t433 * (-mrSges(5,1) * t392 + mrSges(5,3) * t387 + Ifges(5,4) * t412 + Ifges(5,2) * t411 + Ifges(5,6) * qJDD(4) + pkin(4) * t443 + pkin(8) * t453 + qJD(4) * t406 + t432 * t328 + t436 * t329) - pkin(7) * t461 + mrSges(4,2) * t398 + mrSges(4,3) * t445 + (qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (t439 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t463) * qJ(3); t333; t447; t442 + pkin(4) * t332 + Ifges(5,3) * qJDD(4) + (t437 * t405 + t433 * t406) * qJD(1) + mrSges(5,1) * t386 - mrSges(5,2) * t387 + Ifges(5,6) * t411 + Ifges(5,5) * t412; t442; t441;];
tauJ  = t1;

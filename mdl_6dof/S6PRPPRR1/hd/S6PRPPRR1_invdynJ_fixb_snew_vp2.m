% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-05-04 21:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:40:14
% EndTime: 2019-05-04 21:40:16
% DurationCPUTime: 1.40s
% Computational Cost: add. (11914->208), mult. (24220->269), div. (0->0), fcn. (17870->14), ass. (0->107)
t427 = cos(pkin(12));
t421 = t427 ^ 2;
t466 = 0.2e1 * t427;
t438 = qJD(2) ^ 2;
t423 = sin(pkin(12));
t432 = sin(qJ(5));
t435 = cos(qJ(5));
t444 = t423 * t432 - t427 * t435;
t403 = t444 * qJD(2);
t445 = t423 * t435 + t427 * t432;
t404 = t445 * qJD(2);
t454 = t404 * qJD(5);
t392 = -t444 * qJDD(2) - t454;
t465 = pkin(4) * t427;
t464 = mrSges(5,2) * t423;
t425 = sin(pkin(10));
t429 = cos(pkin(10));
t410 = g(1) * t425 - g(2) * t429;
t430 = cos(pkin(6));
t463 = t410 * t430;
t462 = t421 * t438;
t426 = sin(pkin(6));
t433 = sin(qJ(2));
t461 = t426 * t433;
t436 = cos(qJ(2));
t460 = t426 * t436;
t411 = -g(1) * t429 - g(2) * t425;
t422 = -g(3) + qJDD(1);
t385 = -t411 * t433 + t422 * t460 + t436 * t463;
t380 = qJDD(2) * pkin(2) + t385;
t386 = t436 * t411 + t422 * t461 + t433 * t463;
t381 = -pkin(2) * t438 + t386;
t424 = sin(pkin(11));
t428 = cos(pkin(11));
t366 = t424 * t380 + t428 * t381;
t364 = -pkin(3) * t438 + qJDD(2) * qJ(4) + t366;
t447 = -t410 * t426 + t430 * t422;
t399 = qJDD(3) + t447;
t453 = qJD(2) * qJD(4);
t457 = t427 * t399 - 0.2e1 * t423 * t453;
t359 = -t364 * t423 + t457;
t443 = mrSges(5,3) * qJDD(2) + t438 * (-mrSges(5,1) * t427 + t464);
t357 = (-pkin(8) * qJDD(2) + t438 * t465 - t364) * t423 + t457;
t360 = t427 * t364 + t423 * t399 + t453 * t466;
t452 = qJDD(2) * t427;
t358 = -pkin(4) * t462 + pkin(8) * t452 + t360;
t353 = t432 * t357 + t435 * t358;
t388 = mrSges(6,1) * t403 + mrSges(6,2) * t404;
t401 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t404;
t391 = pkin(5) * t403 - pkin(9) * t404;
t437 = qJD(5) ^ 2;
t351 = -pkin(5) * t437 + qJDD(5) * pkin(9) - t391 * t403 + t353;
t420 = t423 ^ 2;
t365 = t428 * t380 - t424 * t381;
t446 = qJDD(4) - t365;
t361 = (-pkin(3) - t465) * qJDD(2) + (-qJ(4) + (-t420 - t421) * pkin(8)) * t438 + t446;
t455 = t403 * qJD(5);
t393 = t445 * qJDD(2) - t455;
t354 = (-t393 + t455) * pkin(9) + (-t392 + t454) * pkin(5) + t361;
t431 = sin(qJ(6));
t434 = cos(qJ(6));
t348 = -t351 * t431 + t354 * t434;
t396 = qJD(5) * t434 - t404 * t431;
t373 = qJD(6) * t396 + qJDD(5) * t431 + t393 * t434;
t397 = qJD(5) * t431 + t404 * t434;
t374 = -mrSges(7,1) * t396 + mrSges(7,2) * t397;
t402 = qJD(6) + t403;
t375 = -mrSges(7,2) * t402 + mrSges(7,3) * t396;
t390 = qJDD(6) - t392;
t346 = m(7) * t348 + mrSges(7,1) * t390 - mrSges(7,3) * t373 - t374 * t397 + t375 * t402;
t349 = t351 * t434 + t354 * t431;
t372 = -qJD(6) * t397 + qJDD(5) * t434 - t393 * t431;
t376 = mrSges(7,1) * t402 - mrSges(7,3) * t397;
t347 = m(7) * t349 - mrSges(7,2) * t390 + mrSges(7,3) * t372 + t374 * t396 - t376 * t402;
t448 = -t346 * t431 + t434 * t347;
t336 = m(6) * t353 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t392 - qJD(5) * t401 - t388 * t403 + t448;
t352 = t357 * t435 - t358 * t432;
t400 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t403;
t350 = -qJDD(5) * pkin(5) - pkin(9) * t437 + t391 * t404 - t352;
t442 = -m(7) * t350 + t372 * mrSges(7,1) - mrSges(7,2) * t373 + t396 * t375 - t376 * t397;
t342 = m(6) * t352 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t393 + qJD(5) * t400 - t388 * t404 + t442;
t458 = t432 * t336 + t435 * t342;
t330 = m(5) * t359 - t443 * t423 + t458;
t449 = t435 * t336 - t432 * t342;
t331 = m(5) * t360 + t443 * t427 + t449;
t450 = -t330 * t423 + t427 * t331;
t324 = m(4) * t366 - mrSges(4,1) * t438 - qJDD(2) * mrSges(4,2) + t450;
t363 = -qJDD(2) * pkin(3) - t438 * qJ(4) + t446;
t338 = t434 * t346 + t431 * t347;
t441 = m(6) * t361 - t392 * mrSges(6,1) + t393 * mrSges(6,2) + t403 * t400 + t404 * t401 + t338;
t440 = -m(5) * t363 + mrSges(5,1) * t452 - t441 + (t420 * t438 + t462) * mrSges(5,3);
t333 = (mrSges(4,1) - t464) * qJDD(2) + t440 - t438 * mrSges(4,2) + m(4) * t365;
t459 = t424 * t324 + t428 * t333;
t451 = m(4) * t399 + t427 * t330 + t423 * t331;
t368 = Ifges(7,4) * t397 + Ifges(7,2) * t396 + Ifges(7,6) * t402;
t369 = Ifges(7,1) * t397 + Ifges(7,4) * t396 + Ifges(7,5) * t402;
t439 = mrSges(7,1) * t348 - mrSges(7,2) * t349 + Ifges(7,5) * t373 + Ifges(7,6) * t372 + Ifges(7,3) * t390 + t368 * t397 - t369 * t396;
t384 = Ifges(6,1) * t404 - Ifges(6,4) * t403 + Ifges(6,5) * qJD(5);
t383 = Ifges(6,4) * t404 - Ifges(6,2) * t403 + Ifges(6,6) * qJD(5);
t382 = Ifges(6,5) * t404 - Ifges(6,6) * t403 + Ifges(6,3) * qJD(5);
t367 = Ifges(7,5) * t397 + Ifges(7,6) * t396 + Ifges(7,3) * t402;
t340 = mrSges(7,2) * t350 - mrSges(7,3) * t348 + Ifges(7,1) * t373 + Ifges(7,4) * t372 + Ifges(7,5) * t390 + t367 * t396 - t368 * t402;
t339 = -mrSges(7,1) * t350 + mrSges(7,3) * t349 + Ifges(7,4) * t373 + Ifges(7,2) * t372 + Ifges(7,6) * t390 - t367 * t397 + t369 * t402;
t337 = qJDD(2) * t464 - t440;
t326 = -mrSges(6,1) * t361 + mrSges(6,3) * t353 + Ifges(6,4) * t393 + Ifges(6,2) * t392 + Ifges(6,6) * qJDD(5) - pkin(5) * t338 + qJD(5) * t384 - t382 * t404 - t439;
t325 = mrSges(6,2) * t361 - mrSges(6,3) * t352 + Ifges(6,1) * t393 + Ifges(6,4) * t392 + Ifges(6,5) * qJDD(5) - pkin(9) * t338 - qJD(5) * t383 - t339 * t431 + t340 * t434 - t382 * t403;
t1 = [m(2) * t422 + (m(3) * t386 - mrSges(3,1) * t438 - qJDD(2) * mrSges(3,2) + t324 * t428 - t333 * t424) * t461 + (m(3) * t385 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t438 + t459) * t460 + t430 * (m(3) * t447 + t451); mrSges(3,1) * t385 - mrSges(3,2) * t386 + mrSges(4,1) * t365 - mrSges(4,2) * t366 + t423 * (mrSges(5,2) * t363 - mrSges(5,3) * t359 - pkin(8) * t458 + t435 * t325 - t432 * t326) + t427 * (-mrSges(5,1) * t363 + mrSges(5,3) * t360 - pkin(4) * t441 + pkin(8) * t449 + t432 * t325 + t435 * t326) - pkin(3) * t337 + qJ(4) * t450 + pkin(2) * t459 + (Ifges(5,2) * t421 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t423 + Ifges(5,4) * t466) * t423) * qJDD(2); t451; t337; mrSges(6,1) * t352 - mrSges(6,2) * t353 + Ifges(6,5) * t393 + Ifges(6,6) * t392 + Ifges(6,3) * qJDD(5) + pkin(5) * t442 + pkin(9) * t448 + t434 * t339 + t431 * t340 + t404 * t383 + t403 * t384; t439;];
tauJ  = t1;

% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:00:12
% EndTime: 2019-12-31 20:00:14
% DurationCPUTime: 1.62s
% Computational Cost: add. (10293->251), mult. (23186->313), div. (0->0), fcn. (15227->8), ass. (0->100)
t440 = -2 * qJD(3);
t439 = Ifges(5,1) + Ifges(6,1);
t432 = Ifges(5,4) - Ifges(6,5);
t431 = -Ifges(5,5) - Ifges(6,4);
t438 = Ifges(5,2) + Ifges(6,3);
t430 = Ifges(5,6) - Ifges(6,6);
t437 = -Ifges(5,3) - Ifges(6,2);
t403 = sin(qJ(2));
t405 = cos(qJ(2));
t420 = qJD(1) * qJD(2);
t391 = t403 * qJDD(1) + t405 * t420;
t408 = qJD(1) ^ 2;
t404 = sin(qJ(1));
t406 = cos(qJ(1));
t413 = -t406 * g(1) - t404 * g(2);
t388 = -t408 * pkin(1) + qJDD(1) * pkin(6) + t413;
t428 = t403 * t388;
t434 = pkin(2) * t408;
t350 = qJDD(2) * pkin(2) - t391 * qJ(3) - t428 + (qJ(3) * t420 + t403 * t434 - g(3)) * t405;
t374 = -t403 * g(3) + t405 * t388;
t392 = t405 * qJDD(1) - t403 * t420;
t423 = qJD(1) * t403;
t393 = qJD(2) * pkin(2) - qJ(3) * t423;
t399 = t405 ^ 2;
t351 = t392 * qJ(3) - qJD(2) * t393 - t399 * t434 + t374;
t400 = sin(pkin(8));
t401 = cos(pkin(8));
t383 = (t405 * t400 + t403 * t401) * qJD(1);
t331 = t401 * t350 - t400 * t351 + t383 * t440;
t382 = (t403 * t400 - t405 * t401) * qJD(1);
t370 = t401 * t391 + t400 * t392;
t402 = sin(qJ(4));
t435 = cos(qJ(4));
t371 = -t435 * qJD(2) + t402 * t383;
t344 = -t371 * qJD(4) + t402 * qJDD(2) + t435 * t370;
t372 = t402 * qJD(2) + t435 * t383;
t353 = t371 * mrSges(6,1) - t372 * mrSges(6,3);
t332 = t400 * t350 + t401 * t351 + t382 * t440;
t365 = t382 * pkin(3) - t383 * pkin(7);
t407 = qJD(2) ^ 2;
t328 = -t407 * pkin(3) + qJDD(2) * pkin(7) - t382 * t365 + t332;
t418 = t404 * g(1) - t406 * g(2);
t411 = -qJDD(1) * pkin(1) - t418;
t356 = -t392 * pkin(2) + qJDD(3) + t393 * t423 + (-qJ(3) * t399 - pkin(6)) * t408 + t411;
t369 = -t400 * t391 + t401 * t392;
t330 = (qJD(2) * t382 - t370) * pkin(7) + (qJD(2) * t383 - t369) * pkin(3) + t356;
t324 = -t402 * t328 + t435 * t330;
t352 = t371 * pkin(4) - t372 * qJ(5);
t368 = qJDD(4) - t369;
t381 = qJD(4) + t382;
t380 = t381 ^ 2;
t322 = -t368 * pkin(4) - t380 * qJ(5) + t372 * t352 + qJDD(5) - t324;
t357 = -t371 * mrSges(6,2) + t381 * mrSges(6,3);
t414 = -m(6) * t322 + t368 * mrSges(6,1) + t381 * t357;
t318 = t344 * mrSges(6,2) + t372 * t353 - t414;
t325 = t435 * t328 + t402 * t330;
t321 = -t380 * pkin(4) + t368 * qJ(5) + 0.2e1 * qJD(5) * t381 - t371 * t352 + t325;
t343 = t372 * qJD(4) - t435 * qJDD(2) + t402 * t370;
t360 = -t381 * mrSges(6,1) + t372 * mrSges(6,2);
t419 = m(6) * t321 + t368 * mrSges(6,3) + t381 * t360;
t425 = t432 * t371 - t439 * t372 + t431 * t381;
t426 = t438 * t371 - t432 * t372 - t430 * t381;
t436 = -t430 * t343 - t431 * t344 - t437 * t368 - t425 * t371 - t426 * t372 + mrSges(5,1) * t324 - mrSges(6,1) * t322 - mrSges(5,2) * t325 + mrSges(6,3) * t321 - pkin(4) * t318 + qJ(5) * (-t343 * mrSges(6,2) - t371 * t353 + t419);
t433 = -mrSges(5,3) - mrSges(6,2);
t364 = t382 * mrSges(4,1) + t383 * mrSges(4,2);
t376 = qJD(2) * mrSges(4,1) - t383 * mrSges(4,3);
t359 = t381 * mrSges(5,1) - t372 * mrSges(5,3);
t424 = -t371 * mrSges(5,1) - t372 * mrSges(5,2) - t353;
t314 = m(5) * t325 - t368 * mrSges(5,2) + t433 * t343 - t381 * t359 + t424 * t371 + t419;
t358 = -t381 * mrSges(5,2) - t371 * mrSges(5,3);
t316 = m(5) * t324 + t368 * mrSges(5,1) + t433 * t344 + t381 * t358 + t424 * t372 + t414;
t416 = t435 * t314 - t402 * t316;
t306 = m(4) * t332 - qJDD(2) * mrSges(4,2) + t369 * mrSges(4,3) - qJD(2) * t376 - t382 * t364 + t416;
t375 = -qJD(2) * mrSges(4,2) - t382 * mrSges(4,3);
t327 = -qJDD(2) * pkin(3) - t407 * pkin(7) + t383 * t365 - t331;
t323 = -0.2e1 * qJD(5) * t372 + (t371 * t381 - t344) * qJ(5) + (t372 * t381 + t343) * pkin(4) + t327;
t319 = m(6) * t323 + t343 * mrSges(6,1) - t344 * mrSges(6,3) + t371 * t357 - t372 * t360;
t409 = -m(5) * t327 - t343 * mrSges(5,1) - t344 * mrSges(5,2) - t371 * t358 - t372 * t359 - t319;
t311 = m(4) * t331 + qJDD(2) * mrSges(4,1) - t370 * mrSges(4,3) + qJD(2) * t375 - t383 * t364 + t409;
t302 = t400 * t306 + t401 * t311;
t309 = t402 * t314 + t435 * t316;
t427 = t430 * t371 + t431 * t372 + t437 * t381;
t422 = qJD(1) * t405;
t417 = t401 * t306 - t400 * t311;
t308 = m(4) * t356 - t369 * mrSges(4,1) + t370 * mrSges(4,2) + t382 * t375 + t383 * t376 + t309;
t395 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t422;
t394 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t423;
t390 = (-t405 * mrSges(3,1) + t403 * mrSges(3,2)) * qJD(1);
t387 = -t408 * pkin(6) + t411;
t386 = Ifges(3,5) * qJD(2) + (t403 * Ifges(3,1) + t405 * Ifges(3,4)) * qJD(1);
t385 = Ifges(3,6) * qJD(2) + (t403 * Ifges(3,4) + t405 * Ifges(3,2)) * qJD(1);
t373 = -t405 * g(3) - t428;
t363 = Ifges(4,1) * t383 - Ifges(4,4) * t382 + Ifges(4,5) * qJD(2);
t362 = Ifges(4,4) * t383 - Ifges(4,2) * t382 + Ifges(4,6) * qJD(2);
t361 = Ifges(4,5) * t383 - Ifges(4,6) * t382 + Ifges(4,3) * qJD(2);
t307 = mrSges(5,2) * t327 + mrSges(6,2) * t322 - mrSges(5,3) * t324 - mrSges(6,3) * t323 - qJ(5) * t319 - t432 * t343 + t439 * t344 - t431 * t368 + t427 * t371 + t426 * t381;
t303 = -mrSges(5,1) * t327 - mrSges(6,1) * t323 + mrSges(6,2) * t321 + mrSges(5,3) * t325 - pkin(4) * t319 - t438 * t343 + t432 * t344 + t430 * t368 + t427 * t372 - t425 * t381;
t301 = -mrSges(4,1) * t356 + mrSges(4,3) * t332 + Ifges(4,4) * t370 + Ifges(4,2) * t369 + Ifges(4,6) * qJDD(2) - pkin(3) * t309 + qJD(2) * t363 - t383 * t361 - t436;
t300 = mrSges(4,2) * t356 - mrSges(4,3) * t331 + Ifges(4,1) * t370 + Ifges(4,4) * t369 + Ifges(4,5) * qJDD(2) - pkin(7) * t309 - qJD(2) * t362 - t402 * t303 + t435 * t307 - t382 * t361;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t418 - mrSges(2,2) * t413 + t403 * (mrSges(3,2) * t387 - mrSges(3,3) * t373 + Ifges(3,1) * t391 + Ifges(3,4) * t392 + Ifges(3,5) * qJDD(2) - qJ(3) * t302 - qJD(2) * t385 + t401 * t300 - t400 * t301) + t405 * (-mrSges(3,1) * t387 + mrSges(3,3) * t374 + Ifges(3,4) * t391 + Ifges(3,2) * t392 + Ifges(3,6) * qJDD(2) - pkin(2) * t308 + qJ(3) * t417 + qJD(2) * t386 + t400 * t300 + t401 * t301) + pkin(1) * (-m(3) * t387 + t392 * mrSges(3,1) - t391 * mrSges(3,2) + (-t394 * t403 + t395 * t405) * qJD(1) - t308) + pkin(6) * (t405 * (m(3) * t374 - qJDD(2) * mrSges(3,2) + t392 * mrSges(3,3) - qJD(2) * t394 + t390 * t422 + t417) - t403 * (m(3) * t373 + qJDD(2) * mrSges(3,1) - t391 * mrSges(3,3) + qJD(2) * t395 - t390 * t423 + t302)); Ifges(3,5) * t391 + Ifges(3,6) * t392 + mrSges(3,1) * t373 - mrSges(3,2) * t374 + Ifges(4,5) * t370 + Ifges(4,6) * t369 + t383 * t362 + t382 * t363 + mrSges(4,1) * t331 - mrSges(4,2) * t332 + t402 * t307 + t435 * t303 + pkin(3) * t409 + pkin(7) * t416 + pkin(2) * t302 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t403 * t385 - t405 * t386) * qJD(1); t308; t436; t318;];
tauJ = t1;

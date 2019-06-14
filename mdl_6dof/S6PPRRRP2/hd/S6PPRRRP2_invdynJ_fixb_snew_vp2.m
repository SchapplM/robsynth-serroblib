% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:32:26
% EndTime: 2019-05-04 20:32:27
% DurationCPUTime: 1.40s
% Computational Cost: add. (10020->204), mult. (17897->257), div. (0->0), fcn. (13636->14), ass. (0->98)
t446 = Ifges(6,1) + Ifges(7,1);
t435 = Ifges(6,4) - Ifges(7,5);
t441 = -Ifges(6,5) - Ifges(7,4);
t445 = Ifges(6,2) + Ifges(7,3);
t433 = Ifges(6,6) - Ifges(7,6);
t405 = sin(qJ(5));
t406 = sin(qJ(4));
t426 = qJD(3) * t406;
t437 = cos(qJ(5));
t382 = -qJD(4) * t437 + t405 * t426;
t408 = cos(qJ(4));
t424 = qJD(3) * qJD(4);
t421 = t408 * t424;
t386 = qJDD(3) * t406 + t421;
t361 = -t382 * qJD(5) + t405 * qJDD(4) + t386 * t437;
t383 = t405 * qJD(4) + t426 * t437;
t365 = mrSges(7,1) * t382 - mrSges(7,3) * t383;
t398 = sin(pkin(11));
t402 = cos(pkin(11));
t389 = -g(1) * t402 - g(2) * t398;
t397 = sin(pkin(12));
t401 = cos(pkin(12));
t388 = g(1) * t398 - g(2) * t402;
t396 = -g(3) + qJDD(1);
t400 = sin(pkin(6));
t404 = cos(pkin(6));
t415 = t388 * t404 + t396 * t400;
t349 = t401 * t389 + t397 * t415;
t407 = sin(qJ(3));
t409 = cos(qJ(3));
t348 = -t397 * t389 + t401 * t415;
t372 = -t388 * t400 + t396 * t404 + qJDD(2);
t399 = sin(pkin(7));
t403 = cos(pkin(7));
t443 = t348 * t403 + t372 * t399;
t342 = t409 * t349 + t443 * t407;
t411 = qJD(3) ^ 2;
t340 = -pkin(3) * t411 + qJDD(3) * pkin(9) + t342;
t344 = -t348 * t399 + t372 * t403;
t334 = t408 * t340 + t406 * t344;
t385 = (-pkin(4) * t408 - pkin(10) * t406) * qJD(3);
t410 = qJD(4) ^ 2;
t425 = qJD(3) * t408;
t332 = -pkin(4) * t410 + qJDD(4) * pkin(10) + t385 * t425 + t334;
t341 = -t407 * t349 + t409 * t443;
t339 = -qJDD(3) * pkin(3) - t411 * pkin(9) - t341;
t422 = t406 * t424;
t387 = qJDD(3) * t408 - t422;
t336 = (-t386 - t421) * pkin(10) + (-t387 + t422) * pkin(4) + t339;
t327 = -t405 * t332 + t336 * t437;
t364 = pkin(5) * t382 - qJ(6) * t383;
t381 = qJDD(5) - t387;
t393 = qJD(5) - t425;
t392 = t393 ^ 2;
t326 = -t381 * pkin(5) - t392 * qJ(6) + t383 * t364 + qJDD(6) - t327;
t371 = -mrSges(7,2) * t382 + mrSges(7,3) * t393;
t418 = -m(7) * t326 + t381 * mrSges(7,1) + t393 * t371;
t322 = t361 * mrSges(7,2) + t383 * t365 - t418;
t328 = t437 * t332 + t405 * t336;
t325 = -pkin(5) * t392 + qJ(6) * t381 + 0.2e1 * qJD(6) * t393 - t364 * t382 + t328;
t360 = qJD(5) * t383 - qJDD(4) * t437 + t386 * t405;
t370 = -mrSges(7,1) * t393 + mrSges(7,2) * t383;
t423 = m(7) * t325 + t381 * mrSges(7,3) + t393 * t370;
t429 = t445 * t382 - t435 * t383 - t433 * t393;
t438 = t435 * t382 - t446 * t383 + t441 * t393;
t440 = -Ifges(6,3) - Ifges(7,2);
t444 = -t441 * t361 - t438 * t382 - t433 * t360 - t440 * t381 + mrSges(6,1) * t327 - mrSges(7,1) * t326 - mrSges(6,2) * t328 + mrSges(7,3) * t325 - pkin(5) * t322 + qJ(6) * (-t360 * mrSges(7,2) - t382 * t365 + t423) - t429 * t383;
t436 = -mrSges(6,3) - mrSges(7,2);
t430 = t433 * t382 + t441 * t383 + t440 * t393;
t427 = -mrSges(6,1) * t382 - mrSges(6,2) * t383 - t365;
t384 = (-mrSges(5,1) * t408 + mrSges(5,2) * t406) * qJD(3);
t390 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t426;
t369 = mrSges(6,1) * t393 - mrSges(6,3) * t383;
t319 = m(6) * t328 - t381 * mrSges(6,2) + t360 * t436 - t393 * t369 + t382 * t427 + t423;
t368 = -mrSges(6,2) * t393 - mrSges(6,3) * t382;
t320 = m(6) * t327 + t381 * mrSges(6,1) + t361 * t436 + t393 * t368 + t383 * t427 + t418;
t419 = t437 * t319 - t320 * t405;
t313 = m(5) * t334 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t387 - qJD(4) * t390 + t384 * t425 + t419;
t333 = -t406 * t340 + t408 * t344;
t391 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t425;
t331 = -qJDD(4) * pkin(4) - t410 * pkin(10) + t385 * t426 - t333;
t329 = -0.2e1 * qJD(6) * t383 + (t382 * t393 - t361) * qJ(6) + (t383 * t393 + t360) * pkin(5) + t331;
t323 = m(7) * t329 + mrSges(7,1) * t360 - t361 * mrSges(7,3) - t383 * t370 + t371 * t382;
t412 = -m(6) * t331 - t360 * mrSges(6,1) - mrSges(6,2) * t361 - t382 * t368 - t369 * t383 - t323;
t317 = m(5) * t333 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t386 + qJD(4) * t391 - t384 * t426 + t412;
t420 = t408 * t313 - t317 * t406;
t309 = m(4) * t342 - mrSges(4,1) * t411 - qJDD(3) * mrSges(4,2) + t420;
t316 = t405 * t319 + t437 * t320;
t413 = -m(5) * t339 + t387 * mrSges(5,1) - t386 * mrSges(5,2) - t390 * t426 + t391 * t425 - t316;
t311 = m(4) * t341 + qJDD(3) * mrSges(4,1) - t411 * mrSges(4,2) + t413;
t417 = t309 * t407 + t311 * t409;
t376 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t406 + Ifges(5,4) * t408) * qJD(3);
t375 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t406 + Ifges(5,2) * t408) * qJD(3);
t315 = mrSges(6,2) * t331 + mrSges(7,2) * t326 - mrSges(6,3) * t327 - mrSges(7,3) * t329 - qJ(6) * t323 - t435 * t360 + t446 * t361 - t441 * t381 + t430 * t382 + t429 * t393;
t314 = -mrSges(6,1) * t331 - mrSges(7,1) * t329 + mrSges(7,2) * t325 + mrSges(6,3) * t328 - pkin(5) * t323 - t445 * t360 + t435 * t361 + t433 * t381 + t430 * t383 - t438 * t393;
t310 = m(4) * t344 + t313 * t406 + t317 * t408;
t308 = m(3) * t372 + t403 * t310 + t417 * t399;
t1 = [m(2) * t396 + t404 * t308 + (t397 * (m(3) * t349 + t309 * t409 - t311 * t407) + t401 * (m(3) * t348 - t399 * t310 + t403 * t417)) * t400; t308; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t341 - mrSges(4,2) * t342 + t406 * (mrSges(5,2) * t339 - mrSges(5,3) * t333 + Ifges(5,1) * t386 + Ifges(5,4) * t387 + Ifges(5,5) * qJDD(4) - pkin(10) * t316 - qJD(4) * t375 - t405 * t314 + t315 * t437) + t408 * (-mrSges(5,1) * t339 + mrSges(5,3) * t334 + Ifges(5,4) * t386 + Ifges(5,2) * t387 + Ifges(5,6) * qJDD(4) - pkin(4) * t316 + qJD(4) * t376 - t444) + pkin(3) * t413 + pkin(9) * t420; Ifges(5,5) * t386 + Ifges(5,6) * t387 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t333 - mrSges(5,2) * t334 + t405 * t315 + t437 * t314 + pkin(4) * t412 + pkin(10) * t419 + (t375 * t406 - t376 * t408) * qJD(3); t444; t322;];
tauJ  = t1;

% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRRRP1
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
% Datum: 2019-05-04 20:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:27:23
% EndTime: 2019-05-04 20:27:25
% DurationCPUTime: 1.44s
% Computational Cost: add. (10120->206), mult. (18195->257), div. (0->0), fcn. (13875->14), ass. (0->99)
t447 = Ifges(6,1) + Ifges(7,1);
t438 = Ifges(6,4) + Ifges(7,4);
t437 = Ifges(6,5) + Ifges(7,5);
t446 = Ifges(6,2) + Ifges(7,2);
t436 = Ifges(6,6) + Ifges(7,6);
t406 = sin(qJ(5));
t409 = cos(qJ(5));
t407 = sin(qJ(4));
t429 = qJD(3) * t407;
t386 = qJD(4) * t409 - t406 * t429;
t410 = cos(qJ(4));
t427 = qJD(3) * qJD(4);
t423 = t410 * t427;
t390 = qJDD(3) * t407 + t423;
t366 = qJD(5) * t386 + qJDD(4) * t406 + t390 * t409;
t387 = qJD(4) * t406 + t409 * t429;
t368 = -mrSges(7,1) * t386 + mrSges(7,2) * t387;
t399 = sin(pkin(11));
t403 = cos(pkin(11));
t393 = -g(1) * t403 - g(2) * t399;
t398 = sin(pkin(12));
t402 = cos(pkin(12));
t392 = g(1) * t399 - g(2) * t403;
t397 = -g(3) + qJDD(1);
t401 = sin(pkin(6));
t405 = cos(pkin(6));
t417 = t392 * t405 + t397 * t401;
t353 = t393 * t402 + t398 * t417;
t408 = sin(qJ(3));
t411 = cos(qJ(3));
t352 = -t393 * t398 + t402 * t417;
t376 = -t392 * t401 + t397 * t405 + qJDD(2);
t400 = sin(pkin(7));
t404 = cos(pkin(7));
t444 = t352 * t404 + t376 * t400;
t345 = t411 * t353 + t444 * t408;
t413 = qJD(3) ^ 2;
t343 = -pkin(3) * t413 + qJDD(3) * pkin(9) + t345;
t347 = -t352 * t400 + t376 * t404;
t336 = t410 * t343 + t407 * t347;
t389 = (-pkin(4) * t410 - pkin(10) * t407) * qJD(3);
t412 = qJD(4) ^ 2;
t428 = qJD(3) * t410;
t334 = -pkin(4) * t412 + qJDD(4) * pkin(10) + t389 * t428 + t336;
t344 = -t408 * t353 + t411 * t444;
t342 = -qJDD(3) * pkin(3) - pkin(9) * t413 - t344;
t424 = t407 * t427;
t391 = qJDD(3) * t410 - t424;
t339 = (-t390 - t423) * pkin(10) + (-t391 + t424) * pkin(4) + t342;
t329 = -t334 * t406 + t409 * t339;
t385 = qJDD(5) - t391;
t396 = qJD(5) - t428;
t325 = -0.2e1 * qJD(6) * t387 + (t386 * t396 - t366) * qJ(6) + (t386 * t387 + t385) * pkin(5) + t329;
t371 = -mrSges(7,2) * t396 + mrSges(7,3) * t386;
t426 = m(7) * t325 + t385 * mrSges(7,1) + t396 * t371;
t323 = -mrSges(7,3) * t366 - t368 * t387 + t426;
t330 = t409 * t334 + t406 * t339;
t365 = -qJD(5) * t387 + qJDD(4) * t409 - t390 * t406;
t373 = pkin(5) * t396 - qJ(6) * t387;
t384 = t386 ^ 2;
t327 = -pkin(5) * t384 + qJ(6) * t365 + 0.2e1 * qJD(6) * t386 - t373 * t396 + t330;
t432 = -t386 * t446 - t387 * t438 - t396 * t436;
t440 = t438 * t386 + t387 * t447 + t437 * t396;
t442 = Ifges(6,3) + Ifges(7,3);
t445 = mrSges(6,1) * t329 + mrSges(7,1) * t325 - mrSges(6,2) * t330 - mrSges(7,2) * t327 + pkin(5) * t323 + t436 * t365 + t437 * t366 + t442 * t385 - t440 * t386 - t432 * t387;
t439 = -mrSges(6,2) - mrSges(7,2);
t433 = -t436 * t386 - t437 * t387 - t442 * t396;
t374 = mrSges(7,1) * t396 - mrSges(7,3) * t387;
t430 = -mrSges(6,1) * t396 + mrSges(6,3) * t387 - t374;
t425 = m(7) * t327 + t365 * mrSges(7,3) + t386 * t368;
t388 = (-mrSges(5,1) * t410 + mrSges(5,2) * t407) * qJD(3);
t394 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t429;
t369 = -mrSges(6,1) * t386 + mrSges(6,2) * t387;
t372 = -mrSges(6,2) * t396 + mrSges(6,3) * t386;
t319 = m(6) * t329 + mrSges(6,1) * t385 + t372 * t396 + (-t368 - t369) * t387 + (-mrSges(6,3) - mrSges(7,3)) * t366 + t426;
t321 = m(6) * t330 + mrSges(6,3) * t365 + t369 * t386 + t385 * t439 + t430 * t396 + t425;
t421 = -t319 * t406 + t409 * t321;
t316 = m(5) * t336 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t391 - qJD(4) * t394 + t388 * t428 + t421;
t335 = -t407 * t343 + t347 * t410;
t395 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t428;
t333 = -qJDD(4) * pkin(4) - pkin(10) * t412 + t389 * t429 - t335;
t331 = -pkin(5) * t365 - qJ(6) * t384 + t373 * t387 + qJDD(6) + t333;
t420 = -m(7) * t331 + t365 * mrSges(7,1) + t386 * t371;
t414 = -m(6) * t333 + t365 * mrSges(6,1) + t366 * t439 + t386 * t372 + t430 * t387 + t420;
t322 = m(5) * t335 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t390 + qJD(4) * t395 - t388 * t429 + t414;
t422 = t410 * t316 - t322 * t407;
t311 = m(4) * t345 - mrSges(4,1) * t413 - qJDD(3) * mrSges(4,2) + t422;
t318 = t319 * t409 + t321 * t406;
t415 = -m(5) * t342 + t391 * mrSges(5,1) - mrSges(5,2) * t390 - t394 * t429 + t395 * t428 - t318;
t314 = m(4) * t344 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t413 + t415;
t419 = t311 * t408 + t314 * t411;
t380 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t407 + Ifges(5,4) * t410) * qJD(3);
t379 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t407 + Ifges(5,2) * t410) * qJD(3);
t328 = mrSges(7,2) * t366 + t374 * t387 - t420;
t317 = mrSges(6,2) * t333 + mrSges(7,2) * t331 - mrSges(6,3) * t329 - mrSges(7,3) * t325 - qJ(6) * t323 + t438 * t365 + t366 * t447 + t437 * t385 - t433 * t386 + t432 * t396;
t313 = -mrSges(6,1) * t333 + mrSges(6,3) * t330 - mrSges(7,1) * t331 + mrSges(7,3) * t327 - pkin(5) * t328 + qJ(6) * t425 + (-qJ(6) * t374 + t440) * t396 + t433 * t387 + (-mrSges(7,2) * qJ(6) + t436) * t385 + t438 * t366 + t446 * t365;
t312 = m(4) * t347 + t316 * t407 + t322 * t410;
t310 = m(3) * t376 + t312 * t404 + t419 * t400;
t1 = [m(2) * t397 + t405 * t310 + (t398 * (m(3) * t353 + t311 * t411 - t314 * t408) + t402 * (m(3) * t352 - t312 * t400 + t404 * t419)) * t401; t310; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t344 - mrSges(4,2) * t345 + t407 * (mrSges(5,2) * t342 - mrSges(5,3) * t335 + Ifges(5,1) * t390 + Ifges(5,4) * t391 + Ifges(5,5) * qJDD(4) - pkin(10) * t318 - qJD(4) * t379 - t313 * t406 + t317 * t409) + t410 * (-mrSges(5,1) * t342 + mrSges(5,3) * t336 + Ifges(5,4) * t390 + Ifges(5,2) * t391 + Ifges(5,6) * qJDD(4) - pkin(4) * t318 + qJD(4) * t380 - t445) + pkin(3) * t415 + pkin(9) * t422; Ifges(5,5) * t390 + Ifges(5,6) * t391 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t335 - mrSges(5,2) * t336 + t406 * t317 + t409 * t313 + pkin(4) * t414 + pkin(10) * t421 + (t379 * t407 - t380 * t410) * qJD(3); t445; t328;];
tauJ  = t1;

% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:09
% EndTime: 2019-12-31 18:32:11
% DurationCPUTime: 1.31s
% Computational Cost: add. (6857->229), mult. (16491->280), div. (0->0), fcn. (11174->8), ass. (0->104)
t449 = Ifges(4,1) + Ifges(5,2);
t442 = Ifges(4,4) + Ifges(5,6);
t441 = Ifges(4,5) - Ifges(5,4);
t448 = -Ifges(4,2) - Ifges(5,3);
t440 = Ifges(4,6) - Ifges(5,5);
t447 = Ifges(4,3) + Ifges(5,1);
t410 = qJD(1) ^ 2;
t446 = -2 * qJD(4);
t445 = cos(qJ(3));
t444 = pkin(2) * t410;
t443 = mrSges(4,1) - mrSges(5,2);
t439 = pkin(6) * qJDD(1);
t406 = sin(qJ(1));
t408 = cos(qJ(1));
t423 = -g(1) * t408 - g(2) * t406;
t391 = -pkin(1) * t410 + qJDD(1) * qJ(2) + t423;
t402 = sin(pkin(8));
t403 = cos(pkin(8));
t429 = qJD(1) * qJD(2);
t426 = -t403 * g(3) - 0.2e1 * t402 * t429;
t354 = (t403 * t444 - t391 - t439) * t402 + t426;
t374 = -g(3) * t402 + (t391 + 0.2e1 * t429) * t403;
t400 = t403 ^ 2;
t355 = -t400 * t444 + t403 * t439 + t374;
t405 = sin(qJ(3));
t335 = t354 * t445 - t405 * t355;
t428 = t403 * t445;
t432 = qJD(1) * t402;
t389 = -qJD(1) * t428 + t405 * t432;
t419 = t402 * t445 + t403 * t405;
t390 = t419 * qJD(1);
t365 = mrSges(4,1) * t389 + mrSges(4,2) * t390;
t430 = t389 * qJD(3);
t372 = qJDD(1) * t419 - t430;
t431 = qJD(3) * t390;
t371 = t431 + (t402 * t405 - t428) * qJDD(1);
t382 = pkin(4) * t390 - qJD(3) * pkin(7);
t388 = t389 ^ 2;
t427 = t406 * g(1) - t408 * g(2);
t422 = qJDD(2) - t427;
t433 = -t402 ^ 2 - t400;
t370 = (-pkin(2) * t403 - pkin(1)) * qJDD(1) + (pkin(6) * t433 - qJ(2)) * t410 + t422;
t411 = pkin(3) * t431 + t390 * t446 + (-t372 + t430) * qJ(4) + t370;
t326 = -pkin(4) * t388 - t382 * t390 + (pkin(3) + pkin(7)) * t371 + t411;
t364 = pkin(3) * t389 - qJ(4) * t390;
t409 = qJD(3) ^ 2;
t333 = -qJDD(3) * pkin(3) - t409 * qJ(4) + t390 * t364 + qJDD(4) - t335;
t327 = (t389 * t390 - qJDD(3)) * pkin(7) + (t372 + t430) * pkin(4) + t333;
t404 = sin(qJ(5));
t407 = cos(qJ(5));
t324 = -t326 * t404 + t327 * t407;
t375 = -qJD(3) * t404 + t389 * t407;
t345 = qJD(5) * t375 + qJDD(3) * t407 + t371 * t404;
t376 = qJD(3) * t407 + t389 * t404;
t346 = -mrSges(6,1) * t375 + mrSges(6,2) * t376;
t386 = qJD(5) + t390;
t350 = -mrSges(6,2) * t386 + mrSges(6,3) * t375;
t369 = qJDD(5) + t372;
t321 = m(6) * t324 + mrSges(6,1) * t369 - mrSges(6,3) * t345 - t346 * t376 + t350 * t386;
t325 = t326 * t407 + t327 * t404;
t344 = -qJD(5) * t376 - qJDD(3) * t404 + t371 * t407;
t351 = mrSges(6,1) * t386 - mrSges(6,3) * t376;
t322 = m(6) * t325 - mrSges(6,2) * t369 + mrSges(6,3) * t344 + t346 * t375 - t351 * t386;
t314 = t407 * t321 + t404 * t322;
t366 = -mrSges(5,2) * t389 - mrSges(5,3) * t390;
t416 = -m(5) * t333 - t372 * mrSges(5,1) - t390 * t366 - t314;
t380 = mrSges(5,1) * t389 - qJD(3) * mrSges(5,3);
t434 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t389 - t380;
t311 = m(4) * t335 - t372 * mrSges(4,3) + qJD(3) * t434 + qJDD(3) * t443 - t390 * t365 + t416;
t336 = t405 * t354 + t445 * t355;
t379 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t390;
t415 = -t409 * pkin(3) + qJDD(3) * qJ(4) - t389 * t364 + t336;
t332 = qJD(3) * t446 - t415;
t381 = mrSges(5,1) * t390 + qJD(3) * mrSges(5,2);
t329 = -t371 * pkin(4) - t388 * pkin(7) + ((2 * qJD(4)) + t382) * qJD(3) + t415;
t417 = -m(6) * t329 + mrSges(6,1) * t344 - t345 * mrSges(6,2) + t350 * t375 - t376 * t351;
t413 = -m(5) * t332 + qJDD(3) * mrSges(5,3) + qJD(3) * t381 - t417;
t319 = m(4) * t336 - qJDD(3) * mrSges(4,2) - qJD(3) * t379 + (-t365 - t366) * t389 + (-mrSges(4,3) - mrSges(5,1)) * t371 + t413;
t438 = t311 * t445 + t319 * t405;
t437 = -qJD(3) * t447 + t389 * t440 - t390 * t441;
t436 = qJD(3) * t440 + t389 * t448 + t390 * t442;
t435 = qJD(3) * t441 - t389 * t442 + t390 * t449;
t425 = -t405 * t311 + t319 * t445;
t424 = -t404 * t321 + t322 * t407;
t421 = -mrSges(3,1) * t403 + mrSges(3,2) * t402;
t420 = mrSges(3,3) * qJDD(1) + t410 * t421;
t331 = t371 * pkin(3) + t411;
t418 = m(5) * t331 - t372 * mrSges(5,3) - t390 * t381 + t424;
t338 = Ifges(6,4) * t376 + Ifges(6,2) * t375 + Ifges(6,6) * t386;
t339 = Ifges(6,1) * t376 + Ifges(6,4) * t375 + Ifges(6,5) * t386;
t414 = mrSges(6,1) * t324 - mrSges(6,2) * t325 + Ifges(6,5) * t345 + Ifges(6,6) * t344 + Ifges(6,3) * t369 + t376 * t338 - t375 * t339;
t412 = m(4) * t370 + t372 * mrSges(4,2) + t371 * t443 + t390 * t379 + t389 * t434 + t418;
t393 = (Ifges(3,5) * t402 + Ifges(3,6) * t403) * qJD(1);
t387 = -qJDD(1) * pkin(1) - t410 * qJ(2) + t422;
t373 = -t402 * t391 + t426;
t337 = Ifges(6,5) * t376 + Ifges(6,6) * t375 + Ifges(6,3) * t386;
t316 = mrSges(6,2) * t329 - mrSges(6,3) * t324 + Ifges(6,1) * t345 + Ifges(6,4) * t344 + Ifges(6,5) * t369 + t337 * t375 - t338 * t386;
t315 = -mrSges(6,1) * t329 + mrSges(6,3) * t325 + Ifges(6,4) * t345 + Ifges(6,2) * t344 + Ifges(6,6) * t369 - t337 * t376 + t339 * t386;
t313 = qJDD(3) * mrSges(5,2) + qJD(3) * t380 - t416;
t312 = -t371 * mrSges(5,2) - t389 * t380 + t418;
t309 = mrSges(3,3) * t410 * t433 + m(3) * t387 + qJDD(1) * t421 + t412;
t308 = mrSges(5,1) * t333 + mrSges(4,2) * t370 - mrSges(4,3) * t335 - mrSges(5,3) * t331 + pkin(4) * t314 - qJ(4) * t312 - t436 * qJD(3) + t441 * qJDD(3) - t442 * t371 + t372 * t449 + t437 * t389 + t414;
t307 = -mrSges(4,1) * t370 - mrSges(5,1) * t332 + mrSges(5,2) * t331 + mrSges(4,3) * t336 - pkin(3) * t312 - pkin(4) * t417 - pkin(7) * t424 + t435 * qJD(3) + t440 * qJDD(3) - t407 * t315 - t404 * t316 + t371 * t448 + t442 * t372 + t437 * t390;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t427 - mrSges(2,2) * t423 + t402 * (t403 * qJD(1) * t393 + mrSges(3,2) * t387 - mrSges(3,3) * t373 + t445 * t308 - t405 * t307 - pkin(6) * t438 + (Ifges(3,1) * t402 + Ifges(3,4) * t403) * qJDD(1)) + t403 * (-t393 * t432 - mrSges(3,1) * t387 + mrSges(3,3) * t374 + t405 * t308 + t445 * t307 - pkin(2) * t412 + pkin(6) * t425 + (Ifges(3,4) * t402 + Ifges(3,2) * t403) * qJDD(1)) - pkin(1) * t309 + qJ(2) * ((m(3) * t374 + t403 * t420 + t425) * t403 + (-m(3) * t373 + t402 * t420 - t438) * t402); t309; mrSges(4,1) * t335 - mrSges(4,2) * t336 + mrSges(5,2) * t333 - mrSges(5,3) * t332 + t407 * t316 - t404 * t315 - pkin(7) * t314 - pkin(3) * t313 + qJ(4) * t413 + t436 * t390 + (-qJ(4) * t366 + t435) * t389 + t441 * t372 + (-mrSges(5,1) * qJ(4) - t440) * t371 + t447 * qJDD(3); t313; t414;];
tauJ = t1;

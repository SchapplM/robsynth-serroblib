% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:52
% EndTime: 2019-12-31 18:48:54
% DurationCPUTime: 1.72s
% Computational Cost: add. (10437->221), mult. (25386->273), div. (0->0), fcn. (18395->8), ass. (0->97)
t448 = Ifges(5,1) + Ifges(6,1);
t438 = Ifges(5,4) - Ifges(6,5);
t446 = Ifges(6,4) + Ifges(5,5);
t447 = Ifges(5,2) + Ifges(6,3);
t445 = Ifges(5,6) - Ifges(6,6);
t444 = Ifges(5,3) + Ifges(6,2);
t407 = sin(pkin(8));
t408 = cos(pkin(8));
t410 = sin(qJ(3));
t412 = cos(qJ(3));
t419 = -t407 * t410 + t408 * t412;
t390 = t419 * qJD(1);
t420 = t407 * t412 + t408 * t410;
t391 = t420 * qJD(1);
t409 = sin(qJ(4));
t441 = cos(qJ(4));
t373 = -t441 * t390 + t409 * t391;
t374 = t409 * t390 + t441 * t391;
t406 = qJD(3) + qJD(4);
t443 = t447 * t373 - t438 * t374 - t445 * t406;
t442 = -t438 * t373 + t448 * t374 + t446 * t406;
t414 = qJD(1) ^ 2;
t440 = pkin(2) * t414;
t439 = -mrSges(5,3) - mrSges(6,2);
t437 = pkin(6) * qJDD(1);
t411 = sin(qJ(1));
t413 = cos(qJ(1));
t423 = -t413 * g(1) - t411 * g(2);
t392 = -t414 * pkin(1) + qJDD(1) * qJ(2) + t423;
t430 = qJD(1) * qJD(2);
t427 = -t408 * g(3) - 0.2e1 * t407 * t430;
t368 = (t408 * t440 - t392 - t437) * t407 + t427;
t382 = -t407 * g(3) + (t392 + 0.2e1 * t430) * t408;
t405 = t408 ^ 2;
t369 = -t405 * t440 + t408 * t437 + t382;
t346 = t412 * t368 - t410 * t369;
t431 = t390 * qJD(3);
t380 = t420 * qJDD(1) + t431;
t331 = (-t380 + t431) * pkin(7) + (t390 * t391 + qJDD(3)) * pkin(3) + t346;
t347 = t410 * t368 + t412 * t369;
t379 = -t391 * qJD(3) + t419 * qJDD(1);
t385 = qJD(3) * pkin(3) - t391 * pkin(7);
t389 = t390 ^ 2;
t333 = -t389 * pkin(3) + t379 * pkin(7) - qJD(3) * t385 + t347;
t329 = t409 * t331 + t441 * t333;
t344 = t374 * qJD(4) - t441 * t379 + t409 * t380;
t364 = t406 * mrSges(5,1) - t374 * mrSges(5,3);
t403 = qJDD(3) + qJDD(4);
t356 = t373 * pkin(4) - t374 * qJ(5);
t402 = t406 ^ 2;
t323 = -t402 * pkin(4) + t403 * qJ(5) + 0.2e1 * qJD(5) * t406 - t373 * t356 + t329;
t365 = -t406 * mrSges(6,1) + t374 * mrSges(6,2);
t429 = m(6) * t323 + t403 * mrSges(6,3) + t406 * t365;
t357 = t373 * mrSges(6,1) - t374 * mrSges(6,3);
t434 = -t373 * mrSges(5,1) - t374 * mrSges(5,2) - t357;
t314 = m(5) * t329 - t403 * mrSges(5,2) + t439 * t344 - t406 * t364 + t434 * t373 + t429;
t328 = t441 * t331 - t409 * t333;
t345 = -t373 * qJD(4) + t409 * t379 + t441 * t380;
t363 = -t406 * mrSges(5,2) - t373 * mrSges(5,3);
t324 = -t403 * pkin(4) - t402 * qJ(5) + t374 * t356 + qJDD(5) - t328;
t366 = -t373 * mrSges(6,2) + t406 * mrSges(6,3);
t424 = -m(6) * t324 + t403 * mrSges(6,1) + t406 * t366;
t316 = m(5) * t328 + t403 * mrSges(5,1) + t439 * t345 + t406 * t363 + t434 * t374 + t424;
t310 = t409 * t314 + t441 * t316;
t377 = -t390 * mrSges(4,1) + t391 * mrSges(4,2);
t383 = -qJD(3) * mrSges(4,2) + t390 * mrSges(4,3);
t308 = m(4) * t346 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t380 + qJD(3) * t383 - t377 * t391 + t310;
t384 = qJD(3) * mrSges(4,1) - t391 * mrSges(4,3);
t425 = t441 * t314 - t409 * t316;
t309 = m(4) * t347 - qJDD(3) * mrSges(4,2) + t379 * mrSges(4,3) - qJD(3) * t384 + t390 * t377 + t425;
t436 = t412 * t308 + t410 * t309;
t435 = t445 * t373 - t446 * t374 - t444 * t406;
t433 = -t407 ^ 2 - t405;
t428 = t411 * g(1) - t413 * g(2);
t426 = -t410 * t308 + t412 * t309;
t422 = qJDD(2) - t428;
t421 = -t408 * mrSges(3,1) + t407 * mrSges(3,2);
t418 = mrSges(3,3) * qJDD(1) + t414 * t421;
t378 = (-pkin(2) * t408 - pkin(1)) * qJDD(1) + (t433 * pkin(6) - qJ(2)) * t414 + t422;
t337 = -t379 * pkin(3) - t389 * pkin(7) + t391 * t385 + t378;
t326 = -0.2e1 * qJD(5) * t374 + (t373 * t406 - t345) * qJ(5) + (t374 * t406 + t344) * pkin(4) + t337;
t317 = m(6) * t326 + t344 * mrSges(6,1) - t345 * mrSges(6,3) - t374 * t365 + t373 * t366;
t417 = m(5) * t337 + t344 * mrSges(5,1) + t345 * mrSges(5,2) + t373 * t363 + t374 * t364 + t317;
t320 = t345 * mrSges(6,2) + t374 * t357 - t424;
t416 = mrSges(5,1) * t328 - mrSges(6,1) * t324 - mrSges(5,2) * t329 + mrSges(6,3) * t323 - pkin(4) * t320 + qJ(5) * t429 + t444 * t403 - t443 * t374 + (-qJ(5) * t357 + t442) * t373 + t446 * t345 + (-qJ(5) * mrSges(6,2) - t445) * t344;
t415 = m(4) * t378 - t379 * mrSges(4,1) + t380 * mrSges(4,2) - t390 * t383 + t391 * t384 + t417;
t388 = -qJDD(1) * pkin(1) - t414 * qJ(2) + t422;
t381 = -t407 * t392 + t427;
t372 = Ifges(4,1) * t391 + Ifges(4,4) * t390 + Ifges(4,5) * qJD(3);
t371 = Ifges(4,4) * t391 + Ifges(4,2) * t390 + Ifges(4,6) * qJD(3);
t370 = Ifges(4,5) * t391 + Ifges(4,6) * t390 + Ifges(4,3) * qJD(3);
t311 = t433 * t414 * mrSges(3,3) + m(3) * t388 + t421 * qJDD(1) + t415;
t304 = mrSges(5,2) * t337 + mrSges(6,2) * t324 - mrSges(5,3) * t328 - mrSges(6,3) * t326 - qJ(5) * t317 - t438 * t344 + t448 * t345 + t435 * t373 + t446 * t403 + t443 * t406;
t303 = -mrSges(5,1) * t337 - mrSges(6,1) * t326 + mrSges(6,2) * t323 + mrSges(5,3) * t329 - pkin(4) * t317 - t447 * t344 + t438 * t345 + t435 * t374 + t445 * t403 + t442 * t406;
t302 = mrSges(4,2) * t378 - mrSges(4,3) * t346 + Ifges(4,1) * t380 + Ifges(4,4) * t379 + Ifges(4,5) * qJDD(3) - pkin(7) * t310 - qJD(3) * t371 - t409 * t303 + t441 * t304 + t390 * t370;
t301 = -mrSges(4,1) * t378 + mrSges(4,3) * t347 + Ifges(4,4) * t380 + Ifges(4,2) * t379 + Ifges(4,6) * qJDD(3) - pkin(3) * t417 + pkin(7) * t425 + qJD(3) * t372 + t441 * t303 + t409 * t304 - t391 * t370;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t428 - mrSges(2,2) * t423 + t407 * (mrSges(3,2) * t388 - mrSges(3,3) * t381 + t412 * t302 - t410 * t301 - pkin(6) * t436 + (Ifges(3,1) * t407 + Ifges(3,4) * t408) * qJDD(1)) + t408 * (-mrSges(3,1) * t388 + mrSges(3,3) * t382 + t410 * t302 + t412 * t301 - pkin(2) * t415 + pkin(6) * t426 + (Ifges(3,4) * t407 + Ifges(3,2) * t408) * qJDD(1)) - pkin(1) * t311 + qJ(2) * ((m(3) * t382 + t418 * t408 + t426) * t408 + (-m(3) * t381 + t418 * t407 - t436) * t407); t311; mrSges(4,1) * t346 - mrSges(4,2) * t347 + Ifges(4,5) * t380 + Ifges(4,6) * t379 + Ifges(4,3) * qJDD(3) + pkin(3) * t310 + t391 * t371 - t390 * t372 + t416; t416; t320;];
tauJ = t1;

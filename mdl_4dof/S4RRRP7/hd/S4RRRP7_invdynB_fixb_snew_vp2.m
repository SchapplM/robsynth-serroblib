% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:16
% EndTime: 2019-12-31 17:20:19
% DurationCPUTime: 1.23s
% Computational Cost: add. (7983->214), mult. (15535->261), div. (0->0), fcn. (8838->6), ass. (0->84)
t431 = Ifges(4,1) + Ifges(5,1);
t425 = Ifges(4,4) - Ifges(5,5);
t430 = -Ifges(4,5) - Ifges(5,4);
t429 = Ifges(4,2) + Ifges(5,3);
t423 = Ifges(4,6) - Ifges(5,6);
t428 = -Ifges(4,3) - Ifges(5,2);
t427 = cos(qJ(3));
t426 = -mrSges(4,3) - mrSges(5,2);
t401 = sin(qJ(1));
t403 = cos(qJ(1));
t392 = -t403 * g(1) - t401 * g(2);
t405 = qJD(1) ^ 2;
t377 = -t405 * pkin(1) + qJDD(1) * pkin(5) + t392;
t400 = sin(qJ(2));
t402 = cos(qJ(2));
t369 = -t400 * g(3) + t402 * t377;
t385 = (-mrSges(3,1) * t402 + mrSges(3,2) * t400) * qJD(1);
t415 = qJD(1) * qJD(2);
t413 = t400 * t415;
t388 = t402 * qJDD(1) - t413;
t417 = qJD(1) * t400;
t389 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t417;
t391 = t401 * g(1) - t403 * g(2);
t376 = -qJDD(1) * pkin(1) - t405 * pkin(5) - t391;
t412 = t402 * t415;
t387 = t400 * qJDD(1) + t412;
t343 = (-t387 - t412) * pkin(6) + (-t388 + t413) * pkin(2) + t376;
t386 = (-pkin(2) * t402 - pkin(6) * t400) * qJD(1);
t404 = qJD(2) ^ 2;
t416 = t402 * qJD(1);
t346 = -t404 * pkin(2) + qJDD(2) * pkin(6) + t386 * t416 + t369;
t399 = sin(qJ(3));
t341 = t399 * t343 + t427 * t346;
t384 = t399 * qJD(2) + t427 * t417;
t357 = t384 * qJD(3) - t427 * qJDD(2) + t399 * t387;
t394 = qJD(3) - t416;
t365 = t394 * mrSges(4,1) - t384 * mrSges(4,3);
t382 = qJDD(3) - t388;
t383 = -t427 * qJD(2) + t399 * t417;
t361 = t383 * pkin(3) - t384 * qJ(4);
t393 = t394 ^ 2;
t337 = -t393 * pkin(3) + t382 * qJ(4) + 0.2e1 * qJD(4) * t394 - t383 * t361 + t341;
t366 = -t394 * mrSges(5,1) + t384 * mrSges(5,2);
t414 = m(5) * t337 + t382 * mrSges(5,3) + t394 * t366;
t362 = t383 * mrSges(5,1) - t384 * mrSges(5,3);
t418 = -t383 * mrSges(4,1) - t384 * mrSges(4,2) - t362;
t333 = m(4) * t341 - t382 * mrSges(4,2) + t426 * t357 - t394 * t365 + t418 * t383 + t414;
t340 = t427 * t343 - t399 * t346;
t358 = -t383 * qJD(3) + t399 * qJDD(2) + t427 * t387;
t364 = -t394 * mrSges(4,2) - t383 * mrSges(4,3);
t339 = -t382 * pkin(3) - t393 * qJ(4) + t384 * t361 + qJDD(4) - t340;
t367 = -t383 * mrSges(5,2) + t394 * mrSges(5,3);
t408 = -m(5) * t339 + t382 * mrSges(5,1) + t394 * t367;
t334 = m(4) * t340 + t382 * mrSges(4,1) + t426 * t358 + t394 * t364 + t418 * t384 + t408;
t409 = t427 * t333 - t399 * t334;
t328 = m(3) * t369 - qJDD(2) * mrSges(3,2) + t388 * mrSges(3,3) - qJD(2) * t389 + t385 * t416 + t409;
t368 = -t402 * g(3) - t400 * t377;
t390 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t416;
t345 = -qJDD(2) * pkin(2) - t404 * pkin(6) + t386 * t417 - t368;
t338 = -0.2e1 * qJD(4) * t384 + (t383 * t394 - t358) * qJ(4) + (t384 * t394 + t357) * pkin(3) + t345;
t335 = m(5) * t338 + t357 * mrSges(5,1) - t358 * mrSges(5,3) - t384 * t366 + t383 * t367;
t406 = -m(4) * t345 - t357 * mrSges(4,1) - t358 * mrSges(4,2) - t383 * t364 - t384 * t365 - t335;
t331 = m(3) * t368 + qJDD(2) * mrSges(3,1) - t387 * mrSges(3,3) + qJD(2) * t390 - t385 * t417 + t406;
t410 = t402 * t328 - t400 * t331;
t320 = m(2) * t392 - t405 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t410;
t329 = t399 * t333 + t427 * t334;
t407 = -m(3) * t376 + t388 * mrSges(3,1) - t387 * mrSges(3,2) - t389 * t417 + t390 * t416 - t329;
t325 = m(2) * t391 + qJDD(1) * mrSges(2,1) - t405 * mrSges(2,2) + t407;
t422 = t401 * t320 + t403 * t325;
t321 = t400 * t328 + t402 * t331;
t421 = t429 * t383 - t425 * t384 - t423 * t394;
t420 = t423 * t383 + t430 * t384 + t428 * t394;
t419 = -t425 * t383 + t431 * t384 - t430 * t394;
t411 = t403 * t320 - t401 * t325;
t375 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t400 + Ifges(3,4) * t402) * qJD(1);
t374 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t400 + Ifges(3,2) * t402) * qJD(1);
t373 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t400 + Ifges(3,6) * t402) * qJD(1);
t323 = mrSges(4,2) * t345 + mrSges(5,2) * t339 - mrSges(4,3) * t340 - mrSges(5,3) * t338 - qJ(4) * t335 - t425 * t357 + t431 * t358 - t382 * t430 + t420 * t383 + t421 * t394;
t322 = -mrSges(4,1) * t345 - mrSges(5,1) * t338 + mrSges(5,2) * t337 + mrSges(4,3) * t341 - pkin(3) * t335 - t429 * t357 + t425 * t358 + t423 * t382 + t420 * t384 + t419 * t394;
t317 = Ifges(3,4) * t387 + Ifges(3,2) * t388 + Ifges(3,6) * qJDD(2) - t373 * t417 + qJD(2) * t375 - mrSges(3,1) * t376 + mrSges(3,3) * t369 - mrSges(4,1) * t340 + mrSges(4,2) * t341 + mrSges(5,1) * t339 - mrSges(5,3) * t337 - pkin(3) * t408 - qJ(4) * t414 - pkin(2) * t329 + (pkin(3) * t362 + t421) * t384 + (qJ(4) * t362 - t419) * t383 + t428 * t382 + (pkin(3) * mrSges(5,2) + t430) * t358 + (qJ(4) * mrSges(5,2) + t423) * t357;
t316 = mrSges(3,2) * t376 - mrSges(3,3) * t368 + Ifges(3,1) * t387 + Ifges(3,4) * t388 + Ifges(3,5) * qJDD(2) - pkin(6) * t329 - qJD(2) * t374 - t399 * t322 + t427 * t323 + t373 * t416;
t315 = Ifges(2,6) * qJDD(1) + t405 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t392 - Ifges(3,5) * t387 - Ifges(3,6) * t388 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t368 + mrSges(3,2) * t369 - t399 * t323 - t427 * t322 - pkin(2) * t406 - pkin(6) * t409 - pkin(1) * t321 + (-t400 * t374 + t402 * t375) * qJD(1);
t314 = -mrSges(2,2) * g(3) - mrSges(2,3) * t391 + Ifges(2,5) * qJDD(1) - t405 * Ifges(2,6) - pkin(5) * t321 + t402 * t316 - t400 * t317;
t1 = [-m(1) * g(1) + t411; -m(1) * g(2) + t422; (-m(1) - m(2)) * g(3) + t321; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t422 + t403 * t314 - t401 * t315; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t411 + t401 * t314 + t403 * t315; -mrSges(1,1) * g(2) + mrSges(2,1) * t391 + mrSges(1,2) * g(1) - mrSges(2,2) * t392 + Ifges(2,3) * qJDD(1) + pkin(1) * t407 + pkin(5) * t410 + t400 * t316 + t402 * t317;];
tauB = t1;

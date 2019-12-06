% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:00
% EndTime: 2019-12-05 15:05:02
% DurationCPUTime: 1.35s
% Computational Cost: add. (14156->162), mult. (21681->210), div. (0->0), fcn. (14524->10), ass. (0->76)
t398 = sin(pkin(7));
t401 = cos(pkin(7));
t390 = -t401 * g(1) - t398 * g(2);
t395 = -g(3) + qJDD(1);
t397 = sin(pkin(8));
t400 = cos(pkin(8));
t377 = t400 * t390 + t397 * t395;
t389 = t398 * g(1) - t401 * g(2);
t388 = qJDD(2) - t389;
t403 = sin(qJ(3));
t405 = cos(qJ(3));
t372 = -t403 * t377 + t405 * t388;
t370 = qJDD(3) * pkin(3) + t372;
t373 = t405 * t377 + t403 * t388;
t406 = qJD(3) ^ 2;
t371 = -t406 * pkin(3) + t373;
t396 = sin(pkin(9));
t399 = cos(pkin(9));
t367 = t396 * t370 + t399 * t371;
t365 = -t406 * pkin(4) + qJDD(3) * pkin(6) + t367;
t376 = t397 * t390 - t400 * t395;
t375 = qJDD(4) + t376;
t402 = sin(qJ(5));
t404 = cos(qJ(5));
t362 = -t402 * t365 + t404 * t375;
t363 = t404 * t365 + t402 * t375;
t379 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t402 + Ifges(6,2) * t404) * qJD(3);
t380 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t402 + Ifges(6,4) * t404) * qJD(3);
t415 = qJD(3) * qJD(5);
t386 = t402 * qJDD(3) + t404 * t415;
t387 = t404 * qJDD(3) - t402 * t415;
t421 = mrSges(6,1) * t362 - mrSges(6,2) * t363 + Ifges(6,5) * t386 + Ifges(6,6) * t387 + Ifges(6,3) * qJDD(5) + (t379 * t402 - t380 * t404) * qJD(3);
t385 = (-mrSges(6,1) * t404 + mrSges(6,2) * t402) * qJD(3);
t416 = qJD(3) * t404;
t392 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t416;
t417 = qJD(3) * t402;
t359 = m(6) * t362 + qJDD(5) * mrSges(6,1) - t386 * mrSges(6,3) + qJD(5) * t392 - t385 * t417;
t391 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t417;
t360 = m(6) * t363 - qJDD(5) * mrSges(6,2) + t387 * mrSges(6,3) - qJD(5) * t391 + t385 * t416;
t351 = -t402 * t359 + t404 * t360;
t346 = m(5) * t367 - t406 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t351;
t366 = t399 * t370 - t396 * t371;
t364 = -qJDD(3) * pkin(4) - t406 * pkin(6) - t366;
t361 = -m(6) * t364 + t387 * mrSges(6,1) - t386 * mrSges(6,2) - t391 * t417 + t392 * t416;
t355 = m(5) * t366 + qJDD(3) * mrSges(5,1) - t406 * mrSges(5,2) + t361;
t343 = t396 * t346 + t399 * t355;
t378 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t402 + Ifges(6,6) * t404) * qJD(3);
t352 = -mrSges(6,1) * t364 + mrSges(6,3) * t363 + Ifges(6,4) * t386 + Ifges(6,2) * t387 + Ifges(6,6) * qJDD(5) + qJD(5) * t380 - t378 * t417;
t353 = mrSges(6,2) * t364 - mrSges(6,3) * t362 + Ifges(6,1) * t386 + Ifges(6,4) * t387 + Ifges(6,5) * qJDD(5) - qJD(5) * t379 + t378 * t416;
t420 = mrSges(4,1) * t372 + mrSges(5,1) * t366 - mrSges(4,2) * t373 - mrSges(5,2) * t367 + pkin(3) * t343 + pkin(4) * t361 + pkin(6) * t351 + t404 * t352 + t402 * t353 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t341 = m(4) * t372 + qJDD(3) * mrSges(4,1) - t406 * mrSges(4,2) + t343;
t410 = t399 * t346 - t396 * t355;
t342 = m(4) * t373 - t406 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t410;
t411 = -t403 * t341 + t405 * t342;
t335 = m(3) * t377 + t411;
t350 = t404 * t359 + t402 * t360;
t349 = m(5) * t375 + t350;
t348 = (-m(3) - m(4)) * t376 - t349;
t412 = t400 * t335 - t397 * t348;
t329 = m(2) * t390 + t412;
t337 = t405 * t341 + t403 * t342;
t336 = m(3) * t388 + t337;
t334 = m(2) * t389 - t336;
t418 = t398 * t329 + t401 * t334;
t330 = t397 * t335 + t400 * t348;
t414 = m(2) * t395 + t330;
t413 = t401 * t329 - t398 * t334;
t339 = -mrSges(5,1) * t375 + mrSges(5,3) * t367 + t406 * Ifges(5,5) + Ifges(5,6) * qJDD(3) - pkin(4) * t350 - t421;
t338 = mrSges(5,2) * t375 - mrSges(5,3) * t366 + Ifges(5,5) * qJDD(3) - t406 * Ifges(5,6) - pkin(6) * t350 - t402 * t352 + t404 * t353;
t326 = mrSges(4,2) * t376 - mrSges(4,3) * t372 + Ifges(4,5) * qJDD(3) - t406 * Ifges(4,6) - qJ(4) * t343 + t399 * t338 - t396 * t339;
t325 = -mrSges(4,1) * t376 + mrSges(4,3) * t373 + t406 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t349 + qJ(4) * t410 + t396 * t338 + t399 * t339;
t324 = -mrSges(3,1) * t388 + mrSges(3,3) * t377 - pkin(2) * t337 - t420;
t323 = mrSges(3,2) * t388 + mrSges(3,3) * t376 - pkin(5) * t337 - t403 * t325 + t405 * t326;
t322 = -mrSges(2,1) * t395 + mrSges(2,3) * t390 + mrSges(3,1) * t376 + mrSges(3,2) * t377 - t403 * t326 - t405 * t325 - pkin(2) * (-m(4) * t376 - t349) - pkin(5) * t411 - pkin(1) * t330;
t321 = mrSges(2,2) * t395 - mrSges(2,3) * t389 - qJ(2) * t330 + t400 * t323 - t397 * t324;
t1 = [-m(1) * g(1) + t413; -m(1) * g(2) + t418; -m(1) * g(3) + t414; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t418 + t401 * t321 - t398 * t322; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t413 + t398 * t321 + t401 * t322; -mrSges(1,1) * g(2) + mrSges(2,1) * t389 + mrSges(1,2) * g(1) - mrSges(2,2) * t390 - pkin(1) * t336 + qJ(2) * t412 + t397 * t323 + t400 * t324; t414; t336; t420; t349; t421;];
tauJB = t1;

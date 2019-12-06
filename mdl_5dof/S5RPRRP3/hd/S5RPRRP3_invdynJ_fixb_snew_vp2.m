% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP3
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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:23
% EndTime: 2019-12-05 18:03:26
% DurationCPUTime: 1.05s
% Computational Cost: add. (4142->197), mult. (8267->242), div. (0->0), fcn. (4901->8), ass. (0->82)
t399 = Ifges(5,4) + Ifges(6,4);
t405 = Ifges(5,2) + Ifges(6,2);
t402 = Ifges(5,6) + Ifges(6,6);
t404 = Ifges(5,1) + Ifges(6,1);
t403 = Ifges(5,5) + Ifges(6,5);
t401 = Ifges(5,3) + Ifges(6,3);
t376 = sin(qJ(4));
t377 = sin(qJ(3));
t379 = cos(qJ(4));
t380 = cos(qJ(3));
t349 = (-t376 * t377 + t379 * t380) * qJD(1);
t350 = (t376 * t380 + t377 * t379) * qJD(1);
t371 = qJD(3) + qJD(4);
t400 = t405 * t349 + t399 * t350 + t402 * t371;
t378 = sin(qJ(1));
t381 = cos(qJ(1));
t396 = t381 * g(2) + t378 * g(3);
t354 = qJDD(1) * pkin(1) + t396;
t382 = qJD(1) ^ 2;
t389 = t378 * g(2) - g(3) * t381;
t356 = -pkin(1) * t382 + t389;
t374 = sin(pkin(8));
t375 = cos(pkin(8));
t335 = t374 * t354 + t375 * t356;
t331 = -pkin(2) * t382 + qJDD(1) * pkin(6) + t335;
t373 = -g(1) + qJDD(2);
t316 = -t331 * t377 + t380 * t373;
t393 = qJD(1) * qJD(3);
t390 = t380 * t393;
t357 = qJDD(1) * t377 + t390;
t302 = (-t357 + t390) * pkin(7) + (t377 * t380 * t382 + qJDD(3)) * pkin(3) + t316;
t317 = t380 * t331 + t377 * t373;
t358 = qJDD(1) * t380 - t377 * t393;
t395 = qJD(1) * t377;
t361 = qJD(3) * pkin(3) - pkin(7) * t395;
t372 = t380 ^ 2;
t303 = -pkin(3) * t372 * t382 + pkin(7) * t358 - qJD(3) * t361 + t317;
t297 = t379 * t302 - t303 * t376;
t319 = qJD(4) * t349 + t357 * t379 + t358 * t376;
t332 = -mrSges(6,1) * t349 + mrSges(6,2) * t350;
t333 = -mrSges(5,1) * t349 + mrSges(5,2) * t350;
t338 = -mrSges(5,2) * t371 + mrSges(5,3) * t349;
t370 = qJDD(3) + qJDD(4);
t291 = -0.2e1 * qJD(5) * t350 + (t349 * t371 - t319) * qJ(5) + (t349 * t350 + t370) * pkin(4) + t297;
t337 = -mrSges(6,2) * t371 + mrSges(6,3) * t349;
t392 = m(6) * t291 + t370 * mrSges(6,1) + t371 * t337;
t282 = m(5) * t297 + mrSges(5,1) * t370 + t338 * t371 + (-t332 - t333) * t350 + (-mrSges(5,3) - mrSges(6,3)) * t319 + t392;
t298 = t376 * t302 + t379 * t303;
t318 = -qJD(4) * t350 - t357 * t376 + t358 * t379;
t340 = mrSges(6,1) * t371 - mrSges(6,3) * t350;
t341 = mrSges(5,1) * t371 - mrSges(5,3) * t350;
t339 = pkin(4) * t371 - qJ(5) * t350;
t342 = t349 ^ 2;
t293 = -pkin(4) * t342 + qJ(5) * t318 + 0.2e1 * qJD(5) * t349 - t339 * t371 + t298;
t391 = m(6) * t293 + t318 * mrSges(6,3) + t349 * t332;
t285 = m(5) * t298 + mrSges(5,3) * t318 + t333 * t349 + (-t340 - t341) * t371 + (-mrSges(5,2) - mrSges(6,2)) * t370 + t391;
t280 = t379 * t282 + t376 * t285;
t398 = -t402 * t349 - t403 * t350 - t401 * t371;
t397 = -t399 * t349 - t404 * t350 - t403 * t371;
t394 = qJD(1) * t380;
t355 = (-mrSges(4,1) * t380 + mrSges(4,2) * t377) * qJD(1);
t360 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t394;
t278 = m(4) * t316 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t357 + qJD(3) * t360 - t355 * t395 + t280;
t359 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t395;
t387 = -t282 * t376 + t379 * t285;
t279 = m(4) * t317 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t358 - qJD(3) * t359 + t355 * t394 + t387;
t388 = -t278 * t377 + t380 * t279;
t334 = t354 * t375 - t374 * t356;
t386 = -qJDD(1) * pkin(2) - t334;
t304 = -pkin(3) * t358 + t361 * t395 + (-pkin(7) * t372 - pkin(6)) * t382 + t386;
t295 = -pkin(4) * t318 - qJ(5) * t342 + t339 * t350 + qJDD(5) + t304;
t288 = m(6) * t295 - t318 * mrSges(6,1) + t319 * mrSges(6,2) - t349 * t337 + t350 * t340;
t385 = m(5) * t304 - t318 * mrSges(5,1) + mrSges(5,2) * t319 - t349 * t338 + t341 * t350 + t288;
t287 = -mrSges(6,3) * t319 - t332 * t350 + t392;
t384 = mrSges(5,1) * t297 + mrSges(6,1) * t291 - mrSges(5,2) * t298 - mrSges(6,2) * t293 + pkin(4) * t287 + t402 * t318 + t403 * t319 + t397 * t349 + t400 * t350 + t401 * t370;
t330 = -pkin(6) * t382 + t386;
t383 = -m(4) * t330 + t358 * mrSges(4,1) - mrSges(4,2) * t357 - t359 * t395 + t360 * t394 - t385;
t348 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t377 + Ifges(4,4) * t380) * qJD(1);
t347 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t377 + Ifges(4,2) * t380) * qJD(1);
t276 = mrSges(5,2) * t304 + mrSges(6,2) * t295 - mrSges(5,3) * t297 - mrSges(6,3) * t291 - qJ(5) * t287 + t399 * t318 + t404 * t319 - t398 * t349 + t403 * t370 - t400 * t371;
t275 = -mrSges(5,1) * t304 + mrSges(5,3) * t298 - mrSges(6,1) * t295 + mrSges(6,3) * t293 - pkin(4) * t288 + qJ(5) * t391 + (-qJ(5) * t340 - t397) * t371 + (-mrSges(6,2) * qJ(5) + t402) * t370 + t398 * t350 + t399 * t319 + t405 * t318;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t396 - mrSges(2,2) * t389 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t334 - mrSges(3,2) * t335 + t377 * (mrSges(4,2) * t330 - mrSges(4,3) * t316 + Ifges(4,1) * t357 + Ifges(4,4) * t358 + Ifges(4,5) * qJDD(3) - pkin(7) * t280 - qJD(3) * t347 - t376 * t275 + t379 * t276) + t380 * (-mrSges(4,1) * t330 + mrSges(4,3) * t317 + Ifges(4,4) * t357 + Ifges(4,2) * t358 + Ifges(4,6) * qJDD(3) - pkin(3) * t385 + pkin(7) * t387 + qJD(3) * t348 + t379 * t275 + t376 * t276) + pkin(2) * t383 + pkin(6) * t388 + pkin(1) * (t374 * (m(3) * t335 - mrSges(3,1) * t382 - qJDD(1) * mrSges(3,2) + t388) + t375 * (m(3) * t334 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t382 + t383)); m(3) * t373 + t278 * t380 + t279 * t377; pkin(3) * t280 + t384 + mrSges(4,1) * t316 - mrSges(4,2) * t317 + Ifges(4,5) * t357 + Ifges(4,6) * t358 + Ifges(4,3) * qJDD(3) + (t347 * t377 - t348 * t380) * qJD(1); t384; t288;];
tauJ = t1;

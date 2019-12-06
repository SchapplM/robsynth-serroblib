% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:26
% EndTime: 2019-12-05 17:59:28
% DurationCPUTime: 1.06s
% Computational Cost: add. (3629->199), mult. (7237->237), div. (0->0), fcn. (4179->6), ass. (0->78)
t395 = Ifges(5,4) + Ifges(6,4);
t403 = Ifges(5,2) + Ifges(6,2);
t400 = Ifges(5,6) + Ifges(6,6);
t402 = Ifges(5,1) + Ifges(6,1);
t401 = Ifges(5,5) + Ifges(6,5);
t399 = Ifges(5,3) + Ifges(6,3);
t370 = sin(qJ(4));
t371 = sin(qJ(3));
t373 = cos(qJ(4));
t374 = cos(qJ(3));
t351 = (-t374 * t370 - t371 * t373) * qJD(1);
t352 = (-t371 * t370 + t374 * t373) * qJD(1);
t367 = qJD(3) + qJD(4);
t398 = t403 * t351 + t395 * t352 + t400 * t367;
t372 = sin(qJ(1));
t375 = cos(qJ(1));
t382 = -g(1) * t375 - g(2) * t372;
t379 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t382;
t390 = qJD(1) * qJD(3);
t355 = -qJDD(1) * t371 - t374 * t390;
t391 = t374 * qJD(1);
t359 = qJD(3) * pkin(3) - pkin(7) * t391;
t369 = t371 ^ 2;
t376 = qJD(1) ^ 2;
t396 = -pkin(1) - pkin(6);
t312 = -pkin(3) * t355 + t359 * t391 + (-pkin(7) * t369 + t396) * t376 + t379;
t385 = t371 * t390;
t356 = qJDD(1) * t374 - t385;
t322 = qJD(4) * t351 + t355 * t370 + t356 * t373;
t342 = mrSges(5,1) * t367 - mrSges(5,3) * t352;
t397 = m(5) * t312 + t322 * mrSges(5,2) + t352 * t342;
t384 = g(1) * t372 - t375 * g(2);
t378 = -qJ(2) * t376 + qJDD(2) - t384;
t344 = t396 * qJDD(1) + t378;
t333 = t371 * g(3) + t374 * t344;
t307 = (-t356 - t385) * pkin(7) + (-t371 * t374 * t376 + qJDD(3)) * pkin(3) + t333;
t334 = -g(3) * t374 + t371 * t344;
t308 = -pkin(3) * t369 * t376 + pkin(7) * t355 - qJD(3) * t359 + t334;
t302 = t373 * t307 - t308 * t370;
t331 = -mrSges(6,1) * t351 + mrSges(6,2) * t352;
t332 = -mrSges(5,1) * t351 + mrSges(5,2) * t352;
t339 = -mrSges(5,2) * t367 + mrSges(5,3) * t351;
t366 = qJDD(3) + qJDD(4);
t296 = -0.2e1 * qJD(5) * t352 + (t351 * t367 - t322) * qJ(5) + (t351 * t352 + t366) * pkin(4) + t302;
t338 = -mrSges(6,2) * t367 + mrSges(6,3) * t351;
t388 = m(6) * t296 + t366 * mrSges(6,1) + t367 * t338;
t287 = m(5) * t302 + mrSges(5,1) * t366 + t339 * t367 + (-t331 - t332) * t352 + (-mrSges(5,3) - mrSges(6,3)) * t322 + t388;
t303 = t370 * t307 + t373 * t308;
t321 = -qJD(4) * t352 + t355 * t373 - t356 * t370;
t341 = mrSges(6,1) * t367 - mrSges(6,3) * t352;
t340 = pkin(4) * t367 - qJ(5) * t352;
t347 = t351 ^ 2;
t298 = -pkin(4) * t347 + qJ(5) * t321 + 0.2e1 * qJD(5) * t351 - t340 * t367 + t303;
t387 = m(6) * t298 + t321 * mrSges(6,3) + t351 * t331;
t290 = m(5) * t303 + mrSges(5,3) * t321 + t332 * t351 + (-t341 - t342) * t367 + (-mrSges(5,2) - mrSges(6,2)) * t366 + t387;
t285 = t373 * t287 + t370 * t290;
t394 = -t400 * t351 - t401 * t352 - t399 * t367;
t393 = -t395 * t351 - t402 * t352 - t401 * t367;
t392 = qJD(1) * t371;
t300 = -pkin(4) * t321 - qJ(5) * t347 + t340 * t352 + qJDD(5) + t312;
t386 = m(6) * t300 + t322 * mrSges(6,2) + t352 * t341;
t383 = -t287 * t370 + t373 * t290;
t354 = (t371 * mrSges(4,1) + t374 * mrSges(4,2)) * qJD(1);
t357 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t392;
t358 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t391;
t381 = (m(4) * t333 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t356 + qJD(3) * t357 - t354 * t391 + t285) * t374 + t371 * (m(4) * t334 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t355 - qJD(3) * t358 - t354 * t392 + t383);
t293 = -mrSges(6,1) * t321 - t338 * t351 + t386;
t292 = -mrSges(6,3) * t322 - t331 * t352 + t388;
t377 = mrSges(5,1) * t302 + mrSges(6,1) * t296 - mrSges(5,2) * t303 - mrSges(6,2) * t298 + pkin(4) * t292 + t400 * t321 + t401 * t322 + t393 * t351 + t398 * t352 + t399 * t366;
t350 = Ifges(4,5) * qJD(3) + (t374 * Ifges(4,1) - t371 * Ifges(4,4)) * qJD(1);
t349 = Ifges(4,6) * qJD(3) + (t374 * Ifges(4,4) - t371 * Ifges(4,2)) * qJD(1);
t346 = -qJDD(1) * pkin(1) + t378;
t345 = pkin(1) * t376 - t379;
t343 = t396 * t376 + t379;
t282 = mrSges(5,2) * t312 + mrSges(6,2) * t300 - mrSges(5,3) * t302 - mrSges(6,3) * t296 - qJ(5) * t292 + t395 * t321 + t402 * t322 - t394 * t351 + t401 * t366 - t398 * t367;
t281 = -mrSges(5,1) * t312 + mrSges(5,3) * t303 - mrSges(6,1) * t300 + mrSges(6,3) * t298 - pkin(4) * t293 + qJ(5) * t387 + (-qJ(5) * t341 - t393) * t367 + (-qJ(5) * mrSges(6,2) + t400) * t366 + t394 * t352 + t395 * t322 + t403 * t321;
t280 = m(3) * t346 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t376 + t381;
t1 = [mrSges(2,1) * t384 - mrSges(2,2) * t382 + mrSges(3,2) * t346 - mrSges(3,3) * t345 + t374 * (mrSges(4,2) * t343 - mrSges(4,3) * t333 + Ifges(4,1) * t356 + Ifges(4,4) * t355 + Ifges(4,5) * qJDD(3) - pkin(7) * t285 - qJD(3) * t349 - t281 * t370 + t282 * t373) - t371 * (Ifges(4,4) * t356 + Ifges(4,2) * t355 + Ifges(4,6) * qJDD(3) + qJD(3) * t350 - mrSges(4,1) * t343 + mrSges(4,3) * t334 + t370 * t282 + t373 * t281 - pkin(3) * ((-t338 - t339) * t351 + (-mrSges(5,1) - mrSges(6,1)) * t321 + t386 + t397) + pkin(7) * t383) - pkin(6) * t381 - pkin(1) * t280 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t345 + m(4) * t343 - mrSges(4,1) * t355 - mrSges(5,1) * t321 + mrSges(3,2) * t376 + mrSges(4,2) * t356 - t339 * t351 + t293 + qJDD(1) * mrSges(3,3) + (t357 * t371 + t358 * t374) * qJD(1) + t397) * qJ(2); t280; t377 + (t374 * t349 + t371 * t350) * qJD(1) + Ifges(4,3) * qJDD(3) + pkin(3) * t285 + mrSges(4,1) * t333 - mrSges(4,2) * t334 + Ifges(4,6) * t355 + Ifges(4,5) * t356; t377; t293;];
tauJ = t1;

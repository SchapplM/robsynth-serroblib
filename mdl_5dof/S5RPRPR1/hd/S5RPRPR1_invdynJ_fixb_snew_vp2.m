% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:19
% EndTime: 2019-12-05 17:47:21
% DurationCPUTime: 1.16s
% Computational Cost: add. (7884->215), mult. (17161->278), div. (0->0), fcn. (10810->8), ass. (0->85)
t363 = sin(qJ(1));
t366 = cos(qJ(1));
t374 = -g(1) * t366 - g(2) * t363;
t371 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t374;
t383 = -pkin(1) - pkin(6);
t367 = qJD(1) ^ 2;
t377 = t363 * g(1) - t366 * g(2);
t370 = -qJ(2) * t367 + qJDD(2) - t377;
t333 = t383 * qJDD(1) + t370;
t362 = sin(qJ(3));
t365 = cos(qJ(3));
t325 = t362 * g(3) + t365 * t333;
t380 = qJD(1) * qJD(3);
t378 = t362 * t380;
t347 = qJDD(1) * t365 - t378;
t309 = (-t347 - t378) * qJ(4) + (-t362 * t365 * t367 + qJDD(3)) * pkin(3) + t325;
t326 = -g(3) * t365 + t362 * t333;
t346 = -qJDD(1) * t362 - t365 * t380;
t381 = t365 * qJD(1);
t349 = qJD(3) * pkin(3) - qJ(4) * t381;
t358 = t362 ^ 2;
t310 = -pkin(3) * t358 * t367 + qJ(4) * t346 - qJD(3) * t349 + t326;
t359 = sin(pkin(8));
t360 = cos(pkin(8));
t340 = (-t362 * t359 + t365 * t360) * qJD(1);
t292 = -0.2e1 * qJD(4) * t340 + t360 * t309 - t310 * t359;
t324 = t346 * t359 + t347 * t360;
t339 = (-t365 * t359 - t362 * t360) * qJD(1);
t289 = (qJD(3) * t339 - t324) * pkin(7) + (t339 * t340 + qJDD(3)) * pkin(4) + t292;
t293 = 0.2e1 * qJD(4) * t339 + t359 * t309 + t360 * t310;
t323 = t346 * t360 - t347 * t359;
t332 = qJD(3) * pkin(4) - pkin(7) * t340;
t338 = t339 ^ 2;
t290 = -pkin(4) * t338 + pkin(7) * t323 - qJD(3) * t332 + t293;
t361 = sin(qJ(5));
t364 = cos(qJ(5));
t287 = t289 * t364 - t290 * t361;
t318 = t339 * t364 - t340 * t361;
t300 = qJD(5) * t318 + t323 * t361 + t324 * t364;
t319 = t339 * t361 + t340 * t364;
t305 = -mrSges(6,1) * t318 + mrSges(6,2) * t319;
t356 = qJD(3) + qJD(5);
t313 = -mrSges(6,2) * t356 + mrSges(6,3) * t318;
t355 = qJDD(3) + qJDD(5);
t283 = m(6) * t287 + mrSges(6,1) * t355 - mrSges(6,3) * t300 - t305 * t319 + t313 * t356;
t288 = t289 * t361 + t290 * t364;
t299 = -qJD(5) * t319 + t323 * t364 - t324 * t361;
t314 = mrSges(6,1) * t356 - mrSges(6,3) * t319;
t284 = m(6) * t288 - mrSges(6,2) * t355 + mrSges(6,3) * t299 + t305 * t318 - t314 * t356;
t277 = t364 * t283 + t361 * t284;
t321 = -mrSges(5,1) * t339 + mrSges(5,2) * t340;
t330 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t339;
t275 = m(5) * t292 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t324 + qJD(3) * t330 - t321 * t340 + t277;
t331 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t340;
t375 = -t361 * t283 + t364 * t284;
t276 = m(5) * t293 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t323 - qJD(3) * t331 + t321 * t339 + t375;
t271 = t360 * t275 + t359 * t276;
t382 = qJD(1) * t362;
t376 = -t275 * t359 + t360 * t276;
t345 = (t362 * mrSges(4,1) + t365 * mrSges(4,2)) * qJD(1);
t348 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t382;
t350 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t381;
t373 = (m(4) * t325 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t347 + qJD(3) * t348 - t345 * t381 + t271) * t365 + (m(4) * t326 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t346 - qJD(3) * t350 - t345 * t382 + t376) * t362;
t312 = -pkin(3) * t346 + qJDD(4) + t349 * t381 + (-qJ(4) * t358 + t383) * t367 + t371;
t294 = -pkin(4) * t323 - pkin(7) * t338 + t332 * t340 + t312;
t369 = m(6) * t294 - mrSges(6,1) * t299 + t300 * mrSges(6,2) - t313 * t318 + t319 * t314;
t302 = Ifges(6,4) * t319 + Ifges(6,2) * t318 + Ifges(6,6) * t356;
t303 = Ifges(6,1) * t319 + Ifges(6,4) * t318 + Ifges(6,5) * t356;
t368 = mrSges(6,1) * t287 - mrSges(6,2) * t288 + Ifges(6,5) * t300 + Ifges(6,6) * t299 + Ifges(6,3) * t355 + t319 * t302 - t318 * t303;
t285 = m(5) * t312 - mrSges(5,1) * t323 + t324 * mrSges(5,2) - t330 * t339 + t340 * t331 + t369;
t343 = (Ifges(4,5) * qJD(3)) + (t365 * Ifges(4,1) - t362 * Ifges(4,4)) * qJD(1);
t342 = (Ifges(4,6) * qJD(3)) + (t365 * Ifges(4,4) - t362 * Ifges(4,2)) * qJD(1);
t337 = -qJDD(1) * pkin(1) + t370;
t334 = pkin(1) * t367 - t371;
t329 = t383 * t367 + t371;
t317 = Ifges(5,1) * t340 + Ifges(5,4) * t339 + (Ifges(5,5) * qJD(3));
t316 = Ifges(5,4) * t340 + Ifges(5,2) * t339 + (Ifges(5,6) * qJD(3));
t315 = Ifges(5,5) * t340 + Ifges(5,6) * t339 + (Ifges(5,3) * qJD(3));
t301 = Ifges(6,5) * t319 + Ifges(6,6) * t318 + Ifges(6,3) * t356;
t279 = mrSges(6,2) * t294 - mrSges(6,3) * t287 + Ifges(6,1) * t300 + Ifges(6,4) * t299 + Ifges(6,5) * t355 + t301 * t318 - t302 * t356;
t278 = -mrSges(6,1) * t294 + mrSges(6,3) * t288 + Ifges(6,4) * t300 + Ifges(6,2) * t299 + Ifges(6,6) * t355 - t301 * t319 + t303 * t356;
t268 = mrSges(5,2) * t312 - mrSges(5,3) * t292 + Ifges(5,1) * t324 + Ifges(5,4) * t323 + Ifges(5,5) * qJDD(3) - pkin(7) * t277 - qJD(3) * t316 - t278 * t361 + t279 * t364 + t315 * t339;
t267 = -mrSges(5,1) * t312 + mrSges(5,3) * t293 + Ifges(5,4) * t324 + Ifges(5,2) * t323 + Ifges(5,6) * qJDD(3) - pkin(4) * t369 + pkin(7) * t375 + qJD(3) * t317 + t364 * t278 + t361 * t279 - t340 * t315;
t266 = m(3) * t337 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t367) + t373;
t1 = [mrSges(2,1) * t377 - mrSges(2,2) * t374 + mrSges(3,2) * t337 - mrSges(3,3) * t334 + t365 * (mrSges(4,2) * t329 - mrSges(4,3) * t325 + Ifges(4,1) * t347 + Ifges(4,4) * t346 + Ifges(4,5) * qJDD(3) - qJ(4) * t271 - qJD(3) * t342 - t267 * t359 + t268 * t360) - t362 * (-mrSges(4,1) * t329 + mrSges(4,3) * t326 + Ifges(4,4) * t347 + Ifges(4,2) * t346 + Ifges(4,6) * qJDD(3) - pkin(3) * t285 + qJ(4) * t376 + qJD(3) * t343 + t360 * t267 + t359 * t268) - pkin(6) * t373 - pkin(1) * t266 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t334 + m(4) * t329 - mrSges(4,1) * t346 + mrSges(3,2) * t367 + mrSges(4,2) * t347 + t285 + qJDD(1) * mrSges(3,3) + (t348 * t362 + t350 * t365) * qJD(1)) * qJ(2); t266; (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t365 * t342 + t362 * t343) * qJD(1) + t368 + Ifges(4,6) * t346 + Ifges(4,5) * t347 - t339 * t317 + t340 * t316 + Ifges(5,6) * t323 + Ifges(5,5) * t324 + mrSges(4,1) * t325 - mrSges(4,2) * t326 + mrSges(5,1) * t292 - mrSges(5,2) * t293 + pkin(4) * t277 + pkin(3) * t271; t285; t368;];
tauJ = t1;

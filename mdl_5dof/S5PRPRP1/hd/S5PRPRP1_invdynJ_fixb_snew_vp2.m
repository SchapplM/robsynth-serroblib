% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:28
% EndTime: 2019-12-05 15:28:29
% DurationCPUTime: 0.60s
% Computational Cost: add. (2191->154), mult. (4918->188), div. (0->0), fcn. (3179->8), ass. (0->78)
t379 = Ifges(5,1) + Ifges(6,1);
t371 = Ifges(5,4) - Ifges(6,5);
t370 = Ifges(5,5) + Ifges(6,4);
t378 = -Ifges(5,2) - Ifges(6,3);
t369 = Ifges(5,6) - Ifges(6,6);
t377 = Ifges(5,3) + Ifges(6,2);
t344 = qJD(2) ^ 2;
t336 = sin(pkin(8));
t338 = cos(pkin(8));
t340 = sin(qJ(4));
t374 = cos(qJ(4));
t376 = t336 * t340 - t338 * t374;
t332 = t338 ^ 2;
t375 = 0.2e1 * t338;
t373 = pkin(3) * t344;
t372 = -mrSges(5,3) - mrSges(6,2);
t368 = pkin(6) * qJDD(2);
t337 = sin(pkin(7));
t339 = cos(pkin(7));
t323 = g(1) * t337 - g(2) * t339;
t324 = -g(1) * t339 - g(2) * t337;
t341 = sin(qJ(2));
t342 = cos(qJ(2));
t361 = t341 * t323 + t342 * t324;
t309 = -pkin(2) * t344 + qJDD(2) * qJ(3) + t361;
t335 = -g(3) + qJDD(1);
t356 = qJD(2) * qJD(3);
t360 = t338 * t335 - 0.2e1 * t336 * t356;
t289 = (t338 * t373 - t309 - t368) * t336 + t360;
t293 = t338 * t309 + t336 * t335 + t356 * t375;
t290 = -t332 * t373 + t338 * t368 + t293;
t286 = t340 * t289 + t374 * t290;
t346 = t336 * t374 + t338 * t340;
t317 = t346 * qJD(2);
t358 = qJD(4) * t317;
t307 = qJDD(2) * t376 + t358;
t313 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t317;
t316 = t376 * qJD(2);
t302 = pkin(4) * t316 - qJ(5) * t317;
t343 = qJD(4) ^ 2;
t281 = -pkin(4) * t343 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t302 * t316 + t286;
t314 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t317;
t355 = m(6) * t281 + qJDD(4) * mrSges(6,3) + qJD(4) * t314;
t303 = mrSges(6,1) * t316 - mrSges(6,3) * t317;
t362 = -mrSges(5,1) * t316 - mrSges(5,2) * t317 - t303;
t276 = m(5) * t286 - qJDD(4) * mrSges(5,2) - qJD(4) * t313 + t307 * t372 + t316 * t362 + t355;
t285 = t289 * t374 - t340 * t290;
t357 = t316 * qJD(4);
t308 = qJDD(2) * t346 - t357;
t312 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t316;
t282 = -qJDD(4) * pkin(4) - t343 * qJ(5) + t317 * t302 + qJDD(5) - t285;
t315 = -mrSges(6,2) * t316 + qJD(4) * mrSges(6,3);
t350 = -m(6) * t282 + qJDD(4) * mrSges(6,1) + qJD(4) * t315;
t277 = m(5) * t285 + qJDD(4) * mrSges(5,1) + qJD(4) * t312 + t308 * t372 + t317 * t362 + t350;
t366 = t340 * t276 + t374 * t277;
t365 = -qJD(4) * t377 + t316 * t369 - t317 * t370;
t364 = qJD(4) * t369 + t316 * t378 + t317 * t371;
t363 = qJD(4) * t370 - t316 * t371 + t317 * t379;
t359 = -t336 ^ 2 - t332;
t351 = t342 * t323 - t341 * t324;
t348 = qJDD(3) - t351;
t291 = (-pkin(3) * t338 - pkin(2)) * qJDD(2) + (pkin(6) * t359 - qJ(3)) * t344 + t348;
t284 = -0.2e1 * qJD(5) * t317 + (-t308 + t357) * qJ(5) + (t307 + t358) * pkin(4) + t291;
t354 = m(6) * t284 + t307 * mrSges(6,1) + t316 * t315;
t352 = t374 * t276 - t340 * t277;
t349 = -mrSges(4,1) * t338 + mrSges(4,2) * t336;
t347 = mrSges(4,3) * qJDD(2) + t344 * t349;
t345 = m(5) * t291 + t307 * mrSges(5,1) + t316 * t312 + (t313 - t314) * t317 + (mrSges(5,2) - mrSges(6,3)) * t308 + t354;
t306 = -qJDD(2) * pkin(2) - t344 * qJ(3) + t348;
t292 = -t309 * t336 + t360;
t279 = t308 * mrSges(6,2) + t317 * t303 - t350;
t278 = -t308 * mrSges(6,3) - t317 * t314 + t354;
t272 = mrSges(4,3) * t344 * t359 + m(4) * t306 + qJDD(2) * t349 + t345;
t271 = m(4) * t293 + t338 * t347 + t352;
t270 = m(4) * t292 - t336 * t347 + t366;
t269 = mrSges(5,2) * t291 + mrSges(6,2) * t282 - mrSges(5,3) * t285 - mrSges(6,3) * t284 - qJ(5) * t278 - t364 * qJD(4) + t370 * qJDD(4) - t371 * t307 + t308 * t379 + t365 * t316;
t268 = -mrSges(5,1) * t291 - mrSges(6,1) * t284 + mrSges(6,2) * t281 + mrSges(5,3) * t286 - pkin(4) * t278 + t363 * qJD(4) + t369 * qJDD(4) + t307 * t378 + t371 * t308 + t365 * t317;
t1 = [t338 * t270 + t336 * t271 + (m(2) + m(3)) * t335; mrSges(3,1) * t351 - mrSges(3,2) * t361 + t336 * (mrSges(4,2) * t306 - mrSges(4,3) * t292 - pkin(6) * t366 - t340 * t268 + t269 * t374) + t338 * (-mrSges(4,1) * t306 + mrSges(4,3) * t293 - pkin(3) * t345 + pkin(6) * t352 + t268 * t374 + t340 * t269) - pkin(2) * t272 + qJ(3) * (-t270 * t336 + t271 * t338) + (Ifges(4,2) * t332 + Ifges(3,3) + (Ifges(4,1) * t336 + Ifges(4,4) * t375) * t336) * qJDD(2); t272; mrSges(5,1) * t285 - mrSges(5,2) * t286 - mrSges(6,1) * t282 + mrSges(6,3) * t281 - pkin(4) * t279 + qJ(5) * t355 + t364 * t317 + (-qJ(5) * t303 + t363) * t316 + t370 * t308 + (-mrSges(6,2) * qJ(5) - t369) * t307 + t377 * qJDD(4); t279;];
tauJ = t1;

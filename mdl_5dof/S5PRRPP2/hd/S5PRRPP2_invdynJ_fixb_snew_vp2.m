% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:56
% EndTime: 2019-12-05 16:08:57
% DurationCPUTime: 0.77s
% Computational Cost: add. (2776->188), mult. (5960->234), div. (0->0), fcn. (3613->8), ass. (0->78)
t367 = Ifges(5,1) + Ifges(6,1);
t362 = Ifges(5,4) - Ifges(6,5);
t361 = Ifges(5,5) + Ifges(6,4);
t366 = -Ifges(5,2) - Ifges(6,3);
t365 = -Ifges(6,2) - Ifges(5,3);
t360 = Ifges(5,6) - Ifges(6,6);
t364 = -2 * qJD(4);
t363 = -mrSges(5,3) - mrSges(6,2);
t359 = cos(pkin(8));
t336 = sin(pkin(7));
t337 = cos(pkin(7));
t325 = -t337 * g(1) - t336 * g(2);
t334 = -g(3) + qJDD(1);
t339 = sin(qJ(2));
t341 = cos(qJ(2));
t306 = t341 * t325 + t339 * t334;
t343 = qJD(2) ^ 2;
t298 = -t343 * pkin(2) + qJDD(2) * pkin(6) + t306;
t324 = -t336 * g(1) + t337 * g(2);
t338 = sin(qJ(3));
t340 = cos(qJ(3));
t277 = -t338 * t298 + t340 * t324;
t352 = qJD(2) * qJD(3);
t350 = t340 * t352;
t322 = t338 * qJDD(2) + t350;
t274 = (-t322 + t350) * qJ(4) + (t338 * t340 * t343 + qJDD(3)) * pkin(3) + t277;
t278 = t340 * t298 + t338 * t324;
t323 = t340 * qJDD(2) - t338 * t352;
t353 = t338 * qJD(2);
t326 = qJD(3) * pkin(3) - qJ(4) * t353;
t333 = t340 ^ 2;
t275 = -t333 * t343 * pkin(3) + t323 * qJ(4) - qJD(3) * t326 + t278;
t335 = sin(pkin(8));
t354 = qJD(2) * t340;
t308 = t335 * t353 - t359 * t354;
t271 = t335 * t274 + t359 * t275 + t308 * t364;
t294 = t335 * t322 - t359 * t323;
t309 = (t335 * t340 + t359 * t338) * qJD(2);
t302 = qJD(3) * mrSges(5,1) - t309 * mrSges(5,3);
t288 = t308 * pkin(4) - t309 * qJ(5);
t342 = qJD(3) ^ 2;
t266 = -t342 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t308 * t288 + t271;
t303 = -qJD(3) * mrSges(6,1) + t309 * mrSges(6,2);
t351 = m(6) * t266 + qJDD(3) * mrSges(6,3) + qJD(3) * t303;
t289 = t308 * mrSges(6,1) - t309 * mrSges(6,3);
t355 = -t308 * mrSges(5,1) - t309 * mrSges(5,2) - t289;
t260 = m(5) * t271 - qJDD(3) * mrSges(5,2) - qJD(3) * t302 + t363 * t294 + t355 * t308 + t351;
t345 = t359 * t274 - t335 * t275;
t270 = t309 * t364 + t345;
t295 = t359 * t322 + t335 * t323;
t301 = -qJD(3) * mrSges(5,2) - t308 * mrSges(5,3);
t267 = -qJDD(3) * pkin(4) - t342 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t288) * t309 - t345;
t304 = -t308 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t347 = -m(6) * t267 + qJDD(3) * mrSges(6,1) + qJD(3) * t304;
t261 = m(5) * t270 + qJDD(3) * mrSges(5,1) + qJD(3) * t301 + t363 * t295 + t355 * t309 + t347;
t256 = t335 * t260 + t359 * t261;
t358 = t365 * qJD(3) + t360 * t308 - t361 * t309;
t357 = t360 * qJD(3) + t366 * t308 + t362 * t309;
t356 = t361 * qJD(3) - t362 * t308 + t367 * t309;
t321 = (-t340 * mrSges(4,1) + t338 * mrSges(4,2)) * qJD(2);
t327 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t353;
t328 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t354;
t348 = t359 * t260 - t335 * t261;
t349 = -t338 * (m(4) * t277 + qJDD(3) * mrSges(4,1) - t322 * mrSges(4,3) + qJD(3) * t328 - t321 * t353 + t256) + t340 * (m(4) * t278 - qJDD(3) * mrSges(4,2) + t323 * mrSges(4,3) - qJD(3) * t327 + t321 * t354 + t348);
t305 = -t339 * t325 + t341 * t334;
t346 = -qJDD(2) * pkin(2) - t305;
t276 = -t323 * pkin(3) + qJDD(4) + t326 * t353 + (-qJ(4) * t333 - pkin(6)) * t343 + t346;
t269 = -0.2e1 * qJD(5) * t309 + (qJD(3) * t308 - t295) * qJ(5) + (qJD(3) * t309 + t294) * pkin(4) + t276;
t264 = m(6) * t269 + t294 * mrSges(6,1) - t295 * mrSges(6,3) - t309 * t303 + t308 * t304;
t262 = m(5) * t276 + t294 * mrSges(5,1) + t295 * mrSges(5,2) + t308 * t301 + t309 * t302 + t264;
t297 = -t343 * pkin(6) + t346;
t344 = -m(4) * t297 + t323 * mrSges(4,1) - t322 * mrSges(4,2) - t327 * t353 + t328 * t354 - t262;
t312 = Ifges(4,5) * qJD(3) + (t338 * Ifges(4,1) + t340 * Ifges(4,4)) * qJD(2);
t311 = Ifges(4,6) * qJD(3) + (t338 * Ifges(4,4) + t340 * Ifges(4,2)) * qJD(2);
t263 = t295 * mrSges(6,2) + t309 * t289 - t347;
t255 = mrSges(5,2) * t276 + mrSges(6,2) * t267 - mrSges(5,3) * t270 - mrSges(6,3) * t269 - qJ(5) * t264 - t357 * qJD(3) + t361 * qJDD(3) - t362 * t294 + t367 * t295 + t358 * t308;
t254 = -mrSges(5,1) * t276 - mrSges(6,1) * t269 + mrSges(6,2) * t266 + mrSges(5,3) * t271 - pkin(4) * t264 + t356 * qJD(3) + t360 * qJDD(3) + t366 * t294 + t362 * t295 + t358 * t309;
t1 = [m(2) * t334 + t339 * (m(3) * t306 - t343 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t349) + t341 * (m(3) * t305 + qJDD(2) * mrSges(3,1) - t343 * mrSges(3,2) + t344); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t305 - mrSges(3,2) * t306 + t338 * (mrSges(4,2) * t297 - mrSges(4,3) * t277 + Ifges(4,1) * t322 + Ifges(4,4) * t323 + Ifges(4,5) * qJDD(3) - qJ(4) * t256 - qJD(3) * t311 - t335 * t254 + t359 * t255) + t340 * (-mrSges(4,1) * t297 + mrSges(4,3) * t278 + Ifges(4,4) * t322 + Ifges(4,2) * t323 + Ifges(4,6) * qJDD(3) - pkin(3) * t262 + qJ(4) * t348 + qJD(3) * t312 + t359 * t254 + t335 * t255) + pkin(2) * t344 + pkin(6) * t349; Ifges(4,5) * t322 + Ifges(4,6) * t323 + mrSges(4,1) * t277 - mrSges(4,2) * t278 + mrSges(5,1) * t270 - mrSges(5,2) * t271 - mrSges(6,1) * t267 + mrSges(6,3) * t266 - pkin(4) * t263 + qJ(5) * t351 + pkin(3) * t256 + t357 * t309 + (-qJ(5) * t289 + t356) * t308 + t361 * t295 + (-qJ(5) * mrSges(6,2) - t360) * t294 + (t338 * t311 - t340 * t312) * qJD(2) + (Ifges(4,3) - t365) * qJDD(3); t262; t263;];
tauJ = t1;

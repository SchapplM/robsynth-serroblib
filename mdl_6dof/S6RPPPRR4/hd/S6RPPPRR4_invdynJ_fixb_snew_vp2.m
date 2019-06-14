% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-05-05 13:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:45:56
% EndTime: 2019-05-05 13:45:57
% DurationCPUTime: 0.69s
% Computational Cost: add. (4011->181), mult. (6802->221), div. (0->0), fcn. (2919->8), ass. (0->78)
t350 = qJD(1) ^ 2;
t345 = sin(qJ(1));
t348 = cos(qJ(1));
t360 = -g(1) * t348 - g(2) * t345;
t356 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t360;
t370 = -pkin(1) - pkin(2);
t313 = t350 * t370 + t356;
t362 = g(1) * t345 - t348 * g(2);
t354 = -qJ(2) * t350 + qJDD(2) - t362;
t314 = qJDD(1) * t370 + t354;
t341 = sin(pkin(9));
t342 = cos(pkin(9));
t299 = t342 * t313 + t341 * t314;
t371 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t299;
t294 = (-pkin(3) - pkin(7)) * t350 + t371;
t344 = sin(qJ(5));
t347 = cos(qJ(5));
t366 = qJD(1) * qJD(5);
t363 = t347 * t366;
t326 = qJDD(1) * t344 + t363;
t364 = t344 * t366;
t327 = -qJDD(1) * t347 + t364;
t287 = (-t327 - t364) * pkin(8) + (-t326 - t363) * pkin(5) + t294;
t298 = -t341 * t313 + t314 * t342;
t297 = qJDD(1) * pkin(3) - qJ(4) * t350 + qJDD(4) - t298;
t295 = qJDD(1) * pkin(7) + t297;
t337 = g(3) + qJDD(3);
t291 = t344 * t295 + t347 * t337;
t325 = (-t344 * pkin(5) + t347 * pkin(8)) * qJD(1);
t349 = qJD(5) ^ 2;
t368 = qJD(1) * t344;
t289 = -pkin(5) * t349 + qJDD(5) * pkin(8) + t325 * t368 + t291;
t343 = sin(qJ(6));
t346 = cos(qJ(6));
t285 = t287 * t346 - t289 * t343;
t367 = qJD(1) * t347;
t322 = qJD(5) * t346 + t343 * t367;
t306 = qJD(6) * t322 + qJDD(5) * t343 + t327 * t346;
t323 = qJD(5) * t343 - t346 * t367;
t307 = -mrSges(7,1) * t322 + mrSges(7,2) * t323;
t330 = qJD(6) - t368;
t311 = -mrSges(7,2) * t330 + mrSges(7,3) * t322;
t321 = qJDD(6) - t326;
t283 = m(7) * t285 + mrSges(7,1) * t321 - mrSges(7,3) * t306 - t307 * t323 + t311 * t330;
t286 = t287 * t343 + t289 * t346;
t305 = -qJD(6) * t323 + qJDD(5) * t346 - t327 * t343;
t312 = mrSges(7,1) * t330 - mrSges(7,3) * t323;
t284 = m(7) * t286 - mrSges(7,2) * t321 + mrSges(7,3) * t305 + t307 * t322 - t312 * t330;
t276 = t346 * t283 + t343 * t284;
t296 = pkin(3) * t350 - t371;
t328 = -(qJD(5) * mrSges(6,2)) + mrSges(6,3) * t368;
t329 = (qJD(5) * mrSges(6,1)) + mrSges(6,3) * t367;
t372 = -m(5) * t296 + m(6) * t294 - mrSges(6,1) * t326 + t350 * mrSges(5,2) + t327 * mrSges(6,2) + t276 + (-t328 * t344 - t329 * t347) * qJD(1);
t369 = t337 * t344;
t361 = -t283 * t343 + t346 * t284;
t324 = (-t344 * mrSges(6,1) - t347 * mrSges(6,2)) * qJD(1);
t275 = m(6) * t291 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t326 - qJD(5) * t329 + t324 * t368 + t361;
t290 = t295 * t347 - t369;
t288 = -qJDD(5) * pkin(5) - pkin(8) * t349 + t369 + (-qJD(1) * t325 - t295) * t347;
t353 = -m(7) * t288 + t305 * mrSges(7,1) - mrSges(7,2) * t306 + t322 * t311 - t312 * t323;
t279 = m(6) * t290 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t327 + qJD(5) * t328 + t324 * t367 + t353;
t357 = t275 * t344 + t347 * t279;
t273 = m(5) * t297 - qJDD(1) * mrSges(5,2) - t350 * mrSges(5,3) + t357;
t272 = m(4) * t298 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t350 - t273;
t274 = m(4) * t299 - mrSges(4,1) * t350 + (mrSges(4,2) - mrSges(5,3)) * qJDD(1) + t372;
t358 = t272 * t342 + t274 * t341;
t301 = Ifges(7,4) * t323 + Ifges(7,2) * t322 + Ifges(7,6) * t330;
t302 = Ifges(7,1) * t323 + Ifges(7,4) * t322 + Ifges(7,5) * t330;
t351 = mrSges(7,1) * t285 - mrSges(7,2) * t286 + Ifges(7,5) * t306 + Ifges(7,6) * t305 + Ifges(7,3) * t321 + t301 * t323 - t322 * t302;
t319 = (Ifges(6,5) * qJD(5)) + (-t347 * Ifges(6,1) + t344 * Ifges(6,4)) * qJD(1);
t318 = (Ifges(6,6) * qJD(5)) + (-t347 * Ifges(6,4) + t344 * Ifges(6,2)) * qJD(1);
t316 = -qJDD(1) * pkin(1) + t354;
t315 = -pkin(1) * t350 + t356;
t300 = Ifges(7,5) * t323 + Ifges(7,6) * t322 + Ifges(7,3) * t330;
t278 = mrSges(7,2) * t288 - mrSges(7,3) * t285 + Ifges(7,1) * t306 + Ifges(7,4) * t305 + Ifges(7,5) * t321 + t300 * t322 - t301 * t330;
t277 = -mrSges(7,1) * t288 + mrSges(7,3) * t286 + Ifges(7,4) * t306 + Ifges(7,2) * t305 + Ifges(7,6) * t321 - t300 * t323 + t302 * t330;
t271 = m(3) * t316 - qJDD(1) * mrSges(3,1) - mrSges(3,3) * t350 + t358;
t1 = [qJ(2) * (m(3) * t315 - mrSges(3,1) * t350 - t272 * t341 + t274 * t342) - pkin(1) * t271 - mrSges(2,2) * t360 + mrSges(2,1) * t362 - pkin(2) * t358 + mrSges(3,3) * t315 - mrSges(3,1) * t316 + pkin(3) * t273 + pkin(7) * t357 - mrSges(4,1) * t298 + mrSges(4,2) * t299 - t347 * (mrSges(6,2) * t294 - mrSges(6,3) * t290 + Ifges(6,1) * t327 + Ifges(6,4) * t326 + Ifges(6,5) * qJDD(5) - pkin(8) * t276 - qJD(5) * t318 - t277 * t343 + t346 * t278) + t344 * (-mrSges(6,1) * t294 + mrSges(6,3) * t291 + Ifges(6,4) * t327 + Ifges(6,2) * t326 + Ifges(6,6) * qJDD(5) - pkin(5) * t276 + qJD(5) * t319 - t351) - mrSges(5,2) * t297 + mrSges(5,3) * t296 + (qJ(2) * mrSges(3,3) + Ifges(5,1) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1) + (qJDD(1) * mrSges(5,3) - t372) * qJ(4); t271; t275 * t347 - t279 * t344 + (m(4) + m(5)) * t337; t273; Ifges(6,5) * t327 + Ifges(6,6) * t326 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t290 - mrSges(6,2) * t291 + t343 * t278 + t346 * t277 + pkin(5) * t353 + pkin(8) * t361 + (-t347 * t318 - t344 * t319) * qJD(1); t351;];
tauJ  = t1;

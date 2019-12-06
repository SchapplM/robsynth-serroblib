% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPP1
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:25
% EndTime: 2019-12-05 16:06:27
% DurationCPUTime: 0.78s
% Computational Cost: add. (2682->188), mult. (5805->231), div. (0->0), fcn. (3539->8), ass. (0->81)
t370 = Ifges(5,1) + Ifges(6,1);
t364 = Ifges(5,4) - Ifges(6,5);
t363 = Ifges(5,5) + Ifges(6,4);
t369 = -Ifges(5,2) - Ifges(6,3);
t368 = -Ifges(6,2) - Ifges(5,3);
t362 = Ifges(5,6) - Ifges(6,6);
t336 = sin(qJ(3));
t338 = cos(qJ(3));
t351 = qJD(2) * qJD(3);
t348 = t338 * t351;
t319 = qJDD(2) * t336 + t348;
t320 = qJDD(2) * t338 - t336 * t351;
t333 = sin(pkin(8));
t360 = cos(pkin(8));
t300 = t360 * t319 + t333 * t320;
t309 = (t333 * t338 + t360 * t336) * qJD(2);
t305 = -qJD(3) * mrSges(6,1) + mrSges(6,2) * t309;
t367 = -mrSges(6,3) * t300 - t305 * t309;
t366 = -2 * qJD(4);
t365 = -mrSges(5,3) - mrSges(6,2);
t341 = qJD(2) ^ 2;
t334 = sin(pkin(7));
t335 = cos(pkin(7));
t321 = g(1) * t334 - g(2) * t335;
t322 = -g(1) * t335 - g(2) * t334;
t337 = sin(qJ(2));
t339 = cos(qJ(2));
t354 = t337 * t321 + t339 * t322;
t298 = -pkin(2) * t341 + qJDD(2) * pkin(6) + t354;
t332 = -g(3) + qJDD(1);
t283 = -t298 * t336 + t338 * t332;
t280 = (-t319 + t348) * qJ(4) + (t336 * t338 * t341 + qJDD(3)) * pkin(3) + t283;
t284 = t338 * t298 + t336 * t332;
t352 = t336 * qJD(2);
t323 = qJD(3) * pkin(3) - qJ(4) * t352;
t331 = t338 ^ 2;
t281 = -pkin(3) * t331 * t341 + qJ(4) * t320 - qJD(3) * t323 + t284;
t353 = qJD(2) * t338;
t308 = t333 * t352 - t360 * t353;
t277 = t333 * t280 + t360 * t281 + t308 * t366;
t299 = t319 * t333 - t360 * t320;
t304 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t309;
t292 = pkin(4) * t308 - qJ(5) * t309;
t340 = qJD(3) ^ 2;
t272 = -pkin(4) * t340 + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t292 * t308 + t277;
t350 = m(6) * t272 + qJDD(3) * mrSges(6,3) + qJD(3) * t305;
t293 = mrSges(6,1) * t308 - mrSges(6,3) * t309;
t355 = -mrSges(5,1) * t308 - mrSges(5,2) * t309 - t293;
t266 = m(5) * t277 - qJDD(3) * mrSges(5,2) - qJD(3) * t304 + t365 * t299 + t355 * t308 + t350;
t343 = t360 * t280 - t333 * t281;
t276 = t309 * t366 + t343;
t303 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t308;
t273 = -qJDD(3) * pkin(4) - t340 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t292) * t309 - t343;
t306 = -mrSges(6,2) * t308 + qJD(3) * mrSges(6,3);
t345 = -m(6) * t273 + qJDD(3) * mrSges(6,1) + qJD(3) * t306;
t267 = m(5) * t276 + qJDD(3) * mrSges(5,1) + qJD(3) * t303 + t365 * t300 + t355 * t309 + t345;
t262 = t333 * t266 + t360 * t267;
t358 = t368 * qJD(3) + t362 * t308 - t363 * t309;
t357 = t362 * qJD(3) + t369 * t308 + t364 * t309;
t356 = t363 * qJD(3) - t364 * t308 + t370 * t309;
t346 = t321 * t339 - t337 * t322;
t344 = -qJDD(2) * pkin(2) - t346;
t282 = -pkin(3) * t320 + qJDD(4) + t323 * t352 + (-qJ(4) * t331 - pkin(6)) * t341 + t344;
t275 = -0.2e1 * qJD(5) * t309 + (qJD(3) * t308 - t300) * qJ(5) + (qJD(3) * t309 + t299) * pkin(4) + t282;
t349 = m(6) * t275 + t299 * mrSges(6,1) + t308 * t306;
t347 = t360 * t266 - t267 * t333;
t342 = m(5) * t282 + mrSges(5,1) * t299 + t303 * t308 + t349;
t325 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t353;
t324 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t352;
t318 = (-t338 * mrSges(4,1) + t336 * mrSges(4,2)) * qJD(2);
t312 = Ifges(4,5) * qJD(3) + (t336 * Ifges(4,1) + t338 * Ifges(4,4)) * qJD(2);
t311 = Ifges(4,6) * qJD(3) + (t336 * Ifges(4,4) + t338 * Ifges(4,2)) * qJD(2);
t297 = -pkin(6) * t341 + t344;
t270 = t349 + t367;
t269 = mrSges(6,2) * t300 + t293 * t309 - t345;
t268 = (t304 - t305) * t309 + (mrSges(5,2) - mrSges(6,3)) * t300 + t342;
t261 = mrSges(5,2) * t282 + mrSges(6,2) * t273 - mrSges(5,3) * t276 - mrSges(6,3) * t275 - qJ(5) * t270 - t357 * qJD(3) + t363 * qJDD(3) - t364 * t299 + t370 * t300 + t358 * t308;
t260 = -mrSges(5,1) * t282 - mrSges(6,1) * t275 + mrSges(6,2) * t272 + mrSges(5,3) * t277 - pkin(4) * t270 + t356 * qJD(3) + t362 * qJDD(3) + t369 * t299 + t364 * t300 + t358 * t309;
t259 = m(4) * t284 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t320 - qJD(3) * t324 + t318 * t353 + t347;
t258 = m(4) * t283 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t319 + qJD(3) * t325 - t318 * t352 + t262;
t1 = [t258 * t338 + t259 * t336 + (m(2) + m(3)) * t332; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t346 - mrSges(3,2) * t354 + t336 * (mrSges(4,2) * t297 - mrSges(4,3) * t283 + Ifges(4,1) * t319 + Ifges(4,4) * t320 + Ifges(4,5) * qJDD(3) - qJ(4) * t262 - qJD(3) * t311 - t333 * t260 + t360 * t261) + t338 * (-mrSges(4,1) * t297 + mrSges(4,3) * t284 + Ifges(4,4) * t319 + Ifges(4,2) * t320 + Ifges(4,6) * qJDD(3) - pkin(3) * t268 + qJ(4) * t347 + qJD(3) * t312 + t360 * t260 + t333 * t261) + pkin(6) * (-t258 * t336 + t259 * t338) + (-m(4) * t297 + mrSges(4,1) * t320 - mrSges(4,2) * t319 - mrSges(5,2) * t300 - t304 * t309 - t342 + (-t324 * t336 + t325 * t338) * qJD(2) - t367) * pkin(2); Ifges(4,5) * t319 + Ifges(4,6) * t320 + mrSges(4,1) * t283 - mrSges(4,2) * t284 + mrSges(5,1) * t276 - mrSges(5,2) * t277 - mrSges(6,1) * t273 + mrSges(6,3) * t272 - pkin(4) * t269 + qJ(5) * t350 + pkin(3) * t262 + t357 * t309 + (-qJ(5) * t293 + t356) * t308 + t363 * t300 + (-qJ(5) * mrSges(6,2) - t362) * t299 + (t336 * t311 - t338 * t312) * qJD(2) + (Ifges(4,3) - t368) * qJDD(3); t268; t269;];
tauJ = t1;

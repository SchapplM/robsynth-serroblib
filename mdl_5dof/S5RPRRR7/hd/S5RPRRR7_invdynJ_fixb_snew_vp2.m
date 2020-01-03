% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:30
% EndTime: 2019-12-31 19:03:32
% DurationCPUTime: 1.07s
% Computational Cost: add. (9300->219), mult. (17871->279), div. (0->0), fcn. (11226->10), ass. (0->93)
t347 = sin(qJ(1));
t351 = cos(qJ(1));
t363 = t347 * g(1) - t351 * g(2);
t328 = qJDD(1) * pkin(1) + t363;
t353 = qJD(1) ^ 2;
t359 = -t351 * g(1) - t347 * g(2);
t330 = -t353 * pkin(1) + t359;
t342 = sin(pkin(9));
t343 = cos(pkin(9));
t307 = t343 * t328 - t342 * t330;
t296 = -qJDD(1) * pkin(2) - t353 * pkin(6) - t307;
t346 = sin(qJ(3));
t350 = cos(qJ(3));
t365 = qJD(1) * qJD(3);
t364 = t350 * t365;
t332 = t346 * qJDD(1) + t364;
t339 = t346 * t365;
t333 = t350 * qJDD(1) - t339;
t284 = (-t332 - t364) * pkin(7) + (-t333 + t339) * pkin(3) + t296;
t308 = t342 * t328 + t343 * t330;
t297 = -t353 * pkin(2) + qJDD(1) * pkin(6) + t308;
t341 = -g(3) + qJDD(2);
t291 = t350 * t297 + t346 * t341;
t331 = (-t350 * pkin(3) - t346 * pkin(7)) * qJD(1);
t352 = qJD(3) ^ 2;
t366 = t350 * qJD(1);
t288 = -t352 * pkin(3) + qJDD(3) * pkin(7) + t331 * t366 + t291;
t345 = sin(qJ(4));
t349 = cos(qJ(4));
t271 = t349 * t284 - t345 * t288;
t367 = t346 * qJD(1);
t326 = t349 * qJD(3) - t345 * t367;
t304 = t326 * qJD(4) + t345 * qJDD(3) + t349 * t332;
t325 = qJDD(4) - t333;
t327 = t345 * qJD(3) + t349 * t367;
t337 = qJD(4) - t366;
t268 = (t326 * t337 - t304) * pkin(8) + (t326 * t327 + t325) * pkin(4) + t271;
t272 = t345 * t284 + t349 * t288;
t303 = -t327 * qJD(4) + t349 * qJDD(3) - t345 * t332;
t312 = t337 * pkin(4) - t327 * pkin(8);
t324 = t326 ^ 2;
t269 = -t324 * pkin(4) + t303 * pkin(8) - t337 * t312 + t272;
t344 = sin(qJ(5));
t348 = cos(qJ(5));
t266 = t348 * t268 - t344 * t269;
t305 = t348 * t326 - t344 * t327;
t278 = t305 * qJD(5) + t344 * t303 + t348 * t304;
t306 = t344 * t326 + t348 * t327;
t289 = -t305 * mrSges(6,1) + t306 * mrSges(6,2);
t336 = qJD(5) + t337;
t292 = -t336 * mrSges(6,2) + t305 * mrSges(6,3);
t321 = qJDD(5) + t325;
t263 = m(6) * t266 + t321 * mrSges(6,1) - t278 * mrSges(6,3) - t306 * t289 + t336 * t292;
t267 = t344 * t268 + t348 * t269;
t277 = -t306 * qJD(5) + t348 * t303 - t344 * t304;
t293 = t336 * mrSges(6,1) - t306 * mrSges(6,3);
t264 = m(6) * t267 - t321 * mrSges(6,2) + t277 * mrSges(6,3) + t305 * t289 - t336 * t293;
t256 = t348 * t263 + t344 * t264;
t329 = (-t350 * mrSges(4,1) + t346 * mrSges(4,2)) * qJD(1);
t334 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t367;
t309 = -t326 * mrSges(5,1) + t327 * mrSges(5,2);
t310 = -t337 * mrSges(5,2) + t326 * mrSges(5,3);
t254 = m(5) * t271 + t325 * mrSges(5,1) - t304 * mrSges(5,3) - t327 * t309 + t337 * t310 + t256;
t311 = t337 * mrSges(5,1) - t327 * mrSges(5,3);
t360 = -t344 * t263 + t348 * t264;
t255 = m(5) * t272 - t325 * mrSges(5,2) + t303 * mrSges(5,3) + t326 * t309 - t337 * t311 + t360;
t361 = -t345 * t254 + t349 * t255;
t251 = m(4) * t291 - qJDD(3) * mrSges(4,2) + t333 * mrSges(4,3) - qJD(3) * t334 + t329 * t366 + t361;
t290 = -t346 * t297 + t350 * t341;
t335 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t366;
t287 = -qJDD(3) * pkin(3) - t352 * pkin(7) + t331 * t367 - t290;
t270 = -t303 * pkin(4) - t324 * pkin(8) + t327 * t312 + t287;
t358 = m(6) * t270 - t277 * mrSges(6,1) + t278 * mrSges(6,2) - t305 * t292 + t306 * t293;
t355 = -m(5) * t287 + t303 * mrSges(5,1) - t304 * mrSges(5,2) + t326 * t310 - t327 * t311 - t358;
t259 = m(4) * t290 + qJDD(3) * mrSges(4,1) - t332 * mrSges(4,3) + qJD(3) * t335 - t329 * t367 + t355;
t362 = t350 * t251 - t346 * t259;
t252 = t349 * t254 + t345 * t255;
t282 = Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t336;
t283 = Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t336;
t357 = -mrSges(6,1) * t266 + mrSges(6,2) * t267 - Ifges(6,5) * t278 - Ifges(6,6) * t277 - Ifges(6,3) * t321 - t306 * t282 + t305 * t283;
t356 = -m(4) * t296 + t333 * mrSges(4,1) - t332 * mrSges(4,2) - t334 * t367 + t335 * t366 - t252;
t299 = Ifges(5,4) * t327 + Ifges(5,2) * t326 + Ifges(5,6) * t337;
t300 = Ifges(5,1) * t327 + Ifges(5,4) * t326 + Ifges(5,5) * t337;
t354 = mrSges(5,1) * t271 - mrSges(5,2) * t272 + Ifges(5,5) * t304 + Ifges(5,6) * t303 + Ifges(5,3) * t325 + pkin(4) * t256 + t327 * t299 - t326 * t300 - t357;
t320 = Ifges(4,5) * qJD(3) + (t346 * Ifges(4,1) + t350 * Ifges(4,4)) * qJD(1);
t319 = Ifges(4,6) * qJD(3) + (t346 * Ifges(4,4) + t350 * Ifges(4,2)) * qJD(1);
t298 = Ifges(5,5) * t327 + Ifges(5,6) * t326 + Ifges(5,3) * t337;
t281 = Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t336;
t258 = mrSges(6,2) * t270 - mrSges(6,3) * t266 + Ifges(6,1) * t278 + Ifges(6,4) * t277 + Ifges(6,5) * t321 + t305 * t281 - t336 * t282;
t257 = -mrSges(6,1) * t270 + mrSges(6,3) * t267 + Ifges(6,4) * t278 + Ifges(6,2) * t277 + Ifges(6,6) * t321 - t306 * t281 + t336 * t283;
t249 = mrSges(5,2) * t287 - mrSges(5,3) * t271 + Ifges(5,1) * t304 + Ifges(5,4) * t303 + Ifges(5,5) * t325 - pkin(8) * t256 - t344 * t257 + t348 * t258 + t326 * t298 - t337 * t299;
t248 = -mrSges(5,1) * t287 + mrSges(5,3) * t272 + Ifges(5,4) * t304 + Ifges(5,2) * t303 + Ifges(5,6) * t325 - pkin(4) * t358 + pkin(8) * t360 + t348 * t257 + t344 * t258 - t327 * t298 + t337 * t300;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t363 - mrSges(2,2) * t359 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t307 - mrSges(3,2) * t308 + t346 * (mrSges(4,2) * t296 - mrSges(4,3) * t290 + Ifges(4,1) * t332 + Ifges(4,4) * t333 + Ifges(4,5) * qJDD(3) - pkin(7) * t252 - qJD(3) * t319 - t345 * t248 + t349 * t249) + t350 * (-mrSges(4,1) * t296 + mrSges(4,3) * t291 + Ifges(4,4) * t332 + Ifges(4,2) * t333 + Ifges(4,6) * qJDD(3) - pkin(3) * t252 + qJD(3) * t320 - t354) + pkin(2) * t356 + pkin(6) * t362 + pkin(1) * (t342 * (m(3) * t308 - t353 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t362) + t343 * (m(3) * t307 + qJDD(1) * mrSges(3,1) - t353 * mrSges(3,2) + t356)); m(3) * t341 + t346 * t251 + t350 * t259; Ifges(4,5) * t332 + Ifges(4,6) * t333 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t290 - mrSges(4,2) * t291 + t345 * t249 + t349 * t248 + pkin(3) * t355 + pkin(7) * t361 + (t346 * t319 - t350 * t320) * qJD(1); t354; -t357;];
tauJ = t1;

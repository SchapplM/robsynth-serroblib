% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:50
% EndTime: 2024-09-28 18:07:51
% DurationCPUTime: 0.38s
% Computational Cost: add. (8625->118), mult. (11764->168), div. (0->0), fcn. (9011->14), ass. (0->74)
t309 = sin(pkin(11));
t312 = cos(pkin(11));
t296 = g(1) * t309 - g(2) * t312;
t308 = -g(3) + qJDD(1);
t311 = sin(pkin(5));
t314 = cos(pkin(5));
t340 = t296 * t314 + t308 * t311;
t310 = sin(pkin(6));
t339 = pkin(10) * t310;
t307 = qJD(2) + qJD(3);
t304 = qJD(4) + t307;
t337 = t304 * t310;
t315 = sin(qJ(5));
t335 = t310 * t315;
t319 = cos(qJ(5));
t334 = t310 * t319;
t297 = -g(1) * t312 - g(2) * t309;
t318 = sin(qJ(2));
t322 = cos(qJ(2));
t280 = -t297 * t318 + t340 * t322;
t278 = qJDD(2) * pkin(2) + t280;
t281 = t322 * t297 + t340 * t318;
t323 = qJD(2) ^ 2;
t279 = -pkin(2) * t323 + t281;
t317 = sin(qJ(3));
t321 = cos(qJ(3));
t273 = t321 * t278 - t279 * t317;
t305 = t307 ^ 2;
t306 = qJDD(2) + qJDD(3);
t270 = pkin(3) * t306 + t273;
t274 = t317 * t278 + t321 * t279;
t271 = -pkin(3) * t305 + t274;
t316 = sin(qJ(4));
t320 = cos(qJ(4));
t265 = t320 * t270 - t271 * t316;
t302 = t304 ^ 2;
t303 = qJDD(4) + t306;
t266 = t316 * t270 + t320 * t271;
t263 = -pkin(4) * t302 + t303 * t339 + t266;
t262 = pkin(4) * t303 + t302 * t339 + t265;
t290 = -t296 * t311 + t308 * t314;
t313 = cos(pkin(6));
t326 = t262 * t313 + t290 * t310;
t258 = -t315 * t263 + t326 * t319;
t295 = t304 * t313 + qJD(5);
t329 = t304 * t334;
t286 = -mrSges(6,2) * t295 + mrSges(6,3) * t329;
t287 = (-mrSges(6,1) * t319 + mrSges(6,2) * t315) * t337;
t331 = qJD(5) * t304;
t288 = (t303 * t315 + t319 * t331) * t310;
t294 = t303 * t313 + qJDD(5);
t330 = t304 * t335;
t256 = m(6) * t258 + mrSges(6,1) * t294 - mrSges(6,3) * t288 + t286 * t295 - t287 * t330;
t259 = t319 * t263 + t326 * t315;
t285 = mrSges(6,1) * t295 - mrSges(6,3) * t330;
t289 = (t303 * t319 - t315 * t331) * t310;
t257 = m(6) * t259 - mrSges(6,2) * t294 + mrSges(6,3) * t289 - t285 * t295 + t287 * t329;
t261 = -t262 * t310 + t290 * t313;
t260 = m(6) * t261 - mrSges(6,1) * t289 + mrSges(6,2) * t288 + (t285 * t315 - t286 * t319) * t337;
t327 = -t260 * t310 + (t256 * t319 + t257 * t315) * t313;
t244 = m(5) * t265 + mrSges(5,1) * t303 - mrSges(5,2) * t302 + t327;
t328 = -t256 * t315 + t319 * t257;
t247 = m(5) * t266 - mrSges(5,1) * t302 - mrSges(5,2) * t303 + t328;
t332 = t320 * t244 + t316 * t247;
t240 = m(4) * t273 + mrSges(4,1) * t306 - mrSges(4,2) * t305 + t332;
t241 = m(4) * t274 - mrSges(4,1) * t305 - mrSges(4,2) * t306 - t244 * t316 + t247 * t320;
t333 = t321 * t240 + t317 * t241;
t283 = Ifges(6,6) * t295 + (Ifges(6,4) * t315 + Ifges(6,2) * t319) * t337;
t284 = Ifges(6,5) * t295 + (Ifges(6,1) * t315 + Ifges(6,4) * t319) * t337;
t250 = mrSges(6,1) * t258 - mrSges(6,2) * t259 + Ifges(6,5) * t288 + Ifges(6,6) * t289 + Ifges(6,3) * t294 + (t283 * t315 - t284 * t319) * t337;
t282 = Ifges(6,3) * t295 + (Ifges(6,5) * t315 + Ifges(6,6) * t319) * t337;
t325 = -mrSges(5,2) * t266 + pkin(4) * t327 + t313 * t250 + t328 * t339 + (mrSges(6,2) * t261 - mrSges(6,3) * t258 + Ifges(6,1) * t288 + Ifges(6,4) * t289 + Ifges(6,5) * t294 + t282 * t329 - t283 * t295) * t335 + (-mrSges(6,1) * t261 + mrSges(6,3) * t259 + Ifges(6,4) * t288 + Ifges(6,2) * t289 + Ifges(6,6) * t294 - t282 * t330 + t284 * t295) * t334 + mrSges(5,1) * t265 + Ifges(5,3) * t303;
t324 = mrSges(4,1) * t273 - mrSges(4,2) * t274 + Ifges(4,3) * t306 + pkin(3) * t332 + t325;
t1 = [m(2) * t308 + (t318 * (m(3) * t281 - mrSges(3,1) * t323 - qJDD(2) * mrSges(3,2) - t240 * t317 + t241 * t321) + t322 * (m(3) * t280 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t323 + t333)) * t311 + (t256 * t334 + t257 * t335 + t260 * t313 + (m(3) + m(4) + m(5)) * t290) * t314; mrSges(3,1) * t280 - mrSges(3,2) * t281 + Ifges(3,3) * qJDD(2) + pkin(2) * t333 + t324; t324; t325; t250;];
tauJ = t1;

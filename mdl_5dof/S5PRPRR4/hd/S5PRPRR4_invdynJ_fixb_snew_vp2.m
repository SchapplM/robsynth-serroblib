% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:50:07
% EndTime: 2019-12-05 15:50:09
% DurationCPUTime: 0.59s
% Computational Cost: add. (3389->162), mult. (6131->213), div. (0->0), fcn. (4190->12), ass. (0->80)
t307 = sin(pkin(9));
t310 = cos(pkin(9));
t297 = g(1) * t307 - g(2) * t310;
t311 = cos(pkin(5));
t336 = t297 * t311;
t308 = sin(pkin(5));
t314 = sin(qJ(2));
t335 = t308 * t314;
t317 = cos(qJ(2));
t334 = t308 * t317;
t305 = -g(3) + qJDD(1);
t323 = -t297 * t308 + t311 * t305;
t281 = qJDD(3) + t323;
t316 = cos(qJ(4));
t333 = t316 * t281;
t298 = -g(1) * t310 - g(2) * t307;
t267 = -t298 * t314 + t305 * t334 + t317 * t336;
t265 = qJDD(2) * pkin(2) + t267;
t268 = t317 * t298 + t305 * t335 + t314 * t336;
t319 = qJD(2) ^ 2;
t266 = -pkin(2) * t319 + t268;
t306 = sin(pkin(10));
t309 = cos(pkin(10));
t261 = t306 * t265 + t309 * t266;
t259 = -pkin(3) * t319 + qJDD(2) * pkin(7) + t261;
t313 = sin(qJ(4));
t256 = t316 * t259 + t313 * t281;
t293 = (-mrSges(5,1) * t316 + mrSges(5,2) * t313) * qJD(2);
t329 = qJD(2) * qJD(4);
t327 = t313 * t329;
t296 = qJDD(2) * t316 - t327;
t331 = qJD(2) * t313;
t299 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t331;
t294 = (-pkin(4) * t316 - pkin(8) * t313) * qJD(2);
t318 = qJD(4) ^ 2;
t330 = t316 * qJD(2);
t253 = -pkin(4) * t318 + qJDD(4) * pkin(8) + t294 * t330 + t256;
t260 = t309 * t265 - t306 * t266;
t258 = -qJDD(2) * pkin(3) - t319 * pkin(7) - t260;
t326 = t316 * t329;
t295 = qJDD(2) * t313 + t326;
t254 = (-t295 - t326) * pkin(8) + (-t296 + t327) * pkin(4) + t258;
t312 = sin(qJ(5));
t315 = cos(qJ(5));
t250 = -t253 * t312 + t254 * t315;
t291 = qJD(4) * t315 - t312 * t331;
t275 = qJD(5) * t291 + qJDD(4) * t312 + t295 * t315;
t292 = qJD(4) * t312 + t315 * t331;
t276 = -mrSges(6,1) * t291 + mrSges(6,2) * t292;
t303 = qJD(5) - t330;
t279 = -mrSges(6,2) * t303 + mrSges(6,3) * t291;
t289 = qJDD(5) - t296;
t248 = m(6) * t250 + mrSges(6,1) * t289 - mrSges(6,3) * t275 - t276 * t292 + t279 * t303;
t251 = t253 * t315 + t254 * t312;
t274 = -qJD(5) * t292 + qJDD(4) * t315 - t295 * t312;
t280 = mrSges(6,1) * t303 - mrSges(6,3) * t292;
t249 = m(6) * t251 - mrSges(6,2) * t289 + mrSges(6,3) * t274 + t276 * t291 - t280 * t303;
t324 = -t248 * t312 + t315 * t249;
t241 = m(5) * t256 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t296 - qJD(4) * t299 + t293 * t330 + t324;
t255 = -t259 * t313 + t333;
t300 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t330;
t252 = -qJDD(4) * pkin(4) - t318 * pkin(8) - t333 + (qJD(2) * t294 + t259) * t313;
t322 = -m(6) * t252 + t274 * mrSges(6,1) - mrSges(6,2) * t275 + t291 * t279 - t280 * t292;
t246 = m(5) * t255 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t295 + qJD(4) * t300 - t293 * t331 + t322;
t325 = t316 * t241 - t246 * t313;
t236 = m(4) * t261 - mrSges(4,1) * t319 - qJDD(2) * mrSges(4,2) + t325;
t242 = t248 * t315 + t249 * t312;
t321 = -m(5) * t258 + t296 * mrSges(5,1) - mrSges(5,2) * t295 - t299 * t331 + t300 * t330 - t242;
t238 = m(4) * t260 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t319 + t321;
t332 = t306 * t236 + t309 * t238;
t328 = m(4) * t281 + t313 * t241 + t316 * t246;
t270 = Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t303;
t271 = Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t303;
t320 = mrSges(6,1) * t250 - mrSges(6,2) * t251 + Ifges(6,5) * t275 + Ifges(6,6) * t274 + Ifges(6,3) * t289 + t270 * t292 - t271 * t291;
t286 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t313 + Ifges(5,4) * t316) * qJD(2);
t285 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t313 + Ifges(5,2) * t316) * qJD(2);
t269 = Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * t303;
t244 = mrSges(6,2) * t252 - mrSges(6,3) * t250 + Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t289 + t269 * t291 - t270 * t303;
t243 = -mrSges(6,1) * t252 + mrSges(6,3) * t251 + Ifges(6,4) * t275 + Ifges(6,2) * t274 + Ifges(6,6) * t289 - t269 * t292 + t271 * t303;
t1 = [m(2) * t305 + (m(3) * t268 - mrSges(3,1) * t319 - qJDD(2) * mrSges(3,2) + t236 * t309 - t238 * t306) * t335 + (m(3) * t267 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t319 + t332) * t334 + t311 * (m(3) * t323 + t328); mrSges(3,1) * t267 - mrSges(3,2) * t268 + mrSges(4,1) * t260 - mrSges(4,2) * t261 + t313 * (mrSges(5,2) * t258 - mrSges(5,3) * t255 + Ifges(5,1) * t295 + Ifges(5,4) * t296 + Ifges(5,5) * qJDD(4) - pkin(8) * t242 - qJD(4) * t285 - t243 * t312 + t244 * t315) + t316 * (-mrSges(5,1) * t258 + mrSges(5,3) * t256 + Ifges(5,4) * t295 + Ifges(5,2) * t296 + Ifges(5,6) * qJDD(4) - pkin(4) * t242 + qJD(4) * t286 - t320) + pkin(3) * t321 + pkin(7) * t325 + pkin(2) * t332 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2); t328; Ifges(5,5) * t295 + Ifges(5,6) * t296 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t255 - mrSges(5,2) * t256 + t312 * t244 + t315 * t243 + pkin(4) * t322 + pkin(8) * t324 + (t285 * t313 - t286 * t316) * qJD(2); t320;];
tauJ = t1;

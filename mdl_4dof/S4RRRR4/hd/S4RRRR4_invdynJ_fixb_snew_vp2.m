% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:48
% EndTime: 2019-12-31 17:25:50
% DurationCPUTime: 0.88s
% Computational Cost: add. (5542->202), mult. (11325->265), div. (0->0), fcn. (7293->8), ass. (0->83)
t301 = sin(qJ(3));
t302 = sin(qJ(2));
t305 = cos(qJ(3));
t306 = cos(qJ(2));
t283 = (t302 * t301 - t306 * t305) * qJD(1);
t308 = qJD(1) ^ 2;
t323 = pkin(2) * t308;
t303 = sin(qJ(1));
t307 = cos(qJ(1));
t315 = -t307 * g(1) - t303 * g(2);
t286 = -t308 * pkin(1) + qJDD(1) * pkin(5) + t315;
t322 = t302 * t286;
t319 = qJD(1) * qJD(2);
t289 = t302 * qJDD(1) + t306 * t319;
t257 = qJDD(2) * pkin(2) - t289 * pkin(6) - t322 + (pkin(6) * t319 + t302 * t323 - g(3)) * t306;
t275 = -t302 * g(3) + t306 * t286;
t290 = t306 * qJDD(1) - t302 * t319;
t321 = qJD(1) * t302;
t293 = qJD(2) * pkin(2) - pkin(6) * t321;
t299 = t306 ^ 2;
t258 = t290 * pkin(6) - qJD(2) * t293 - t299 * t323 + t275;
t245 = t301 * t257 + t305 * t258;
t284 = (t306 * t301 + t302 * t305) * qJD(1);
t262 = -t284 * qJD(3) - t301 * t289 + t305 * t290;
t270 = t283 * mrSges(4,1) + t284 * mrSges(4,2);
t298 = qJD(2) + qJD(3);
t277 = t298 * mrSges(4,1) - t284 * mrSges(4,3);
t297 = qJDD(2) + qJDD(3);
t263 = -t283 * qJD(3) + t305 * t289 + t301 * t290;
t318 = t303 * g(1) - t307 * g(2);
t313 = -qJDD(1) * pkin(1) - t318;
t264 = -t290 * pkin(2) + t293 * t321 + (-pkin(6) * t299 - pkin(5)) * t308 + t313;
t240 = (t283 * t298 - t263) * pkin(7) + (t284 * t298 - t262) * pkin(3) + t264;
t271 = t283 * pkin(3) - t284 * pkin(7);
t296 = t298 ^ 2;
t242 = -t296 * pkin(3) + t297 * pkin(7) - t283 * t271 + t245;
t300 = sin(qJ(4));
t304 = cos(qJ(4));
t238 = t304 * t240 - t300 * t242;
t272 = -t300 * t284 + t304 * t298;
t248 = t272 * qJD(4) + t304 * t263 + t300 * t297;
t273 = t304 * t284 + t300 * t298;
t255 = -t272 * mrSges(5,1) + t273 * mrSges(5,2);
t261 = qJDD(4) - t262;
t279 = qJD(4) + t283;
t265 = -t279 * mrSges(5,2) + t272 * mrSges(5,3);
t235 = m(5) * t238 + t261 * mrSges(5,1) - t248 * mrSges(5,3) - t273 * t255 + t279 * t265;
t239 = t300 * t240 + t304 * t242;
t247 = -t273 * qJD(4) - t300 * t263 + t304 * t297;
t266 = t279 * mrSges(5,1) - t273 * mrSges(5,3);
t236 = m(5) * t239 - t261 * mrSges(5,2) + t247 * mrSges(5,3) + t272 * t255 - t279 * t266;
t316 = -t300 * t235 + t304 * t236;
t223 = m(4) * t245 - t297 * mrSges(4,2) + t262 * mrSges(4,3) - t283 * t270 - t298 * t277 + t316;
t244 = t305 * t257 - t301 * t258;
t276 = -t298 * mrSges(4,2) - t283 * mrSges(4,3);
t241 = -t297 * pkin(3) - t296 * pkin(7) + t284 * t271 - t244;
t312 = -m(5) * t241 + t247 * mrSges(5,1) - t248 * mrSges(5,2) + t272 * t265 - t273 * t266;
t231 = m(4) * t244 + t297 * mrSges(4,1) - t263 * mrSges(4,3) - t284 * t270 + t298 * t276 + t312;
t220 = t301 * t223 + t305 * t231;
t225 = t304 * t235 + t300 * t236;
t320 = qJD(1) * t306;
t317 = t305 * t223 - t301 * t231;
t249 = Ifges(5,5) * t273 + Ifges(5,6) * t272 + Ifges(5,3) * t279;
t251 = Ifges(5,1) * t273 + Ifges(5,4) * t272 + Ifges(5,5) * t279;
t228 = -mrSges(5,1) * t241 + mrSges(5,3) * t239 + Ifges(5,4) * t248 + Ifges(5,2) * t247 + Ifges(5,6) * t261 - t273 * t249 + t279 * t251;
t250 = Ifges(5,4) * t273 + Ifges(5,2) * t272 + Ifges(5,6) * t279;
t229 = mrSges(5,2) * t241 - mrSges(5,3) * t238 + Ifges(5,1) * t248 + Ifges(5,4) * t247 + Ifges(5,5) * t261 + t272 * t249 - t279 * t250;
t268 = Ifges(4,4) * t284 - Ifges(4,2) * t283 + Ifges(4,6) * t298;
t269 = Ifges(4,1) * t284 - Ifges(4,4) * t283 + Ifges(4,5) * t298;
t311 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t263 + Ifges(4,6) * t262 + Ifges(4,3) * t297 + pkin(3) * t312 + pkin(7) * t316 + t304 * t228 + t300 * t229 + t284 * t268 + t283 * t269;
t310 = m(4) * t264 - t262 * mrSges(4,1) + t263 * mrSges(4,2) + t283 * t276 + t284 * t277 + t225;
t309 = mrSges(5,1) * t238 - mrSges(5,2) * t239 + Ifges(5,5) * t248 + Ifges(5,6) * t247 + Ifges(5,3) * t261 + t273 * t250 - t272 * t251;
t292 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t320;
t291 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t321;
t288 = (-t306 * mrSges(3,1) + t302 * mrSges(3,2)) * qJD(1);
t285 = -t308 * pkin(5) + t313;
t282 = Ifges(3,5) * qJD(2) + (t302 * Ifges(3,1) + t306 * Ifges(3,4)) * qJD(1);
t281 = Ifges(3,6) * qJD(2) + (t302 * Ifges(3,4) + t306 * Ifges(3,2)) * qJD(1);
t274 = -t306 * g(3) - t322;
t267 = Ifges(4,5) * t284 - Ifges(4,6) * t283 + Ifges(4,3) * t298;
t219 = -mrSges(4,1) * t264 + mrSges(4,3) * t245 + Ifges(4,4) * t263 + Ifges(4,2) * t262 + Ifges(4,6) * t297 - pkin(3) * t225 - t284 * t267 + t298 * t269 - t309;
t218 = mrSges(4,2) * t264 - mrSges(4,3) * t244 + Ifges(4,1) * t263 + Ifges(4,4) * t262 + Ifges(4,5) * t297 - pkin(7) * t225 - t300 * t228 + t304 * t229 - t283 * t267 - t298 * t268;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t318 - mrSges(2,2) * t315 + t302 * (mrSges(3,2) * t285 - mrSges(3,3) * t274 + Ifges(3,1) * t289 + Ifges(3,4) * t290 + Ifges(3,5) * qJDD(2) - pkin(6) * t220 - qJD(2) * t281 + t305 * t218 - t301 * t219) + t306 * (-mrSges(3,1) * t285 + mrSges(3,3) * t275 + Ifges(3,4) * t289 + Ifges(3,2) * t290 + Ifges(3,6) * qJDD(2) - pkin(2) * t310 + pkin(6) * t317 + qJD(2) * t282 + t301 * t218 + t305 * t219) + pkin(1) * (-m(3) * t285 + t290 * mrSges(3,1) - t289 * mrSges(3,2) + (-t291 * t302 + t292 * t306) * qJD(1) - t310) + pkin(5) * (t306 * (m(3) * t275 - qJDD(2) * mrSges(3,2) + t290 * mrSges(3,3) - qJD(2) * t291 + t288 * t320 + t317) - t302 * (m(3) * t274 + qJDD(2) * mrSges(3,1) - t289 * mrSges(3,3) + qJD(2) * t292 - t288 * t321 + t220)); (t302 * t281 - t306 * t282) * qJD(1) + t311 + Ifges(3,3) * qJDD(2) + pkin(2) * t220 + mrSges(3,1) * t274 - mrSges(3,2) * t275 + Ifges(3,5) * t289 + Ifges(3,6) * t290; t311; t309;];
tauJ = t1;

% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:51:54
% EndTime: 2019-05-05 03:52:19
% DurationCPUTime: 14.44s
% Computational Cost: add. (258799->336), mult. (531642->418), div. (0->0), fcn. (373593->12), ass. (0->133)
t295 = sin(pkin(10));
t298 = cos(pkin(10));
t285 = g(1) * t295 - g(2) * t298;
t286 = -g(1) * t298 - g(2) * t295;
t293 = -g(3) + qJDD(1);
t304 = cos(qJ(2));
t299 = cos(pkin(6));
t302 = sin(qJ(2));
t328 = t299 * t302;
t296 = sin(pkin(6));
t329 = t296 * t302;
t241 = t285 * t328 + t304 * t286 + t293 * t329;
t306 = qJD(2) ^ 2;
t235 = -pkin(2) * t306 + qJDD(2) * pkin(8) + t241;
t258 = -t285 * t296 + t293 * t299;
t301 = sin(qJ(3));
t303 = cos(qJ(3));
t226 = t303 * t235 + t301 * t258;
t280 = (-pkin(3) * t303 - qJ(4) * t301) * qJD(2);
t305 = qJD(3) ^ 2;
t324 = qJD(2) * t303;
t200 = -pkin(3) * t305 + qJDD(3) * qJ(4) + t280 * t324 + t226;
t240 = -t302 * t286 + (t285 * t299 + t293 * t296) * t304;
t234 = -qJDD(2) * pkin(2) - t306 * pkin(8) - t240;
t323 = qJD(2) * qJD(3);
t321 = t303 * t323;
t282 = qJDD(2) * t301 + t321;
t292 = t301 * t323;
t283 = qJDD(2) * t303 - t292;
t213 = (-t282 - t321) * qJ(4) + (-t283 + t292) * pkin(3) + t234;
t294 = sin(pkin(11));
t297 = cos(pkin(11));
t325 = qJD(2) * t301;
t275 = qJD(3) * t294 + t297 * t325;
t193 = -0.2e1 * qJD(4) * t275 - t294 * t200 + t297 * t213;
t256 = qJDD(3) * t294 + t282 * t297;
t274 = qJD(3) * t297 - t294 * t325;
t190 = (-t274 * t324 - t256) * pkin(9) + (t274 * t275 - t283) * pkin(4) + t193;
t194 = 0.2e1 * qJD(4) * t274 + t297 * t200 + t294 * t213;
t255 = qJDD(3) * t297 - t282 * t294;
t257 = -pkin(4) * t324 - pkin(9) * t275;
t273 = t274 ^ 2;
t192 = -pkin(4) * t273 + pkin(9) * t255 + t257 * t324 + t194;
t300 = sin(qJ(5));
t334 = cos(qJ(5));
t186 = t300 * t190 + t334 * t192;
t248 = t300 * t274 + t334 * t275;
t214 = t248 * qJD(5) - t334 * t255 + t300 * t256;
t291 = qJD(5) - t324;
t237 = mrSges(6,1) * t291 - mrSges(6,3) * t248;
t247 = -t334 * t274 + t300 * t275;
t277 = qJDD(5) - t283;
t227 = pkin(5) * t247 - qJ(6) * t248;
t290 = t291 ^ 2;
t181 = -pkin(5) * t290 + qJ(6) * t277 + 0.2e1 * qJD(6) * t291 - t227 * t247 + t186;
t238 = -mrSges(7,1) * t291 + mrSges(7,2) * t248;
t322 = m(7) * t181 + t277 * mrSges(7,3) + t291 * t238;
t228 = mrSges(7,1) * t247 - mrSges(7,3) * t248;
t326 = -mrSges(6,1) * t247 - mrSges(6,2) * t248 - t228;
t332 = -mrSges(6,3) - mrSges(7,2);
t168 = m(6) * t186 - t277 * mrSges(6,2) + t332 * t214 - t291 * t237 + t326 * t247 + t322;
t185 = t334 * t190 - t300 * t192;
t215 = -t247 * qJD(5) + t300 * t255 + t334 * t256;
t236 = -mrSges(6,2) * t291 - mrSges(6,3) * t247;
t183 = -t277 * pkin(5) - t290 * qJ(6) + t248 * t227 + qJDD(6) - t185;
t239 = -mrSges(7,2) * t247 + mrSges(7,3) * t291;
t318 = -m(7) * t183 + t277 * mrSges(7,1) + t291 * t239;
t170 = m(6) * t185 + t277 * mrSges(6,1) + t332 * t215 + t291 * t236 + t326 * t248 + t318;
t163 = t300 * t168 + t334 * t170;
t249 = -mrSges(5,1) * t274 + mrSges(5,2) * t275;
t253 = mrSges(5,2) * t324 + mrSges(5,3) * t274;
t161 = m(5) * t193 - mrSges(5,1) * t283 - mrSges(5,3) * t256 - t249 * t275 - t253 * t324 + t163;
t254 = -mrSges(5,1) * t324 - mrSges(5,3) * t275;
t319 = t334 * t168 - t170 * t300;
t162 = m(5) * t194 + mrSges(5,2) * t283 + mrSges(5,3) * t255 + t249 * t274 + t254 * t324 + t319;
t159 = -t161 * t294 + t297 * t162;
t281 = (-mrSges(4,1) * t303 + mrSges(4,2) * t301) * qJD(2);
t287 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t325;
t157 = m(4) * t226 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t283 - qJD(3) * t287 + t281 * t324 + t159;
t225 = -t301 * t235 + t303 * t258;
t199 = -qJDD(3) * pkin(3) - t305 * qJ(4) + t280 * t325 + qJDD(4) - t225;
t195 = -t255 * pkin(4) - t273 * pkin(9) + t275 * t257 + t199;
t188 = -0.2e1 * qJD(6) * t248 + (t247 * t291 - t215) * qJ(6) + (t248 * t291 + t214) * pkin(5) + t195;
t178 = m(7) * t188 + t214 * mrSges(7,1) - t215 * mrSges(7,3) - t248 * t238 + t247 * t239;
t311 = m(6) * t195 + t214 * mrSges(6,1) + t215 * mrSges(6,2) + t247 * t236 + t248 * t237 + t178;
t173 = -m(5) * t199 + t255 * mrSges(5,1) - t256 * mrSges(5,2) + t274 * t253 - t275 * t254 - t311;
t288 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t324;
t172 = m(4) * t225 + qJDD(3) * mrSges(4,1) - t282 * mrSges(4,3) + qJD(3) * t288 - t281 * t325 + t173;
t152 = t301 * t157 + t303 * t172;
t220 = Ifges(7,1) * t248 + Ifges(7,4) * t291 + Ifges(7,5) * t247;
t221 = Ifges(6,1) * t248 - Ifges(6,4) * t247 + Ifges(6,5) * t291;
t317 = -mrSges(7,1) * t188 + mrSges(7,2) * t181;
t218 = Ifges(7,4) * t248 + Ifges(7,2) * t291 + Ifges(7,6) * t247;
t327 = -Ifges(6,5) * t248 + Ifges(6,6) * t247 - Ifges(6,3) * t291 - t218;
t164 = -mrSges(6,1) * t195 + mrSges(6,3) * t186 - pkin(5) * t178 + (t220 + t221) * t291 + (Ifges(6,6) - Ifges(7,6)) * t277 + t327 * t248 + (Ifges(6,4) - Ifges(7,5)) * t215 + (-Ifges(6,2) - Ifges(7,3)) * t214 + t317;
t219 = Ifges(6,4) * t248 - Ifges(6,2) * t247 + Ifges(6,6) * t291;
t216 = Ifges(7,5) * t248 + Ifges(7,6) * t291 + Ifges(7,3) * t247;
t314 = mrSges(7,2) * t183 - mrSges(7,3) * t188 + Ifges(7,1) * t215 + Ifges(7,4) * t277 + Ifges(7,5) * t214 + t291 * t216;
t165 = mrSges(6,2) * t195 - mrSges(6,3) * t185 + Ifges(6,1) * t215 - Ifges(6,4) * t214 + Ifges(6,5) * t277 - qJ(6) * t178 - t291 * t219 + t327 * t247 + t314;
t242 = Ifges(5,5) * t275 + Ifges(5,6) * t274 - Ifges(5,3) * t324;
t244 = Ifges(5,1) * t275 + Ifges(5,4) * t274 - Ifges(5,5) * t324;
t146 = -mrSges(5,1) * t199 + mrSges(5,3) * t194 + Ifges(5,4) * t256 + Ifges(5,2) * t255 - Ifges(5,6) * t283 - pkin(4) * t311 + pkin(9) * t319 + t334 * t164 + t300 * t165 - t275 * t242 - t244 * t324;
t243 = Ifges(5,4) * t275 + Ifges(5,2) * t274 - Ifges(5,6) * t324;
t147 = mrSges(5,2) * t199 - mrSges(5,3) * t193 + Ifges(5,1) * t256 + Ifges(5,4) * t255 - Ifges(5,5) * t283 - pkin(9) * t163 - t300 * t164 + t334 * t165 + t274 * t242 + t243 * t324;
t265 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t301 + Ifges(4,2) * t303) * qJD(2);
t266 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t301 + Ifges(4,4) * t303) * qJD(2);
t335 = mrSges(4,1) * t225 - mrSges(4,2) * t226 + Ifges(4,5) * t282 + Ifges(4,6) * t283 + Ifges(4,3) * qJDD(3) + pkin(3) * t173 + qJ(4) * t159 + t297 * t146 + t294 * t147 + (t265 * t301 - t266 * t303) * qJD(2);
t136 = -mrSges(3,1) * t258 + mrSges(3,3) * t241 + t306 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t152 - t335;
t320 = t303 * t157 - t172 * t301;
t150 = m(3) * t241 - mrSges(3,1) * t306 - qJDD(2) * mrSges(3,2) + t320;
t158 = t161 * t297 + t162 * t294;
t310 = -m(4) * t234 + t283 * mrSges(4,1) - mrSges(4,2) * t282 - t287 * t325 + t288 * t324 - t158;
t154 = m(3) * t240 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t306 + t310;
t144 = t304 * t150 - t154 * t302;
t336 = pkin(7) * t144 + t136 * t304;
t330 = t154 * t304;
t151 = m(3) * t258 + t152;
t141 = t150 * t328 - t151 * t296 + t299 * t330;
t264 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t301 + Ifges(4,6) * t303) * qJD(2);
t137 = mrSges(4,2) * t234 - mrSges(4,3) * t225 + Ifges(4,1) * t282 + Ifges(4,4) * t283 + Ifges(4,5) * qJDD(3) - qJ(4) * t158 - qJD(3) * t265 - t146 * t294 + t147 * t297 + t264 * t324;
t312 = mrSges(7,1) * t183 - mrSges(7,3) * t181 - Ifges(7,4) * t215 - Ifges(7,2) * t277 - Ifges(7,6) * t214 + t248 * t216 - t247 * t220;
t309 = mrSges(6,2) * t186 - t247 * t221 - qJ(6) * (-t214 * mrSges(7,2) - t247 * t228 + t322) - pkin(5) * (-t215 * mrSges(7,2) - t248 * t228 + t318) - mrSges(6,1) * t185 - t248 * t219 + Ifges(6,6) * t214 - Ifges(6,5) * t215 - Ifges(6,3) * t277 + t312;
t307 = mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t256 + Ifges(5,6) * t255 + pkin(4) * t163 + t275 * t243 - t274 * t244 - t309;
t145 = Ifges(4,4) * t282 + qJD(3) * t266 + mrSges(4,3) * t226 - mrSges(4,1) * t234 - t307 + Ifges(4,6) * qJDD(3) - pkin(3) * t158 - t264 * t325 + (Ifges(5,3) + Ifges(4,2)) * t283;
t132 = mrSges(3,1) * t240 - mrSges(3,2) * t241 + Ifges(3,3) * qJDD(2) + pkin(2) * t310 + pkin(8) * t320 + t301 * t137 + t303 * t145;
t134 = mrSges(3,2) * t258 - mrSges(3,3) * t240 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t306 - pkin(8) * t152 + t137 * t303 - t145 * t301;
t313 = mrSges(2,1) * t285 - mrSges(2,2) * t286 + pkin(1) * t141 + t299 * t132 + t134 * t329 + t336 * t296;
t142 = m(2) * t286 + t144;
t140 = t299 * t151 + (t150 * t302 + t330) * t296;
t138 = m(2) * t285 + t141;
t130 = mrSges(2,2) * t293 - mrSges(2,3) * t285 + t304 * t134 - t302 * t136 + (-t140 * t296 - t141 * t299) * pkin(7);
t129 = -mrSges(2,1) * t293 + mrSges(2,3) * t286 - pkin(1) * t140 - t296 * t132 + (t134 * t302 + t336) * t299;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t298 * t130 - t295 * t129 - qJ(1) * (t138 * t298 + t142 * t295), t130, t134, t137, t147, t165, -t218 * t247 + t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t295 * t130 + t298 * t129 + qJ(1) * (-t138 * t295 + t142 * t298), t129, t136, t145, t146, t164, -t312; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t313, t313, t132, t335, -Ifges(5,3) * t283 + t307, -t309, Ifges(7,5) * t215 + Ifges(7,6) * t277 + Ifges(7,3) * t214 + t248 * t218 - t291 * t220 - t317;];
m_new  = t1;

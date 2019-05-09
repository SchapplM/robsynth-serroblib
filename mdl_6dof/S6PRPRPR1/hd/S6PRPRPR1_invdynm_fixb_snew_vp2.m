% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:07:48
% EndTime: 2019-05-04 22:08:09
% DurationCPUTime: 15.52s
% Computational Cost: add. (295574->298), mult. (589482->389), div. (0->0), fcn. (416602->14), ass. (0->133)
t275 = sin(pkin(12));
t279 = cos(pkin(12));
t284 = sin(qJ(4));
t287 = cos(qJ(4));
t248 = (t275 * t284 - t279 * t287) * qJD(2);
t277 = sin(pkin(10));
t281 = cos(pkin(10));
t263 = g(1) * t277 - g(2) * t281;
t264 = -g(1) * t281 - g(2) * t277;
t274 = -g(3) + qJDD(1);
t285 = sin(qJ(2));
t282 = cos(pkin(6));
t288 = cos(qJ(2));
t312 = t282 * t288;
t278 = sin(pkin(6));
t314 = t278 * t288;
t226 = t263 * t312 - t264 * t285 + t274 * t314;
t221 = qJDD(2) * pkin(2) + t226;
t313 = t282 * t285;
t315 = t278 * t285;
t227 = t263 * t313 + t288 * t264 + t274 * t315;
t290 = qJD(2) ^ 2;
t222 = -pkin(2) * t290 + t227;
t276 = sin(pkin(11));
t280 = cos(pkin(11));
t206 = t276 * t221 + t280 * t222;
t203 = -pkin(3) * t290 + qJDD(2) * pkin(8) + t206;
t245 = -t263 * t278 + t282 * t274;
t242 = qJDD(3) + t245;
t199 = -t284 * t203 + t287 * t242;
t309 = qJD(2) * qJD(4);
t307 = t287 * t309;
t260 = qJDD(2) * t284 + t307;
t196 = (-t260 + t307) * qJ(5) + (t284 * t287 * t290 + qJDD(4)) * pkin(4) + t199;
t200 = t287 * t203 + t284 * t242;
t261 = qJDD(2) * t287 - t284 * t309;
t311 = qJD(2) * t284;
t265 = qJD(4) * pkin(4) - qJ(5) * t311;
t273 = t287 ^ 2;
t197 = -pkin(4) * t273 * t290 + qJ(5) * t261 - qJD(4) * t265 + t200;
t317 = 2 * qJD(5);
t192 = t275 * t196 + t279 * t197 - t248 * t317;
t249 = (t275 * t287 + t279 * t284) * qJD(2);
t229 = mrSges(6,1) * t248 + mrSges(6,2) * t249;
t235 = -t260 * t275 + t261 * t279;
t244 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t249;
t230 = pkin(5) * t248 - pkin(9) * t249;
t289 = qJD(4) ^ 2;
t189 = -pkin(5) * t289 + qJDD(4) * pkin(9) - t230 * t248 + t192;
t205 = t280 * t221 - t276 * t222;
t299 = -qJDD(2) * pkin(3) - t205;
t198 = -t261 * pkin(4) + qJDD(5) + t265 * t311 + (-qJ(5) * t273 - pkin(8)) * t290 + t299;
t236 = t260 * t279 + t261 * t275;
t193 = (qJD(4) * t248 - t236) * pkin(9) + (qJD(4) * t249 - t235) * pkin(5) + t198;
t283 = sin(qJ(6));
t286 = cos(qJ(6));
t186 = -t189 * t283 + t193 * t286;
t239 = qJD(4) * t286 - t249 * t283;
t213 = qJD(6) * t239 + qJDD(4) * t283 + t236 * t286;
t240 = qJD(4) * t283 + t249 * t286;
t215 = -mrSges(7,1) * t239 + mrSges(7,2) * t240;
t247 = qJD(6) + t248;
t216 = -mrSges(7,2) * t247 + mrSges(7,3) * t239;
t234 = qJDD(6) - t235;
t182 = m(7) * t186 + mrSges(7,1) * t234 - mrSges(7,3) * t213 - t215 * t240 + t216 * t247;
t187 = t189 * t286 + t193 * t283;
t212 = -qJD(6) * t240 + qJDD(4) * t286 - t236 * t283;
t217 = mrSges(7,1) * t247 - mrSges(7,3) * t240;
t183 = m(7) * t187 - mrSges(7,2) * t234 + mrSges(7,3) * t212 + t215 * t239 - t217 * t247;
t303 = -t182 * t283 + t286 * t183;
t168 = m(6) * t192 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t235 - qJD(4) * t244 - t229 * t248 + t303;
t302 = -t279 * t196 + t275 * t197;
t191 = -0.2e1 * qJD(5) * t249 - t302;
t243 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t248;
t188 = -qJDD(4) * pkin(5) - t289 * pkin(9) + (t317 + t230) * t249 + t302;
t296 = -m(7) * t188 + t212 * mrSges(7,1) - mrSges(7,2) * t213 + t239 * t216 - t217 * t240;
t178 = m(6) * t191 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t236 + qJD(4) * t243 - t229 * t249 + t296;
t163 = t275 * t168 + t279 * t178;
t253 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t284 + Ifges(5,2) * t287) * qJD(2);
t254 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t284 + Ifges(5,4) * t287) * qJD(2);
t207 = Ifges(7,5) * t240 + Ifges(7,6) * t239 + Ifges(7,3) * t247;
t209 = Ifges(7,1) * t240 + Ifges(7,4) * t239 + Ifges(7,5) * t247;
t175 = -mrSges(7,1) * t188 + mrSges(7,3) * t187 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t234 - t207 * t240 + t209 * t247;
t208 = Ifges(7,4) * t240 + Ifges(7,2) * t239 + Ifges(7,6) * t247;
t176 = mrSges(7,2) * t188 - mrSges(7,3) * t186 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t234 + t207 * t239 - t208 * t247;
t224 = Ifges(6,4) * t249 - Ifges(6,2) * t248 + Ifges(6,6) * qJD(4);
t225 = Ifges(6,1) * t249 - Ifges(6,4) * t248 + Ifges(6,5) * qJD(4);
t294 = -mrSges(6,1) * t191 + mrSges(6,2) * t192 - Ifges(6,5) * t236 - Ifges(6,6) * t235 - Ifges(6,3) * qJDD(4) - pkin(5) * t296 - pkin(9) * t303 - t286 * t175 - t283 * t176 - t249 * t224 - t248 * t225;
t318 = mrSges(5,1) * t199 - mrSges(5,2) * t200 + Ifges(5,5) * t260 + Ifges(5,6) * t261 + Ifges(5,3) * qJDD(4) + pkin(4) * t163 + (t253 * t284 - t254 * t287) * qJD(2) - t294;
t259 = (-mrSges(5,1) * t287 + mrSges(5,2) * t284) * qJD(2);
t310 = qJD(2) * t287;
t267 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t310;
t161 = m(5) * t199 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t260 + qJD(4) * t267 - t259 * t311 + t163;
t266 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t311;
t304 = t279 * t168 - t178 * t275;
t162 = m(5) * t200 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t261 - qJD(4) * t266 + t259 * t310 + t304;
t305 = -t161 * t284 + t287 * t162;
t152 = m(4) * t206 - mrSges(4,1) * t290 - qJDD(2) * mrSges(4,2) + t305;
t202 = -t290 * pkin(8) + t299;
t171 = t286 * t182 + t283 * t183;
t295 = m(6) * t198 - t235 * mrSges(6,1) + mrSges(6,2) * t236 + t248 * t243 + t244 * t249 + t171;
t292 = -m(5) * t202 + t261 * mrSges(5,1) - mrSges(5,2) * t260 - t266 * t311 + t267 * t310 - t295;
t165 = m(4) * t205 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t290 + t292;
t149 = t276 * t152 + t280 * t165;
t147 = m(3) * t226 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t290 + t149;
t306 = t280 * t152 - t165 * t276;
t148 = m(3) * t227 - mrSges(3,1) * t290 - qJDD(2) * mrSges(3,2) + t306;
t138 = -t147 * t285 + t288 * t148;
t316 = pkin(7) * t138;
t155 = t287 * t161 + t284 * t162;
t308 = m(4) * t242 + t155;
t153 = m(3) * t245 + t308;
t135 = t147 * t312 + t148 * t313 - t153 * t278;
t223 = Ifges(6,5) * t249 - Ifges(6,6) * t248 + Ifges(6,3) * qJD(4);
t156 = mrSges(6,2) * t198 - mrSges(6,3) * t191 + Ifges(6,1) * t236 + Ifges(6,4) * t235 + Ifges(6,5) * qJDD(4) - pkin(9) * t171 - qJD(4) * t224 - t175 * t283 + t176 * t286 - t223 * t248;
t293 = mrSges(7,1) * t186 - mrSges(7,2) * t187 + Ifges(7,5) * t213 + Ifges(7,6) * t212 + Ifges(7,3) * t234 + t208 * t240 - t209 * t239;
t157 = -mrSges(6,1) * t198 + mrSges(6,3) * t192 + Ifges(6,4) * t236 + Ifges(6,2) * t235 + Ifges(6,6) * qJDD(4) - pkin(5) * t171 + qJD(4) * t225 - t223 * t249 - t293;
t252 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t284 + Ifges(5,6) * t287) * qJD(2);
t141 = -mrSges(5,1) * t202 + mrSges(5,3) * t200 + Ifges(5,4) * t260 + Ifges(5,2) * t261 + Ifges(5,6) * qJDD(4) - pkin(4) * t295 + qJ(5) * t304 + qJD(4) * t254 + t275 * t156 + t279 * t157 - t252 * t311;
t143 = mrSges(5,2) * t202 - mrSges(5,3) * t199 + Ifges(5,1) * t260 + Ifges(5,4) * t261 + Ifges(5,5) * qJDD(4) - qJ(5) * t163 - qJD(4) * t253 + t156 * t279 - t157 * t275 + t252 * t310;
t131 = mrSges(4,2) * t242 - mrSges(4,3) * t205 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t290 - pkin(8) * t155 - t141 * t284 + t143 * t287;
t139 = -mrSges(4,1) * t242 + mrSges(4,3) * t206 + t290 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t155 - t318;
t126 = -mrSges(3,1) * t245 + mrSges(3,3) * t227 + t290 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t308 + qJ(3) * t306 + t276 * t131 + t280 * t139;
t128 = mrSges(3,2) * t245 - mrSges(3,3) * t226 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t290 - qJ(3) * t149 + t131 * t280 - t139 * t276;
t297 = mrSges(4,1) * t205 - mrSges(4,2) * t206 + Ifges(4,3) * qJDD(2) + pkin(3) * t292 + pkin(8) * t305 + t287 * t141 + t284 * t143;
t130 = mrSges(3,1) * t226 - mrSges(3,2) * t227 + Ifges(3,3) * qJDD(2) + pkin(2) * t149 + t297;
t298 = mrSges(2,1) * t263 - mrSges(2,2) * t264 + pkin(1) * t135 + t126 * t314 + t128 * t315 + t282 * t130 + t278 * t316;
t136 = m(2) * t264 + t138;
t134 = t282 * t153 + (t147 * t288 + t148 * t285) * t278;
t132 = m(2) * t263 + t135;
t124 = mrSges(2,2) * t274 - mrSges(2,3) * t263 - t285 * t126 + t288 * t128 + (-t134 * t278 - t135 * t282) * pkin(7);
t123 = -mrSges(2,1) * t274 + mrSges(2,3) * t264 - pkin(1) * t134 - t278 * t130 + (t126 * t288 + t128 * t285 + t316) * t282;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t281 * t124 - t277 * t123 - qJ(1) * (t132 * t281 + t136 * t277), t124, t128, t131, t143, t156, t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t277 * t124 + t281 * t123 + qJ(1) * (-t132 * t277 + t136 * t281), t123, t126, t139, t141, t157, t175; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t298, t298, t130, t297, t318, -t294, t293;];
m_new  = t1;

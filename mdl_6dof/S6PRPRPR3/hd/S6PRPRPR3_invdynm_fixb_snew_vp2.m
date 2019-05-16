% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 22:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:33:07
% EndTime: 2019-05-04 22:33:16
% DurationCPUTime: 6.37s
% Computational Cost: add. (114462->305), mult. (213629->374), div. (0->0), fcn. (134281->12), ass. (0->133)
t276 = sin(pkin(10));
t279 = cos(pkin(10));
t256 = g(1) * t276 - g(2) * t279;
t257 = -g(1) * t279 - g(2) * t276;
t274 = -g(3) + qJDD(1);
t283 = sin(qJ(2));
t280 = cos(pkin(6));
t286 = cos(qJ(2));
t316 = t280 * t286;
t277 = sin(pkin(6));
t318 = t277 * t286;
t208 = t256 * t316 - t257 * t283 + t274 * t318;
t206 = qJDD(2) * pkin(2) + t208;
t317 = t280 * t283;
t319 = t277 * t283;
t209 = t256 * t317 + t286 * t257 + t274 * t319;
t288 = qJD(2) ^ 2;
t207 = -pkin(2) * t288 + t209;
t275 = sin(pkin(11));
t278 = cos(pkin(11));
t201 = t275 * t206 + t278 * t207;
t198 = -pkin(3) * t288 + qJDD(2) * pkin(8) + t201;
t282 = sin(qJ(4));
t195 = t282 * t198;
t225 = -t256 * t277 + t280 * t274;
t224 = qJDD(3) + t225;
t285 = cos(qJ(4));
t315 = t285 * t224;
t193 = -t195 + t315;
t194 = t285 * t198 + t282 * t224;
t233 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t282 + Ifges(5,4) * t285) * qJD(2);
t250 = (mrSges(6,2) * t285 - mrSges(6,3) * t282) * qJD(2);
t309 = qJD(2) * qJD(4);
t307 = t285 * t309;
t252 = qJDD(2) * t282 + t307;
t306 = t282 * t309;
t253 = qJDD(2) * t285 - t306;
t310 = qJD(2) * t285;
t260 = -mrSges(6,1) * t310 - qJD(4) * mrSges(6,3);
t249 = (-pkin(4) * t285 - qJ(5) * t282) * qJD(2);
t287 = qJD(4) ^ 2;
t311 = qJD(2) * t282;
t303 = -t287 * qJ(5) + t249 * t311 + qJDD(5) + t195;
t322 = pkin(9) * t288;
t324 = -pkin(4) - pkin(9);
t187 = t252 * pkin(5) + t324 * qJDD(4) + (-pkin(5) * t309 - t282 * t322 - t224) * t285 + t303;
t264 = pkin(5) * t311 - qJD(4) * pkin(9);
t273 = t285 ^ 2;
t200 = t278 * t206 - t275 * t207;
t302 = -qJDD(2) * pkin(3) - t200;
t325 = -2 * qJD(5);
t292 = pkin(4) * t306 + t311 * t325 + (-t252 - t307) * qJ(5) + t302;
t188 = -t264 * t311 + (-pkin(5) * t273 - pkin(8)) * t288 + t324 * t253 + t292;
t281 = sin(qJ(6));
t284 = cos(qJ(6));
t182 = t187 * t284 - t188 * t281;
t247 = -qJD(4) * t281 - t284 * t310;
t218 = qJD(6) * t247 + qJDD(4) * t284 - t253 * t281;
t248 = qJD(4) * t284 - t281 * t310;
t219 = -mrSges(7,1) * t247 + mrSges(7,2) * t248;
t266 = qJD(6) + t311;
t222 = -mrSges(7,2) * t266 + mrSges(7,3) * t247;
t245 = qJDD(6) + t252;
t179 = m(7) * t182 + mrSges(7,1) * t245 - mrSges(7,3) * t218 - t219 * t248 + t222 * t266;
t183 = t187 * t281 + t188 * t284;
t217 = -qJD(6) * t248 - qJDD(4) * t281 - t253 * t284;
t223 = mrSges(7,1) * t266 - mrSges(7,3) * t248;
t180 = m(7) * t183 - mrSges(7,2) * t245 + mrSges(7,3) * t217 + t219 * t247 - t223 * t266;
t168 = t284 * t179 + t281 * t180;
t296 = -t287 * pkin(4) + qJDD(4) * qJ(5) + t249 * t310 + t194;
t186 = -t273 * t322 + t253 * pkin(5) + ((2 * qJD(5)) + t264) * qJD(4) + t296;
t210 = Ifges(7,5) * t248 + Ifges(7,6) * t247 + Ifges(7,3) * t266;
t212 = Ifges(7,1) * t248 + Ifges(7,4) * t247 + Ifges(7,5) * t266;
t171 = -mrSges(7,1) * t186 + mrSges(7,3) * t183 + Ifges(7,4) * t218 + Ifges(7,2) * t217 + Ifges(7,6) * t245 - t210 * t248 + t212 * t266;
t211 = Ifges(7,4) * t248 + Ifges(7,2) * t247 + Ifges(7,6) * t266;
t172 = mrSges(7,2) * t186 - mrSges(7,3) * t182 + Ifges(7,1) * t218 + Ifges(7,4) * t217 + Ifges(7,5) * t245 + t210 * t247 - t211 * t266;
t189 = qJD(4) * t325 - t296;
t191 = -qJDD(4) * pkin(4) + t303 - t315;
t235 = Ifges(6,4) * qJD(4) + (-Ifges(6,2) * t282 - Ifges(6,6) * t285) * qJD(2);
t294 = -mrSges(6,2) * t191 + mrSges(6,3) * t189 - Ifges(6,1) * qJDD(4) + Ifges(6,4) * t252 + Ifges(6,5) * t253 + pkin(9) * t168 + t281 * t171 - t284 * t172 - t235 * t310;
t184 = -m(7) * t186 + t217 * mrSges(7,1) - t218 * mrSges(7,2) + t247 * t222 - t248 * t223;
t261 = mrSges(6,1) * t311 + qJD(4) * mrSges(6,2);
t295 = -m(6) * t189 + qJDD(4) * mrSges(6,3) + qJD(4) * t261 + t250 * t310 - t184;
t300 = -m(6) * t191 - t252 * mrSges(6,1) - t168;
t234 = Ifges(6,5) * qJD(4) + (-Ifges(6,6) * t282 - Ifges(6,3) * t285) * qJD(2);
t312 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t282 + Ifges(5,2) * t285) * qJD(2) - t234;
t327 = (-t285 * t233 + t312 * t282) * qJD(2) + mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t252 + Ifges(5,6) * t253 + Ifges(5,3) * qJDD(4) + pkin(4) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t260 - t250 * t311 + t300) + qJ(5) * (mrSges(6,1) * t253 + t295) - t294;
t251 = (-mrSges(5,1) * t285 + mrSges(5,2) * t282) * qJD(2);
t259 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t310;
t164 = m(5) * t193 - t252 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t259 - t260) * qJD(4) + (-t250 - t251) * t311 + t300;
t258 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t311;
t175 = t251 * t310 + m(5) * t194 - qJDD(4) * mrSges(5,2) - qJD(4) * t258 + (mrSges(5,3) + mrSges(6,1)) * t253 + t295;
t304 = -t164 * t282 + t285 * t175;
t157 = m(4) * t201 - mrSges(4,1) * t288 - qJDD(2) * mrSges(4,2) + t304;
t321 = t288 * pkin(8);
t197 = t302 - t321;
t169 = -t281 * t179 + t284 * t180;
t192 = -t253 * pkin(4) + t292 - t321;
t301 = -m(6) * t192 - t253 * mrSges(6,2) + t261 * t311 - t169;
t290 = -m(5) * t197 + t259 * t310 + t253 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t252 + (-t258 * t282 - t260 * t285) * qJD(2) + t301;
t162 = m(4) * t200 + qJDD(2) * mrSges(4,1) - t288 * mrSges(4,2) + t290;
t154 = t275 * t157 + t278 * t162;
t152 = m(3) * t208 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t288 + t154;
t305 = t278 * t157 - t162 * t275;
t153 = m(3) * t209 - mrSges(3,1) * t288 - qJDD(2) * mrSges(3,2) + t305;
t143 = -t152 * t283 + t286 * t153;
t323 = pkin(7) * t143;
t320 = Ifges(5,4) + Ifges(6,6);
t160 = t285 * t164 + t282 * t175;
t236 = Ifges(6,1) * qJD(4) + (-Ifges(6,4) * t282 - Ifges(6,5) * t285) * qJD(2);
t313 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t282 + Ifges(5,6) * t285) * qJD(2) + t236;
t308 = m(4) * t224 + t160;
t158 = m(3) * t225 + t308;
t140 = t152 * t316 + t153 * t317 - t158 * t277;
t166 = -t252 * mrSges(6,3) + t260 * t310 - t301;
t293 = -mrSges(6,1) * t189 + mrSges(6,2) * t192 - pkin(5) * t184 - pkin(9) * t169 - t284 * t171 - t281 * t172;
t146 = -mrSges(5,1) * t197 + mrSges(5,3) * t194 - pkin(4) * t166 + (Ifges(5,2) + Ifges(6,3)) * t253 + t320 * t252 + (Ifges(5,6) - Ifges(6,5)) * qJDD(4) + (t233 - t235) * qJD(4) - t313 * t311 + t293;
t297 = mrSges(7,1) * t182 - mrSges(7,2) * t183 + Ifges(7,5) * t218 + Ifges(7,6) * t217 + Ifges(7,3) * t245 + t248 * t211 - t247 * t212;
t291 = mrSges(6,1) * t191 - mrSges(6,3) * t192 + pkin(5) * t168 + t297;
t148 = t320 * t253 + (Ifges(5,1) + Ifges(6,2)) * t252 + (Ifges(5,5) - Ifges(6,4)) * qJDD(4) - t312 * qJD(4) + t313 * t310 - mrSges(5,3) * t193 + mrSges(5,2) * t197 - qJ(5) * t166 + t291;
t136 = mrSges(4,2) * t224 - mrSges(4,3) * t200 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t288 - pkin(8) * t160 - t146 * t282 + t148 * t285;
t144 = -mrSges(4,1) * t224 + mrSges(4,3) * t201 + t288 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t160 - t327;
t131 = -mrSges(3,1) * t225 + mrSges(3,3) * t209 + t288 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t308 + qJ(3) * t305 + t275 * t136 + t278 * t144;
t133 = mrSges(3,2) * t225 - mrSges(3,3) * t208 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t288 - qJ(3) * t154 + t136 * t278 - t144 * t275;
t298 = mrSges(4,1) * t200 - mrSges(4,2) * t201 + Ifges(4,3) * qJDD(2) + pkin(3) * t290 + pkin(8) * t304 + t285 * t146 + t282 * t148;
t135 = mrSges(3,1) * t208 - mrSges(3,2) * t209 + Ifges(3,3) * qJDD(2) + pkin(2) * t154 + t298;
t299 = mrSges(2,1) * t256 - mrSges(2,2) * t257 + pkin(1) * t140 + t131 * t318 + t133 * t319 + t280 * t135 + t277 * t323;
t141 = m(2) * t257 + t143;
t139 = t280 * t158 + (t152 * t286 + t153 * t283) * t277;
t137 = m(2) * t256 + t140;
t129 = mrSges(2,2) * t274 - mrSges(2,3) * t256 - t283 * t131 + t286 * t133 + (-t139 * t277 - t140 * t280) * pkin(7);
t128 = -mrSges(2,1) * t274 + mrSges(2,3) * t257 - pkin(1) * t139 - t277 * t135 + (t131 * t286 + t133 * t283 + t323) * t280;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t279 * t129 - t276 * t128 - qJ(1) * (t137 * t279 + t141 * t276), t129, t133, t136, t148, -t234 * t311 - t294, t172; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t276 * t129 + t279 * t128 + qJ(1) * (-t137 * t276 + t141 * t279), t128, t131, t144, t146, Ifges(6,4) * qJDD(4) - Ifges(6,2) * t252 - Ifges(6,6) * t253 - qJD(4) * t234 - t236 * t310 - t291, t171; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t299, t299, t135, t298, t327, Ifges(6,5) * qJDD(4) - Ifges(6,6) * t252 - Ifges(6,3) * t253 + qJD(4) * t235 + t236 * t311 - t293, t297;];
m_new  = t1;

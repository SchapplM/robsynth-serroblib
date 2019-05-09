% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 20:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:01:32
% EndTime: 2019-05-05 20:01:59
% DurationCPUTime: 14.70s
% Computational Cost: add. (272221->341), mult. (575440->417), div. (0->0), fcn. (384750->10), ass. (0->137)
t301 = sin(qJ(1));
t305 = cos(qJ(1));
t283 = -t305 * g(1) - t301 * g(2);
t335 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t283;
t334 = (-pkin(1) - pkin(7));
t333 = mrSges(2,1) - mrSges(3,2);
t332 = -Ifges(3,4) + Ifges(2,5);
t331 = (Ifges(3,5) - Ifges(2,6));
t307 = qJD(1) ^ 2;
t245 = (t307 * t334) - t335;
t300 = sin(qJ(3));
t304 = cos(qJ(3));
t329 = qJD(1) * qJD(3);
t326 = t304 * t329;
t277 = qJDD(1) * t300 + t326;
t327 = t300 * t329;
t278 = qJDD(1) * t304 - t327;
t223 = (-t278 + t327) * qJ(4) + (t277 + t326) * pkin(3) + t245;
t282 = t301 * g(1) - g(2) * t305;
t320 = -t307 * qJ(2) + qJDD(2) - t282;
t251 = qJDD(1) * t334 + t320;
t241 = -g(3) * t304 + t251 * t300;
t275 = (pkin(3) * t300 - qJ(4) * t304) * qJD(1);
t286 = t300 * qJD(1);
t306 = qJD(3) ^ 2;
t226 = -pkin(3) * t306 + qJDD(3) * qJ(4) - t275 * t286 + t241;
t296 = sin(pkin(10));
t297 = cos(pkin(10));
t330 = qJD(1) * t304;
t270 = qJD(3) * t296 + t297 * t330;
t203 = -0.2e1 * qJD(4) * t270 + t223 * t297 - t296 * t226;
t249 = qJDD(3) * t296 + t278 * t297;
t269 = qJD(3) * t297 - t296 * t330;
t194 = (t269 * t286 - t249) * pkin(8) + (t269 * t270 + t277) * pkin(4) + t203;
t204 = 0.2e1 * qJD(4) * t269 + t223 * t296 + t226 * t297;
t248 = qJDD(3) * t297 - t278 * t296;
t250 = pkin(4) * t286 - pkin(8) * t270;
t268 = t269 ^ 2;
t196 = -pkin(4) * t268 + pkin(8) * t248 - t250 * t286 + t204;
t299 = sin(qJ(5));
t303 = cos(qJ(5));
t181 = t194 * t303 - t299 * t196;
t236 = t269 * t303 - t270 * t299;
t213 = qJD(5) * t236 + t248 * t299 + t249 * t303;
t237 = t269 * t299 + t270 * t303;
t274 = qJDD(5) + t277;
t285 = t286 + qJD(5);
t178 = (t236 * t285 - t213) * pkin(9) + (t236 * t237 + t274) * pkin(5) + t181;
t182 = t194 * t299 + t196 * t303;
t212 = -qJD(5) * t237 + t248 * t303 - t249 * t299;
t229 = pkin(5) * t285 - pkin(9) * t237;
t235 = t236 ^ 2;
t179 = -pkin(5) * t235 + pkin(9) * t212 - t229 * t285 + t182;
t298 = sin(qJ(6));
t302 = cos(qJ(6));
t176 = t178 * t302 - t179 * t298;
t218 = t236 * t302 - t237 * t298;
t190 = qJD(6) * t218 + t212 * t298 + t213 * t302;
t219 = t236 * t298 + t237 * t302;
t202 = -mrSges(7,1) * t218 + mrSges(7,2) * t219;
t284 = qJD(6) + t285;
t207 = -mrSges(7,2) * t284 + mrSges(7,3) * t218;
t265 = qJDD(6) + t274;
t171 = m(7) * t176 + mrSges(7,1) * t265 - mrSges(7,3) * t190 - t202 * t219 + t207 * t284;
t177 = t178 * t298 + t179 * t302;
t189 = -qJD(6) * t219 + t212 * t302 - t213 * t298;
t208 = mrSges(7,1) * t284 - mrSges(7,3) * t219;
t172 = m(7) * t177 - mrSges(7,2) * t265 + mrSges(7,3) * t189 + t202 * t218 - t208 * t284;
t164 = t171 * t302 + t172 * t298;
t220 = -mrSges(6,1) * t236 + mrSges(6,2) * t237;
t227 = -mrSges(6,2) * t285 + mrSges(6,3) * t236;
t161 = m(6) * t181 + mrSges(6,1) * t274 - mrSges(6,3) * t213 - t220 * t237 + t227 * t285 + t164;
t228 = mrSges(6,1) * t285 - mrSges(6,3) * t237;
t323 = -t171 * t298 + t172 * t302;
t162 = m(6) * t182 - mrSges(6,2) * t274 + mrSges(6,3) * t212 + t220 * t236 - t228 * t285 + t323;
t157 = t161 * t303 + t162 * t299;
t238 = -mrSges(5,1) * t269 + mrSges(5,2) * t270;
t246 = -mrSges(5,2) * t286 + mrSges(5,3) * t269;
t155 = m(5) * t203 + mrSges(5,1) * t277 - mrSges(5,3) * t249 - t238 * t270 + t246 * t286 + t157;
t247 = mrSges(5,1) * t286 - mrSges(5,3) * t270;
t324 = -t161 * t299 + t162 * t303;
t156 = m(5) * t204 - mrSges(5,2) * t277 + mrSges(5,3) * t248 + t238 * t269 - t247 * t286 + t324;
t148 = t155 * t297 + t156 * t296;
t276 = (mrSges(4,1) * t300 + mrSges(4,2) * t304) * qJD(1);
t281 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t325 = -t155 * t296 + t156 * t297;
t146 = m(4) * t241 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t277 - qJD(3) * t281 - t276 * t286 + t325;
t240 = g(3) * t300 + t251 * t304;
t280 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t286;
t225 = -qJDD(3) * pkin(3) - qJ(4) * t306 + t275 * t330 + qJDD(4) - t240;
t205 = -pkin(4) * t248 - pkin(8) * t268 + t250 * t270 + t225;
t184 = -pkin(5) * t212 - pkin(9) * t235 + t229 * t237 + t205;
t322 = m(7) * t184 - mrSges(7,1) * t189 + mrSges(7,2) * t190 - t207 * t218 + t208 * t219;
t313 = m(6) * t205 - mrSges(6,1) * t212 + mrSges(6,2) * t213 - t227 * t236 + t228 * t237 + t322;
t309 = -m(5) * t225 + mrSges(5,1) * t248 - mrSges(5,2) * t249 + t246 * t269 - t247 * t270 - t313;
t167 = m(4) * t240 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t278 + qJD(3) * t280 - t276 * t330 + t309;
t141 = t146 * t304 - t167 * t300;
t140 = t300 * t146 + t304 * t167;
t256 = -qJDD(1) * pkin(1) + t320;
t319 = -m(3) * t256 + (mrSges(3,3) * t307) - t140;
t144 = -m(4) * t245 - mrSges(4,1) * t277 - mrSges(4,2) * t278 - t280 * t286 - t281 * t330 - t148;
t198 = Ifges(7,4) * t219 + Ifges(7,2) * t218 + Ifges(7,6) * t284;
t199 = Ifges(7,1) * t219 + Ifges(7,4) * t218 + Ifges(7,5) * t284;
t318 = -mrSges(7,1) * t176 + mrSges(7,2) * t177 - Ifges(7,5) * t190 - Ifges(7,6) * t189 - Ifges(7,3) * t265 - t198 * t219 + t218 * t199;
t197 = Ifges(7,5) * t219 + Ifges(7,6) * t218 + Ifges(7,3) * t284;
t165 = -mrSges(7,1) * t184 + mrSges(7,3) * t177 + Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t265 - t197 * t219 + t199 * t284;
t166 = mrSges(7,2) * t184 - mrSges(7,3) * t176 + Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t265 + t197 * t218 - t198 * t284;
t214 = Ifges(6,5) * t237 + Ifges(6,6) * t236 + Ifges(6,3) * t285;
t216 = Ifges(6,1) * t237 + Ifges(6,4) * t236 + Ifges(6,5) * t285;
t150 = -mrSges(6,1) * t205 + mrSges(6,3) * t182 + Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t274 - pkin(5) * t322 + pkin(9) * t323 + t302 * t165 + t298 * t166 - t237 * t214 + t285 * t216;
t215 = Ifges(6,4) * t237 + Ifges(6,2) * t236 + Ifges(6,6) * t285;
t151 = mrSges(6,2) * t205 - mrSges(6,3) * t181 + Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t274 - pkin(9) * t164 - t165 * t298 + t166 * t302 + t214 * t236 - t215 * t285;
t230 = Ifges(5,5) * t270 + Ifges(5,6) * t269 + Ifges(5,3) * t286;
t232 = Ifges(5,1) * t270 + Ifges(5,4) * t269 + Ifges(5,5) * t286;
t134 = -mrSges(5,1) * t225 + mrSges(5,3) * t204 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t277 - pkin(4) * t313 + pkin(8) * t324 + t303 * t150 + t299 * t151 - t270 * t230 + t232 * t286;
t231 = Ifges(5,4) * t270 + Ifges(5,2) * t269 + Ifges(5,6) * t286;
t136 = mrSges(5,2) * t225 - mrSges(5,3) * t203 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t277 - pkin(8) * t157 - t150 * t299 + t151 * t303 + t230 * t269 - t231 * t286;
t261 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t304 - Ifges(4,6) * t300) * qJD(1);
t262 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t304 - Ifges(4,2) * t300) * qJD(1);
t131 = mrSges(4,2) * t245 - mrSges(4,3) * t240 + Ifges(4,1) * t278 - Ifges(4,4) * t277 + Ifges(4,5) * qJDD(3) - qJ(4) * t148 - qJD(3) * t262 - t134 * t296 + t136 * t297 - t261 * t286;
t263 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t304 - Ifges(4,4) * t300) * qJD(1);
t310 = -mrSges(6,1) * t181 + mrSges(6,2) * t182 - Ifges(6,5) * t213 - Ifges(6,6) * t212 - Ifges(6,3) * t274 - pkin(5) * t164 - t215 * t237 + t236 * t216 + t318;
t308 = mrSges(5,1) * t203 - mrSges(5,2) * t204 + Ifges(5,5) * t249 + Ifges(5,6) * t248 + pkin(4) * t157 + t270 * t231 - t269 * t232 - t310;
t132 = Ifges(4,6) * qJDD(3) + (-Ifges(5,3) - Ifges(4,2)) * t277 + Ifges(4,4) * t278 + qJD(3) * t263 + mrSges(4,3) * t241 - mrSges(4,1) * t245 - t308 - pkin(3) * t148 - t261 * t330;
t254 = t307 * pkin(1) + t335;
t317 = mrSges(3,2) * t256 - mrSges(3,3) * t254 + Ifges(3,1) * qJDD(1) - pkin(7) * t140 + t131 * t304 - t132 * t300;
t316 = -mrSges(3,1) * t254 - pkin(2) * t144 - pkin(7) * t141 - t300 * t131 - t304 * t132;
t315 = mrSges(4,1) * t240 - mrSges(4,2) * t241 + Ifges(4,5) * t278 - Ifges(4,6) * t277 + Ifges(4,3) * qJDD(3) + pkin(3) * t309 + qJ(4) * t325 + t134 * t297 + t136 * t296 + t262 * t330 + t263 * t286;
t314 = -m(3) * t254 + mrSges(3,2) * t307 + qJDD(1) * mrSges(3,3) - t144;
t312 = -mrSges(2,2) * t283 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t319) + qJ(2) * t314 + mrSges(2,1) * t282 + Ifges(2,3) * qJDD(1) + t317;
t311 = mrSges(3,1) * t256 + pkin(2) * t140 + t315;
t142 = m(2) * t283 - mrSges(2,1) * t307 - qJDD(1) * mrSges(2,2) + t314;
t139 = -m(3) * g(3) + t141;
t137 = m(2) * t282 - t307 * mrSges(2,2) + qJDD(1) * t333 + t319;
t129 = t311 - mrSges(2,3) * t282 - qJ(2) * t139 + (t331 * t307) + t332 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t128 = mrSges(2,3) * t283 - pkin(1) * t139 + g(3) * t333 - qJDD(1) * t331 + t307 * t332 + t316;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t305 * t129 - t301 * t128 - pkin(6) * (t137 * t305 + t142 * t301), t129, t317, t131, t136, t151, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t129 + t305 * t128 + pkin(6) * (-t137 * t301 + t142 * t305), t128, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t307 * Ifges(3,5)) - t311, t132, t134, t150, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t312, t312, mrSges(3,2) * g(3) + t307 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t316, t315, Ifges(5,3) * t277 + t308, -t310, -t318;];
m_new  = t1;

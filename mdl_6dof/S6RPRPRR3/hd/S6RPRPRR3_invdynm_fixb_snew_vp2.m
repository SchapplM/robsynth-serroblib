% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:33:06
% EndTime: 2019-05-05 18:33:39
% DurationCPUTime: 22.96s
% Computational Cost: add. (431584->341), mult. (902463->428), div. (0->0), fcn. (615358->12), ass. (0->137)
t298 = sin(qJ(1));
t302 = cos(qJ(1));
t281 = t298 * g(1) - g(2) * t302;
t272 = qJDD(1) * pkin(1) + t281;
t282 = -g(1) * t302 - g(2) * t298;
t304 = qJD(1) ^ 2;
t275 = -pkin(1) * t304 + t282;
t292 = sin(pkin(10));
t294 = cos(pkin(10));
t245 = t294 * t272 - t292 * t275;
t234 = -qJDD(1) * pkin(2) - t304 * pkin(7) - t245;
t297 = sin(qJ(3));
t301 = cos(qJ(3));
t320 = qJD(1) * qJD(3);
t319 = t301 * t320;
t276 = qJDD(1) * t297 + t319;
t286 = t297 * t320;
t277 = qJDD(1) * t301 - t286;
t220 = (-t276 - t319) * qJ(4) + (-t277 + t286) * pkin(3) + t234;
t246 = t292 * t272 + t294 * t275;
t235 = -pkin(2) * t304 + qJDD(1) * pkin(7) + t246;
t290 = -g(3) + qJDD(2);
t228 = t301 * t235 + t297 * t290;
t273 = (-pkin(3) * t301 - qJ(4) * t297) * qJD(1);
t303 = qJD(3) ^ 2;
t321 = qJD(1) * t301;
t226 = -pkin(3) * t303 + qJDD(3) * qJ(4) + t273 * t321 + t228;
t291 = sin(pkin(11));
t293 = cos(pkin(11));
t322 = qJD(1) * t297;
t269 = qJD(3) * t291 + t293 * t322;
t197 = -0.2e1 * qJD(4) * t269 + t293 * t220 - t291 * t226;
t252 = qJDD(3) * t291 + t276 * t293;
t268 = qJD(3) * t293 - t291 * t322;
t191 = (-t268 * t321 - t252) * pkin(8) + (t268 * t269 - t277) * pkin(4) + t197;
t198 = 0.2e1 * qJD(4) * t268 + t291 * t220 + t293 * t226;
t251 = qJDD(3) * t293 - t276 * t291;
t253 = -pkin(4) * t321 - pkin(8) * t269;
t267 = t268 ^ 2;
t196 = -pkin(4) * t267 + pkin(8) * t251 + t253 * t321 + t198;
t296 = sin(qJ(5));
t300 = cos(qJ(5));
t181 = t300 * t191 - t296 * t196;
t242 = t268 * t300 - t269 * t296;
t213 = qJD(5) * t242 + t251 * t296 + t252 * t300;
t243 = t268 * t296 + t269 * t300;
t271 = qJDD(5) - t277;
t284 = qJD(5) - t321;
t178 = (t242 * t284 - t213) * pkin(9) + (t242 * t243 + t271) * pkin(5) + t181;
t182 = t296 * t191 + t300 * t196;
t212 = -qJD(5) * t243 + t251 * t300 - t252 * t296;
t231 = pkin(5) * t284 - pkin(9) * t243;
t241 = t242 ^ 2;
t179 = -pkin(5) * t241 + pkin(9) * t212 - t231 * t284 + t182;
t295 = sin(qJ(6));
t299 = cos(qJ(6));
t177 = t178 * t295 + t179 * t299;
t227 = -t297 * t235 + t290 * t301;
t224 = -qJDD(3) * pkin(3) - qJ(4) * t303 + t273 * t322 + qJDD(4) - t227;
t205 = -pkin(4) * t251 - pkin(8) * t267 + t269 * t253 + t224;
t184 = -pkin(5) * t212 - pkin(9) * t241 + t231 * t243 + t205;
t223 = t242 * t295 + t243 * t299;
t192 = -qJD(6) * t223 + t212 * t299 - t213 * t295;
t222 = t242 * t299 - t243 * t295;
t193 = qJD(6) * t222 + t212 * t295 + t213 * t299;
t283 = qJD(6) + t284;
t199 = Ifges(7,5) * t223 + Ifges(7,6) * t222 + Ifges(7,3) * t283;
t201 = Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t283;
t265 = qJDD(6) + t271;
t165 = -mrSges(7,1) * t184 + mrSges(7,3) * t177 + Ifges(7,4) * t193 + Ifges(7,2) * t192 + Ifges(7,6) * t265 - t199 * t223 + t201 * t283;
t176 = t178 * t299 - t179 * t295;
t200 = Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t283;
t166 = mrSges(7,2) * t184 - mrSges(7,3) * t176 + Ifges(7,1) * t193 + Ifges(7,4) * t192 + Ifges(7,5) * t265 + t199 * t222 - t200 * t283;
t216 = Ifges(6,5) * t243 + Ifges(6,6) * t242 + Ifges(6,3) * t284;
t218 = Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t284;
t207 = -mrSges(7,2) * t283 + mrSges(7,3) * t222;
t208 = mrSges(7,1) * t283 - mrSges(7,3) * t223;
t314 = m(7) * t184 - t192 * mrSges(7,1) + t193 * mrSges(7,2) - t222 * t207 + t223 * t208;
t204 = -mrSges(7,1) * t222 + mrSges(7,2) * t223;
t170 = m(7) * t176 + mrSges(7,1) * t265 - mrSges(7,3) * t193 - t204 * t223 + t207 * t283;
t171 = m(7) * t177 - mrSges(7,2) * t265 + mrSges(7,3) * t192 + t204 * t222 - t208 * t283;
t315 = -t170 * t295 + t299 * t171;
t152 = -mrSges(6,1) * t205 + mrSges(6,3) * t182 + Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t271 - pkin(5) * t314 + pkin(9) * t315 + t299 * t165 + t295 * t166 - t243 * t216 + t284 * t218;
t164 = t299 * t170 + t295 * t171;
t217 = Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t284;
t153 = mrSges(6,2) * t205 - mrSges(6,3) * t181 + Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t271 - pkin(9) * t164 - t165 * t295 + t166 * t299 + t216 * t242 - t217 * t284;
t236 = Ifges(5,5) * t269 + Ifges(5,6) * t268 - Ifges(5,3) * t321;
t238 = Ifges(5,1) * t269 + Ifges(5,4) * t268 - Ifges(5,5) * t321;
t229 = -mrSges(6,2) * t284 + mrSges(6,3) * t242;
t230 = mrSges(6,1) * t284 - mrSges(6,3) * t243;
t309 = m(6) * t205 - t212 * mrSges(6,1) + t213 * mrSges(6,2) - t242 * t229 + t243 * t230 + t314;
t225 = -mrSges(6,1) * t242 + mrSges(6,2) * t243;
t161 = m(6) * t181 + mrSges(6,1) * t271 - mrSges(6,3) * t213 - t225 * t243 + t229 * t284 + t164;
t162 = m(6) * t182 - mrSges(6,2) * t271 + mrSges(6,3) * t212 + t225 * t242 - t230 * t284 + t315;
t316 = -t161 * t296 + t300 * t162;
t137 = -mrSges(5,1) * t224 + mrSges(5,3) * t198 + Ifges(5,4) * t252 + Ifges(5,2) * t251 - Ifges(5,6) * t277 - pkin(4) * t309 + pkin(8) * t316 + t300 * t152 + t296 * t153 - t269 * t236 - t238 * t321;
t157 = t300 * t161 + t296 * t162;
t237 = Ifges(5,4) * t269 + Ifges(5,2) * t268 - Ifges(5,6) * t321;
t138 = mrSges(5,2) * t224 - mrSges(5,3) * t197 + Ifges(5,1) * t252 + Ifges(5,4) * t251 - Ifges(5,5) * t277 - pkin(8) * t157 - t152 * t296 + t153 * t300 + t236 * t268 + t237 * t321;
t247 = -mrSges(5,1) * t268 + mrSges(5,2) * t269;
t249 = mrSges(5,2) * t321 + mrSges(5,3) * t268;
t155 = m(5) * t197 - mrSges(5,1) * t277 - mrSges(5,3) * t252 - t247 * t269 - t249 * t321 + t157;
t250 = -mrSges(5,1) * t321 - mrSges(5,3) * t269;
t156 = m(5) * t198 + mrSges(5,2) * t277 + mrSges(5,3) * t251 + t247 * t268 + t250 * t321 + t316;
t151 = -t155 * t291 + t293 * t156;
t174 = -m(5) * t224 + t251 * mrSges(5,1) - t252 * mrSges(5,2) + t268 * t249 - t269 * t250 - t309;
t262 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t297 + Ifges(4,2) * t301) * qJD(1);
t263 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t297 + Ifges(4,4) * t301) * qJD(1);
t323 = mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t276 + Ifges(4,6) * t277 + Ifges(4,3) * qJDD(3) + pkin(3) * t174 + qJ(4) * t151 + t293 * t137 + t291 * t138 + (t262 * t297 - t263 * t301) * qJD(1);
t274 = (-mrSges(4,1) * t301 + mrSges(4,2) * t297) * qJD(1);
t279 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t322;
t149 = m(4) * t228 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t277 - qJD(3) * t279 + t274 * t321 + t151;
t280 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t321;
t173 = m(4) * t227 + qJDD(3) * mrSges(4,1) - t276 * mrSges(4,3) + qJD(3) * t280 - t274 * t322 + t174;
t317 = t301 * t149 - t173 * t297;
t141 = m(3) * t246 - mrSges(3,1) * t304 - qJDD(1) * mrSges(3,2) + t317;
t150 = t155 * t293 + t156 * t291;
t308 = -m(4) * t234 + t277 * mrSges(4,1) - mrSges(4,2) * t276 - t279 * t322 + t280 * t321 - t150;
t145 = m(3) * t245 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t304 + t308;
t134 = t292 * t141 + t294 * t145;
t143 = t297 * t149 + t301 * t173;
t318 = t294 * t141 - t145 * t292;
t261 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t297 + Ifges(4,6) * t301) * qJD(1);
t130 = mrSges(4,2) * t234 - mrSges(4,3) * t227 + Ifges(4,1) * t276 + Ifges(4,4) * t277 + Ifges(4,5) * qJDD(3) - qJ(4) * t150 - qJD(3) * t262 - t137 * t291 + t138 * t293 + t261 * t321;
t311 = -mrSges(7,1) * t176 + mrSges(7,2) * t177 - Ifges(7,5) * t193 - Ifges(7,6) * t192 - Ifges(7,3) * t265 - t223 * t200 + t222 * t201;
t307 = -mrSges(6,1) * t181 + mrSges(6,2) * t182 - Ifges(6,5) * t213 - Ifges(6,6) * t212 - Ifges(6,3) * t271 - pkin(5) * t164 - t243 * t217 + t242 * t218 + t311;
t305 = mrSges(5,1) * t197 - mrSges(5,2) * t198 + Ifges(5,5) * t252 + Ifges(5,6) * t251 + pkin(4) * t157 + t269 * t237 - t268 * t238 - t307;
t136 = Ifges(4,6) * qJDD(3) + (Ifges(4,2) + Ifges(5,3)) * t277 + Ifges(4,4) * t276 + qJD(3) * t263 + mrSges(4,3) * t228 - mrSges(4,1) * t234 - pkin(3) * t150 - t261 * t322 - t305;
t312 = mrSges(3,1) * t245 - mrSges(3,2) * t246 + Ifges(3,3) * qJDD(1) + pkin(2) * t308 + pkin(7) * t317 + t297 * t130 + t301 * t136;
t310 = mrSges(2,1) * t281 - mrSges(2,2) * t282 + Ifges(2,3) * qJDD(1) + pkin(1) * t134 + t312;
t132 = m(2) * t282 - mrSges(2,1) * t304 - qJDD(1) * mrSges(2,2) + t318;
t131 = m(2) * t281 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t304 + t134;
t128 = -mrSges(3,1) * t290 + mrSges(3,3) * t246 + t304 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t143 - t323;
t127 = mrSges(3,2) * t290 - mrSges(3,3) * t245 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t304 - pkin(7) * t143 + t130 * t301 - t136 * t297;
t126 = -mrSges(2,2) * g(3) - mrSges(2,3) * t281 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t304 - qJ(2) * t134 + t127 * t294 - t128 * t292;
t125 = Ifges(2,6) * qJDD(1) + t304 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t282 + t292 * t127 + t294 * t128 - pkin(1) * (m(3) * t290 + t143) + qJ(2) * t318;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t302 * t126 - t298 * t125 - pkin(6) * (t131 * t302 + t132 * t298), t126, t127, t130, t138, t153, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t126 + t302 * t125 + pkin(6) * (-t131 * t298 + t132 * t302), t125, t128, t136, t137, t152, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t310, t310, t312, t323, -Ifges(5,3) * t277 + t305, -t307, -t311;];
m_new  = t1;

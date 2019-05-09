% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:57:00
% EndTime: 2019-05-06 02:57:41
% DurationCPUTime: 24.57s
% Computational Cost: add. (482044->341), mult. (950945->426), div. (0->0), fcn. (656287->12), ass. (0->139)
t297 = sin(qJ(1));
t302 = cos(qJ(1));
t281 = t297 * g(1) - t302 * g(2);
t271 = qJDD(1) * pkin(1) + t281;
t282 = -t302 * g(1) - t297 * g(2);
t304 = qJD(1) ^ 2;
t273 = -t304 * pkin(1) + t282;
t291 = sin(pkin(11));
t292 = cos(pkin(11));
t248 = t292 * t271 - t291 * t273;
t234 = -qJDD(1) * pkin(2) - t304 * pkin(7) - t248;
t296 = sin(qJ(3));
t301 = cos(qJ(3));
t320 = qJD(1) * qJD(3);
t319 = t301 * t320;
t275 = t296 * qJDD(1) + t319;
t286 = t296 * t320;
t276 = t301 * qJDD(1) - t286;
t219 = (-t275 - t319) * pkin(8) + (-t276 + t286) * pkin(3) + t234;
t249 = t291 * t271 + t292 * t273;
t235 = -t304 * pkin(2) + qJDD(1) * pkin(7) + t249;
t290 = -g(3) + qJDD(2);
t228 = t301 * t235 + t296 * t290;
t274 = (-pkin(3) * t301 - pkin(8) * t296) * qJD(1);
t303 = qJD(3) ^ 2;
t321 = t301 * qJD(1);
t225 = -t303 * pkin(3) + qJDD(3) * pkin(8) + t274 * t321 + t228;
t295 = sin(qJ(4));
t300 = cos(qJ(4));
t203 = t300 * t219 - t295 * t225;
t322 = qJD(1) * t296;
t269 = t300 * qJD(3) - t295 * t322;
t243 = t269 * qJD(4) + t295 * qJDD(3) + t300 * t275;
t268 = qJDD(4) - t276;
t270 = t295 * qJD(3) + t300 * t322;
t284 = qJD(4) - t321;
t193 = (t269 * t284 - t243) * pkin(9) + (t269 * t270 + t268) * pkin(4) + t203;
t204 = t295 * t219 + t300 * t225;
t242 = -t270 * qJD(4) + t300 * qJDD(3) - t295 * t275;
t253 = t284 * pkin(4) - t270 * pkin(9);
t267 = t269 ^ 2;
t195 = -t267 * pkin(4) + t242 * pkin(9) - t284 * t253 + t204;
t294 = sin(qJ(5));
t299 = cos(qJ(5));
t181 = t299 * t193 - t294 * t195;
t246 = t299 * t269 - t294 * t270;
t211 = t246 * qJD(5) + t294 * t242 + t299 * t243;
t247 = t294 * t269 + t299 * t270;
t264 = qJDD(5) + t268;
t283 = qJD(5) + t284;
t178 = (t246 * t283 - t211) * pkin(10) + (t246 * t247 + t264) * pkin(5) + t181;
t182 = t294 * t193 + t299 * t195;
t210 = -t247 * qJD(5) + t299 * t242 - t294 * t243;
t231 = t283 * pkin(5) - t247 * pkin(10);
t245 = t246 ^ 2;
t179 = -t245 * pkin(5) + t210 * pkin(10) - t283 * t231 + t182;
t293 = sin(qJ(6));
t298 = cos(qJ(6));
t177 = t293 * t178 + t298 * t179;
t227 = -t296 * t235 + t301 * t290;
t224 = -qJDD(3) * pkin(3) - t303 * pkin(8) + t274 * t322 - t227;
t202 = -t242 * pkin(4) - t267 * pkin(9) + t270 * t253 + t224;
t184 = -t210 * pkin(5) - t245 * pkin(10) + t247 * t231 + t202;
t223 = t293 * t246 + t298 * t247;
t189 = -t223 * qJD(6) + t298 * t210 - t293 * t211;
t222 = t298 * t246 - t293 * t247;
t190 = t222 * qJD(6) + t293 * t210 + t298 * t211;
t278 = qJD(6) + t283;
t197 = Ifges(7,5) * t223 + Ifges(7,6) * t222 + Ifges(7,3) * t278;
t199 = Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t278;
t257 = qJDD(6) + t264;
t165 = -mrSges(7,1) * t184 + mrSges(7,3) * t177 + Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t257 - t223 * t197 + t278 * t199;
t176 = t298 * t178 - t293 * t179;
t198 = Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t278;
t166 = mrSges(7,2) * t184 - mrSges(7,3) * t176 + Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t257 + t222 * t197 - t278 * t198;
t216 = Ifges(6,5) * t247 + Ifges(6,6) * t246 + Ifges(6,3) * t283;
t218 = Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t283;
t212 = -t278 * mrSges(7,2) + t222 * mrSges(7,3);
t213 = t278 * mrSges(7,1) - t223 * mrSges(7,3);
t314 = m(7) * t184 - t189 * mrSges(7,1) + t190 * mrSges(7,2) - t222 * t212 + t223 * t213;
t205 = -t222 * mrSges(7,1) + t223 * mrSges(7,2);
t172 = m(7) * t176 + t257 * mrSges(7,1) - t190 * mrSges(7,3) - t223 * t205 + t278 * t212;
t173 = m(7) * t177 - t257 * mrSges(7,2) + t189 * mrSges(7,3) + t222 * t205 - t278 * t213;
t315 = -t293 * t172 + t298 * t173;
t152 = -mrSges(6,1) * t202 + mrSges(6,3) * t182 + Ifges(6,4) * t211 + Ifges(6,2) * t210 + Ifges(6,6) * t264 - pkin(5) * t314 + pkin(10) * t315 + t298 * t165 + t293 * t166 - t247 * t216 + t283 * t218;
t164 = t298 * t172 + t293 * t173;
t217 = Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t283;
t153 = mrSges(6,2) * t202 - mrSges(6,3) * t181 + Ifges(6,1) * t211 + Ifges(6,4) * t210 + Ifges(6,5) * t264 - pkin(10) * t164 - t293 * t165 + t298 * t166 + t246 * t216 - t283 * t217;
t236 = Ifges(5,5) * t270 + Ifges(5,6) * t269 + Ifges(5,3) * t284;
t238 = Ifges(5,1) * t270 + Ifges(5,4) * t269 + Ifges(5,5) * t284;
t229 = -t283 * mrSges(6,2) + t246 * mrSges(6,3);
t230 = t283 * mrSges(6,1) - t247 * mrSges(6,3);
t309 = m(6) * t202 - t210 * mrSges(6,1) + t211 * mrSges(6,2) - t246 * t229 + t247 * t230 + t314;
t226 = -t246 * mrSges(6,1) + t247 * mrSges(6,2);
t161 = m(6) * t181 + t264 * mrSges(6,1) - t211 * mrSges(6,3) - t247 * t226 + t283 * t229 + t164;
t162 = m(6) * t182 - t264 * mrSges(6,2) + t210 * mrSges(6,3) + t246 * t226 - t283 * t230 + t315;
t316 = -t294 * t161 + t299 * t162;
t137 = -mrSges(5,1) * t224 + mrSges(5,3) * t204 + Ifges(5,4) * t243 + Ifges(5,2) * t242 + Ifges(5,6) * t268 - pkin(4) * t309 + pkin(9) * t316 + t299 * t152 + t294 * t153 - t270 * t236 + t284 * t238;
t157 = t299 * t161 + t294 * t162;
t237 = Ifges(5,4) * t270 + Ifges(5,2) * t269 + Ifges(5,6) * t284;
t138 = mrSges(5,2) * t224 - mrSges(5,3) * t203 + Ifges(5,1) * t243 + Ifges(5,4) * t242 + Ifges(5,5) * t268 - pkin(9) * t157 - t294 * t152 + t299 * t153 + t269 * t236 - t284 * t237;
t250 = -t269 * mrSges(5,1) + t270 * mrSges(5,2);
t251 = -t284 * mrSges(5,2) + t269 * mrSges(5,3);
t155 = m(5) * t203 + t268 * mrSges(5,1) - t243 * mrSges(5,3) - t270 * t250 + t284 * t251 + t157;
t252 = t284 * mrSges(5,1) - t270 * mrSges(5,3);
t156 = m(5) * t204 - t268 * mrSges(5,2) + t242 * mrSges(5,3) + t269 * t250 - t284 * t252 + t316;
t151 = -t295 * t155 + t300 * t156;
t174 = -m(5) * t224 + t242 * mrSges(5,1) - t243 * mrSges(5,2) + t269 * t251 - t270 * t252 - t309;
t262 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t296 + Ifges(4,2) * t301) * qJD(1);
t263 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t296 + Ifges(4,4) * t301) * qJD(1);
t323 = mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t275 + Ifges(4,6) * t276 + Ifges(4,3) * qJDD(3) + pkin(3) * t174 + pkin(8) * t151 + t300 * t137 + t295 * t138 + (t296 * t262 - t301 * t263) * qJD(1);
t272 = (-mrSges(4,1) * t301 + mrSges(4,2) * t296) * qJD(1);
t279 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t322;
t149 = m(4) * t228 - qJDD(3) * mrSges(4,2) + t276 * mrSges(4,3) - qJD(3) * t279 + t272 * t321 + t151;
t280 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t321;
t168 = m(4) * t227 + qJDD(3) * mrSges(4,1) - t275 * mrSges(4,3) + qJD(3) * t280 - t272 * t322 + t174;
t317 = t301 * t149 - t296 * t168;
t141 = m(3) * t249 - t304 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t317;
t150 = t300 * t155 + t295 * t156;
t308 = -m(4) * t234 + t276 * mrSges(4,1) - t275 * mrSges(4,2) - t279 * t322 + t280 * t321 - t150;
t145 = m(3) * t248 + qJDD(1) * mrSges(3,1) - t304 * mrSges(3,2) + t308;
t134 = t291 * t141 + t292 * t145;
t143 = t296 * t149 + t301 * t168;
t318 = t292 * t141 - t291 * t145;
t261 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t296 + Ifges(4,6) * t301) * qJD(1);
t130 = mrSges(4,2) * t234 - mrSges(4,3) * t227 + Ifges(4,1) * t275 + Ifges(4,4) * t276 + Ifges(4,5) * qJDD(3) - pkin(8) * t150 - qJD(3) * t262 - t295 * t137 + t300 * t138 + t261 * t321;
t311 = -mrSges(7,1) * t176 + mrSges(7,2) * t177 - Ifges(7,5) * t190 - Ifges(7,6) * t189 - Ifges(7,3) * t257 - t223 * t198 + t222 * t199;
t307 = -mrSges(6,1) * t181 + mrSges(6,2) * t182 - Ifges(6,5) * t211 - Ifges(6,6) * t210 - Ifges(6,3) * t264 - pkin(5) * t164 - t247 * t217 + t246 * t218 + t311;
t305 = mrSges(5,1) * t203 - mrSges(5,2) * t204 + Ifges(5,5) * t243 + Ifges(5,6) * t242 + Ifges(5,3) * t268 + pkin(4) * t157 + t270 * t237 - t269 * t238 - t307;
t136 = -mrSges(4,1) * t234 + mrSges(4,3) * t228 + Ifges(4,4) * t275 + Ifges(4,2) * t276 + Ifges(4,6) * qJDD(3) - pkin(3) * t150 + qJD(3) * t263 - t261 * t322 - t305;
t312 = mrSges(3,1) * t248 - mrSges(3,2) * t249 + Ifges(3,3) * qJDD(1) + pkin(2) * t308 + pkin(7) * t317 + t296 * t130 + t301 * t136;
t310 = mrSges(2,1) * t281 - mrSges(2,2) * t282 + Ifges(2,3) * qJDD(1) + pkin(1) * t134 + t312;
t132 = m(2) * t282 - t304 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t318;
t131 = m(2) * t281 + qJDD(1) * mrSges(2,1) - t304 * mrSges(2,2) + t134;
t128 = -mrSges(3,1) * t290 + mrSges(3,3) * t249 + t304 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t143 - t323;
t127 = mrSges(3,2) * t290 - mrSges(3,3) * t248 + Ifges(3,5) * qJDD(1) - t304 * Ifges(3,6) - pkin(7) * t143 + t301 * t130 - t296 * t136;
t126 = -mrSges(2,2) * g(3) - mrSges(2,3) * t281 + Ifges(2,5) * qJDD(1) - t304 * Ifges(2,6) - qJ(2) * t134 + t292 * t127 - t291 * t128;
t125 = Ifges(2,6) * qJDD(1) + t304 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t282 + t291 * t127 + t292 * t128 - pkin(1) * (m(3) * t290 + t143) + qJ(2) * t318;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t302 * t126 - t297 * t125 - pkin(6) * (t302 * t131 + t297 * t132), t126, t127, t130, t138, t153, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t297 * t126 + t302 * t125 + pkin(6) * (-t297 * t131 + t302 * t132), t125, t128, t136, t137, t152, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t310, t310, t312, t323, t305, -t307, -t311;];
m_new  = t1;

% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:04:10
% EndTime: 2019-05-05 22:04:48
% DurationCPUTime: 23.73s
% Computational Cost: add. (457197->340), mult. (921292->428), div. (0->0), fcn. (627319->12), ass. (0->137)
t298 = sin(qJ(1));
t302 = cos(qJ(1));
t281 = t298 * g(1) - t302 * g(2);
t272 = qJDD(1) * pkin(1) + t281;
t282 = -t302 * g(1) - t298 * g(2);
t304 = qJD(1) ^ 2;
t274 = -t304 * pkin(1) + t282;
t292 = sin(pkin(10));
t294 = cos(pkin(10));
t250 = t294 * t272 - t292 * t274;
t234 = -qJDD(1) * pkin(2) - t304 * pkin(7) - t250;
t297 = sin(qJ(3));
t301 = cos(qJ(3));
t320 = qJD(1) * qJD(3);
t319 = t301 * t320;
t276 = t297 * qJDD(1) + t319;
t286 = t297 * t320;
t277 = t301 * qJDD(1) - t286;
t217 = (-t276 - t319) * pkin(8) + (-t277 + t286) * pkin(3) + t234;
t251 = t292 * t272 + t294 * t274;
t235 = -t304 * pkin(2) + qJDD(1) * pkin(7) + t251;
t290 = -g(3) + qJDD(2);
t228 = t301 * t235 + t297 * t290;
t275 = (-pkin(3) * t301 - pkin(8) * t297) * qJD(1);
t303 = qJD(3) ^ 2;
t321 = t301 * qJD(1);
t225 = -t303 * pkin(3) + qJDD(3) * pkin(8) + t275 * t321 + t228;
t296 = sin(qJ(4));
t300 = cos(qJ(4));
t203 = t300 * t217 - t296 * t225;
t322 = qJD(1) * t297;
t270 = t300 * qJD(3) - t296 * t322;
t245 = t270 * qJD(4) + t296 * qJDD(3) + t300 * t276;
t269 = qJDD(4) - t277;
t271 = t296 * qJD(3) + t300 * t322;
t284 = qJD(4) - t321;
t187 = (t270 * t284 - t245) * qJ(5) + (t270 * t271 + t269) * pkin(4) + t203;
t204 = t296 * t217 + t300 * t225;
t244 = -t271 * qJD(4) + t300 * qJDD(3) - t296 * t276;
t254 = t284 * pkin(4) - t271 * qJ(5);
t268 = t270 ^ 2;
t195 = -t268 * pkin(4) + t244 * qJ(5) - t284 * t254 + t204;
t291 = sin(pkin(11));
t293 = cos(pkin(11));
t249 = t291 * t270 + t293 * t271;
t181 = -0.2e1 * qJD(5) * t249 + t293 * t187 - t291 * t195;
t219 = t291 * t244 + t293 * t245;
t248 = t293 * t270 - t291 * t271;
t178 = (t248 * t284 - t219) * pkin(9) + (t248 * t249 + t269) * pkin(5) + t181;
t182 = 0.2e1 * qJD(5) * t248 + t291 * t187 + t293 * t195;
t218 = t293 * t244 - t291 * t245;
t231 = t284 * pkin(5) - t249 * pkin(9);
t246 = t248 ^ 2;
t179 = -t246 * pkin(5) + t218 * pkin(9) - t284 * t231 + t182;
t295 = sin(qJ(6));
t299 = cos(qJ(6));
t177 = t295 * t178 + t299 * t179;
t227 = -t297 * t235 + t301 * t290;
t224 = -qJDD(3) * pkin(3) - t303 * pkin(8) + t275 * t322 - t227;
t200 = -t244 * pkin(4) - t268 * qJ(5) + t271 * t254 + qJDD(5) + t224;
t184 = -t218 * pkin(5) - t246 * pkin(9) + t249 * t231 + t200;
t223 = t295 * t248 + t299 * t249;
t193 = -t223 * qJD(6) + t299 * t218 - t295 * t219;
t222 = t299 * t248 - t295 * t249;
t194 = t222 * qJD(6) + t295 * t218 + t299 * t219;
t283 = qJD(6) + t284;
t197 = Ifges(7,5) * t223 + Ifges(7,6) * t222 + Ifges(7,3) * t283;
t199 = Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t283;
t265 = qJDD(6) + t269;
t165 = -mrSges(7,1) * t184 + mrSges(7,3) * t177 + Ifges(7,4) * t194 + Ifges(7,2) * t193 + Ifges(7,6) * t265 - t223 * t197 + t283 * t199;
t176 = t299 * t178 - t295 * t179;
t198 = Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t283;
t166 = mrSges(7,2) * t184 - mrSges(7,3) * t176 + Ifges(7,1) * t194 + Ifges(7,4) * t193 + Ifges(7,5) * t265 + t222 * t197 - t283 * t198;
t214 = Ifges(6,5) * t249 + Ifges(6,6) * t248 + Ifges(6,3) * t284;
t216 = Ifges(6,1) * t249 + Ifges(6,4) * t248 + Ifges(6,5) * t284;
t207 = -t283 * mrSges(7,2) + t222 * mrSges(7,3);
t208 = t283 * mrSges(7,1) - t223 * mrSges(7,3);
t314 = m(7) * t184 - t193 * mrSges(7,1) + t194 * mrSges(7,2) - t222 * t207 + t223 * t208;
t205 = -t222 * mrSges(7,1) + t223 * mrSges(7,2);
t170 = m(7) * t176 + t265 * mrSges(7,1) - t194 * mrSges(7,3) - t223 * t205 + t283 * t207;
t171 = m(7) * t177 - t265 * mrSges(7,2) + t193 * mrSges(7,3) + t222 * t205 - t283 * t208;
t315 = -t295 * t170 + t299 * t171;
t152 = -mrSges(6,1) * t200 + mrSges(6,3) * t182 + Ifges(6,4) * t219 + Ifges(6,2) * t218 + Ifges(6,6) * t269 - pkin(5) * t314 + pkin(9) * t315 + t299 * t165 + t295 * t166 - t249 * t214 + t284 * t216;
t164 = t299 * t170 + t295 * t171;
t215 = Ifges(6,4) * t249 + Ifges(6,2) * t248 + Ifges(6,6) * t284;
t153 = mrSges(6,2) * t200 - mrSges(6,3) * t181 + Ifges(6,1) * t219 + Ifges(6,4) * t218 + Ifges(6,5) * t269 - pkin(9) * t164 - t295 * t165 + t299 * t166 + t248 * t214 - t284 * t215;
t236 = Ifges(5,5) * t271 + Ifges(5,6) * t270 + Ifges(5,3) * t284;
t238 = Ifges(5,1) * t271 + Ifges(5,4) * t270 + Ifges(5,5) * t284;
t229 = -t284 * mrSges(6,2) + t248 * mrSges(6,3);
t230 = t284 * mrSges(6,1) - t249 * mrSges(6,3);
t309 = m(6) * t200 - t218 * mrSges(6,1) + t219 * mrSges(6,2) - t248 * t229 + t249 * t230 + t314;
t226 = -t248 * mrSges(6,1) + t249 * mrSges(6,2);
t161 = m(6) * t181 + t269 * mrSges(6,1) - t219 * mrSges(6,3) - t249 * t226 + t284 * t229 + t164;
t162 = m(6) * t182 - t269 * mrSges(6,2) + t218 * mrSges(6,3) + t248 * t226 - t284 * t230 + t315;
t316 = -t291 * t161 + t293 * t162;
t137 = -mrSges(5,1) * t224 + mrSges(5,3) * t204 + Ifges(5,4) * t245 + Ifges(5,2) * t244 + Ifges(5,6) * t269 - pkin(4) * t309 + qJ(5) * t316 + t293 * t152 + t291 * t153 - t271 * t236 + t284 * t238;
t157 = t293 * t161 + t291 * t162;
t237 = Ifges(5,4) * t271 + Ifges(5,2) * t270 + Ifges(5,6) * t284;
t138 = mrSges(5,2) * t224 - mrSges(5,3) * t203 + Ifges(5,1) * t245 + Ifges(5,4) * t244 + Ifges(5,5) * t269 - qJ(5) * t157 - t291 * t152 + t293 * t153 + t270 * t236 - t284 * t237;
t252 = -t270 * mrSges(5,1) + t271 * mrSges(5,2);
t253 = -t284 * mrSges(5,2) + t270 * mrSges(5,3);
t155 = m(5) * t203 + t269 * mrSges(5,1) - t245 * mrSges(5,3) - t271 * t252 + t284 * t253 + t157;
t255 = t284 * mrSges(5,1) - t271 * mrSges(5,3);
t156 = m(5) * t204 - t269 * mrSges(5,2) + t244 * mrSges(5,3) + t270 * t252 - t284 * t255 + t316;
t151 = -t296 * t155 + t300 * t156;
t174 = -m(5) * t224 + t244 * mrSges(5,1) - t245 * mrSges(5,2) + t270 * t253 - t271 * t255 - t309;
t262 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t297 + Ifges(4,2) * t301) * qJD(1);
t263 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t297 + Ifges(4,4) * t301) * qJD(1);
t323 = mrSges(4,1) * t227 - mrSges(4,2) * t228 + Ifges(4,5) * t276 + Ifges(4,6) * t277 + Ifges(4,3) * qJDD(3) + pkin(3) * t174 + pkin(8) * t151 + t300 * t137 + t296 * t138 + (t297 * t262 - t301 * t263) * qJD(1);
t273 = (-mrSges(4,1) * t301 + mrSges(4,2) * t297) * qJD(1);
t279 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t322;
t149 = m(4) * t228 - qJDD(3) * mrSges(4,2) + t277 * mrSges(4,3) - qJD(3) * t279 + t273 * t321 + t151;
t280 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t321;
t173 = m(4) * t227 + qJDD(3) * mrSges(4,1) - t276 * mrSges(4,3) + qJD(3) * t280 - t273 * t322 + t174;
t317 = t301 * t149 - t297 * t173;
t141 = m(3) * t251 - t304 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t317;
t150 = t300 * t155 + t296 * t156;
t308 = -m(4) * t234 + t277 * mrSges(4,1) - t276 * mrSges(4,2) - t279 * t322 + t280 * t321 - t150;
t145 = m(3) * t250 + qJDD(1) * mrSges(3,1) - t304 * mrSges(3,2) + t308;
t134 = t292 * t141 + t294 * t145;
t143 = t297 * t149 + t301 * t173;
t318 = t294 * t141 - t292 * t145;
t261 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t297 + Ifges(4,6) * t301) * qJD(1);
t130 = mrSges(4,2) * t234 - mrSges(4,3) * t227 + Ifges(4,1) * t276 + Ifges(4,4) * t277 + Ifges(4,5) * qJDD(3) - pkin(8) * t150 - qJD(3) * t262 - t296 * t137 + t300 * t138 + t261 * t321;
t311 = -mrSges(7,1) * t176 + mrSges(7,2) * t177 - Ifges(7,5) * t194 - Ifges(7,6) * t193 - Ifges(7,3) * t265 - t223 * t198 + t222 * t199;
t307 = -mrSges(6,1) * t181 + mrSges(6,2) * t182 - Ifges(6,5) * t219 - Ifges(6,6) * t218 - Ifges(6,3) * t269 - pkin(5) * t164 - t249 * t215 + t248 * t216 + t311;
t305 = mrSges(5,1) * t203 - mrSges(5,2) * t204 + Ifges(5,5) * t245 + Ifges(5,6) * t244 + Ifges(5,3) * t269 + pkin(4) * t157 + t271 * t237 - t270 * t238 - t307;
t136 = -mrSges(4,1) * t234 + mrSges(4,3) * t228 + Ifges(4,4) * t276 + Ifges(4,2) * t277 + Ifges(4,6) * qJDD(3) - pkin(3) * t150 + qJD(3) * t263 - t261 * t322 - t305;
t312 = mrSges(3,1) * t250 - mrSges(3,2) * t251 + Ifges(3,3) * qJDD(1) + pkin(2) * t308 + pkin(7) * t317 + t297 * t130 + t301 * t136;
t310 = mrSges(2,1) * t281 - mrSges(2,2) * t282 + Ifges(2,3) * qJDD(1) + pkin(1) * t134 + t312;
t132 = m(2) * t282 - t304 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t318;
t131 = m(2) * t281 + qJDD(1) * mrSges(2,1) - t304 * mrSges(2,2) + t134;
t128 = -mrSges(3,1) * t290 + mrSges(3,3) * t251 + t304 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t143 - t323;
t127 = mrSges(3,2) * t290 - mrSges(3,3) * t250 + Ifges(3,5) * qJDD(1) - t304 * Ifges(3,6) - pkin(7) * t143 + t301 * t130 - t297 * t136;
t126 = -mrSges(2,2) * g(3) - mrSges(2,3) * t281 + Ifges(2,5) * qJDD(1) - t304 * Ifges(2,6) - qJ(2) * t134 + t294 * t127 - t292 * t128;
t125 = Ifges(2,6) * qJDD(1) + t304 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t282 + t292 * t127 + t294 * t128 - pkin(1) * (m(3) * t290 + t143) + qJ(2) * t318;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t302 * t126 - t298 * t125 - pkin(6) * (t302 * t131 + t298 * t132), t126, t127, t130, t138, t153, t166; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t126 + t302 * t125 + pkin(6) * (-t298 * t131 + t302 * t132), t125, t128, t136, t137, t152, t165; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t310, t310, t312, t323, t305, -t307, -t311;];
m_new  = t1;

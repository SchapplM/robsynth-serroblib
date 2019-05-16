% Calculate vector of cutting torques with Newton-Euler for
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-05-04 20:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:03:25
% EndTime: 2019-05-04 20:04:06
% DurationCPUTime: 32.14s
% Computational Cost: add. (619390->291), mult. (1173204->384), div. (0->0), fcn. (915565->16), ass. (0->134)
t262 = sin(pkin(11));
t267 = cos(pkin(11));
t254 = -g(1) * t267 - g(2) * t262;
t261 = sin(pkin(12));
t266 = cos(pkin(12));
t253 = g(1) * t262 - g(2) * t267;
t259 = -g(3) + qJDD(1);
t264 = sin(pkin(6));
t269 = cos(pkin(6));
t285 = t253 * t269 + t259 * t264;
t218 = -t261 * t254 + t266 * t285;
t219 = t266 * t254 + t261 * t285;
t229 = -t253 * t264 + t259 * t269 + qJDD(2);
t275 = cos(qJ(3));
t268 = cos(pkin(7));
t272 = sin(qJ(3));
t294 = t268 * t272;
t263 = sin(pkin(7));
t295 = t263 * t272;
t199 = t218 * t294 + t219 * t275 + t229 * t295;
t277 = qJD(3) ^ 2;
t197 = -pkin(3) * t277 + qJDD(3) * pkin(9) + t199;
t212 = -t218 * t263 + t229 * t268;
t271 = sin(qJ(4));
t274 = cos(qJ(4));
t190 = t197 * t274 + t212 * t271;
t248 = (-pkin(4) * t274 - qJ(5) * t271) * qJD(3);
t276 = qJD(4) ^ 2;
t292 = t274 * qJD(3);
t188 = -pkin(4) * t276 + qJDD(4) * qJ(5) + t248 * t292 + t190;
t198 = -t272 * t219 + (t218 * t268 + t229 * t263) * t275;
t196 = -qJDD(3) * pkin(3) - t277 * pkin(9) - t198;
t291 = qJD(3) * qJD(4);
t290 = t274 * t291;
t250 = qJDD(3) * t271 + t290;
t258 = t271 * t291;
t251 = qJDD(3) * t274 - t258;
t193 = (-t250 - t290) * qJ(5) + (-t251 + t258) * pkin(4) + t196;
t260 = sin(pkin(13));
t265 = cos(pkin(13));
t293 = qJD(3) * t271;
t245 = qJD(4) * t260 + t265 * t293;
t183 = -0.2e1 * qJD(5) * t245 - t260 * t188 + t193 * t265;
t233 = qJDD(4) * t260 + t250 * t265;
t244 = qJD(4) * t265 - t260 * t293;
t181 = (-t244 * t292 - t233) * pkin(10) + (t244 * t245 - t251) * pkin(5) + t183;
t184 = 0.2e1 * qJD(5) * t244 + t188 * t265 + t193 * t260;
t232 = qJDD(4) * t265 - t250 * t260;
t234 = -pkin(5) * t292 - pkin(10) * t245;
t243 = t244 ^ 2;
t182 = -pkin(5) * t243 + pkin(10) * t232 + t234 * t292 + t184;
t270 = sin(qJ(6));
t273 = cos(qJ(6));
t178 = t181 * t273 - t182 * t270;
t224 = t244 * t273 - t245 * t270;
t205 = qJD(6) * t224 + t232 * t270 + t233 * t273;
t225 = t244 * t270 + t245 * t273;
t211 = -mrSges(7,1) * t224 + mrSges(7,2) * t225;
t257 = qJD(6) - t292;
t214 = -mrSges(7,2) * t257 + mrSges(7,3) * t224;
t247 = qJDD(6) - t251;
t173 = m(7) * t178 + mrSges(7,1) * t247 - mrSges(7,3) * t205 - t211 * t225 + t214 * t257;
t179 = t181 * t270 + t182 * t273;
t204 = -qJD(6) * t225 + t232 * t273 - t233 * t270;
t215 = mrSges(7,1) * t257 - mrSges(7,3) * t225;
t174 = m(7) * t179 - mrSges(7,2) * t247 + mrSges(7,3) * t204 + t211 * t224 - t215 * t257;
t167 = t173 * t273 + t174 * t270;
t226 = -mrSges(6,1) * t244 + mrSges(6,2) * t245;
t230 = mrSges(6,2) * t292 + mrSges(6,3) * t244;
t165 = m(6) * t183 - mrSges(6,1) * t251 - mrSges(6,3) * t233 - t226 * t245 - t230 * t292 + t167;
t231 = -mrSges(6,1) * t292 - mrSges(6,3) * t245;
t288 = -t173 * t270 + t174 * t273;
t166 = m(6) * t184 + mrSges(6,2) * t251 + mrSges(6,3) * t232 + t226 * t244 + t231 * t292 + t288;
t163 = -t165 * t260 + t166 * t265;
t249 = (-mrSges(5,1) * t274 + mrSges(5,2) * t271) * qJD(3);
t255 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t293;
t161 = m(5) * t190 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t251 - qJD(4) * t255 + t249 * t292 + t163;
t189 = -t197 * t271 + t212 * t274;
t187 = -qJDD(4) * pkin(4) - qJ(5) * t276 + t248 * t293 + qJDD(5) - t189;
t185 = -pkin(5) * t232 - pkin(10) * t243 + t234 * t245 + t187;
t282 = m(7) * t185 - mrSges(7,1) * t204 + mrSges(7,2) * t205 - t214 * t224 + t215 * t225;
t180 = -m(6) * t187 + mrSges(6,1) * t232 - mrSges(6,2) * t233 + t230 * t244 - t231 * t245 - t282;
t256 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t292;
t176 = m(5) * t189 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t250 + qJD(4) * t256 - t249 * t293 + t180;
t289 = t161 * t274 - t176 * t271;
t150 = m(4) * t199 - mrSges(4,1) * t277 - qJDD(3) * mrSges(4,2) + t289;
t153 = t161 * t271 + t176 * t274;
t152 = m(4) * t212 + t153;
t162 = t165 * t265 + t166 * t260;
t280 = -m(5) * t196 + mrSges(5,1) * t251 - mrSges(5,2) * t250 - t255 * t293 + t256 * t292 - t162;
t158 = m(4) * t198 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t277 + t280;
t296 = t158 * t275;
t140 = t150 * t294 - t152 * t263 + t268 * t296;
t137 = m(3) * t218 + t140;
t145 = t150 * t275 - t158 * t272;
t144 = m(3) * t219 + t145;
t304 = t137 * t266 + t144 * t261;
t206 = Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t257;
t208 = Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t257;
t168 = -mrSges(7,1) * t185 + mrSges(7,3) * t179 + Ifges(7,4) * t205 + Ifges(7,2) * t204 + Ifges(7,6) * t247 - t206 * t225 + t208 * t257;
t207 = Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t257;
t169 = mrSges(7,2) * t185 - mrSges(7,3) * t178 + Ifges(7,1) * t205 + Ifges(7,4) * t204 + Ifges(7,5) * t247 + t206 * t224 - t207 * t257;
t220 = Ifges(6,5) * t245 + Ifges(6,6) * t244 - Ifges(6,3) * t292;
t222 = Ifges(6,1) * t245 + Ifges(6,4) * t244 - Ifges(6,5) * t292;
t154 = -mrSges(6,1) * t187 + mrSges(6,3) * t184 + Ifges(6,4) * t233 + Ifges(6,2) * t232 - Ifges(6,6) * t251 - pkin(5) * t282 + pkin(10) * t288 + t273 * t168 + t270 * t169 - t245 * t220 - t222 * t292;
t221 = Ifges(6,4) * t245 + Ifges(6,2) * t244 - Ifges(6,6) * t292;
t155 = mrSges(6,2) * t187 - mrSges(6,3) * t183 + Ifges(6,1) * t233 + Ifges(6,4) * t232 - Ifges(6,5) * t251 - pkin(10) * t167 - t168 * t270 + t169 * t273 + t220 * t244 + t221 * t292;
t238 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t271 + Ifges(5,6) * t274) * qJD(3);
t239 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t271 + Ifges(5,2) * t274) * qJD(3);
t141 = mrSges(5,2) * t196 - mrSges(5,3) * t189 + Ifges(5,1) * t250 + Ifges(5,4) * t251 + Ifges(5,5) * qJDD(4) - qJ(5) * t162 - qJD(4) * t239 - t154 * t260 + t155 * t265 + t238 * t292;
t240 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t271 + Ifges(5,4) * t274) * qJD(3);
t281 = -mrSges(7,1) * t178 + mrSges(7,2) * t179 - Ifges(7,5) * t205 - Ifges(7,6) * t204 - Ifges(7,3) * t247 - t207 * t225 + t224 * t208;
t278 = -mrSges(6,1) * t183 + mrSges(6,2) * t184 - Ifges(6,5) * t233 - Ifges(6,6) * t232 - pkin(5) * t167 - t245 * t221 + t244 * t222 + t281;
t146 = Ifges(5,6) * qJDD(4) + t278 + (Ifges(5,2) + Ifges(6,3)) * t251 + Ifges(5,4) * t250 + qJD(4) * t240 + mrSges(5,3) * t190 - mrSges(5,1) * t196 - pkin(4) * t162 - t238 * t293;
t130 = mrSges(4,1) * t198 - mrSges(4,2) * t199 + Ifges(4,3) * qJDD(3) + pkin(3) * t280 + pkin(9) * t289 + t271 * t141 + t274 * t146;
t139 = t150 * t295 + t152 * t268 + t263 * t296;
t131 = mrSges(4,2) * t212 - mrSges(4,3) * t198 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t277 - pkin(9) * t153 + t141 * t274 - t146 * t271;
t302 = mrSges(5,1) * t189 - mrSges(5,2) * t190 + Ifges(5,5) * t250 + Ifges(5,6) * t251 + Ifges(5,3) * qJDD(4) + pkin(4) * t180 + qJ(5) * t163 + t265 * t154 + t260 * t155 + (t239 * t271 - t240 * t274) * qJD(3);
t135 = -mrSges(4,1) * t212 + mrSges(4,3) * t199 + t277 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t153 - t302;
t284 = pkin(8) * t145 + t131 * t272 + t135 * t275;
t123 = -mrSges(3,1) * t229 + mrSges(3,3) * t219 - pkin(2) * t139 - t263 * t130 + t268 * t284;
t125 = mrSges(3,2) * t229 - mrSges(3,3) * t218 + t275 * t131 - t272 * t135 + (-t139 * t263 - t140 * t268) * pkin(8);
t134 = -t137 * t261 + t144 * t266;
t303 = qJ(2) * t134 + t123 * t266 + t125 * t261;
t138 = m(3) * t229 + t139;
t129 = -t264 * t138 + t269 * t304;
t121 = mrSges(3,1) * t218 - mrSges(3,2) * t219 + pkin(2) * t140 + t268 * t130 + t263 * t284;
t283 = mrSges(2,1) * t253 - mrSges(2,2) * t254 + pkin(1) * t129 + t269 * t121 + t264 * t303;
t132 = m(2) * t254 + t134;
t128 = t269 * t138 + t264 * t304;
t126 = m(2) * t253 + t129;
t119 = mrSges(2,2) * t259 - mrSges(2,3) * t253 - t261 * t123 + t266 * t125 + (-t128 * t264 - t129 * t269) * qJ(2);
t118 = -mrSges(2,1) * t259 + mrSges(2,3) * t254 - pkin(1) * t128 - t264 * t121 + t269 * t303;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t267 * t119 - t262 * t118 - qJ(1) * (t126 * t267 + t132 * t262), t119, t125, t131, t141, t155, t169; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t262 * t119 + t267 * t118 + qJ(1) * (-t126 * t262 + t132 * t267), t118, t123, t135, t146, t154, t168; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t283, t283, t121, t130, t302, -Ifges(6,3) * t251 - t278, -t281;];
m_new  = t1;

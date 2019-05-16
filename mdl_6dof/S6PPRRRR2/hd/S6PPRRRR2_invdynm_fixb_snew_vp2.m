% Calculate vector of cutting torques with Newton-Euler for
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:54:48
% EndTime: 2019-05-04 20:55:21
% DurationCPUTime: 32.30s
% Computational Cost: add. (658262->291), mult. (1199988->382), div. (0->0), fcn. (940363->16), ass. (0->136)
t261 = sin(pkin(12));
t265 = cos(pkin(12));
t253 = -t265 * g(1) - t261 * g(2);
t260 = sin(pkin(13));
t264 = cos(pkin(13));
t252 = t261 * g(1) - t265 * g(2);
t259 = -g(3) + qJDD(1);
t263 = sin(pkin(6));
t267 = cos(pkin(6));
t285 = t252 * t267 + t259 * t263;
t218 = -t260 * t253 + t264 * t285;
t219 = t264 * t253 + t260 * t285;
t234 = -t263 * t252 + t267 * t259 + qJDD(2);
t275 = cos(qJ(3));
t266 = cos(pkin(7));
t271 = sin(qJ(3));
t294 = t266 * t271;
t262 = sin(pkin(7));
t295 = t262 * t271;
t199 = t218 * t294 + t275 * t219 + t234 * t295;
t277 = qJD(3) ^ 2;
t197 = -t277 * pkin(3) + qJDD(3) * pkin(9) + t199;
t211 = -t262 * t218 + t266 * t234;
t270 = sin(qJ(4));
t274 = cos(qJ(4));
t190 = t274 * t197 + t270 * t211;
t248 = (-pkin(4) * t274 - pkin(10) * t270) * qJD(3);
t276 = qJD(4) ^ 2;
t292 = t274 * qJD(3);
t188 = -t276 * pkin(4) + qJDD(4) * pkin(10) + t248 * t292 + t190;
t198 = -t271 * t219 + (t218 * t266 + t234 * t262) * t275;
t196 = -qJDD(3) * pkin(3) - t277 * pkin(9) - t198;
t291 = qJD(3) * qJD(4);
t290 = t274 * t291;
t249 = t270 * qJDD(3) + t290;
t258 = t270 * t291;
t250 = t274 * qJDD(3) - t258;
t193 = (-t249 - t290) * pkin(10) + (-t250 + t258) * pkin(4) + t196;
t269 = sin(qJ(5));
t273 = cos(qJ(5));
t183 = -t269 * t188 + t273 * t193;
t293 = qJD(3) * t270;
t245 = t273 * qJD(4) - t269 * t293;
t226 = t245 * qJD(5) + t269 * qJDD(4) + t273 * t249;
t244 = qJDD(5) - t250;
t246 = t269 * qJD(4) + t273 * t293;
t257 = qJD(5) - t292;
t181 = (t245 * t257 - t226) * pkin(11) + (t245 * t246 + t244) * pkin(5) + t183;
t184 = t273 * t188 + t269 * t193;
t225 = -t246 * qJD(5) + t273 * qJDD(4) - t269 * t249;
t233 = t257 * pkin(5) - t246 * pkin(11);
t243 = t245 ^ 2;
t182 = -t243 * pkin(5) + t225 * pkin(11) - t257 * t233 + t184;
t268 = sin(qJ(6));
t272 = cos(qJ(6));
t179 = t272 * t181 - t268 * t182;
t227 = t272 * t245 - t268 * t246;
t205 = t227 * qJD(6) + t268 * t225 + t272 * t226;
t228 = t268 * t245 + t272 * t246;
t212 = -t227 * mrSges(7,1) + t228 * mrSges(7,2);
t256 = qJD(6) + t257;
t214 = -t256 * mrSges(7,2) + t227 * mrSges(7,3);
t240 = qJDD(6) + t244;
t175 = m(7) * t179 + t240 * mrSges(7,1) - t205 * mrSges(7,3) - t228 * t212 + t256 * t214;
t180 = t268 * t181 + t272 * t182;
t204 = -t228 * qJD(6) + t272 * t225 - t268 * t226;
t215 = t256 * mrSges(7,1) - t228 * mrSges(7,3);
t176 = m(7) * t180 - t240 * mrSges(7,2) + t204 * mrSges(7,3) + t227 * t212 - t256 * t215;
t167 = t272 * t175 + t268 * t176;
t229 = -t245 * mrSges(6,1) + t246 * mrSges(6,2);
t231 = -t257 * mrSges(6,2) + t245 * mrSges(6,3);
t165 = m(6) * t183 + t244 * mrSges(6,1) - t226 * mrSges(6,3) - t246 * t229 + t257 * t231 + t167;
t232 = t257 * mrSges(6,1) - t246 * mrSges(6,3);
t288 = -t268 * t175 + t272 * t176;
t166 = m(6) * t184 - t244 * mrSges(6,2) + t225 * mrSges(6,3) + t245 * t229 - t257 * t232 + t288;
t163 = -t269 * t165 + t273 * t166;
t247 = (-mrSges(5,1) * t274 + mrSges(5,2) * t270) * qJD(3);
t254 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t293;
t161 = m(5) * t190 - qJDD(4) * mrSges(5,2) + t250 * mrSges(5,3) - qJD(4) * t254 + t247 * t292 + t163;
t189 = -t270 * t197 + t274 * t211;
t187 = -qJDD(4) * pkin(4) - t276 * pkin(10) + t248 * t293 - t189;
t185 = -t225 * pkin(5) - t243 * pkin(11) + t246 * t233 + t187;
t282 = m(7) * t185 - t204 * mrSges(7,1) + t205 * mrSges(7,2) - t227 * t214 + t228 * t215;
t177 = -m(6) * t187 + t225 * mrSges(6,1) - t226 * mrSges(6,2) + t245 * t231 - t246 * t232 - t282;
t255 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t292;
t171 = m(5) * t189 + qJDD(4) * mrSges(5,1) - t249 * mrSges(5,3) + qJD(4) * t255 - t247 * t293 + t177;
t289 = t274 * t161 - t270 * t171;
t150 = m(4) * t199 - t277 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t289;
t153 = t270 * t161 + t274 * t171;
t152 = m(4) * t211 + t153;
t162 = t273 * t165 + t269 * t166;
t280 = -m(5) * t196 + t250 * mrSges(5,1) - t249 * mrSges(5,2) - t254 * t293 + t255 * t292 - t162;
t158 = m(4) * t198 + qJDD(3) * mrSges(4,1) - t277 * mrSges(4,2) + t280;
t296 = t158 * t275;
t140 = t150 * t294 - t262 * t152 + t266 * t296;
t137 = m(3) * t218 + t140;
t145 = t275 * t150 - t271 * t158;
t144 = m(3) * t219 + t145;
t304 = t137 * t264 + t144 * t260;
t206 = Ifges(7,5) * t228 + Ifges(7,6) * t227 + Ifges(7,3) * t256;
t208 = Ifges(7,1) * t228 + Ifges(7,4) * t227 + Ifges(7,5) * t256;
t168 = -mrSges(7,1) * t185 + mrSges(7,3) * t180 + Ifges(7,4) * t205 + Ifges(7,2) * t204 + Ifges(7,6) * t240 - t228 * t206 + t256 * t208;
t207 = Ifges(7,4) * t228 + Ifges(7,2) * t227 + Ifges(7,6) * t256;
t169 = mrSges(7,2) * t185 - mrSges(7,3) * t179 + Ifges(7,1) * t205 + Ifges(7,4) * t204 + Ifges(7,5) * t240 + t227 * t206 - t256 * t207;
t220 = Ifges(6,5) * t246 + Ifges(6,6) * t245 + Ifges(6,3) * t257;
t222 = Ifges(6,1) * t246 + Ifges(6,4) * t245 + Ifges(6,5) * t257;
t154 = -mrSges(6,1) * t187 + mrSges(6,3) * t184 + Ifges(6,4) * t226 + Ifges(6,2) * t225 + Ifges(6,6) * t244 - pkin(5) * t282 + pkin(11) * t288 + t272 * t168 + t268 * t169 - t246 * t220 + t257 * t222;
t221 = Ifges(6,4) * t246 + Ifges(6,2) * t245 + Ifges(6,6) * t257;
t155 = mrSges(6,2) * t187 - mrSges(6,3) * t183 + Ifges(6,1) * t226 + Ifges(6,4) * t225 + Ifges(6,5) * t244 - pkin(11) * t167 - t268 * t168 + t272 * t169 + t245 * t220 - t257 * t221;
t237 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t270 + Ifges(5,6) * t274) * qJD(3);
t238 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t270 + Ifges(5,2) * t274) * qJD(3);
t141 = mrSges(5,2) * t196 - mrSges(5,3) * t189 + Ifges(5,1) * t249 + Ifges(5,4) * t250 + Ifges(5,5) * qJDD(4) - pkin(10) * t162 - qJD(4) * t238 - t269 * t154 + t273 * t155 + t237 * t292;
t239 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t270 + Ifges(5,4) * t274) * qJD(3);
t281 = -mrSges(7,1) * t179 + mrSges(7,2) * t180 - Ifges(7,5) * t205 - Ifges(7,6) * t204 - Ifges(7,3) * t240 - t228 * t207 + t227 * t208;
t278 = mrSges(6,1) * t183 - mrSges(6,2) * t184 + Ifges(6,5) * t226 + Ifges(6,6) * t225 + Ifges(6,3) * t244 + pkin(5) * t167 + t246 * t221 - t245 * t222 - t281;
t146 = -mrSges(5,1) * t196 + mrSges(5,3) * t190 + Ifges(5,4) * t249 + Ifges(5,2) * t250 + Ifges(5,6) * qJDD(4) - pkin(4) * t162 + qJD(4) * t239 - t237 * t293 - t278;
t130 = mrSges(4,1) * t198 - mrSges(4,2) * t199 + Ifges(4,3) * qJDD(3) + pkin(3) * t280 + pkin(9) * t289 + t270 * t141 + t274 * t146;
t139 = t150 * t295 + t266 * t152 + t262 * t296;
t131 = mrSges(4,2) * t211 - mrSges(4,3) * t198 + Ifges(4,5) * qJDD(3) - t277 * Ifges(4,6) - pkin(9) * t153 + t274 * t141 - t270 * t146;
t302 = mrSges(5,1) * t189 - mrSges(5,2) * t190 + Ifges(5,5) * t249 + Ifges(5,6) * t250 + Ifges(5,3) * qJDD(4) + pkin(4) * t177 + pkin(10) * t163 + t273 * t154 + t269 * t155 + (t270 * t238 - t274 * t239) * qJD(3);
t135 = -mrSges(4,1) * t211 + mrSges(4,3) * t199 + t277 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t153 - t302;
t284 = pkin(8) * t145 + t131 * t271 + t135 * t275;
t123 = -mrSges(3,1) * t234 + mrSges(3,3) * t219 - pkin(2) * t139 - t262 * t130 + t266 * t284;
t125 = mrSges(3,2) * t234 - mrSges(3,3) * t218 + t275 * t131 - t271 * t135 + (-t139 * t262 - t140 * t266) * pkin(8);
t134 = -t260 * t137 + t264 * t144;
t303 = qJ(2) * t134 + t123 * t264 + t125 * t260;
t138 = m(3) * t234 + t139;
t129 = -t263 * t138 + t304 * t267;
t121 = mrSges(3,1) * t218 - mrSges(3,2) * t219 + pkin(2) * t140 + t266 * t130 + t262 * t284;
t283 = mrSges(2,1) * t252 - mrSges(2,2) * t253 + pkin(1) * t129 + t267 * t121 + t303 * t263;
t132 = m(2) * t253 + t134;
t128 = t267 * t138 + t304 * t263;
t126 = m(2) * t252 + t129;
t119 = mrSges(2,2) * t259 - mrSges(2,3) * t252 - t260 * t123 + t264 * t125 + (-t128 * t263 - t129 * t267) * qJ(2);
t118 = -mrSges(2,1) * t259 + mrSges(2,3) * t253 - pkin(1) * t128 - t263 * t121 + t303 * t267;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t265 * t119 - t261 * t118 - qJ(1) * (t265 * t126 + t261 * t132), t119, t125, t131, t141, t155, t169; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t261 * t119 + t265 * t118 + qJ(1) * (-t261 * t126 + t265 * t132), t118, t123, t135, t146, t154, t168; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t283, t283, t121, t130, t302, t278, -t281;];
m_new  = t1;

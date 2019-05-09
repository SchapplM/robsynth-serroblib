% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-04 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:07:53
% EndTime: 2019-05-04 23:08:09
% DurationCPUTime: 9.31s
% Computational Cost: add. (166280->296), mult. (329524->371), div. (0->0), fcn. (217051->12), ass. (0->129)
t274 = sin(pkin(10));
t277 = cos(pkin(10));
t259 = t274 * g(1) - t277 * g(2);
t260 = -t277 * g(1) - t274 * g(2);
t270 = -g(3) + qJDD(1);
t275 = sin(pkin(6));
t278 = cos(pkin(6));
t281 = sin(qJ(2));
t284 = cos(qJ(2));
t214 = -t281 * t260 + (t259 * t278 + t270 * t275) * t284;
t286 = qJD(2) ^ 2;
t293 = -t286 * qJ(3) + qJDD(3) - t214;
t317 = -pkin(2) - pkin(8);
t208 = t317 * qJDD(2) + t293;
t232 = -t275 * t259 + t278 * t270;
t280 = sin(qJ(4));
t283 = cos(qJ(4));
t203 = t280 * t208 + t283 * t232;
t255 = (mrSges(5,1) * t280 + mrSges(5,2) * t283) * qJD(2);
t306 = qJD(2) * qJD(4);
t303 = t283 * t306;
t256 = t280 * qJDD(2) + t303;
t308 = qJD(2) * t283;
t262 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t308;
t254 = (pkin(4) * t280 - qJ(5) * t283) * qJD(2);
t285 = qJD(4) ^ 2;
t307 = t280 * qJD(2);
t188 = -t285 * pkin(4) + qJDD(4) * qJ(5) - t254 * t307 + t203;
t309 = t278 * t281;
t310 = t275 * t281;
t215 = t259 * t309 + t284 * t260 + t270 * t310;
t318 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t215;
t207 = t317 * t286 - t318;
t304 = t280 * t306;
t257 = t283 * qJDD(2) - t304;
t191 = (-t257 + t304) * qJ(5) + (t256 + t303) * pkin(4) + t207;
t273 = sin(pkin(11));
t276 = cos(pkin(11));
t247 = t273 * qJD(4) + t276 * t308;
t182 = -0.2e1 * qJD(5) * t247 - t273 * t188 + t276 * t191;
t230 = t273 * qJDD(4) + t276 * t257;
t246 = t276 * qJD(4) - t273 * t308;
t180 = (t246 * t307 - t230) * pkin(9) + (t246 * t247 + t256) * pkin(5) + t182;
t183 = 0.2e1 * qJD(5) * t246 + t276 * t188 + t273 * t191;
t229 = t276 * qJDD(4) - t273 * t257;
t231 = pkin(5) * t307 - t247 * pkin(9);
t245 = t246 ^ 2;
t181 = -t245 * pkin(5) + t229 * pkin(9) - t231 * t307 + t183;
t279 = sin(qJ(6));
t282 = cos(qJ(6));
t178 = t282 * t180 - t279 * t181;
t220 = t282 * t246 - t279 * t247;
t196 = t220 * qJD(6) + t279 * t229 + t282 * t230;
t221 = t279 * t246 + t282 * t247;
t204 = -t220 * mrSges(7,1) + t221 * mrSges(7,2);
t264 = qJD(6) + t307;
t212 = -t264 * mrSges(7,2) + t220 * mrSges(7,3);
t251 = qJDD(6) + t256;
t173 = m(7) * t178 + t251 * mrSges(7,1) - t196 * mrSges(7,3) - t221 * t204 + t264 * t212;
t179 = t279 * t180 + t282 * t181;
t195 = -t221 * qJD(6) + t282 * t229 - t279 * t230;
t213 = t264 * mrSges(7,1) - t221 * mrSges(7,3);
t174 = m(7) * t179 - t251 * mrSges(7,2) + t195 * mrSges(7,3) + t220 * t204 - t264 * t213;
t166 = t282 * t173 + t279 * t174;
t222 = -t246 * mrSges(6,1) + t247 * mrSges(6,2);
t227 = -mrSges(6,2) * t307 + t246 * mrSges(6,3);
t164 = m(6) * t182 + t256 * mrSges(6,1) - t230 * mrSges(6,3) - t247 * t222 + t227 * t307 + t166;
t228 = mrSges(6,1) * t307 - t247 * mrSges(6,3);
t301 = -t279 * t173 + t282 * t174;
t165 = m(6) * t183 - t256 * mrSges(6,2) + t229 * mrSges(6,3) + t246 * t222 - t228 * t307 + t301;
t302 = -t273 * t164 + t276 * t165;
t156 = m(5) * t203 - qJDD(4) * mrSges(5,2) - t256 * mrSges(5,3) - qJD(4) * t262 - t255 * t307 + t302;
t202 = t283 * t208 - t280 * t232;
t261 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t307;
t187 = -qJDD(4) * pkin(4) - t285 * qJ(5) + t254 * t308 + qJDD(5) - t202;
t184 = -t229 * pkin(5) - t245 * pkin(9) + t247 * t231 + t187;
t296 = m(7) * t184 - t195 * mrSges(7,1) + t196 * mrSges(7,2) - t220 * t212 + t221 * t213;
t288 = -m(6) * t187 + t229 * mrSges(6,1) - t230 * mrSges(6,2) + t246 * t227 - t247 * t228 - t296;
t169 = m(5) * t202 + qJDD(4) * mrSges(5,1) - t257 * mrSges(5,3) + qJD(4) * t261 - t255 * t308 + t288;
t147 = t283 * t156 - t280 * t169;
t145 = m(4) * t232 + t147;
t197 = Ifges(7,5) * t221 + Ifges(7,6) * t220 + Ifges(7,3) * t264;
t199 = Ifges(7,1) * t221 + Ifges(7,4) * t220 + Ifges(7,5) * t264;
t167 = -mrSges(7,1) * t184 + mrSges(7,3) * t179 + Ifges(7,4) * t196 + Ifges(7,2) * t195 + Ifges(7,6) * t251 - t221 * t197 + t264 * t199;
t198 = Ifges(7,4) * t221 + Ifges(7,2) * t220 + Ifges(7,6) * t264;
t168 = mrSges(7,2) * t184 - mrSges(7,3) * t178 + Ifges(7,1) * t196 + Ifges(7,4) * t195 + Ifges(7,5) * t251 + t220 * t197 - t264 * t198;
t216 = Ifges(6,5) * t247 + Ifges(6,6) * t246 + Ifges(6,3) * t307;
t218 = Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t307;
t149 = -mrSges(6,1) * t187 + mrSges(6,3) * t183 + Ifges(6,4) * t230 + Ifges(6,2) * t229 + Ifges(6,6) * t256 - pkin(5) * t296 + pkin(9) * t301 + t282 * t167 + t279 * t168 - t247 * t216 + t218 * t307;
t217 = Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t307;
t151 = mrSges(6,2) * t187 - mrSges(6,3) * t182 + Ifges(6,1) * t230 + Ifges(6,4) * t229 + Ifges(6,5) * t256 - pkin(9) * t166 - t279 * t167 + t282 * t168 + t246 * t216 - t217 * t307;
t159 = t276 * t164 + t273 * t165;
t239 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t283 - Ifges(5,6) * t280) * qJD(2);
t240 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t283 - Ifges(5,2) * t280) * qJD(2);
t137 = mrSges(5,2) * t207 - mrSges(5,3) * t202 + Ifges(5,1) * t257 - Ifges(5,4) * t256 + Ifges(5,5) * qJDD(4) - qJ(5) * t159 - qJD(4) * t240 - t273 * t149 + t276 * t151 - t239 * t307;
t241 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t283 - Ifges(5,4) * t280) * qJD(2);
t295 = -mrSges(7,1) * t178 + mrSges(7,2) * t179 - Ifges(7,5) * t196 - Ifges(7,6) * t195 - Ifges(7,3) * t251 - t221 * t198 + t220 * t199;
t287 = -mrSges(6,1) * t182 + mrSges(6,2) * t183 - Ifges(6,5) * t230 - Ifges(6,6) * t229 - pkin(5) * t166 - t247 * t217 + t246 * t218 + t295;
t141 = -pkin(4) * t159 + Ifges(5,6) * qJDD(4) - t239 * t308 + (-Ifges(5,2) - Ifges(6,3)) * t256 + mrSges(5,3) * t203 - mrSges(5,1) * t207 + qJD(4) * t241 + Ifges(5,4) * t257 + t287;
t157 = -m(5) * t207 - t256 * mrSges(5,1) - t257 * mrSges(5,2) - t261 * t307 - t262 * t308 - t159;
t209 = t286 * pkin(2) + t318;
t292 = -mrSges(4,1) * t209 - pkin(3) * t157 - pkin(8) * t147 - t280 * t137 - t283 * t141;
t313 = Ifges(4,5) - Ifges(3,6);
t314 = -Ifges(4,4) + Ifges(3,5);
t315 = mrSges(3,1) - mrSges(4,2);
t129 = mrSges(3,3) * t215 - pkin(2) * t145 - t313 * qJDD(2) - t315 * t232 + t314 * t286 + t292;
t146 = t280 * t156 + t283 * t169;
t211 = -qJDD(2) * pkin(2) + t293;
t298 = -m(4) * t211 + t286 * mrSges(4,3) - t146;
t143 = m(3) * t214 - t286 * mrSges(3,2) + t315 * qJDD(2) + t298;
t290 = -m(4) * t209 + t286 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t157;
t154 = m(3) * t215 - t286 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t290;
t140 = -t281 * t143 + t284 * t154;
t319 = pkin(7) * t140 + t129 * t284;
t311 = t143 * t284;
t144 = m(3) * t232 + t145;
t135 = -t275 * t144 + t154 * t309 + t278 * t311;
t294 = mrSges(4,2) * t211 - mrSges(4,3) * t209 + Ifges(4,1) * qJDD(2) - pkin(8) * t146 + t283 * t137 - t280 * t141;
t127 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t214 - mrSges(3,2) * t215 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t298) + qJ(3) * t290 + t294;
t291 = mrSges(5,1) * t202 - mrSges(5,2) * t203 + Ifges(5,5) * t257 - Ifges(5,6) * t256 + Ifges(5,3) * qJDD(4) + pkin(4) * t288 + qJ(5) * t302 + t276 * t149 + t273 * t151 + t240 * t308 + t241 * t307;
t289 = mrSges(4,1) * t211 + pkin(3) * t146 + t291;
t131 = -qJ(3) * t145 + (mrSges(3,2) - mrSges(4,3)) * t232 + t313 * t286 + t314 * qJDD(2) - mrSges(3,3) * t214 + t289;
t297 = mrSges(2,1) * t259 - mrSges(2,2) * t260 + pkin(1) * t135 + t278 * t127 + t131 * t310 + t319 * t275;
t138 = m(2) * t260 + t140;
t134 = t278 * t144 + (t154 * t281 + t311) * t275;
t132 = m(2) * t259 + t135;
t125 = mrSges(2,2) * t270 - mrSges(2,3) * t259 - t281 * t129 + t284 * t131 + (-t134 * t275 - t135 * t278) * pkin(7);
t124 = -mrSges(2,1) * t270 + mrSges(2,3) * t260 - pkin(1) * t134 - t275 * t127 + (t131 * t281 + t319) * t278;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t277 * t125 - t274 * t124 - qJ(1) * (t277 * t132 + t274 * t138), t125, t131, t294, t137, t151, t168; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t274 * t125 + t277 * t124 + qJ(1) * (-t274 * t132 + t277 * t138), t124, t129, mrSges(4,3) * t232 + Ifges(4,4) * qJDD(2) - t286 * Ifges(4,5) - t289, t141, t149, t167; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t297, t297, t127, -mrSges(4,2) * t232 + t286 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t292, t291, Ifges(6,3) * t256 - t287, -t295;];
m_new  = t1;

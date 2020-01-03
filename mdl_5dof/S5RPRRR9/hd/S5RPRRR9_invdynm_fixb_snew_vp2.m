% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:27
% EndTime: 2019-12-31 19:07:41
% DurationCPUTime: 9.17s
% Computational Cost: add. (141559->290), mult. (341067->365), div. (0->0), fcn. (257709->10), ass. (0->125)
t265 = qJD(1) ^ 2;
t260 = sin(qJ(1));
t264 = cos(qJ(1));
t240 = -g(1) * t264 - g(2) * t260;
t233 = -pkin(1) * t265 + qJDD(1) * qJ(2) + t240;
t255 = sin(pkin(9));
t256 = cos(pkin(9));
t289 = qJD(1) * qJD(2);
t287 = -t256 * g(3) - 0.2e1 * t255 * t289;
t294 = pkin(2) * t256;
t206 = (-pkin(6) * qJDD(1) + t265 * t294 - t233) * t255 + t287;
t224 = -t255 * g(3) + (t233 + 0.2e1 * t289) * t256;
t288 = qJDD(1) * t256;
t251 = t256 ^ 2;
t292 = t251 * t265;
t207 = -pkin(2) * t292 + pkin(6) * t288 + t224;
t259 = sin(qJ(3));
t263 = cos(qJ(3));
t187 = t263 * t206 - t259 * t207;
t277 = t255 * t263 + t256 * t259;
t276 = -t255 * t259 + t256 * t263;
t231 = t276 * qJD(1);
t290 = t231 * qJD(3);
t222 = qJDD(1) * t277 + t290;
t232 = t277 * qJD(1);
t165 = (-t222 + t290) * pkin(7) + (t231 * t232 + qJDD(3)) * pkin(3) + t187;
t188 = t259 * t206 + t263 * t207;
t221 = -t232 * qJD(3) + qJDD(1) * t276;
t227 = qJD(3) * pkin(3) - pkin(7) * t232;
t230 = t231 ^ 2;
t170 = -pkin(3) * t230 + pkin(7) * t221 - qJD(3) * t227 + t188;
t258 = sin(qJ(4));
t262 = cos(qJ(4));
t163 = t165 * t258 + t170 * t262;
t213 = t231 * t258 + t232 * t262;
t183 = -qJD(4) * t213 + t221 * t262 - t222 * t258;
t212 = t231 * t262 - t232 * t258;
t196 = -mrSges(5,1) * t212 + mrSges(5,2) * t213;
t252 = qJD(3) + qJD(4);
t204 = mrSges(5,1) * t252 - mrSges(5,3) * t213;
t249 = qJDD(3) + qJDD(4);
t197 = -pkin(4) * t212 - pkin(8) * t213;
t248 = t252 ^ 2;
t159 = -pkin(4) * t248 + pkin(8) * t249 + t197 * t212 + t163;
t250 = t255 ^ 2;
t239 = t260 * g(1) - t264 * g(2);
t282 = qJDD(2) - t239;
t220 = (-pkin(1) - t294) * qJDD(1) + (-qJ(2) + (-t250 - t251) * pkin(6)) * t265 + t282;
t177 = -t221 * pkin(3) - t230 * pkin(7) + t232 * t227 + t220;
t184 = qJD(4) * t212 + t221 * t258 + t222 * t262;
t160 = (-t212 * t252 - t184) * pkin(8) + (t213 * t252 - t183) * pkin(4) + t177;
t257 = sin(qJ(5));
t261 = cos(qJ(5));
t156 = -t159 * t257 + t160 * t261;
t199 = -t213 * t257 + t252 * t261;
t168 = qJD(5) * t199 + t184 * t261 + t249 * t257;
t182 = qJDD(5) - t183;
t200 = t213 * t261 + t252 * t257;
t185 = -mrSges(6,1) * t199 + mrSges(6,2) * t200;
t208 = qJD(5) - t212;
t189 = -mrSges(6,2) * t208 + mrSges(6,3) * t199;
t152 = m(6) * t156 + mrSges(6,1) * t182 - mrSges(6,3) * t168 - t185 * t200 + t189 * t208;
t157 = t159 * t261 + t160 * t257;
t167 = -qJD(5) * t200 - t184 * t257 + t249 * t261;
t190 = mrSges(6,1) * t208 - mrSges(6,3) * t200;
t153 = m(6) * t157 - mrSges(6,2) * t182 + mrSges(6,3) * t167 + t185 * t199 - t190 * t208;
t283 = -t152 * t257 + t153 * t261;
t139 = m(5) * t163 - mrSges(5,2) * t249 + mrSges(5,3) * t183 + t196 * t212 - t204 * t252 + t283;
t162 = t165 * t262 - t170 * t258;
t203 = -mrSges(5,2) * t252 + mrSges(5,3) * t212;
t158 = -pkin(4) * t249 - pkin(8) * t248 + t197 * t213 - t162;
t272 = -m(6) * t158 + mrSges(6,1) * t167 - mrSges(6,2) * t168 + t189 * t199 - t190 * t200;
t148 = m(5) * t162 + mrSges(5,1) * t249 - mrSges(5,3) * t184 - t196 * t213 + t203 * t252 + t272;
t134 = t139 * t258 + t148 * t262;
t216 = -mrSges(4,1) * t231 + mrSges(4,2) * t232;
t225 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t231;
t131 = m(4) * t187 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t222 + qJD(3) * t225 - t216 * t232 + t134;
t226 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t232;
t284 = t139 * t262 - t148 * t258;
t132 = m(4) * t188 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t221 - qJD(3) * t226 + t216 * t231 + t284;
t125 = t131 * t263 + t132 * t259;
t223 = -t255 * t233 + t287;
t210 = Ifges(4,4) * t232 + Ifges(4,2) * t231 + Ifges(4,6) * qJD(3);
t211 = Ifges(4,1) * t232 + Ifges(4,4) * t231 + Ifges(4,5) * qJD(3);
t171 = Ifges(6,5) * t200 + Ifges(6,6) * t199 + Ifges(6,3) * t208;
t173 = Ifges(6,1) * t200 + Ifges(6,4) * t199 + Ifges(6,5) * t208;
t145 = -mrSges(6,1) * t158 + mrSges(6,3) * t157 + Ifges(6,4) * t168 + Ifges(6,2) * t167 + Ifges(6,6) * t182 - t171 * t200 + t173 * t208;
t172 = Ifges(6,4) * t200 + Ifges(6,2) * t199 + Ifges(6,6) * t208;
t146 = mrSges(6,2) * t158 - mrSges(6,3) * t156 + Ifges(6,1) * t168 + Ifges(6,4) * t167 + Ifges(6,5) * t182 + t171 * t199 - t172 * t208;
t192 = Ifges(5,4) * t213 + Ifges(5,2) * t212 + Ifges(5,6) * t252;
t193 = Ifges(5,1) * t213 + Ifges(5,4) * t212 + Ifges(5,5) * t252;
t271 = -mrSges(5,1) * t162 + mrSges(5,2) * t163 - Ifges(5,5) * t184 - Ifges(5,6) * t183 - Ifges(5,3) * t249 - pkin(4) * t272 - pkin(8) * t283 - t145 * t261 - t146 * t257 - t192 * t213 + t212 * t193;
t267 = -mrSges(4,1) * t187 + mrSges(4,2) * t188 - Ifges(4,5) * t222 - Ifges(4,6) * t221 - Ifges(4,3) * qJDD(3) - pkin(3) * t134 - t232 * t210 + t231 * t211 + t271;
t280 = Ifges(3,4) * t255 + Ifges(3,2) * t256;
t281 = Ifges(3,1) * t255 + Ifges(3,4) * t256;
t295 = -mrSges(3,1) * t223 + mrSges(3,2) * t224 - pkin(2) * t125 - (t255 * t280 - t256 * t281) * t265 + t267;
t293 = mrSges(3,2) * t255;
t141 = t152 * t261 + t153 * t257;
t279 = Ifges(3,5) * t255 + Ifges(3,6) * t256;
t291 = t265 * t279;
t275 = mrSges(3,3) * qJDD(1) + t265 * (-mrSges(3,1) * t256 + t293);
t123 = m(3) * t223 - t255 * t275 + t125;
t285 = -t259 * t131 + t132 * t263;
t124 = m(3) * t224 + t256 * t275 + t285;
t286 = -t123 * t255 + t124 * t256;
t274 = m(5) * t177 - mrSges(5,1) * t183 + mrSges(5,2) * t184 - t212 * t203 + t213 * t204 + t141;
t191 = Ifges(5,5) * t213 + Ifges(5,6) * t212 + Ifges(5,3) * t252;
t126 = mrSges(5,2) * t177 - mrSges(5,3) * t162 + Ifges(5,1) * t184 + Ifges(5,4) * t183 + Ifges(5,5) * t249 - pkin(8) * t141 - t145 * t257 + t146 * t261 + t191 * t212 - t192 * t252;
t269 = mrSges(6,1) * t156 - mrSges(6,2) * t157 + Ifges(6,5) * t168 + Ifges(6,6) * t167 + Ifges(6,3) * t182 + t172 * t200 - t173 * t199;
t127 = -mrSges(5,1) * t177 + mrSges(5,3) * t163 + Ifges(5,4) * t184 + Ifges(5,2) * t183 + Ifges(5,6) * t249 - pkin(4) * t141 - t191 * t213 + t193 * t252 - t269;
t209 = Ifges(4,5) * t232 + Ifges(4,6) * t231 + Ifges(4,3) * qJD(3);
t117 = -mrSges(4,1) * t220 + mrSges(4,3) * t188 + Ifges(4,4) * t222 + Ifges(4,2) * t221 + Ifges(4,6) * qJDD(3) - pkin(3) * t274 + pkin(7) * t284 + qJD(3) * t211 + t258 * t126 + t262 * t127 - t232 * t209;
t121 = mrSges(4,2) * t220 - mrSges(4,3) * t187 + Ifges(4,1) * t222 + Ifges(4,4) * t221 + Ifges(4,5) * qJDD(3) - pkin(7) * t134 - qJD(3) * t210 + t126 * t262 - t127 * t258 + t209 * t231;
t229 = -qJDD(1) * pkin(1) - t265 * qJ(2) + t282;
t270 = m(4) * t220 - t221 * mrSges(4,1) + t222 * mrSges(4,2) - t231 * t225 + t232 * t226 + t274;
t113 = -mrSges(3,1) * t229 + mrSges(3,3) * t224 - pkin(2) * t270 + pkin(6) * t285 + qJDD(1) * t280 + t263 * t117 + t259 * t121 - t255 * t291;
t115 = mrSges(3,2) * t229 - mrSges(3,3) * t223 - pkin(6) * t125 + qJDD(1) * t281 - t259 * t117 + t263 * t121 + t256 * t291;
t268 = -m(3) * t229 + mrSges(3,1) * t288 - t270 + (t250 * t265 + t292) * mrSges(3,3);
t273 = -mrSges(2,2) * t240 + qJ(2) * t286 + t256 * t113 + t255 * t115 + pkin(1) * (-qJDD(1) * t293 + t268) + mrSges(2,1) * t239 + Ifges(2,3) * qJDD(1);
t135 = -t265 * mrSges(2,2) + t268 + m(2) * t239 + (mrSges(2,1) - t293) * qJDD(1);
t120 = t123 * t256 + t124 * t255;
t118 = m(2) * t240 - mrSges(2,1) * t265 - qJDD(1) * mrSges(2,2) + t286;
t116 = t265 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t240 - pkin(1) * t120 + (Ifges(2,6) - t279) * qJDD(1) + t295;
t111 = -mrSges(2,2) * g(3) - mrSges(2,3) * t239 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t265 - qJ(2) * t120 - t113 * t255 + t115 * t256;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t264 * t111 - t260 * t116 - pkin(5) * (t118 * t260 + t135 * t264), t111, t115, t121, t126, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t260 * t111 + t264 * t116 + pkin(5) * (t118 * t264 - t135 * t260), t116, t113, t117, t127, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t273, t273, qJDD(1) * t279 - t295, -t267, -t271, t269;];
m_new = t1;

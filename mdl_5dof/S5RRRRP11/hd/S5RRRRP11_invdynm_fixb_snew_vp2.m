% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:37
% EndTime: 2019-12-31 22:14:56
% DurationCPUTime: 8.76s
% Computational Cost: add. (147184->320), mult. (314466->408), div. (0->0), fcn. (239003->10), ass. (0->127)
t261 = sin(pkin(5));
t265 = sin(qJ(2));
t268 = cos(qJ(2));
t284 = qJD(1) * qJD(2);
t248 = (-qJDD(1) * t268 + t265 * t284) * t261;
t286 = qJD(1) * t261;
t246 = (-pkin(2) * t268 - pkin(8) * t265) * t286;
t262 = cos(pkin(5));
t257 = qJD(1) * t262 + qJD(2);
t255 = t257 ^ 2;
t256 = qJDD(1) * t262 + qJDD(2);
t285 = qJD(1) * t268;
t266 = sin(qJ(1));
t269 = cos(qJ(1));
t253 = t266 * g(1) - g(2) * t269;
t270 = qJD(1) ^ 2;
t297 = pkin(7) * t261;
t243 = qJDD(1) * pkin(1) + t270 * t297 + t253;
t254 = -g(1) * t269 - g(2) * t266;
t244 = -pkin(1) * t270 + qJDD(1) * t297 + t254;
t292 = t262 * t265;
t287 = t243 * t292 + t268 * t244;
t195 = -t255 * pkin(2) + t256 * pkin(8) + (-g(3) * t265 + t246 * t285) * t261 + t287;
t247 = (qJDD(1) * t265 + t268 * t284) * t261;
t296 = t262 * g(3);
t196 = t248 * pkin(2) - t247 * pkin(8) - t296 + (-t243 + (pkin(2) * t265 - pkin(8) * t268) * t257 * qJD(1)) * t261;
t264 = sin(qJ(3));
t267 = cos(qJ(3));
t172 = t267 * t195 + t264 * t196;
t282 = t265 * t286;
t236 = t257 * t267 - t264 * t282;
t237 = t257 * t264 + t267 * t282;
t222 = -pkin(3) * t236 - pkin(9) * t237;
t240 = qJDD(3) + t248;
t281 = t261 * t285;
t252 = qJD(3) - t281;
t251 = t252 ^ 2;
t168 = -pkin(3) * t251 + pkin(9) * t240 + t222 * t236 + t172;
t291 = t262 * t268;
t293 = t261 * t268;
t219 = -g(3) * t293 + t243 * t291 - t265 * t244;
t194 = -t256 * pkin(2) - t255 * pkin(8) + t246 * t282 - t219;
t217 = -qJD(3) * t237 - t247 * t264 + t256 * t267;
t218 = qJD(3) * t236 + t247 * t267 + t256 * t264;
t170 = (-t236 * t252 - t218) * pkin(9) + (t237 * t252 - t217) * pkin(3) + t194;
t263 = sin(qJ(4));
t298 = cos(qJ(4));
t164 = -t168 * t263 + t170 * t298;
t165 = t168 * t298 + t170 * t263;
t224 = t237 * t298 + t252 * t263;
t181 = qJD(4) * t224 + t218 * t263 - t240 * t298;
t223 = t237 * t263 - t252 * t298;
t182 = -qJD(4) * t223 + t218 * t298 + t240 * t263;
t234 = qJD(4) - t236;
t183 = Ifges(6,5) * t224 + Ifges(6,6) * t234 + Ifges(6,3) * t223;
t186 = Ifges(5,4) * t224 - Ifges(5,2) * t223 + Ifges(5,6) * t234;
t188 = Ifges(5,1) * t224 - Ifges(5,4) * t223 + Ifges(5,5) * t234;
t200 = mrSges(6,1) * t223 - mrSges(6,3) * t224;
t215 = qJDD(4) - t217;
t199 = pkin(4) * t223 - qJ(5) * t224;
t233 = t234 ^ 2;
t160 = -pkin(4) * t233 + qJ(5) * t215 + 0.2e1 * qJD(5) * t234 - t199 * t223 + t165;
t162 = -pkin(4) * t215 - qJ(5) * t233 + t199 * t224 + qJDD(5) - t164;
t187 = Ifges(6,1) * t224 + Ifges(6,4) * t234 + Ifges(6,5) * t223;
t276 = mrSges(6,1) * t162 - mrSges(6,3) * t160 - Ifges(6,4) * t182 - Ifges(6,2) * t215 - Ifges(6,6) * t181 - t187 * t223;
t203 = -mrSges(6,2) * t223 + mrSges(6,3) * t234;
t279 = -m(6) * t162 + t215 * mrSges(6,1) + t234 * t203;
t206 = -mrSges(6,1) * t234 + mrSges(6,2) * t224;
t283 = m(6) * t160 + t215 * mrSges(6,3) + t234 * t206;
t299 = -(-t186 + t183) * t224 + mrSges(5,1) * t164 - mrSges(5,2) * t165 + Ifges(5,5) * t182 - Ifges(5,6) * t181 + Ifges(5,3) * t215 + pkin(4) * (-t182 * mrSges(6,2) - t224 * t200 + t279) + qJ(5) * (-t181 * mrSges(6,2) - t223 * t200 + t283) + t223 * t188 - t276;
t295 = -mrSges(5,3) - mrSges(6,2);
t294 = t261 * t265;
t205 = mrSges(5,1) * t234 - mrSges(5,3) * t224;
t288 = -mrSges(5,1) * t223 - mrSges(5,2) * t224 - t200;
t152 = m(5) * t165 - t215 * mrSges(5,2) + t181 * t295 - t234 * t205 + t223 * t288 + t283;
t204 = -mrSges(5,2) * t234 - mrSges(5,3) * t223;
t153 = m(5) * t164 + t215 * mrSges(5,1) + t182 * t295 + t234 * t204 + t224 * t288 + t279;
t148 = t152 * t298 - t153 * t263;
t221 = -mrSges(4,1) * t236 + mrSges(4,2) * t237;
t226 = mrSges(4,1) * t252 - mrSges(4,3) * t237;
t145 = m(4) * t172 - mrSges(4,2) * t240 + mrSges(4,3) * t217 + t221 * t236 - t226 * t252 + t148;
t171 = -t264 * t195 + t267 * t196;
t167 = -t240 * pkin(3) - t251 * pkin(9) + t237 * t222 - t171;
t163 = -0.2e1 * qJD(5) * t224 + (t223 * t234 - t182) * qJ(5) + (t224 * t234 + t181) * pkin(4) + t167;
t157 = m(6) * t163 + mrSges(6,1) * t181 - mrSges(6,3) * t182 + t203 * t223 - t224 * t206;
t154 = -m(5) * t167 - mrSges(5,1) * t181 - mrSges(5,2) * t182 - t223 * t204 - t205 * t224 - t157;
t225 = -mrSges(4,2) * t252 + mrSges(4,3) * t236;
t150 = m(4) * t171 + mrSges(4,1) * t240 - mrSges(4,3) * t218 - t221 * t237 + t225 * t252 + t154;
t139 = t145 * t264 + t150 * t267;
t185 = Ifges(6,4) * t224 + Ifges(6,2) * t234 + Ifges(6,6) * t223;
t290 = -Ifges(5,5) * t224 + Ifges(5,6) * t223 - Ifges(5,3) * t234 - t185;
t220 = -g(3) * t294 + t287;
t241 = mrSges(3,1) * t257 - mrSges(3,3) * t282;
t245 = (-mrSges(3,1) * t268 + mrSges(3,2) * t265) * t286;
t280 = t145 * t267 - t150 * t264;
t137 = m(3) * t220 - mrSges(3,2) * t256 - mrSges(3,3) * t248 - t241 * t257 + t245 * t281 + t280;
t242 = -mrSges(3,2) * t257 + mrSges(3,3) * t281;
t147 = t152 * t263 + t153 * t298;
t273 = -m(4) * t194 + t217 * mrSges(4,1) - mrSges(4,2) * t218 + t236 * t225 - t226 * t237 - t147;
t141 = m(3) * t219 + mrSges(3,1) * t256 - mrSges(3,3) * t247 + t242 * t257 - t245 * t282 + t273;
t133 = t137 * t268 - t141 * t265;
t230 = -t261 * t243 - t296;
t138 = m(3) * t230 + t248 * mrSges(3,1) + t247 * mrSges(3,2) + (t241 * t265 - t242 * t268) * t286 + t139;
t129 = t137 * t292 - t138 * t261 + t141 * t291;
t278 = -mrSges(6,1) * t163 + mrSges(6,2) * t160;
t275 = mrSges(6,2) * t162 - mrSges(6,3) * t163 + Ifges(6,1) * t182 + Ifges(6,4) * t215 + Ifges(6,5) * t181 + t183 * t234;
t142 = -mrSges(5,1) * t167 + mrSges(5,3) * t165 - pkin(4) * t157 + (t187 + t188) * t234 + t290 * t224 + (Ifges(5,6) - Ifges(6,6)) * t215 + (Ifges(5,4) - Ifges(6,5)) * t182 + (-Ifges(5,2) - Ifges(6,3)) * t181 + t278;
t146 = mrSges(5,2) * t167 - mrSges(5,3) * t164 + Ifges(5,1) * t182 - Ifges(5,4) * t181 + Ifges(5,5) * t215 - qJ(5) * t157 - t234 * t186 + t223 * t290 + t275;
t211 = Ifges(4,5) * t237 + Ifges(4,6) * t236 + Ifges(4,3) * t252;
t212 = Ifges(4,4) * t237 + Ifges(4,2) * t236 + Ifges(4,6) * t252;
t130 = mrSges(4,2) * t194 - mrSges(4,3) * t171 + Ifges(4,1) * t218 + Ifges(4,4) * t217 + Ifges(4,5) * t240 - pkin(9) * t147 - t142 * t263 + t146 * t298 + t211 * t236 - t212 * t252;
t213 = Ifges(4,1) * t237 + Ifges(4,4) * t236 + Ifges(4,5) * t252;
t134 = -mrSges(4,1) * t194 + mrSges(4,3) * t172 + Ifges(4,4) * t218 + Ifges(4,2) * t217 + Ifges(4,6) * t240 - pkin(3) * t147 - t237 * t211 + t252 * t213 - t299;
t228 = Ifges(3,6) * t257 + (Ifges(3,4) * t265 + Ifges(3,2) * t268) * t286;
t229 = Ifges(3,5) * t257 + (Ifges(3,1) * t265 + Ifges(3,4) * t268) * t286;
t121 = Ifges(3,5) * t247 - Ifges(3,6) * t248 + Ifges(3,3) * t256 + mrSges(3,1) * t219 - mrSges(3,2) * t220 + t264 * t130 + t267 * t134 + pkin(2) * t273 + pkin(8) * t280 + (t228 * t265 - t229 * t268) * t286;
t227 = Ifges(3,3) * t257 + (Ifges(3,5) * t265 + Ifges(3,6) * t268) * t286;
t123 = mrSges(3,2) * t230 - mrSges(3,3) * t219 + Ifges(3,1) * t247 - Ifges(3,4) * t248 + Ifges(3,5) * t256 - pkin(8) * t139 + t130 * t267 - t134 * t264 + t227 * t281 - t228 * t257;
t271 = mrSges(4,1) * t171 - mrSges(4,2) * t172 + Ifges(4,5) * t218 + Ifges(4,6) * t217 + Ifges(4,3) * t240 + pkin(3) * t154 + pkin(9) * t148 + t142 * t298 + t146 * t263 + t212 * t237 - t213 * t236;
t125 = -mrSges(3,1) * t230 + mrSges(3,3) * t220 + Ifges(3,4) * t247 - Ifges(3,2) * t248 + Ifges(3,6) * t256 - pkin(2) * t139 - t227 * t282 + t229 * t257 - t271;
t274 = mrSges(2,1) * t253 - mrSges(2,2) * t254 + Ifges(2,3) * qJDD(1) + pkin(1) * t129 + t121 * t262 + t123 * t294 + t125 * t293 + t133 * t297;
t131 = m(2) * t254 - mrSges(2,1) * t270 - qJDD(1) * mrSges(2,2) + t133;
t128 = t262 * t138 + (t137 * t265 + t141 * t268) * t261;
t126 = m(2) * t253 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t270 + t129;
t119 = -mrSges(2,2) * g(3) - mrSges(2,3) * t253 + Ifges(2,5) * qJDD(1) - t270 * Ifges(2,6) + t268 * t123 - t265 * t125 + (-t128 * t261 - t129 * t262) * pkin(7);
t118 = mrSges(2,1) * g(3) + mrSges(2,3) * t254 + t270 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t128 - t261 * t121 + (pkin(7) * t133 + t123 * t265 + t125 * t268) * t262;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t119 - t266 * t118 - pkin(6) * (t126 * t269 + t131 * t266), t119, t123, t130, t146, -t185 * t223 + t275; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t266 * t119 + t269 * t118 + pkin(6) * (-t126 * t266 + t131 * t269), t118, t125, t134, t142, -t224 * t183 - t276; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t274, t274, t121, t271, t299, Ifges(6,5) * t182 + Ifges(6,6) * t215 + Ifges(6,3) * t181 + t224 * t185 - t234 * t187 - t278;];
m_new = t1;

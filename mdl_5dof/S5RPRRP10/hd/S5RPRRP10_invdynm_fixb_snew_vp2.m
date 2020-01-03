% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:51:00
% EndTime: 2019-12-31 18:51:08
% DurationCPUTime: 5.04s
% Computational Cost: add. (56619->285), mult. (133136->346), div. (0->0), fcn. (93221->8), ass. (0->116)
t261 = qJD(1) ^ 2;
t252 = sin(pkin(8));
t253 = cos(pkin(8));
t255 = sin(qJ(3));
t258 = cos(qJ(3));
t273 = t252 * t255 - t253 * t258;
t232 = t273 * qJD(1);
t274 = t252 * t258 + t253 * t255;
t290 = t232 * qJD(3);
t221 = t274 * qJDD(1) - t290;
t233 = t274 * qJD(1);
t254 = sin(qJ(4));
t257 = cos(qJ(4));
t225 = qJD(3) * t257 - t233 * t254;
t189 = qJD(4) * t225 + qJDD(3) * t254 + t221 * t257;
t226 = qJD(3) * t254 + t233 * t257;
t192 = -mrSges(6,1) * t225 + mrSges(6,2) * t226;
t256 = sin(qJ(1));
t259 = cos(qJ(1));
t241 = -g(1) * t259 - g(2) * t256;
t234 = -pkin(1) * t261 + qJDD(1) * qJ(2) + t241;
t288 = qJD(1) * qJD(2);
t284 = -t253 * g(3) - 0.2e1 * t252 * t288;
t297 = pkin(2) * t253;
t203 = (-pkin(6) * qJDD(1) + t261 * t297 - t234) * t252 + t284;
t223 = -g(3) * t252 + (t234 + 0.2e1 * t288) * t253;
t287 = qJDD(1) * t253;
t249 = t253 ^ 2;
t294 = t249 * t261;
t204 = -pkin(2) * t294 + pkin(6) * t287 + t223;
t168 = t255 * t203 + t258 * t204;
t218 = pkin(3) * t232 - pkin(7) * t233;
t260 = qJD(3) ^ 2;
t160 = -pkin(3) * t260 + qJDD(3) * pkin(7) - t218 * t232 + t168;
t248 = t252 ^ 2;
t240 = t256 * g(1) - t259 * g(2);
t279 = qJDD(2) - t240;
t219 = (-pkin(1) - t297) * qJDD(1) + (-qJ(2) + (-t248 - t249) * pkin(6)) * t261 + t279;
t289 = t233 * qJD(3);
t220 = -t273 * qJDD(1) - t289;
t163 = (-t221 + t290) * pkin(7) + (-t220 + t289) * pkin(3) + t219;
t156 = -t254 * t160 + t257 * t163;
t217 = qJDD(4) - t220;
t230 = qJD(4) + t232;
t150 = -0.2e1 * qJD(5) * t226 + (t225 * t230 - t189) * qJ(5) + (t225 * t226 + t217) * pkin(4) + t156;
t196 = -mrSges(6,2) * t230 + mrSges(6,3) * t225;
t286 = m(6) * t150 + t217 * mrSges(6,1) + t230 * t196;
t147 = -mrSges(6,3) * t189 - t192 * t226 + t286;
t157 = t257 * t160 + t254 * t163;
t174 = Ifges(5,4) * t226 + Ifges(5,2) * t225 + Ifges(5,6) * t230;
t175 = Ifges(6,1) * t226 + Ifges(6,4) * t225 + Ifges(6,5) * t230;
t176 = Ifges(5,1) * t226 + Ifges(5,4) * t225 + Ifges(5,5) * t230;
t188 = -qJD(4) * t226 + qJDD(3) * t257 - t221 * t254;
t198 = pkin(4) * t230 - qJ(5) * t226;
t224 = t225 ^ 2;
t153 = -pkin(4) * t224 + qJ(5) * t188 + 0.2e1 * qJD(5) * t225 - t198 * t230 + t157;
t173 = Ifges(6,4) * t226 + Ifges(6,2) * t225 + Ifges(6,6) * t230;
t270 = -mrSges(6,1) * t150 + mrSges(6,2) * t153 - Ifges(6,5) * t189 - Ifges(6,6) * t188 - Ifges(6,3) * t217 - t226 * t173;
t299 = mrSges(5,1) * t156 - mrSges(5,2) * t157 + Ifges(5,5) * t189 + Ifges(5,6) * t188 + Ifges(5,3) * t217 + pkin(4) * t147 + t226 * t174 - (t176 + t175) * t225 - t270;
t213 = mrSges(4,1) * t232 + mrSges(4,2) * t233;
t228 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t233;
t193 = -mrSges(5,1) * t225 + mrSges(5,2) * t226;
t197 = -mrSges(5,2) * t230 + mrSges(5,3) * t225;
t139 = m(5) * t156 + mrSges(5,1) * t217 + t197 * t230 + (-t192 - t193) * t226 + (-mrSges(5,3) - mrSges(6,3)) * t189 + t286;
t285 = m(6) * t153 + t188 * mrSges(6,3) + t225 * t192;
t199 = mrSges(6,1) * t230 - mrSges(6,3) * t226;
t292 = -mrSges(5,1) * t230 + mrSges(5,3) * t226 - t199;
t296 = -mrSges(5,2) - mrSges(6,2);
t142 = m(5) * t157 + mrSges(5,3) * t188 + t193 * t225 + t296 * t217 + t292 * t230 + t285;
t281 = -t139 * t254 + t257 * t142;
t132 = m(4) * t168 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t220 - qJD(3) * t228 - t213 * t232 + t281;
t167 = t203 * t258 - t255 * t204;
t227 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t232;
t159 = -qJDD(3) * pkin(3) - pkin(7) * t260 + t233 * t218 - t167;
t155 = -pkin(4) * t188 - qJ(5) * t224 + t198 * t226 + qJDD(5) + t159;
t280 = -m(6) * t155 + t188 * mrSges(6,1) + t225 * t196;
t265 = -m(5) * t159 + t188 * mrSges(5,1) + t296 * t189 + t225 * t197 + t292 * t226 + t280;
t144 = m(4) * t167 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t221 + qJD(3) * t227 - t213 * t233 + t265;
t125 = t255 * t132 + t258 * t144;
t222 = -t252 * t234 + t284;
t171 = Ifges(6,5) * t226 + Ifges(6,6) * t225 + Ifges(6,3) * t230;
t172 = Ifges(5,5) * t226 + Ifges(5,6) * t225 + Ifges(5,3) * t230;
t271 = -mrSges(6,1) * t155 + mrSges(6,3) * t153 + Ifges(6,4) * t189 + Ifges(6,2) * t188 + Ifges(6,6) * t217 + t230 * t175;
t127 = Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * t217 + t230 * t176 - mrSges(5,1) * t159 + mrSges(5,3) * t157 - pkin(4) * (mrSges(6,2) * t189 - t280) + qJ(5) * (-mrSges(6,2) * t217 - t199 * t230 + t285) + (-pkin(4) * t199 - t171 - t172) * t226 + t271;
t269 = mrSges(6,2) * t155 - mrSges(6,3) * t150 + Ifges(6,1) * t189 + Ifges(6,4) * t188 + Ifges(6,5) * t217 + t225 * t171;
t134 = mrSges(5,2) * t159 - mrSges(5,3) * t156 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * t217 - qJ(5) * t147 + t172 * t225 + (-t173 - t174) * t230 + t269;
t206 = Ifges(4,4) * t233 - Ifges(4,2) * t232 + Ifges(4,6) * qJD(3);
t207 = Ifges(4,1) * t233 - Ifges(4,4) * t232 + Ifges(4,5) * qJD(3);
t266 = -mrSges(4,1) * t167 + mrSges(4,2) * t168 - Ifges(4,5) * t221 - Ifges(4,6) * t220 - Ifges(4,3) * qJDD(3) - pkin(3) * t265 - pkin(7) * t281 - t257 * t127 - t254 * t134 - t233 * t206 - t232 * t207;
t277 = Ifges(3,4) * t252 + Ifges(3,2) * t253;
t278 = Ifges(3,1) * t252 + Ifges(3,4) * t253;
t298 = -mrSges(3,1) * t222 + mrSges(3,2) * t223 - pkin(2) * t125 - (t252 * t277 - t253 * t278) * t261 + t266;
t295 = mrSges(3,2) * t252;
t136 = t257 * t139 + t254 * t142;
t276 = Ifges(3,5) * t252 + Ifges(3,6) * t253;
t291 = t261 * t276;
t272 = mrSges(3,3) * qJDD(1) + t261 * (-mrSges(3,1) * t253 + t295);
t123 = m(3) * t222 - t272 * t252 + t125;
t282 = t258 * t132 - t255 * t144;
t124 = m(3) * t223 + t272 * t253 + t282;
t283 = -t123 * t252 + t253 * t124;
t205 = Ifges(4,5) * t233 - Ifges(4,6) * t232 + Ifges(4,3) * qJD(3);
t117 = mrSges(4,2) * t219 - mrSges(4,3) * t167 + Ifges(4,1) * t221 + Ifges(4,4) * t220 + Ifges(4,5) * qJDD(3) - pkin(7) * t136 - qJD(3) * t206 - t127 * t254 + t134 * t257 - t205 * t232;
t121 = -mrSges(4,1) * t219 + mrSges(4,3) * t168 + Ifges(4,4) * t221 + Ifges(4,2) * t220 + Ifges(4,6) * qJDD(3) - pkin(3) * t136 + qJD(3) * t207 - t233 * t205 - t299;
t231 = -qJDD(1) * pkin(1) - t261 * qJ(2) + t279;
t267 = m(4) * t219 - t220 * mrSges(4,1) + mrSges(4,2) * t221 + t232 * t227 + t228 * t233 + t136;
t114 = -mrSges(3,1) * t231 + mrSges(3,3) * t223 - pkin(2) * t267 + pkin(6) * t282 + t277 * qJDD(1) + t255 * t117 + t258 * t121 - t252 * t291;
t116 = mrSges(3,2) * t231 - mrSges(3,3) * t222 - pkin(6) * t125 + t278 * qJDD(1) + t117 * t258 - t121 * t255 + t253 * t291;
t264 = -m(3) * t231 + mrSges(3,1) * t287 - t267 + (t248 * t261 + t294) * mrSges(3,3);
t268 = -mrSges(2,2) * t241 + qJ(2) * t283 + t253 * t114 + t252 * t116 + pkin(1) * (-qJDD(1) * t295 + t264) + mrSges(2,1) * t240 + Ifges(2,3) * qJDD(1);
t128 = t264 + m(2) * t240 - mrSges(2,2) * t261 + (mrSges(2,1) - t295) * qJDD(1);
t120 = t123 * t253 + t124 * t252;
t118 = m(2) * t241 - mrSges(2,1) * t261 - qJDD(1) * mrSges(2,2) + t283;
t112 = mrSges(2,1) * g(3) + (Ifges(2,6) - t276) * qJDD(1) + t261 * Ifges(2,5) + mrSges(2,3) * t241 - pkin(1) * t120 + t298;
t111 = -mrSges(2,2) * g(3) - mrSges(2,3) * t240 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t261 - qJ(2) * t120 - t114 * t252 + t116 * t253;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t259 * t111 - t256 * t112 - pkin(5) * (t118 * t256 + t128 * t259), t111, t116, t117, t134, -t173 * t230 + t269; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t256 * t111 + t259 * t112 + pkin(5) * (t118 * t259 - t128 * t256), t112, t114, t121, t127, -t226 * t171 + t271; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t268, t268, t276 * qJDD(1) - t298, -t266, t299, -t225 * t175 - t270;];
m_new = t1;

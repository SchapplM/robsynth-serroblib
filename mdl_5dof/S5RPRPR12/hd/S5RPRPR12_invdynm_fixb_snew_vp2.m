% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:34
% EndTime: 2019-12-31 18:29:47
% DurationCPUTime: 8.22s
% Computational Cost: add. (123098->291), mult. (302369->370), div. (0->0), fcn. (219141->10), ass. (0->127)
t268 = qJD(1) ^ 2;
t259 = sin(pkin(8));
t296 = qJD(1) * t259;
t261 = cos(pkin(8));
t302 = qJD(1) * t261;
t264 = sin(qJ(1));
t266 = cos(qJ(1));
t245 = -t266 * g(1) - t264 * g(2);
t238 = -t268 * pkin(1) + qJDD(1) * qJ(2) + t245;
t293 = qJD(1) * qJD(2);
t289 = -t261 * g(3) - 0.2e1 * t259 * t293;
t299 = pkin(2) * t261;
t202 = (-pkin(6) * qJDD(1) + t268 * t299 - t238) * t259 + t289;
t225 = -t259 * g(3) + (t238 + 0.2e1 * t293) * t261;
t291 = qJDD(1) * t261;
t255 = t261 ^ 2;
t297 = t255 * t268;
t209 = -pkin(2) * t297 + pkin(6) * t291 + t225;
t263 = sin(qJ(3));
t300 = cos(qJ(3));
t187 = t263 * t202 + t300 * t209;
t290 = t261 * t300;
t236 = -qJD(1) * t290 + t263 * t296;
t278 = t300 * t259 + t261 * t263;
t237 = t278 * qJD(1);
t216 = t236 * mrSges(4,1) + t237 * mrSges(4,2);
t292 = qJDD(1) * t259;
t294 = t237 * qJD(3);
t222 = -qJDD(1) * t290 + t263 * t292 + t294;
t232 = qJD(3) * mrSges(4,1) - t237 * mrSges(4,3);
t215 = t236 * pkin(3) - t237 * qJ(4);
t267 = qJD(3) ^ 2;
t176 = -t267 * pkin(3) + qJDD(3) * qJ(4) - t236 * t215 + t187;
t254 = t259 ^ 2;
t244 = t264 * g(1) - t266 * g(2);
t284 = qJDD(2) - t244;
t221 = (-pkin(1) - t299) * qJDD(1) + (-qJ(2) + (-t254 - t255) * pkin(6)) * t268 + t284;
t295 = t236 * qJD(3);
t223 = t278 * qJDD(1) - t295;
t179 = (-t223 + t295) * qJ(4) + (t222 + t294) * pkin(3) + t221;
t258 = sin(pkin(9));
t260 = cos(pkin(9));
t230 = t258 * qJD(3) + t260 * t237;
t165 = -0.2e1 * qJD(4) * t230 - t258 * t176 + t260 * t179;
t208 = t258 * qJDD(3) + t260 * t223;
t229 = t260 * qJD(3) - t258 * t237;
t163 = (t229 * t236 - t208) * pkin(7) + (t229 * t230 + t222) * pkin(4) + t165;
t166 = 0.2e1 * qJD(4) * t229 + t260 * t176 + t258 * t179;
t206 = t236 * pkin(4) - t230 * pkin(7);
t207 = t260 * qJDD(3) - t258 * t223;
t228 = t229 ^ 2;
t164 = -t228 * pkin(4) + t207 * pkin(7) - t236 * t206 + t166;
t262 = sin(qJ(5));
t265 = cos(qJ(5));
t161 = t265 * t163 - t262 * t164;
t194 = t265 * t229 - t262 * t230;
t175 = t194 * qJD(5) + t262 * t207 + t265 * t208;
t195 = t262 * t229 + t265 * t230;
t184 = -t194 * mrSges(6,1) + t195 * mrSges(6,2);
t234 = qJD(5) + t236;
t188 = -t234 * mrSges(6,2) + t194 * mrSges(6,3);
t220 = qJDD(5) + t222;
t156 = m(6) * t161 + t220 * mrSges(6,1) - t175 * mrSges(6,3) - t195 * t184 + t234 * t188;
t162 = t262 * t163 + t265 * t164;
t174 = -t195 * qJD(5) + t265 * t207 - t262 * t208;
t189 = t234 * mrSges(6,1) - t195 * mrSges(6,3);
t157 = m(6) * t162 - t220 * mrSges(6,2) + t174 * mrSges(6,3) + t194 * t184 - t234 * t189;
t148 = t265 * t156 + t262 * t157;
t197 = -t229 * mrSges(5,1) + t230 * mrSges(5,2);
t204 = -t236 * mrSges(5,2) + t229 * mrSges(5,3);
t146 = m(5) * t165 + t222 * mrSges(5,1) - t208 * mrSges(5,3) - t230 * t197 + t236 * t204 + t148;
t205 = t236 * mrSges(5,1) - t230 * mrSges(5,3);
t285 = -t262 * t156 + t265 * t157;
t147 = m(5) * t166 - t222 * mrSges(5,2) + t207 * mrSges(5,3) + t229 * t197 - t236 * t205 + t285;
t286 = -t258 * t146 + t260 * t147;
t139 = m(4) * t187 - qJDD(3) * mrSges(4,2) - t222 * mrSges(4,3) - qJD(3) * t232 - t236 * t216 + t286;
t186 = t300 * t202 - t263 * t209;
t231 = -qJD(3) * mrSges(4,2) - t236 * mrSges(4,3);
t173 = -qJDD(3) * pkin(3) - t267 * qJ(4) + t237 * t215 + qJDD(4) - t186;
t167 = -t207 * pkin(4) - t228 * pkin(7) + t230 * t206 + t173;
t276 = m(6) * t167 - t174 * mrSges(6,1) + t175 * mrSges(6,2) - t194 * t188 + t195 * t189;
t271 = -m(5) * t173 + t207 * mrSges(5,1) - t208 * mrSges(5,2) + t229 * t204 - t230 * t205 - t276;
t152 = m(4) * t186 + qJDD(3) * mrSges(4,1) - t223 * mrSges(4,3) + qJD(3) * t231 - t237 * t216 + t271;
t130 = t263 * t139 + t300 * t152;
t224 = -t259 * t238 + t289;
t180 = Ifges(6,5) * t195 + Ifges(6,6) * t194 + Ifges(6,3) * t234;
t182 = Ifges(6,1) * t195 + Ifges(6,4) * t194 + Ifges(6,5) * t234;
t149 = -mrSges(6,1) * t167 + mrSges(6,3) * t162 + Ifges(6,4) * t175 + Ifges(6,2) * t174 + Ifges(6,6) * t220 - t195 * t180 + t234 * t182;
t181 = Ifges(6,4) * t195 + Ifges(6,2) * t194 + Ifges(6,6) * t234;
t150 = mrSges(6,2) * t167 - mrSges(6,3) * t161 + Ifges(6,1) * t175 + Ifges(6,4) * t174 + Ifges(6,5) * t220 + t194 * t180 - t234 * t181;
t190 = Ifges(5,5) * t230 + Ifges(5,6) * t229 + Ifges(5,3) * t236;
t192 = Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * t236;
t132 = -mrSges(5,1) * t173 + mrSges(5,3) * t166 + Ifges(5,4) * t208 + Ifges(5,2) * t207 + Ifges(5,6) * t222 - pkin(4) * t276 + pkin(7) * t285 + t265 * t149 + t262 * t150 - t230 * t190 + t236 * t192;
t191 = Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * t236;
t134 = mrSges(5,2) * t173 - mrSges(5,3) * t165 + Ifges(5,1) * t208 + Ifges(5,4) * t207 + Ifges(5,5) * t222 - pkin(7) * t148 - t262 * t149 + t265 * t150 + t229 * t190 - t236 * t191;
t211 = Ifges(4,4) * t237 - Ifges(4,2) * t236 + Ifges(4,6) * qJD(3);
t212 = Ifges(4,1) * t237 - Ifges(4,4) * t236 + Ifges(4,5) * qJD(3);
t273 = -mrSges(4,1) * t186 + mrSges(4,2) * t187 - Ifges(4,5) * t223 + Ifges(4,6) * t222 - Ifges(4,3) * qJDD(3) - pkin(3) * t271 - qJ(4) * t286 - t260 * t132 - t258 * t134 - t237 * t211 - t236 * t212;
t282 = Ifges(3,4) * t259 + Ifges(3,2) * t261;
t283 = Ifges(3,1) * t259 + Ifges(3,4) * t261;
t301 = -mrSges(3,1) * t224 + mrSges(3,2) * t225 - pkin(2) * t130 - (t282 * t296 - t283 * t302) * qJD(1) + t273;
t298 = mrSges(3,2) * t259;
t141 = t260 * t146 + t258 * t147;
t279 = mrSges(3,3) * qJDD(1) + t268 * (-mrSges(3,1) * t261 + t298);
t128 = m(3) * t224 - t279 * t259 + t130;
t287 = t300 * t139 - t263 * t152;
t129 = m(3) * t225 + t279 * t261 + t287;
t288 = -t259 * t128 + t261 * t129;
t281 = Ifges(3,5) * t259 + Ifges(3,6) * t261;
t210 = Ifges(4,5) * t237 - Ifges(4,6) * t236 + Ifges(4,3) * qJD(3);
t122 = mrSges(4,2) * t221 - mrSges(4,3) * t186 + Ifges(4,1) * t223 - Ifges(4,4) * t222 + Ifges(4,5) * qJDD(3) - qJ(4) * t141 - qJD(3) * t211 - t258 * t132 + t260 * t134 - t236 * t210;
t275 = -mrSges(6,1) * t161 + mrSges(6,2) * t162 - Ifges(6,5) * t175 - Ifges(6,6) * t174 - Ifges(6,3) * t220 - t195 * t181 + t194 * t182;
t269 = -mrSges(5,1) * t165 + mrSges(5,2) * t166 - Ifges(5,5) * t208 - Ifges(5,6) * t207 - pkin(4) * t148 - t230 * t191 + t229 * t192 + t275;
t126 = t269 + (-Ifges(4,2) - Ifges(5,3)) * t222 + Ifges(4,6) * qJDD(3) - t237 * t210 - mrSges(4,1) * t221 + Ifges(4,4) * t223 + qJD(3) * t212 + mrSges(4,3) * t187 - pkin(3) * t141;
t235 = -qJDD(1) * pkin(1) - t268 * qJ(2) + t284;
t240 = t281 * qJD(1);
t274 = m(4) * t221 + t222 * mrSges(4,1) + t223 * mrSges(4,2) + t236 * t231 + t237 * t232 + t141;
t118 = -mrSges(3,1) * t235 + mrSges(3,3) * t225 - pkin(2) * t274 + pkin(6) * t287 + t282 * qJDD(1) + t263 * t122 + t300 * t126 - t240 * t296;
t121 = mrSges(3,2) * t235 - mrSges(3,3) * t224 - pkin(6) * t130 + t283 * qJDD(1) + t300 * t122 - t263 * t126 + t240 * t302;
t272 = -m(3) * t235 + mrSges(3,1) * t291 - t274 + (t254 * t268 + t297) * mrSges(3,3);
t277 = -mrSges(2,2) * t245 + qJ(2) * t288 + t261 * t118 + t259 * t121 + pkin(1) * (-mrSges(3,2) * t292 + t272) + mrSges(2,1) * t244 + Ifges(2,3) * qJDD(1);
t135 = -t268 * mrSges(2,2) + m(2) * t244 + t272 + (mrSges(2,1) - t298) * qJDD(1);
t125 = t261 * t128 + t259 * t129;
t123 = m(2) * t245 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t288;
t119 = (Ifges(2,6) - t281) * qJDD(1) + t268 * Ifges(2,5) + mrSges(2,3) * t245 + mrSges(2,1) * g(3) - pkin(1) * t125 + t301;
t116 = -mrSges(2,2) * g(3) - mrSges(2,3) * t244 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - qJ(2) * t125 - t259 * t118 + t261 * t121;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t266 * t116 - t264 * t119 - pkin(5) * (t264 * t123 + t266 * t135), t116, t121, t122, t134, t150; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t116 + t266 * t119 + pkin(5) * (t266 * t123 - t264 * t135), t119, t118, t126, t132, t149; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t277, t277, t281 * qJDD(1) - t301, -t273, Ifges(5,3) * t222 - t269, -t275;];
m_new = t1;

% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR10
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:52
% EndTime: 2019-12-31 19:10:06
% DurationCPUTime: 8.42s
% Computational Cost: add. (132632->291), mult. (312967->367), div. (0->0), fcn. (229305->10), ass. (0->128)
t268 = qJD(1) ^ 2;
t257 = sin(pkin(9));
t295 = qJD(1) * t257;
t258 = cos(pkin(9));
t294 = qJD(1) * t258;
t262 = sin(qJ(1));
t266 = cos(qJ(1));
t244 = -t266 * g(1) - t262 * g(2);
t237 = -t268 * pkin(1) + qJDD(1) * qJ(2) + t244;
t292 = qJD(1) * qJD(2);
t289 = -t258 * g(3) - 0.2e1 * t257 * t292;
t298 = pkin(2) * t258;
t206 = (-pkin(6) * qJDD(1) + t268 * t298 - t237) * t257 + t289;
t224 = -t257 * g(3) + (t237 + 0.2e1 * t292) * t258;
t290 = qJDD(1) * t258;
t254 = t258 ^ 2;
t296 = t254 * t268;
t207 = -pkin(2) * t296 + pkin(6) * t290 + t224;
t261 = sin(qJ(3));
t265 = cos(qJ(3));
t185 = t261 * t206 + t265 * t207;
t235 = -t261 * t295 + t265 * t294;
t279 = t257 * t265 + t258 * t261;
t236 = t279 * qJD(1);
t214 = -t235 * mrSges(4,1) + t236 * mrSges(4,2);
t232 = t236 * qJD(3);
t291 = qJDD(1) * t257;
t221 = -t261 * t291 + t265 * t290 - t232;
t229 = qJD(3) * mrSges(4,1) - t236 * mrSges(4,3);
t219 = -t235 * pkin(3) - t236 * pkin(7);
t267 = qJD(3) ^ 2;
t174 = -t267 * pkin(3) + qJDD(3) * pkin(7) + t235 * t219 + t185;
t253 = t257 ^ 2;
t243 = t262 * g(1) - t266 * g(2);
t284 = qJDD(2) - t243;
t220 = (-pkin(1) - t298) * qJDD(1) + (-qJ(2) + (-t253 - t254) * pkin(6)) * t268 + t284;
t293 = t235 * qJD(3);
t222 = t279 * qJDD(1) + t293;
t177 = (-t222 - t293) * pkin(7) + (-t221 + t232) * pkin(3) + t220;
t260 = sin(qJ(4));
t264 = cos(qJ(4));
t164 = -t260 * t174 + t264 * t177;
t226 = t264 * qJD(3) - t260 * t236;
t194 = t226 * qJD(4) + t260 * qJDD(3) + t264 * t222;
t218 = qJDD(4) - t221;
t227 = t260 * qJD(3) + t264 * t236;
t233 = qJD(4) - t235;
t161 = (t226 * t233 - t194) * pkin(8) + (t226 * t227 + t218) * pkin(4) + t164;
t165 = t264 * t174 + t260 * t177;
t193 = -t227 * qJD(4) + t264 * qJDD(3) - t260 * t222;
t205 = t233 * pkin(4) - t227 * pkin(8);
t225 = t226 ^ 2;
t162 = -t225 * pkin(4) + t193 * pkin(8) - t233 * t205 + t165;
t259 = sin(qJ(5));
t263 = cos(qJ(5));
t159 = t263 * t161 - t259 * t162;
t195 = t263 * t226 - t259 * t227;
t170 = t195 * qJD(5) + t259 * t193 + t263 * t194;
t196 = t259 * t226 + t263 * t227;
t182 = -t195 * mrSges(6,1) + t196 * mrSges(6,2);
t231 = qJD(5) + t233;
t186 = -t231 * mrSges(6,2) + t195 * mrSges(6,3);
t213 = qJDD(5) + t218;
t154 = m(6) * t159 + t213 * mrSges(6,1) - t170 * mrSges(6,3) - t196 * t182 + t231 * t186;
t160 = t259 * t161 + t263 * t162;
t169 = -t196 * qJD(5) + t263 * t193 - t259 * t194;
t187 = t231 * mrSges(6,1) - t196 * mrSges(6,3);
t155 = m(6) * t160 - t213 * mrSges(6,2) + t169 * mrSges(6,3) + t195 * t182 - t231 * t187;
t146 = t263 * t154 + t259 * t155;
t198 = -t226 * mrSges(5,1) + t227 * mrSges(5,2);
t201 = -t233 * mrSges(5,2) + t226 * mrSges(5,3);
t144 = m(5) * t164 + t218 * mrSges(5,1) - t194 * mrSges(5,3) - t227 * t198 + t233 * t201 + t146;
t202 = t233 * mrSges(5,1) - t227 * mrSges(5,3);
t285 = -t259 * t154 + t263 * t155;
t145 = m(5) * t165 - t218 * mrSges(5,2) + t193 * mrSges(5,3) + t226 * t198 - t233 * t202 + t285;
t286 = -t260 * t144 + t264 * t145;
t137 = m(4) * t185 - qJDD(3) * mrSges(4,2) + t221 * mrSges(4,3) - qJD(3) * t229 + t235 * t214 + t286;
t184 = t265 * t206 - t261 * t207;
t228 = -qJD(3) * mrSges(4,2) + t235 * mrSges(4,3);
t173 = -qJDD(3) * pkin(3) - t267 * pkin(7) + t236 * t219 - t184;
t163 = -t193 * pkin(4) - t225 * pkin(8) + t227 * t205 + t173;
t276 = m(6) * t163 - t169 * mrSges(6,1) + t170 * mrSges(6,2) - t195 * t186 + t196 * t187;
t271 = -m(5) * t173 + t193 * mrSges(5,1) - t194 * mrSges(5,2) + t226 * t201 - t227 * t202 - t276;
t150 = m(4) * t184 + qJDD(3) * mrSges(4,1) - t222 * mrSges(4,3) + qJD(3) * t228 - t236 * t214 + t271;
t130 = t261 * t137 + t265 * t150;
t223 = -t257 * t237 + t289;
t178 = Ifges(6,5) * t196 + Ifges(6,6) * t195 + Ifges(6,3) * t231;
t180 = Ifges(6,1) * t196 + Ifges(6,4) * t195 + Ifges(6,5) * t231;
t147 = -mrSges(6,1) * t163 + mrSges(6,3) * t160 + Ifges(6,4) * t170 + Ifges(6,2) * t169 + Ifges(6,6) * t213 - t196 * t178 + t231 * t180;
t179 = Ifges(6,4) * t196 + Ifges(6,2) * t195 + Ifges(6,6) * t231;
t148 = mrSges(6,2) * t163 - mrSges(6,3) * t159 + Ifges(6,1) * t170 + Ifges(6,4) * t169 + Ifges(6,5) * t213 + t195 * t178 - t231 * t179;
t188 = Ifges(5,5) * t227 + Ifges(5,6) * t226 + Ifges(5,3) * t233;
t190 = Ifges(5,1) * t227 + Ifges(5,4) * t226 + Ifges(5,5) * t233;
t127 = -mrSges(5,1) * t173 + mrSges(5,3) * t165 + Ifges(5,4) * t194 + Ifges(5,2) * t193 + Ifges(5,6) * t218 - pkin(4) * t276 + pkin(8) * t285 + t263 * t147 + t259 * t148 - t227 * t188 + t233 * t190;
t189 = Ifges(5,4) * t227 + Ifges(5,2) * t226 + Ifges(5,6) * t233;
t132 = mrSges(5,2) * t173 - mrSges(5,3) * t164 + Ifges(5,1) * t194 + Ifges(5,4) * t193 + Ifges(5,5) * t218 - pkin(8) * t146 - t259 * t147 + t263 * t148 + t226 * t188 - t233 * t189;
t209 = Ifges(4,4) * t236 + Ifges(4,2) * t235 + Ifges(4,6) * qJD(3);
t210 = Ifges(4,1) * t236 + Ifges(4,4) * t235 + Ifges(4,5) * qJD(3);
t273 = -mrSges(4,1) * t184 + mrSges(4,2) * t185 - Ifges(4,5) * t222 - Ifges(4,6) * t221 - Ifges(4,3) * qJDD(3) - pkin(3) * t271 - pkin(7) * t286 - t264 * t127 - t260 * t132 - t236 * t209 + t235 * t210;
t282 = Ifges(3,4) * t257 + Ifges(3,2) * t258;
t283 = Ifges(3,1) * t257 + Ifges(3,4) * t258;
t299 = -mrSges(3,1) * t223 + mrSges(3,2) * t224 - pkin(2) * t130 - (t282 * t295 - t283 * t294) * qJD(1) + t273;
t297 = mrSges(3,2) * t257;
t139 = t264 * t144 + t260 * t145;
t278 = mrSges(3,3) * qJDD(1) + t268 * (-mrSges(3,1) * t258 + t297);
t128 = m(3) * t223 - t278 * t257 + t130;
t287 = t265 * t137 - t261 * t150;
t129 = m(3) * t224 + t278 * t258 + t287;
t288 = -t257 * t128 + t258 * t129;
t281 = Ifges(3,5) * t257 + Ifges(3,6) * t258;
t208 = Ifges(4,5) * t236 + Ifges(4,6) * t235 + Ifges(4,3) * qJD(3);
t120 = mrSges(4,2) * t220 - mrSges(4,3) * t184 + Ifges(4,1) * t222 + Ifges(4,4) * t221 + Ifges(4,5) * qJDD(3) - pkin(7) * t139 - qJD(3) * t209 - t260 * t127 + t264 * t132 + t235 * t208;
t275 = -mrSges(6,1) * t159 + mrSges(6,2) * t160 - Ifges(6,5) * t170 - Ifges(6,6) * t169 - Ifges(6,3) * t213 - t196 * t179 + t195 * t180;
t269 = mrSges(5,1) * t164 - mrSges(5,2) * t165 + Ifges(5,5) * t194 + Ifges(5,6) * t193 + Ifges(5,3) * t218 + pkin(4) * t146 + t227 * t189 - t226 * t190 - t275;
t124 = -mrSges(4,1) * t220 + mrSges(4,3) * t185 + Ifges(4,4) * t222 + Ifges(4,2) * t221 + Ifges(4,6) * qJDD(3) - pkin(3) * t139 + qJD(3) * t210 - t236 * t208 - t269;
t234 = -qJDD(1) * pkin(1) - t268 * qJ(2) + t284;
t239 = t281 * qJD(1);
t274 = m(4) * t220 - t221 * mrSges(4,1) + t222 * mrSges(4,2) - t235 * t228 + t236 * t229 + t139;
t116 = -mrSges(3,1) * t234 + mrSges(3,3) * t224 - pkin(2) * t274 + pkin(6) * t287 + t282 * qJDD(1) + t261 * t120 + t265 * t124 - t239 * t295;
t119 = mrSges(3,2) * t234 - mrSges(3,3) * t223 - pkin(6) * t130 + t283 * qJDD(1) + t265 * t120 - t261 * t124 + t239 * t294;
t272 = -m(3) * t234 + mrSges(3,1) * t290 - t274 + (t253 * t268 + t296) * mrSges(3,3);
t277 = -mrSges(2,2) * t244 + qJ(2) * t288 + t258 * t116 + t257 * t119 + pkin(1) * (-mrSges(3,2) * t291 + t272) + mrSges(2,1) * t243 + Ifges(2,3) * qJDD(1);
t133 = -t268 * mrSges(2,2) + m(2) * t243 + t272 + (mrSges(2,1) - t297) * qJDD(1);
t123 = t258 * t128 + t257 * t129;
t121 = m(2) * t244 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t288;
t117 = (Ifges(2,6) - t281) * qJDD(1) + t268 * Ifges(2,5) + mrSges(2,3) * t244 + mrSges(2,1) * g(3) - pkin(1) * t123 + t299;
t114 = -mrSges(2,2) * g(3) - mrSges(2,3) * t243 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - qJ(2) * t123 - t257 * t116 + t258 * t119;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t266 * t114 - t262 * t117 - pkin(5) * (t262 * t121 + t266 * t133), t114, t119, t120, t132, t148; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t262 * t114 + t266 * t117 + pkin(5) * (t266 * t121 - t262 * t133), t117, t116, t124, t127, t147; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t277, t277, t281 * qJDD(1) - t299, -t273, t269, -t275;];
m_new = t1;

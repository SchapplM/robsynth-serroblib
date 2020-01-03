% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:47
% EndTime: 2019-12-31 16:33:51
% DurationCPUTime: 2.26s
% Computational Cost: add. (5297->266), mult. (7789->395), div. (0->0), fcn. (5338->8), ass. (0->184)
t253 = qJD(2) + qJD(3);
t251 = t253 ^ 2;
t263 = cos(qJ(3));
t252 = qJDD(2) + qJDD(3);
t260 = sin(qJ(3));
t288 = t260 * t252;
t218 = t263 * t251 + t288;
t282 = t263 * t252;
t221 = t260 * t251 - t282;
t261 = sin(qJ(2));
t264 = cos(qJ(2));
t178 = t261 * t218 + t264 * t221;
t257 = sin(pkin(7));
t308 = t257 * t178;
t258 = cos(pkin(7));
t307 = t258 * t178;
t236 = t257 * g(1) - t258 * g(2);
t187 = pkin(5) * t218 - t263 * t236;
t304 = pkin(5) * t221 - t260 * t236;
t131 = pkin(4) * t178 + t261 * t187 + t264 * t304;
t267 = t264 * t218 - t261 * t221;
t132 = pkin(4) * t267 + t264 * t187 - t261 * t304;
t237 = t258 * g(1) + t257 * g(2);
t256 = g(3) - qJDD(1);
t210 = -t264 * t237 - t261 * t256;
t301 = qJD(2) ^ 2;
t208 = -t301 * pkin(2) + t210;
t209 = -t261 * t237 + t264 * t256;
t266 = qJDD(2) * pkin(2) - t209;
t152 = t260 * t208 - t263 * t266;
t153 = t263 * t208 + t260 * t266;
t271 = t260 * t152 + t263 * t153;
t121 = t263 * t152 - t260 * t153;
t281 = t264 * t121;
t106 = -t261 * t271 + t281;
t287 = t261 * t121;
t107 = t264 * t271 + t287;
t226 = t258 * t236;
t198 = -t257 * t237 + t226;
t259 = sin(qJ(4));
t254 = t259 ^ 2;
t300 = t254 * t251;
t235 = t264 * qJDD(2) - t261 * t301;
t299 = t257 * t235;
t298 = t257 * t236;
t296 = t257 * t256;
t295 = t258 * t235;
t294 = t258 * t256;
t146 = -t252 * pkin(3) - t251 * pkin(6) + t152;
t293 = t259 * t146;
t262 = cos(qJ(4));
t242 = t259 * t251 * t262;
t231 = qJDD(4) + t242;
t292 = t259 * t231;
t232 = qJDD(4) - t242;
t291 = t259 * t232;
t290 = t259 * t252;
t286 = t262 * t146;
t285 = t262 * t231;
t284 = t262 * t232;
t244 = t262 * t252;
t147 = -t251 * pkin(3) + t252 * pkin(6) + t153;
t140 = t262 * t147 - t259 * t236;
t255 = t262 ^ 2;
t280 = t254 + t255;
t279 = qJD(4) * t253;
t278 = t259 * t279;
t277 = t262 * t279;
t139 = t259 * t147 + t262 * t236;
t112 = t259 * t139 + t262 * t140;
t217 = t280 * t252;
t245 = t255 * t251;
t222 = t245 + t300;
t173 = t260 * t217 + t263 * t222;
t176 = t263 * t217 - t260 * t222;
t141 = t264 * t173 + t261 * t176;
t142 = -t261 * t173 + t264 * t176;
t276 = -pkin(1) * t141 - pkin(2) * t173 - pkin(3) * t222 - pkin(6) * t217 + qJ(1) * t142 - t112;
t275 = pkin(1) * t267 + pkin(2) * t218 + qJ(1) * t178 + t153;
t274 = pkin(1) * t178 + pkin(2) * t221 - qJ(1) * t267 + t152;
t234 = t261 * qJDD(2) + t264 * t301;
t273 = -pkin(1) * t234 + qJ(1) * t235 - t210;
t272 = pkin(1) * t235 + qJ(1) * t234 - t209;
t160 = t261 * t209 + t264 * t210;
t199 = -t258 * t237 - t298;
t269 = t260 * t242;
t268 = t263 * t242;
t206 = pkin(4) * t234 - t264 * t236;
t205 = -pkin(4) * t235 - t261 * t236;
t111 = t262 * t139 - t259 * t140;
t159 = t264 * t209 - t261 * t210;
t265 = qJD(4) ^ 2;
t241 = -t245 - t265;
t240 = t245 - t265;
t239 = -t265 - t300;
t238 = t265 - t300;
t233 = pkin(1) * t236;
t225 = t258 * t234;
t224 = t257 * t234;
t223 = t245 - t300;
t216 = t244 - 0.2e1 * t278;
t215 = t244 - t278;
t214 = t277 + t290;
t213 = 0.2e1 * t277 + t290;
t212 = t280 * t279;
t204 = t260 * qJDD(4) + t263 * t212;
t203 = -t263 * qJDD(4) + t260 * t212;
t197 = t262 * t214 - t254 * t279;
t196 = -t259 * t215 - t255 * t279;
t195 = -t259 * t239 - t284;
t194 = -t259 * t240 - t284;
t193 = -t259 * t238 + t285;
t192 = t262 * t241 - t292;
t191 = t262 * t240 - t291;
t190 = -t262 * t238 - t292;
t189 = t262 * t239 - t291;
t188 = t259 * t241 + t285;
t183 = (-t215 + t278) * t262;
t182 = (-t214 - t277) * t259;
t172 = -t259 * t213 + t262 * t216;
t171 = -t262 * t213 - t259 * t216;
t170 = t258 * t267;
t169 = t257 * t267;
t168 = t263 * t193 + t259 * t288;
t167 = t263 * t191 + t260 * t244;
t166 = t260 * t193 - t259 * t282;
t165 = t260 * t191 - t262 * t282;
t164 = t263 * t197 - t269;
t163 = t263 * t196 + t269;
t162 = t260 * t197 + t268;
t161 = t260 * t196 - t268;
t157 = t263 * t195 + t260 * t213;
t156 = t263 * t192 - t260 * t216;
t155 = t260 * t195 - t263 * t213;
t154 = t260 * t192 + t263 * t216;
t150 = -t261 * t203 + t264 * t204;
t149 = t263 * t172 - t260 * t223;
t148 = t260 * t172 + t263 * t223;
t145 = t258 * t160 - t298;
t144 = t257 * t160 + t226;
t138 = -t261 * t166 + t264 * t168;
t137 = -t261 * t165 + t264 * t167;
t136 = -t261 * t162 + t264 * t164;
t135 = -t261 * t161 + t264 * t163;
t134 = -pkin(6) * t189 + t286;
t133 = -pkin(6) * t188 + t293;
t128 = -t261 * t155 + t264 * t157;
t127 = -t261 * t154 + t264 * t156;
t126 = t264 * t155 + t261 * t157;
t125 = t264 * t154 + t261 * t156;
t124 = -pkin(3) * t189 + t140;
t123 = -pkin(3) * t188 + t139;
t118 = pkin(2) * t236 + pkin(5) * t271;
t117 = -t261 * t148 + t264 * t149;
t116 = t258 * t128 + t257 * t189;
t115 = t258 * t127 + t257 * t188;
t114 = t257 * t128 - t258 * t189;
t113 = t257 * t127 - t258 * t188;
t109 = -pkin(5) * t173 + t263 * t111;
t108 = pkin(5) * t176 + t260 * t111;
t104 = t263 * t112 + t260 * t146;
t103 = t260 * t112 - t263 * t146;
t102 = t258 * t107 - t298;
t101 = t257 * t107 + t226;
t100 = -pkin(5) * t155 - t260 * t124 + t263 * t134;
t99 = -pkin(5) * t154 - t260 * t123 + t263 * t133;
t98 = -pkin(1) * t126 - pkin(2) * t155 + pkin(3) * t213 - pkin(6) * t195 - t293;
t97 = -pkin(1) * t125 - pkin(2) * t154 - pkin(3) * t216 - pkin(6) * t192 + t286;
t96 = -pkin(2) * t189 + pkin(5) * t157 + t263 * t124 + t260 * t134;
t95 = -pkin(2) * t188 + pkin(5) * t156 + t263 * t123 + t260 * t133;
t93 = pkin(1) * t106 + pkin(2) * t121;
t92 = -t261 * t103 + t264 * t104;
t91 = t264 * t103 + t261 * t104;
t90 = pkin(4) * t106 + pkin(5) * t281 - t261 * t118;
t89 = -pkin(4) * t141 - t261 * t108 + t264 * t109;
t88 = -pkin(5) * t103 - (pkin(3) * t260 - pkin(6) * t263) * t111;
t87 = -t111 * t257 + t258 * t92;
t86 = t111 * t258 + t257 * t92;
t85 = -pkin(4) * t126 + t264 * t100 - t261 * t96;
t84 = -pkin(4) * t125 - t261 * t95 + t264 * t99;
t83 = pkin(5) * t104 - (-pkin(3) * t263 - pkin(6) * t260 - pkin(2)) * t111;
t82 = -pkin(1) * t91 - pkin(2) * t103 + pkin(3) * t146 - pkin(6) * t112;
t81 = -pkin(4) * t91 - t261 * t83 + t264 * t88;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, 0, 0, 0, 0, 0, 0, -t225, -t295, 0, t145, 0, 0, 0, 0, 0, 0, -t170, t307, 0, t102, 0, 0, 0, 0, 0, 0, t115, t116, t258 * t142, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, 0, 0, 0, 0, 0, 0, -t224, -t299, 0, t144, 0, 0, 0, 0, 0, 0, -t169, t308, 0, t101, 0, 0, 0, 0, 0, 0, t113, t114, t257 * t142, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, 0, 0, 0, 0, 0, 0, t235, -t234, 0, -t159, 0, 0, 0, 0, 0, 0, -t178, -t267, 0, -t106, 0, 0, 0, 0, 0, 0, t125, t126, t141, t91; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t296, -t294, -t198, -qJ(1) * t198, 0, 0, t295, 0, -t225, t257 * qJDD(2), t258 * t205 + t272 * t257, t258 * t206 + t273 * t257, t258 * t159, -qJ(1) * t144 - (pkin(1) * t257 - pkin(4) * t258) * t159, 0, 0, -t307, 0, -t170, t257 * t252, t258 * t131 - t274 * t257, t258 * t132 - t257 * t275, t258 * t106, -qJ(1) * t101 - t257 * t93 + t258 * t90, t258 * t136 - t257 * t182, t258 * t117 - t257 * t171, t258 * t138 - t257 * t190, t258 * t135 - t257 * t183, t258 * t137 - t257 * t194, t258 * t150, -qJ(1) * t113 - t257 * t97 + t258 * t84, -qJ(1) * t114 - t257 * t98 + t258 * t85, -t257 * t276 + t258 * t89, -qJ(1) * t86 - t257 * t82 + t258 * t81; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t294, -t296, t199, qJ(1) * t199, 0, 0, t299, 0, -t224, -t258 * qJDD(2), t257 * t205 - t272 * t258, t257 * t206 - t273 * t258, t257 * t159, qJ(1) * t145 - (-pkin(1) * t258 - pkin(4) * t257) * t159, 0, 0, -t308, 0, -t169, -t258 * t252, t257 * t131 + t258 * t274, t257 * t132 + t258 * t275, t257 * t106, qJ(1) * t102 + t257 * t90 + t258 * t93, t257 * t136 + t258 * t182, t257 * t117 + t258 * t171, t257 * t138 + t258 * t190, t257 * t135 + t258 * t183, t257 * t137 + t258 * t194, t257 * t150, qJ(1) * t115 + t257 * t84 + t258 * t97, qJ(1) * t116 + t257 * t85 + t258 * t98, t257 * t89 + t258 * t276, qJ(1) * t87 + t257 * t81 + t258 * t82; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t236, t237, 0, 0, 0, 0, t234, 0, t235, 0, -t206, t205, t160, pkin(4) * t160 + t233, 0, 0, t267, 0, -t178, 0, -t132, t131, t107, pkin(4) * t107 + pkin(5) * t287 + t264 * t118 + t233, t264 * t162 + t261 * t164, t264 * t148 + t261 * t149, t264 * t166 + t261 * t168, t264 * t161 + t261 * t163, t264 * t165 + t261 * t167, t264 * t203 + t261 * t204, -pkin(1) * t188 + pkin(4) * t127 + t261 * t99 + t264 * t95, -pkin(1) * t189 + pkin(4) * t128 + t261 * t100 + t264 * t96, pkin(4) * t142 + t264 * t108 + t261 * t109, pkin(1) * t111 + pkin(4) * t92 + t261 * t88 + t264 * t83;];
tauB_reg = t1;

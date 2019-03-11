% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:34
% EndTime: 2019-03-09 03:59:43
% DurationCPUTime: 4.11s
% Computational Cost: add. (4664->426), mult. (9631->563), div. (0->0), fcn. (7128->14), ass. (0->228)
t187 = sin(qJ(3));
t191 = cos(qJ(3));
t285 = sin(pkin(10));
t286 = cos(pkin(10));
t129 = -t187 * t285 + t191 * t286;
t203 = t187 * t286 + t191 * t285;
t307 = t203 * qJD(1);
t317 = qJD(5) + t307;
t124 = t129 * qJD(1);
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t252 = t190 * qJD(3);
t94 = t124 * t186 - t252;
t322 = t317 * t94;
t155 = pkin(3) * t285 + pkin(8);
t177 = qJ(3) + pkin(10);
t163 = sin(t177);
t164 = cos(t177);
t188 = sin(qJ(1));
t192 = cos(qJ(1));
t308 = g(1) * t188 - g(2) * t192;
t202 = -g(3) * t163 + t164 * t308;
t193 = -pkin(1) - pkin(7);
t142 = qJDD(1) * t193 + qJDD(2);
t133 = t191 * t142;
t143 = qJD(1) * t193 + qJD(2);
t245 = t191 * qJDD(1);
t248 = qJD(1) * qJD(4);
t249 = qJD(1) * qJD(3);
t258 = qJD(3) * t187;
t54 = -t191 * t248 - t143 * t258 + qJDD(3) * pkin(3) + t133 + (t187 * t249 - t245) * qJ(4);
t257 = qJD(3) * t191;
t64 = (-qJ(4) * qJD(1) + t143) * t257 + (-qJ(4) * qJDD(1) + t142 - t248) * t187;
t25 = -t285 * t64 + t286 * t54;
t21 = -qJDD(3) * pkin(4) - t25;
t321 = qJD(5) * t155 * t317 + t202 + t21;
t185 = sin(qJ(6));
t189 = cos(qJ(6));
t96 = qJD(3) * t186 + t124 * t190;
t214 = t185 * t94 - t189 * t96;
t40 = t185 * t96 + t189 * t94;
t320 = t214 * t40;
t256 = qJD(5) * t186;
t281 = t307 * t186;
t319 = t256 + t281;
t195 = qJD(1) ^ 2;
t204 = -qJ(2) * t195 - t308;
t226 = t317 * t190;
t227 = qJDD(1) * t285;
t228 = qJDD(1) * t286;
t84 = -qJD(3) * t124 - t187 * t228 - t191 * t227;
t82 = -qJDD(5) + t84;
t290 = t186 * t82;
t318 = t226 * t317 - t290;
t316 = t214 ^ 2 - t40 ^ 2;
t109 = qJD(6) + t317;
t253 = qJD(6) * t189;
t254 = qJD(6) * t185;
t85 = qJD(3) * t307 + t187 * t227 - t191 * t228;
t36 = qJD(5) * t252 + qJDD(3) * t186 - t124 * t256 - t190 * t85;
t231 = -qJDD(3) * t190 - t186 * t85;
t37 = qJD(5) * t96 + t231;
t8 = -t185 * t37 + t189 * t36 - t253 * t94 - t254 * t96;
t315 = t109 * t40 + t8;
t183 = qJ(5) + qJ(6);
t170 = cos(t183);
t274 = t170 * t188;
t169 = sin(t183);
t275 = t169 * t192;
t103 = t163 * t274 + t275;
t273 = t170 * t192;
t276 = t169 * t188;
t105 = t163 * t273 - t276;
t259 = qJD(1) * t191;
t113 = -qJ(4) * t259 + t143 * t191;
t106 = qJD(3) * pkin(3) + t113;
t260 = qJD(1) * t187;
t112 = -qJ(4) * t260 + t143 * t187;
t235 = t286 * t112;
t58 = t106 * t285 + t235;
t48 = qJD(3) * pkin(8) + t58;
t135 = pkin(3) * t260 + qJD(1) * qJ(2) + qJD(4);
t59 = pkin(4) * t307 - pkin(8) * t124 + t135;
t24 = t186 * t59 + t190 * t48;
t15 = -pkin(9) * t94 + t24;
t12 = t15 * t254;
t301 = g(3) * t164;
t99 = t285 * t112;
t57 = t106 * t286 - t99;
t47 = -qJD(3) * pkin(4) - t57;
t32 = pkin(5) * t94 + t47;
t314 = g(1) * t103 - g(2) * t105 + t170 * t301 + t32 * t40 + t12;
t102 = -t163 * t276 + t273;
t104 = t163 * t275 + t274;
t26 = t285 * t54 + t286 * t64;
t22 = qJDD(3) * pkin(8) + t26;
t178 = qJDD(1) * qJ(2);
t179 = qJD(1) * qJD(2);
t239 = t191 * t249;
t246 = t187 * qJDD(1);
t213 = qJDD(4) + t178 + t179 + (t239 + t246) * pkin(3);
t29 = -pkin(4) * t84 + pkin(8) * t85 + t213;
t28 = t190 * t29;
t2 = -pkin(5) * t82 - pkin(9) * t36 - qJD(5) * t24 - t186 * t22 + t28;
t255 = qJD(5) * t190;
t211 = t186 * t29 + t190 * t22 + t255 * t59 - t256 * t48;
t3 = -pkin(9) * t37 + t211;
t243 = -t185 * t3 + t189 * t2;
t23 = -t186 * t48 + t190 * t59;
t14 = -pkin(9) * t96 + t23;
t11 = pkin(5) * t317 + t14;
t291 = t15 * t189;
t5 = t11 * t185 + t291;
t313 = -g(1) * t102 - g(2) * t104 - qJD(6) * t5 + t169 * t301 + t214 * t32 + t243;
t199 = qJD(6) * t214 - t185 * t36 - t189 * t37;
t312 = -t109 * t214 + t199;
t131 = t185 * t186 - t189 * t190;
t80 = -qJDD(6) + t82;
t272 = t185 * t190;
t132 = t186 * t189 + t272;
t306 = qJD(5) + qJD(6);
t88 = t306 * t132;
t311 = -t109 * t88 + t131 * t80;
t75 = t132 * t129;
t218 = g(1) * t192 + g(2) * t188;
t244 = 0.2e1 * t179;
t305 = 0.2e1 * t178 + t244 - t218;
t292 = t132 * t80;
t87 = t306 * t131;
t296 = -t131 * t307 - t87;
t304 = -t109 * t296 + t292;
t300 = g(3) * t187;
t173 = t187 * pkin(3);
t299 = pkin(9) + t155;
t70 = t113 * t286 - t99;
t72 = pkin(3) * t259 + pkin(4) * t124 + pkin(8) * t307;
t298 = t186 * t72 + t190 * t70;
t267 = qJ(2) + t173;
t83 = pkin(4) * t203 - pkin(8) * t129 + t267;
t266 = qJ(4) - t193;
t136 = t266 * t187;
t137 = t266 * t191;
t90 = -t136 * t286 - t137 * t285;
t86 = t190 * t90;
t297 = t186 * t83 + t86;
t67 = t132 * t307;
t295 = t88 + t67;
t294 = t124 * t40;
t293 = t124 * t94;
t74 = t190 * t82;
t289 = t36 * t186;
t288 = t214 * t124;
t287 = t96 * t124;
t284 = pkin(1) * qJDD(1);
t229 = qJD(3) * t285;
t230 = qJD(3) * t286;
t122 = t187 * t229 - t191 * t230;
t282 = t109 * t122;
t123 = -t187 * t230 - t191 * t229;
t280 = t123 * t186;
t279 = t123 * t190;
t278 = t129 * t186;
t277 = t129 * t190;
t271 = t186 * t188;
t270 = t186 * t192;
t269 = t188 * t190;
t268 = t190 * t192;
t265 = pkin(1) * t192 + qJ(2) * t188;
t182 = t191 ^ 2;
t263 = t187 ^ 2 - t182;
t194 = qJD(3) ^ 2;
t262 = -t194 - t195;
t261 = qJD(1) * t135;
t251 = pkin(3) * t257 + qJD(2);
t247 = qJDD(3) * t187;
t242 = t129 * t256;
t241 = t129 * t255;
t238 = qJD(6) * t11 + t3;
t236 = qJD(5) * t299;
t232 = -qJD(5) * t59 - t22;
t225 = qJDD(2) - t284;
t224 = -qJD(5) * t203 - qJD(1);
t223 = -t109 * t67 + t311;
t69 = t113 * t285 + t235;
t222 = pkin(5) * t319 - t69;
t156 = -pkin(3) * t286 - pkin(4);
t221 = -t255 * t48 + t28;
t127 = t299 * t186;
t216 = pkin(9) * t281 + qJD(6) * t127 + t186 * t236 + t298;
t128 = t299 * t190;
t66 = t190 * t72;
t215 = pkin(5) * t124 + qJD(6) * t128 - t186 * t70 + t66 + (pkin(9) * t307 + t236) * t190;
t110 = -t191 * qJD(4) + t258 * t266;
t111 = -qJD(3) * t137 - qJD(4) * t187;
t62 = -t110 * t286 + t111 * t285;
t89 = -t136 * t285 + t137 * t286;
t212 = -t317 * t319 - t74;
t63 = t110 * t285 + t111 * t286;
t71 = -pkin(4) * t122 - pkin(8) * t123 + t251;
t210 = t186 * t71 + t190 * t63 + t255 * t83 - t256 * t90;
t209 = t241 + t280;
t207 = t155 * t82 + t317 * t47;
t206 = t109 * t132;
t205 = 0.2e1 * qJ(2) * t249 + qJDD(3) * t193;
t198 = -t122 * t58 + t123 * t57 + t129 * t25 + t203 * t26 - t308;
t197 = -t193 * t194 + t305;
t184 = -qJ(4) - pkin(7);
t172 = t192 * qJ(2);
t167 = qJDD(3) * t191;
t140 = -pkin(5) * t190 + t156;
t118 = t163 * t268 - t271;
t117 = t163 * t270 + t269;
t116 = t163 * t269 + t270;
t115 = -t163 * t271 + t268;
t78 = t190 * t83;
t76 = t131 * t129;
t61 = t190 * t71;
t55 = pkin(5) * t278 + t89;
t31 = pkin(5) * t209 + t62;
t30 = -pkin(9) * t278 + t297;
t20 = pkin(5) * t203 - pkin(9) * t277 - t186 * t90 + t78;
t17 = t123 * t272 - t185 * t242 - t254 * t278 + (t277 * t306 + t280) * t189;
t16 = -t131 * t123 - t306 * t75;
t10 = pkin(5) * t37 + t21;
t7 = -pkin(9) * t209 + t210;
t6 = -pkin(9) * t279 - pkin(5) * t122 - t186 * t63 + t61 + (-t86 + (pkin(9) * t129 - t83) * t186) * qJD(5);
t4 = t11 * t189 - t15 * t185;
t1 = [qJDD(1), t308, t218, qJDD(2) - t308 - 0.2e1 * t284, t305, -t225 * pkin(1) - g(1) * (-pkin(1) * t188 + t172) - g(2) * t265 + (t244 + t178) * qJ(2), qJDD(1) * t182 - 0.2e1 * t187 * t239, -0.2e1 * t187 * t245 + 0.2e1 * t249 * t263, -t187 * t194 + t167, -t191 * t194 - t247, 0, t187 * t197 + t191 * t205, -t187 * t205 + t191 * t197, t124 * t62 - t307 * t63 + t84 * t90 - t85 * t89 - t198, t26 * t90 + t58 * t63 - t25 * t89 - t57 * t62 + t213 * t267 + t135 * t251 - g(1) * (t192 * t173 + t172 + (-pkin(1) + t184) * t188) - g(2) * (t173 * t188 - t184 * t192 + t265) t96 * t279 + (t190 * t36 - t256 * t96) * t129 (-t186 * t96 - t190 * t94) * t123 + (-t289 - t190 * t37 + (t186 * t94 - t190 * t96) * qJD(5)) * t129, -t82 * t277 - t122 * t96 + t203 * t36 + (-t242 + t279) * t317, t122 * t94 - t203 * t37 - t209 * t317 + t278 * t82, -t122 * t317 - t203 * t82 (-t255 * t90 + t61) * t317 - t78 * t82 + t221 * t203 - t23 * t122 + t62 * t94 + t89 * t37 + t47 * t241 - g(1) * t118 - g(2) * t116 + ((-qJD(5) * t83 - t63) * t317 + t90 * t82 + t232 * t203 + t21 * t129 + t47 * t123) * t186, -t210 * t317 + t297 * t82 - t211 * t203 + t24 * t122 + t62 * t96 + t89 * t36 + t47 * t279 + g(1) * t117 - g(2) * t115 + (t21 * t190 - t256 * t47) * t129, -t16 * t214 - t76 * t8, -t16 * t40 + t17 * t214 - t199 * t76 - t75 * t8, t109 * t16 + t122 * t214 + t203 * t8 + t76 * t80, -t109 * t17 + t122 * t40 + t199 * t203 + t75 * t80, -t203 * t80 - t282 (-t185 * t7 + t189 * t6) * t109 - (-t185 * t30 + t189 * t20) * t80 + t243 * t203 - t4 * t122 + t31 * t40 - t55 * t199 + t10 * t75 + t32 * t17 - g(1) * t105 - g(2) * t103 + ((-t185 * t20 - t189 * t30) * t109 - t5 * t203) * qJD(6), g(1) * t104 - g(2) * t102 - t10 * t76 + t12 * t203 + t5 * t122 + t32 * t16 - t31 * t214 + t55 * t8 + (-(-qJD(6) * t30 + t6) * t109 + t20 * t80 - t2 * t203) * t185 + (-(qJD(6) * t20 + t7) * t109 + t30 * t80 - t238 * t203) * t189; 0, 0, 0, qJDD(1), -t195, t225 + t204, 0, 0, 0, 0, 0, t187 * t262 + t167, t191 * t262 - t247, t122 * t307 - t123 * t124 + t129 * t85 + t203 * t84, t198 - t261, 0, 0, 0, 0, 0, t203 * t290 - t123 * t94 - t129 * t37 + (t122 * t186 + t190 * t224) * t317, t203 * t74 - t123 * t96 - t129 * t36 + (t122 * t190 - t186 * t224) * t317, 0, 0, 0, 0, 0, -t123 * t40 + t129 * t199 + t122 * t206 + t131 * t109 * qJD(1) - (-t109 * t87 - t292) * t203, qJD(1) * t206 + t123 * t214 - t129 * t8 - t131 * t282 - t203 * t311; 0, 0, 0, 0, 0, 0, t191 * t195 * t187, -t263 * t195, t245, -t246, qJDD(3), t191 * t204 + t133 + t300, g(3) * t191 + (-t142 - t204) * t187 (t58 - t69) * t124 - (t57 - t70) * t307 + (t285 * t84 + t286 * t85) * pkin(3), t57 * t69 - t58 * t70 + (t285 * t26 + t286 * t25 + t300 + (-t308 - t261) * t191) * pkin(3), t226 * t96 + t289 (t36 - t322) * t190 + (-t317 * t96 - t37) * t186, -t287 + t318, t212 + t293, -t317 * t124, -t66 * t317 - t23 * t124 + t156 * t37 - t69 * t94 + (t317 * t70 + t207) * t186 - t321 * t190, t24 * t124 + t156 * t36 + t186 * t321 + t190 * t207 + t298 * t317 - t69 * t96, t8 * t132 - t214 * t296, -t131 * t8 + t132 * t199 + t214 * t295 - t296 * t40, t288 - t304, t223 + t294, -t109 * t124 -(-t127 * t189 - t128 * t185) * t80 - t140 * t199 + t10 * t131 - t4 * t124 + t222 * t40 + t295 * t32 + (t185 * t216 - t189 * t215) * t109 - t202 * t170 (-t127 * t185 + t128 * t189) * t80 + t140 * t8 + t10 * t132 + t5 * t124 - t222 * t214 + t296 * t32 + (t185 * t215 + t189 * t216) * t109 + t202 * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124 ^ 2 - t307 ^ 2, t124 * t57 + t307 * t58 + t213 - t218, 0, 0, 0, 0, 0, t212 - t293, -t287 - t318, 0, 0, 0, 0, 0, t223 - t294, t288 + t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t94, -t94 ^ 2 + t96 ^ 2, t36 + t322, -t231 + (-qJD(5) + t317) * t96, -t82, -g(1) * t115 - g(2) * t117 + t317 * t24 - t47 * t96 + (t232 + t301) * t186 + t221, g(1) * t116 - g(2) * t118 + t190 * t301 + t23 * t317 + t47 * t94 - t211, -t320, t316, t315, t312, -t80 -(-t14 * t185 - t291) * t109 + (-t109 * t254 - t189 * t80 - t40 * t96) * pkin(5) + t313 (-t109 * t15 - t2) * t185 + (t109 * t14 - t238) * t189 + (-t109 * t253 + t185 * t80 + t214 * t96) * pkin(5) + t314; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t320, t316, t315, t312, -t80, t109 * t5 + t313, t109 * t4 - t185 * t2 - t189 * t238 + t314;];
tau_reg  = t1;

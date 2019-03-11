% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x34]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:18
% EndTime: 2019-03-09 07:21:26
% DurationCPUTime: 3.62s
% Computational Cost: add. (5372->378), mult. (11492->530), div. (0->0), fcn. (8368->8), ass. (0->213)
t296 = cos(qJ(4));
t234 = t296 * qJD(4);
t194 = -pkin(1) - pkin(7);
t160 = t194 * qJD(1) + qJD(2);
t190 = sin(qJ(3));
t260 = qJD(1) * t190;
t130 = -pkin(8) * t260 + t160 * t190;
t189 = sin(qJ(4));
t121 = t189 * t130;
t193 = cos(qJ(3));
t259 = qJD(1) * t193;
t131 = -pkin(8) * t259 + t193 * t160;
t85 = t296 * t131 - t121;
t313 = pkin(3) * t234 - t85;
t241 = t296 * t193;
t134 = -qJD(1) * t241 + t189 * t260;
t182 = qJD(3) + qJD(4);
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t110 = -t134 * t188 - t192 * t182;
t210 = -t189 * t193 - t296 * t190;
t135 = t210 * qJD(1);
t311 = qJD(5) - t135;
t318 = t110 * t311;
t213 = t134 * t192 - t182 * t188;
t317 = t213 * t311;
t97 = -pkin(4) * t134 - pkin(9) * t135;
t86 = pkin(3) * t259 + t97;
t316 = -t313 * t188 - t192 * t86;
t187 = sin(qJ(6));
t191 = cos(qJ(6));
t214 = t110 * t187 + t191 * t213;
t60 = t191 * t110 - t187 * t213;
t315 = t214 * t60;
t142 = t187 * t188 - t191 * t192;
t299 = qJD(5) + qJD(6);
t281 = (t135 - t299) * t142;
t267 = t187 * t192;
t145 = t188 * t191 + t267;
t103 = t299 * t145;
t314 = -t145 * t135 + t103;
t256 = qJD(4) * t189;
t258 = qJD(3) * t190;
t312 = -t189 * t258 - t190 * t256;
t310 = t214 ^ 2 - t60 ^ 2;
t253 = qJD(6) * t187;
t122 = t296 * t130;
t123 = qJD(3) * pkin(3) + t131;
t77 = t189 * t123 + t122;
t69 = pkin(9) * t182 + t77;
t154 = pkin(3) * t260 + qJD(1) * qJ(2);
t79 = -pkin(4) * t135 + pkin(9) * t134 + t154;
t35 = t188 * t79 + t192 * t69;
t26 = -pkin(10) * t110 + t35;
t24 = t26 * t253;
t76 = t296 * t123 - t121;
t68 = -t182 * pkin(4) - t76;
t45 = t110 * pkin(5) + t68;
t309 = t45 * t60 + t24;
t129 = qJD(6) + t311;
t252 = qJD(6) * t191;
t254 = qJD(5) * t192;
t255 = qJD(5) * t188;
t95 = t135 * t182;
t53 = t134 * t255 + t182 * t254 + t192 * t95;
t54 = -t213 * qJD(5) + t188 * t95;
t13 = -t110 * t252 - t187 * t54 + t191 * t53 + t213 * t253;
t308 = t129 * t60 + t13;
t228 = pkin(8) * qJD(1) - t160;
t124 = t228 * t258;
t257 = qJD(3) * t193;
t125 = t228 * t257;
t39 = t123 * t234 + t189 * t124 - t296 * t125 - t130 * t256;
t183 = qJD(1) * qJD(2);
t249 = qJD(1) * qJD(3);
t233 = t193 * t249;
t150 = pkin(3) * t233 + t183;
t202 = t296 * qJD(3) + t234;
t201 = t202 * t193;
t263 = t312 * qJD(1);
t96 = qJD(1) * t201 + t263;
t44 = pkin(4) * t96 - pkin(9) * t95 + t150;
t42 = t192 * t44;
t199 = -t35 * qJD(5) - t188 * t39 + t42;
t4 = pkin(5) * t96 - pkin(10) * t53 + t199;
t211 = t188 * t44 + t192 * t39 + t79 * t254 - t69 * t255;
t5 = -pkin(10) * t54 + t211;
t244 = -t187 * t5 + t191 * t4;
t34 = -t188 * t69 + t192 * t79;
t25 = pkin(10) * t213 + t34;
t21 = pkin(5) * t311 + t25;
t286 = t191 * t26;
t9 = t187 * t21 + t286;
t307 = -t9 * qJD(6) + t45 * t214 + t244;
t198 = t214 * qJD(6) - t187 * t53 - t191 * t54;
t306 = -t129 * t214 + t198;
t40 = t123 * t256 - t296 * t124 - t189 * t125 + t130 * t234;
t305 = -t40 * t192 + t68 * t255;
t84 = t189 * t131 + t122;
t223 = pkin(3) * t256 - t84;
t304 = t187 * t255 + t188 * t253;
t294 = pkin(8) - t194;
t151 = t294 * t190;
t152 = t294 * t193;
t303 = t189 * t151 - t296 * t152;
t271 = t135 * t188;
t302 = (t255 - t271) * pkin(5);
t301 = -qJD(6) * t192 - t254;
t300 = t188 * t86 - t313 * t192;
t298 = 0.2e1 * t183;
t297 = -pkin(9) - pkin(10);
t295 = t192 * pkin(5);
t174 = pkin(3) * t189 + pkin(9);
t293 = -pkin(10) - t174;
t292 = t188 * t97 + t192 * t76;
t290 = t142 * t96;
t289 = t145 * t96;
t287 = t188 * t96;
t285 = t192 * t96;
t283 = t53 * t188;
t71 = t96 * t210;
t109 = -t296 * t151 - t189 * t152;
t100 = t192 * t109;
t143 = t189 * t190 - t241;
t171 = t190 * pkin(3) + qJ(2);
t98 = -pkin(4) * t210 + pkin(9) * t143 + t171;
t282 = t188 * t98 + t100;
t279 = t302 + t223;
t104 = -t189 * t257 - t190 * t202 - t193 * t256;
t278 = t104 * t188;
t277 = t104 * t192;
t105 = t201 + t312;
t276 = t105 * t129;
t275 = t105 * t182;
t274 = t129 * t134;
t273 = t311 * t134;
t272 = t134 * t135;
t270 = t135 * t192;
t269 = t143 * t188;
t268 = t143 * t192;
t195 = qJD(3) ^ 2;
t266 = t195 * t190;
t265 = t195 * t193;
t196 = qJD(1) ^ 2;
t264 = t196 * qJ(2);
t262 = t190 ^ 2 - t193 ^ 2;
t261 = -t195 - t196;
t161 = pkin(3) * t257 + qJD(2);
t248 = 0.2e1 * qJD(1);
t247 = pkin(10) * t271;
t243 = qJD(5) * t297;
t239 = t143 * t254;
t232 = qJD(6) * t21 + t5;
t230 = -t188 * t76 + t192 * t97;
t229 = qJD(5) * t293;
t227 = t192 * t311;
t226 = -qJD(5) * t210 + qJD(1);
t175 = -t296 * pkin(3) - pkin(4);
t224 = -t35 * t134 + t40 * t188 + t68 * t254;
t222 = -t77 + t302;
t221 = -t134 * pkin(5) - pkin(10) * t270;
t137 = t293 * t188;
t220 = -qJD(6) * t137 - t188 * t229 - t247 + t300;
t181 = t192 * pkin(10);
t138 = t174 * t192 + t181;
t219 = qJD(6) * t138 - t192 * t229 + t221 - t316;
t155 = t297 * t188;
t218 = -qJD(6) * t155 - t188 * t243 - t247 + t292;
t156 = pkin(9) * t192 + t181;
t217 = qJD(6) * t156 - t192 * t243 + t221 + t230;
t216 = -t135 * t68 - t174 * t96;
t212 = t34 * t134 + t305;
t209 = t239 - t278;
t208 = t143 * t255 + t277;
t51 = pkin(4) * t105 - pkin(9) * t104 + t161;
t140 = t294 * t258;
t141 = qJD(3) * t152;
t56 = qJD(4) * t303 + t189 * t140 - t296 * t141;
t207 = -t109 * t255 + t188 * t51 + t192 * t56 + t98 * t254;
t17 = pkin(5) * t54 + t40;
t8 = -t187 * t26 + t191 * t21;
t206 = t8 * t134 + t17 * t142 + t314 * t45;
t205 = -t9 * t134 + t17 * t145 + t281 * t45;
t204 = t134 * t154 - t40;
t203 = t129 * t142;
t197 = -t154 * t135 - t39;
t57 = t109 * qJD(4) - t296 * t140 - t189 * t141;
t176 = -pkin(4) - t295;
t153 = t175 - t295;
t99 = t104 * t182;
t94 = t192 * t98;
t92 = t142 * t143;
t91 = t145 * t143;
t78 = -pkin(5) * t269 - t303;
t72 = t134 ^ 2 - t135 ^ 2;
t67 = -t134 * t182 - t202 * t259 - t263;
t47 = t192 * t51;
t43 = pkin(10) * t269 + t282;
t36 = -pkin(5) * t210 + pkin(10) * t268 - t109 * t188 + t94;
t29 = -t209 * pkin(5) + t57;
t23 = -t134 * t213 + t227 * t311 + t287;
t22 = -t188 * t311 ^ 2 - t110 * t134 + t285;
t20 = -t213 * t227 + t283;
t19 = t104 * t267 + (-t268 * t299 + t278) * t191 + t304 * t143;
t18 = t103 * t143 - t142 * t104;
t12 = -t129 * t314 - t60 * t134 - t290;
t11 = t281 * t129 - t134 * t214 + t289;
t10 = t209 * pkin(10) + t207;
t7 = (t53 - t318) * t192 + (-t54 + t317) * t188;
t6 = -pkin(10) * t277 + pkin(5) * t105 - t188 * t56 + t47 + (-t100 + (-pkin(10) * t143 - t98) * t188) * qJD(5);
t2 = t13 * t145 - t214 * t281;
t1 = -t13 * t142 + t145 * t198 + t214 * t314 - t281 * t60;
t3 = [0, 0, 0, 0, t298, qJ(2) * t298, -0.2e1 * t190 * t233, 0.2e1 * t262 * t249, -t266, -t265, 0, -t194 * t266 + (qJ(2) * t257 + qJD(2) * t190) * t248, -t194 * t265 + (-qJ(2) * t258 + qJD(2) * t193) * t248, -t104 * t134 - t143 * t95, t104 * t135 + t105 * t134 + t143 * t96 + t210 * t95, t99, -t275, 0, t105 * t154 - t135 * t161 - t150 * t210 + t171 * t96 - t182 * t57, t104 * t154 - t134 * t161 - t143 * t150 + t171 * t95 - t182 * t56, -t208 * t213 - t53 * t268 (-t110 * t192 + t188 * t213) * t104 + (t283 + t192 * t54 + (-t110 * t188 - t192 * t213) * qJD(5)) * t143, -t105 * t213 + t208 * t311 - t210 * t53 - t96 * t268, -t110 * t105 + t209 * t311 + t210 * t54 + t96 * t269, t105 * t311 - t71 (-t109 * t254 + t47) * t311 + t94 * t96 - (-t254 * t69 + t42) * t210 + t34 * t105 + t57 * t110 - t303 * t54 - t68 * t239 + ((-qJD(5) * t98 - t56) * t311 - t109 * t96 - (-qJD(5) * t79 - t39) * t210 - t40 * t143 + t68 * t104) * t188, -t35 * t105 + t143 * t305 - t207 * t311 + t210 * t211 - t213 * t57 + t68 * t277 - t282 * t96 - t303 * t53, t13 * t92 - t18 * t214, t13 * t91 - t18 * t60 + t19 * t214 + t198 * t92, -t105 * t214 + t129 * t18 - t13 * t210 + t92 * t96, -t105 * t60 - t129 * t19 - t198 * t210 + t91 * t96, -t71 + t276 (-t187 * t10 + t191 * t6) * t129 + (-t187 * t43 + t191 * t36) * t96 - t244 * t210 + t8 * t105 + t29 * t60 - t78 * t198 - t17 * t91 + t45 * t19 + ((-t187 * t36 - t191 * t43) * t129 + t9 * t210) * qJD(6), -t9 * t105 + t78 * t13 - t24 * t210 + t17 * t92 + t45 * t18 - t29 * t214 + (-(-qJD(6) * t43 + t6) * t129 - t36 * t96 + t4 * t210) * t187 + (-(qJD(6) * t36 + t10) * t129 - t43 * t96 + t232 * t210) * t191; 0, 0, 0, 0, -t196, -t264, 0, 0, 0, 0, 0, t261 * t190, t261 * t193, 0, 0, 0, 0, 0, qJD(1) * t135 + t99, qJD(1) * t134 - t275, 0, 0, 0, 0, 0, t188 * t71 - t104 * t110 + t143 * t54 + (-t105 * t188 - t192 * t226) * t311, t192 * t71 + t104 * t213 + t143 * t53 + (-t105 * t192 + t188 * t226) * t311, 0, 0, 0, 0, 0, -t104 * t60 - t143 * t198 - t145 * t276 + qJD(1) * t203 - ((t191 * t301 + t304) * t129 - t289) * t210, t104 * t214 + t143 * t13 + t105 * t203 + t145 * t129 * qJD(1) - (-(t187 * t301 - t188 * t252 - t191 * t255) * t129 + t290) * t210; 0, 0, 0, 0, 0, 0, t193 * t196 * t190, -t262 * t196, 0, 0, 0, -t193 * t264, t190 * t264, t272, t72, 0, t67, 0, t182 * t84 + (t135 * t259 - t182 * t256) * pkin(3) + t204, t85 * t182 + (t134 * t259 - t182 * t234) * pkin(3) + t197, t20, t7, t23, t22, t273, t175 * t54 + t216 * t188 + t223 * t110 + (-t174 * t254 + t316) * t311 + t212, t175 * t53 + t216 * t192 - t223 * t213 + (t174 * t255 + t300) * t311 + t224, t2, t1, t11, t12, t274 (t137 * t191 - t138 * t187) * t96 - t153 * t198 + t279 * t60 + (t187 * t220 - t191 * t219) * t129 + t206 -(t137 * t187 + t138 * t191) * t96 + t153 * t13 - t279 * t214 + (t187 * t219 + t191 * t220) * t129 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t72, 0, t67, 0, t182 * t77 + t204, t76 * t182 + t197, t20, t7, t23, t22, t273, -pkin(4) * t54 - t230 * t311 - t77 * t110 - t68 * t271 + (-t254 * t311 - t287) * pkin(9) + t212, -pkin(4) * t53 + t292 * t311 + t77 * t213 - t68 * t270 + (t255 * t311 - t285) * pkin(9) + t224, t2, t1, t11, t12, t274 (t155 * t191 - t156 * t187) * t96 - t176 * t198 + t222 * t60 + (t187 * t218 - t191 * t217) * t129 + t206 -(t155 * t187 + t156 * t191) * t96 + t176 * t13 - t222 * t214 + (t187 * t217 + t191 * t218) * t129 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213 * t110, -t110 ^ 2 + t213 ^ 2, t53 + t318, -t54 - t317, t96, t213 * t68 + t311 * t35 + t199, t110 * t68 + t311 * t34 - t211, -t315, t310, t308, t306, t96 -(-t187 * t25 - t286) * t129 + (-t129 * t253 + t191 * t96 + t213 * t60) * pkin(5) + t307 (-t129 * t26 - t4) * t187 + (t129 * t25 - t232) * t191 + (-t129 * t252 - t187 * t96 - t213 * t214) * pkin(5) + t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t315, t310, t308, t306, t96, t129 * t9 + t307, t129 * t8 - t187 * t4 - t191 * t232 + t309;];
tauc_reg  = t3;

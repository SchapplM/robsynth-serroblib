% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:58
% EndTime: 2019-12-31 22:36:15
% DurationCPUTime: 5.84s
% Computational Cost: add. (6738->410), mult. (17706->602), div. (0->0), fcn. (14149->10), ass. (0->206)
t191 = cos(qJ(5));
t255 = qJD(5) * t191;
t186 = cos(pkin(5));
t254 = t186 * qJD(1);
t175 = qJD(2) + t254;
t193 = cos(qJ(3));
t189 = sin(qJ(3));
t190 = sin(qJ(2));
t185 = sin(pkin(5));
t265 = qJD(1) * t185;
t245 = t190 * t265;
t225 = t189 * t245;
t125 = t193 * t175 - t225;
t126 = t189 * t175 + t193 * t245;
t188 = sin(qJ(4));
t192 = cos(qJ(4));
t83 = -t192 * t125 + t188 * t126;
t312 = t191 * t83;
t320 = -t255 - t312;
t194 = cos(qJ(2));
t264 = qJD(1) * t194;
t240 = t185 * t264;
t226 = t189 * t240;
t291 = pkin(8) + pkin(9);
t247 = qJD(3) * t291;
t251 = pkin(1) * t254;
t140 = -pkin(7) * t245 + t194 * t251;
t207 = (pkin(2) * t190 - pkin(8) * t194) * t185;
t141 = qJD(1) * t207;
t267 = t193 * t140 + t189 * t141;
t319 = pkin(9) * t226 - t189 * t247 - t267;
t229 = -t189 * t140 + t193 * t141;
t272 = t193 * t194;
t318 = t193 * t247 + (pkin(3) * t190 - pkin(9) * t272) * t265 + t229;
t166 = -qJD(3) + t240;
t158 = -qJD(4) + t166;
t253 = qJD(1) * qJD(2);
t238 = t185 * t253;
t169 = t190 * t238;
t187 = sin(qJ(5));
t210 = t188 * t125 + t192 * t126;
t256 = qJD(5) * t187;
t257 = qJD(4) * t192;
t258 = qJD(4) * t188;
t224 = t194 * t238;
t259 = qJD(3) * t193;
t97 = -qJD(3) * t225 + t175 * t259 + t193 * t224;
t261 = qJD(2) * t194;
t242 = t189 * t261;
t260 = qJD(3) * t189;
t98 = (t190 * t259 + t242) * t265 + t175 * t260;
t42 = t125 * t257 - t126 * t258 - t188 * t98 + t192 * t97;
t16 = -t158 * t255 + t187 * t169 + t191 * t42 - t210 * t256;
t14 = t16 * t187;
t211 = t187 * t158 - t191 * t210;
t317 = t320 * t211 + t14;
t271 = -qJD(5) - t83;
t43 = qJD(4) * t210 + t188 * t97 + t192 * t98;
t36 = t187 * t43;
t287 = -t255 * t271 + t36;
t316 = t210 * t211 - t271 * t312 + t287;
t310 = t187 * t271;
t38 = t191 * t43;
t64 = t191 * t158 + t187 * t210;
t315 = t210 * t64 - t271 * t310 + t38;
t17 = -qJD(5) * t211 - t191 * t169 + t187 * t42;
t314 = t16 * t191 - t187 * t17 - t211 * t310 + t320 * t64;
t172 = t190 * t251;
t143 = pkin(7) * t240 + t172;
t109 = t175 * pkin(8) + t143;
t138 = (-pkin(2) * t194 - pkin(8) * t190 - pkin(1)) * t185;
t121 = qJD(1) * t138;
t75 = t193 * t109 + t189 * t121;
t57 = t125 * pkin(9) + t75;
t284 = t188 * t57;
t74 = -t189 * t109 + t193 * t121;
t56 = -t126 * pkin(9) + t74;
t53 = -t166 * pkin(3) + t56;
t22 = t192 * t53 - t284;
t20 = t158 * pkin(4) - t22;
t313 = t20 * t83;
t311 = t210 * t83;
t150 = t188 * t189 - t192 * t193;
t300 = qJD(3) + qJD(4);
t101 = t300 * t150;
t111 = t150 * t240;
t270 = -t101 + t111;
t151 = t188 * t193 + t192 * t189;
t269 = (-t240 + t300) * t151;
t309 = t210 ^ 2 - t83 ^ 2;
t51 = pkin(4) * t210 + t83 * pkin(10);
t308 = -t83 * t158 + t42;
t142 = qJD(2) * t207;
t133 = qJD(1) * t142;
t274 = t185 * t190;
t176 = pkin(7) * t274;
t289 = pkin(1) * t194;
t144 = (t186 * t289 - t176) * qJD(2);
t134 = qJD(1) * t144;
t198 = -qJD(3) * t75 + t193 * t133 - t189 * t134;
t30 = pkin(3) * t169 - t97 * pkin(9) + t198;
t205 = -t109 * t260 + t121 * t259 + t189 * t133 + t193 * t134;
t35 = -t98 * pkin(9) + t205;
t235 = -t188 * t30 - t192 * t35 - t53 * t257 + t57 * t258;
t108 = -t175 * pkin(2) - t140;
t87 = -t125 * pkin(3) + t108;
t307 = t87 * t83 + t235;
t179 = t188 * pkin(3) + pkin(10);
t304 = (t126 * pkin(3) + qJD(5) * t179 + t51) * t271;
t303 = t271 * t210;
t167 = t291 * t189;
t168 = t291 * t193;
t114 = -t188 * t167 + t192 * t168;
t302 = qJD(4) * t114 + t319 * t188 + t318 * t192;
t209 = -t192 * t167 - t188 * t168;
t301 = qJD(4) * t209 - t318 * t188 + t319 * t192;
t222 = -t143 + (-t226 + t260) * pkin(3);
t273 = t185 * t194;
t137 = pkin(7) * t273 + (pkin(1) * t190 + pkin(8)) * t186;
t268 = t193 * t137 + t189 * t138;
t281 = t192 * t57;
t23 = t188 * t53 + t281;
t21 = -t158 * pkin(10) + t23;
t31 = t83 * pkin(4) - pkin(10) * t210 + t87;
t218 = t187 * t21 - t191 * t31;
t299 = t20 * t256 + t210 * t218;
t201 = -t23 * qJD(4) - t188 * t35 + t192 * t30;
t5 = -pkin(4) * t169 - t201;
t9 = t187 * t31 + t191 * t21;
t298 = t5 * t187 + t20 * t255 + t9 * t210;
t296 = -t87 * t210 + t201;
t181 = -t193 * pkin(3) - pkin(2);
t99 = t150 * pkin(4) - t151 * pkin(10) + t181;
t295 = t20 * qJD(5) * t151 + (-t269 * pkin(4) + t270 * pkin(10) + qJD(5) * t114 - t222) * t271 + t99 * t43;
t293 = -t158 * t210 - t43;
t147 = t186 * t189 + t193 * t274;
t230 = -t189 * t137 + t193 * t138;
t62 = -pkin(3) * t273 - t147 * pkin(9) + t230;
t146 = -t186 * t193 + t189 * t274;
t69 = -t146 * pkin(9) + t268;
t213 = t188 * t62 + t192 * t69;
t243 = t185 * t261;
t104 = -qJD(3) * t146 + t193 * t243;
t197 = -t268 * qJD(3) + t193 * t142 - t189 * t144;
t263 = qJD(2) * t190;
t244 = t185 * t263;
t44 = pkin(3) * t244 - t104 * pkin(9) + t197;
t103 = qJD(3) * t147 + t185 * t242;
t204 = -t137 * t260 + t138 * t259 + t189 * t142 + t193 * t144;
t47 = -t103 * pkin(9) + t204;
t292 = -qJD(4) * t213 - t188 * t47 + t192 * t44;
t135 = pkin(7) * t224 + qJD(2) * t172;
t73 = t98 * pkin(3) + t135;
t11 = t43 * pkin(4) - t42 * pkin(10) + t73;
t4 = pkin(10) * t169 - t235;
t1 = -t218 * qJD(5) + t187 * t11 + t191 * t4;
t285 = pkin(4) * t245 + t302;
t280 = t125 * t166;
t279 = t126 * t166;
t278 = t151 * t191;
t277 = t166 * t189;
t276 = t166 * t193;
t182 = t185 ^ 2;
t275 = t182 * qJD(1) ^ 2;
t145 = t186 * pkin(1) * t263 + pkin(7) * t243;
t266 = t190 ^ 2 - t194 ^ 2;
t262 = qJD(2) * t193;
t249 = t190 * t275;
t248 = t187 * t273;
t246 = t182 * t264;
t239 = t182 * t253;
t228 = t175 + t254;
t227 = 0.2e1 * t239;
t26 = t188 * t56 + t281;
t223 = pkin(3) * t258 - t26;
t88 = t103 * pkin(3) + t145;
t220 = -0.2e1 * pkin(1) * t239;
t33 = -pkin(10) * t273 + t213;
t94 = t192 * t146 + t188 * t147;
t95 = -t188 * t146 + t192 * t147;
t136 = t176 + (-pkin(2) - t289) * t186;
t96 = t146 * pkin(3) + t136;
t48 = t94 * pkin(4) - t95 * pkin(10) + t96;
t217 = t187 * t48 + t191 * t33;
t216 = -t187 * t33 + t191 * t48;
t214 = -t188 * t69 + t192 * t62;
t79 = t187 * t95 + t191 * t273;
t206 = t188 * t44 + t192 * t47 + t62 * t257 - t69 * t258;
t90 = -t191 * t111 + t187 * t245;
t203 = -t191 * t101 - t151 * t256 - t90;
t2 = -t9 * qJD(5) + t191 * t11 - t187 * t4;
t27 = t192 * t56 - t284;
t199 = -t179 * t43 + t313 - (-pkin(3) * t257 + t27) * t271;
t196 = -t20 * t101 - t114 * t43 + t5 * t151 - (pkin(10) * t245 - qJD(5) * t99 - t301) * t271;
t180 = -t192 * pkin(3) - pkin(4);
t89 = -t187 * t111 - t191 * t245;
t80 = t191 * t95 - t248;
t50 = qJD(4) * t95 + t192 * t103 + t188 * t104;
t49 = -qJD(4) * t94 - t188 * t103 + t192 * t104;
t32 = pkin(4) * t273 - t214;
t25 = -qJD(5) * t248 + t187 * t49 - t191 * t244 + t95 * t255;
t24 = -qJD(5) * t79 + t187 * t244 + t191 * t49;
t12 = t50 * pkin(4) - t49 * pkin(10) + t88;
t7 = -pkin(4) * t244 - t292;
t6 = pkin(10) * t244 + t206;
t3 = [0, 0, 0, t190 * t194 * t227, -t266 * t227, t228 * t243, -t228 * t244, 0, -t135 * t186 - t145 * t175 + t190 * t220, -t134 * t186 - t144 * t175 + t194 * t220, t126 * t104 + t97 * t147, -t126 * t103 + t104 * t125 - t97 * t146 - t147 * t98, -t104 * t166 + (-t194 * t97 + (qJD(1) * t147 + t126) * t263) * t185, t103 * t166 + (t194 * t98 + (-qJD(1) * t146 + t125) * t263) * t185, (-t166 * t185 - t246) * t263, -t197 * t166 - t145 * t125 + t136 * t98 + t135 * t146 + t108 * t103 + (-t198 * t194 + (qJD(1) * t230 + t74) * t263) * t185, t204 * t166 + t145 * t126 + t136 * t97 + t135 * t147 + t108 * t104 + (t205 * t194 + (-qJD(1) * t268 - t75) * t263) * t185, t210 * t49 + t42 * t95, -t210 * t50 - t42 * t94 - t95 * t43 - t49 * t83, -t49 * t158 + (-t194 * t42 + (qJD(1) * t95 + t210) * t263) * t185, t50 * t158 + (t194 * t43 + (-qJD(1) * t94 - t83) * t263) * t185, (-t158 * t185 - t246) * t263, -t292 * t158 + t88 * t83 + t96 * t43 + t73 * t94 + t87 * t50 + (-t201 * t194 + (qJD(1) * t214 + t22) * t263) * t185, t206 * t158 + t88 * t210 + t96 * t42 + t73 * t95 + t87 * t49 + (-t235 * t194 + (-qJD(1) * t213 - t23) * t263) * t185, t16 * t80 - t211 * t24, -t16 * t79 - t80 * t17 + t211 * t25 - t24 * t64, t16 * t94 - t211 * t50 - t24 * t271 + t80 * t43, -t17 * t94 + t25 * t271 - t79 * t43 - t64 * t50, -t271 * t50 + t43 * t94, -(-qJD(5) * t217 + t191 * t12 - t187 * t6) * t271 + t216 * t43 + t2 * t94 - t218 * t50 + t7 * t64 + t32 * t17 + t5 * t79 + t20 * t25, (qJD(5) * t216 + t187 * t12 + t191 * t6) * t271 - t217 * t43 - t1 * t94 - t9 * t50 - t7 * t211 + t32 * t16 + t5 * t80 + t20 * t24; 0, 0, 0, -t194 * t249, t266 * t275, (qJD(2) - t175) * t240, t175 * t245 - t169, 0, pkin(1) * t249 + t143 * t175 - t135, pkin(7) * t169 + t140 * t175 + (-t186 * t253 + t275) * t289, -t126 * t276 + t97 * t189, (t97 - t280) * t193 + (-t98 + t279) * t189, -t166 * t259 + (t166 * t272 + (qJD(2) * t189 - t126) * t190) * t265, t166 * t260 + (-t194 * t277 + (-t125 + t262) * t190) * t265, t166 * t245, -pkin(2) * t98 - t135 * t193 + t229 * t166 + t143 * t125 + (pkin(8) * t276 + t108 * t189) * qJD(3) + (-t74 * t190 + (-pkin(8) * t263 - t108 * t194) * t189) * t265, -pkin(2) * t97 + t135 * t189 - t267 * t166 - t143 * t126 + (-pkin(8) * t277 + t108 * t193) * qJD(3) + (-t108 * t272 + (-pkin(8) * t262 + t75) * t190) * t265, t42 * t151 + t210 * t270, -t42 * t150 - t151 * t43 - t210 * t269 - t270 * t83, -t270 * t158 + (qJD(2) * t151 - t210) * t245, t269 * t158 + (-qJD(2) * t150 + t83) * t245, t158 * t245, t73 * t150 + t181 * t43 + t269 * t87 + t222 * t83 + t302 * t158 + (qJD(2) * t209 - t22) * t245, t73 * t151 + t181 * t42 + t270 * t87 + t222 * t210 + t301 * t158 + (-qJD(2) * t114 + t23) * t245, t16 * t278 - t203 * t211, t90 * t64 - t211 * t89 - (t187 * t211 - t191 * t64) * t101 + (-t14 - t17 * t191 + (t187 * t64 + t191 * t211) * qJD(5)) * t151, t16 * t150 - t203 * t271 - t211 * t269 + t278 * t43, -t151 * t36 - t17 * t150 - t269 * t64 - (t187 * t101 - t151 * t255 + t89) * t271, t43 * t150 - t269 * t271, t2 * t150 - t17 * t209 + t196 * t187 + t295 * t191 - t20 * t89 - t218 * t269 + t285 * t64, -t1 * t150 - t16 * t209 - t295 * t187 + t196 * t191 - t20 * t90 - t211 * t285 - t269 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t125, -t125 ^ 2 + t126 ^ 2, t97 + t280, -t279 - t98, t169, -t108 * t126 - t75 * t166 + t198, -t108 * t125 - t74 * t166 - t205, t311, t309, t308, t293, t169, -t26 * t158 + (-t126 * t83 + t158 * t258 + t169 * t192) * pkin(3) + t296, -t27 * t158 + (-t126 * t210 + t158 * t257 - t169 * t188) * pkin(3) + t307, t317, t314, t316, t315, t303, t180 * t17 + t223 * t64 + (-t5 + t304) * t191 + t199 * t187 + t299, t180 * t16 - t187 * t304 + t191 * t199 - t211 * t223 + t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, t309, t308, t293, t169, -t23 * t158 + t296, -t22 * t158 + t307, t317, t314, t316, t315, t303, -pkin(4) * t17 - t5 * t191 + (-t187 * t22 + t191 * t51) * t271 - t23 * t64 + t187 * t313 - t287 * pkin(10) + t299, -pkin(4) * t16 - (t187 * t51 + t191 * t22) * t271 + t23 * t211 + t20 * t312 + (-t256 * t271 - t38) * pkin(10) + t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211 * t64, t211 ^ 2 - t64 ^ 2, -t271 * t64 + t16, t211 * t271 - t17, t43, t20 * t211 - t271 * t9 + t2, t20 * t64 + t218 * t271 - t1;];
tauc_reg = t3;

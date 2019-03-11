% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:34
% EndTime: 2019-03-09 05:02:43
% DurationCPUTime: 3.92s
% Computational Cost: add. (4672->362), mult. (11223->513), div. (0->0), fcn. (7898->10), ass. (0->205)
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t253 = t200 * qJD(3);
t198 = sin(qJ(3));
t263 = qJD(1) * t198;
t162 = t197 * t263 - t253;
t260 = qJD(3) * t197;
t164 = t200 * t263 + t260;
t192 = sin(pkin(11));
t194 = cos(pkin(11));
t105 = t162 * t192 - t164 * t194;
t196 = sin(qJ(6));
t199 = cos(qJ(6));
t254 = qJD(6) * t196;
t219 = -t162 * t194 - t164 * t192;
t274 = t199 * t219;
t250 = qJD(3) * qJD(4);
t256 = qJD(4) * t198;
t240 = t197 * t256;
t201 = cos(qJ(3));
t241 = t201 * t253;
t307 = -t240 + t241;
t114 = t307 * qJD(1) + t200 * t250;
t258 = qJD(3) * t201;
t242 = t197 * t258;
t255 = qJD(4) * t200;
t308 = t198 * t255 + t242;
t115 = qJD(1) * t308 + t197 * t250;
t69 = -t114 * t192 - t115 * t194;
t70 = t114 * t194 - t115 * t192;
t10 = qJD(6) * t274 + t105 * t254 + t196 * t69 + t199 * t70;
t262 = qJD(1) * t201;
t180 = -qJD(4) + t262;
t176 = -qJD(6) + t180;
t53 = t105 * t196 + t274;
t293 = t176 * t53;
t320 = t10 + t293;
t306 = -t199 * t105 + t196 * t219;
t319 = t306 * t53;
t318 = t306 ^ 2 - t53 ^ 2;
t182 = sin(pkin(10)) * pkin(1) + pkin(7);
t170 = t182 * qJD(1);
t188 = t198 * qJD(2);
t136 = t201 * t170 + t188;
t121 = qJD(3) * pkin(8) + t136;
t184 = -cos(pkin(10)) * pkin(1) - pkin(2);
t152 = -pkin(3) * t201 - pkin(8) * t198 + t184;
t126 = t152 * qJD(1);
t278 = t197 * t126;
t81 = t121 * t200 + t278;
t63 = -qJ(5) * t162 + t81;
t292 = t194 * t63;
t80 = -t121 * t197 + t200 * t126;
t62 = -qJ(5) * t164 + t80;
t55 = -pkin(4) * t180 + t62;
t23 = t192 * t55 + t292;
t305 = pkin(9) * t219;
t16 = t23 + t305;
t15 = t16 * t254;
t251 = qJD(1) * qJD(3);
t185 = t198 * t251;
t303 = qJD(2) * t201 - t198 * t170;
t124 = t303 * qJD(3);
t223 = pkin(3) * t198 - pkin(8) * t201;
t167 = t223 * qJD(3);
t151 = qJD(1) * t167;
t229 = t197 * t124 - t200 * t151;
t204 = -qJD(4) * t81 - t229;
t24 = pkin(4) * t185 - t114 * qJ(5) - t164 * qJD(5) + t204;
t248 = t200 * t124 + t126 * t255 + t197 * t151;
t257 = qJD(4) * t197;
t209 = -t121 * t257 + t248;
t30 = -qJ(5) * t115 - qJD(5) * t162 + t209;
t6 = -t192 * t30 + t194 * t24;
t4 = pkin(5) * t185 - pkin(9) * t70 + t6;
t120 = -qJD(3) * pkin(3) - t303;
t99 = pkin(4) * t162 + qJD(5) + t120;
t56 = -pkin(5) * t219 + t99;
t317 = -t196 * t4 - t56 * t53 + t15;
t11 = qJD(6) * t306 + t196 * t70 - t199 * t69;
t290 = t306 * t176;
t315 = -t11 - t290;
t299 = -qJ(5) - pkin(8);
t232 = qJD(4) * t299;
t243 = t197 * t262;
t252 = t200 * qJD(5);
t166 = t223 * qJD(1);
t269 = t197 * t166 + t200 * t303;
t314 = qJ(5) * t243 + t197 * t232 + t252 - t269;
t273 = t200 * t201;
t214 = pkin(4) * t198 - qJ(5) * t273;
t228 = t200 * t166 - t197 * t303;
t313 = -qJD(1) * t214 - t197 * qJD(5) + t200 * t232 - t228;
t7 = t192 * t24 + t194 * t30;
t5 = pkin(9) * t69 + t7;
t244 = -t196 * t5 + t199 * t4;
t312 = -t56 * t306 + t244;
t311 = pkin(9) * t105;
t156 = t192 * t200 + t194 * t197;
t210 = t156 * t201;
t310 = qJD(1) * t210 - t156 * qJD(4);
t218 = t192 * t197 - t194 * t200;
t309 = t180 * t218;
t295 = -t192 * t314 + t194 * t313;
t294 = t192 * t313 + t194 * t314;
t302 = -t136 + (-t243 + t257) * pkin(4);
t190 = t198 ^ 2;
t217 = qJD(1) * t190 - t180 * t201;
t301 = -t180 * t240 - t217 * t253;
t300 = pkin(4) * t192;
t133 = t156 * t198;
t134 = t218 * t198;
t221 = -t199 * t133 + t134 * t196;
t88 = -t133 * t196 - t134 * t199;
t89 = -qJD(3) * t210 + t218 * t256;
t90 = t156 * t256 + t192 * t242 - t194 * t241;
t27 = qJD(6) * t88 - t196 * t90 - t199 * t89;
t298 = t27 * t176 + t185 * t221;
t165 = t182 * t273;
t259 = qJD(3) * t198;
t279 = t182 * t197;
t267 = t200 * t167 + t259 * t279;
t40 = -t198 * t252 + t214 * qJD(3) + (-t165 + (qJ(5) * t198 - t152) * t197) * qJD(4) + t267;
t268 = t152 * t255 + t197 * t167;
t275 = t198 * t200;
t44 = (-qJ(5) * qJD(4) - qJD(3) * t182) * t275 + (-qJD(5) * t198 + (-qJ(5) * qJD(3) - qJD(4) * t182) * t201) * t197 + t268;
t14 = t192 * t40 + t194 * t44;
t220 = -t156 * t196 - t199 * t218;
t297 = qJD(6) * t220 + t310 * t196 + t309 * t199;
t104 = t156 * t199 - t196 * t218;
t296 = qJD(6) * t104 + t309 * t196 - t310 * t199;
t58 = t192 * t63;
t29 = t194 * t62 - t58;
t140 = t200 * t152;
t91 = -qJ(5) * t275 + t140 + (-pkin(4) - t279) * t201;
t266 = t197 * t152 + t165;
t277 = t197 * t198;
t98 = -qJ(5) * t277 + t266;
t43 = t192 * t91 + t194 * t98;
t22 = t194 * t55 - t58;
t12 = -pkin(5) * t180 + t22 + t311;
t291 = t199 * t12;
t289 = -t310 * pkin(5) + t302;
t288 = t114 * t197;
t287 = t120 * t197;
t286 = t120 * t200;
t125 = qJD(3) * t188 + t170 * t258;
t285 = t125 * t197;
t284 = t125 * t200;
t283 = t162 * t180;
t282 = t162 * t198;
t281 = t164 * t180;
t280 = t180 * t200;
t276 = t197 * t201;
t272 = t201 * t115;
t202 = qJD(3) ^ 2;
t271 = t202 * t198;
t270 = t202 * t201;
t172 = t299 * t197;
t173 = t299 * t200;
t110 = t192 * t172 - t194 * t173;
t265 = pkin(4) * t277 + t198 * t182;
t264 = -t201 ^ 2 + t190;
t171 = qJD(1) * t184;
t246 = t308 * pkin(4) + t182 * t258;
t245 = -pkin(4) * t200 - pkin(3);
t238 = t180 * t255;
t235 = qJD(6) * t12 + t5;
t234 = -t10 * t201 + t259 * t306;
t13 = -t192 * t44 + t194 * t40;
t28 = -t192 * t62 - t292;
t42 = -t192 * t98 + t194 * t91;
t231 = -t114 * t201 + t164 * t259;
t230 = t180 * t182 + t121;
t109 = t194 * t172 + t173 * t192;
t226 = t198 * t238;
t82 = pkin(4) * t115 + t125;
t95 = -pkin(9) * t218 + t110;
t225 = pkin(5) * t263 + t309 * pkin(9) + qJD(6) * t95 - t295;
t94 = -pkin(9) * t156 + t109;
t224 = t310 * pkin(9) + qJD(6) * t94 + t294;
t2 = t196 * t12 + t199 * t16;
t33 = -pkin(5) * t201 + pkin(9) * t134 + t42;
t34 = -pkin(9) * t133 + t43;
t222 = t196 * t33 + t199 * t34;
t215 = 0.2e1 * qJD(3) * t171;
t183 = pkin(4) * t194 + pkin(5);
t213 = t183 * t196 + t199 * t300;
t212 = t183 * t199 - t196 * t300;
t211 = t11 * t201 + t259 * t53;
t208 = t217 * t197;
t26 = qJD(6) * t221 + t196 * t89 - t199 * t90;
t207 = t176 * t26 - t185 * t88;
t203 = qJD(1) ^ 2;
t127 = pkin(5) * t218 + t245;
t100 = pkin(5) * t133 + t265;
t75 = pkin(4) * t164 - pkin(5) * t105;
t57 = -pkin(5) * t89 + t246;
t37 = -pkin(5) * t69 + t82;
t18 = t29 + t311;
t17 = t28 - t305;
t9 = pkin(9) * t89 + t14;
t8 = pkin(5) * t259 + pkin(9) * t90 + t13;
t1 = -t16 * t196 + t291;
t3 = [0, 0, 0, 0, 0.2e1 * t201 * t185, -0.2e1 * t264 * t251, t270, -t271, 0, -t182 * t270 + t198 * t215, t182 * t271 + t201 * t215, t114 * t275 + t307 * t164 (-t162 * t200 - t164 * t197) * t258 + (-t288 - t115 * t200 + (t162 * t197 - t164 * t200) * qJD(4)) * t198, t231 - t301, t226 + t272 + (-t208 - t282) * qJD(3) (-t180 - t262) * t259 -(-t152 * t257 + t267) * t180 + ((t162 * t182 + t287) * qJD(3) + (t200 * t230 + t278) * qJD(4) + t229) * t201 + (t120 * t255 + t182 * t115 + t285 + ((-t182 * t276 + t140) * qJD(1) + t80) * qJD(3)) * t198, t268 * t180 + (-t230 * t257 + (t164 * t182 + t286) * qJD(3) + t248) * t201 + (-t120 * t257 + t182 * t114 + t284 + (-qJD(1) * t266 - t182 * t280 - t81) * qJD(3)) * t198, t105 * t13 - t133 * t7 + t134 * t6 + t14 * t219 + t22 * t90 + t23 * t89 - t42 * t70 + t43 * t69, t22 * t13 + t23 * t14 + t246 * t99 + t265 * t82 + t6 * t42 + t7 * t43, t10 * t88 + t26 * t306, t10 * t221 - t11 * t88 + t26 * t53 - t27 * t306, -t207 + t234, t211 + t298 (-t176 - t262) * t259 -(-t196 * t9 + t199 * t8) * t176 - t244 * t201 - t57 * t53 + t100 * t11 - t37 * t221 + t56 * t27 + (t176 * t222 + t2 * t201) * qJD(6) + ((-t196 * t34 + t199 * t33) * qJD(1) + t1) * t259, t100 * t10 - t15 * t201 + t56 * t26 + t37 * t88 + t57 * t306 + ((-qJD(6) * t34 + t8) * t176 + t4 * t201) * t196 + ((qJD(6) * t33 + t9) * t176 + t235 * t201) * t199 + (-qJD(1) * t222 - t2) * t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271, -t270, 0, 0, 0, 0, 0, t226 - t272 + (-t208 + t282) * qJD(3), t231 + t301, t105 * t89 + t133 * t70 - t134 * t69 - t219 * t90, -t133 * t6 - t134 * t7 - t201 * t82 + t22 * t89 - t23 * t90 + t259 * t99, 0, 0, 0, 0, 0, -t211 + t298, t207 + t234; 0, 0, 0, 0, -t198 * t203 * t201, t264 * t203, 0, 0, 0, qJD(3) * t136 - t171 * t263 - t125, -t171 * t262, -t164 * t280 + t288 (t114 + t283) * t200 + (-t115 + t281) * t197, -t238 + (t180 * t273 + (-t164 + t260) * t198) * qJD(1), t180 * t257 + (-t180 * t276 + (t162 + t253) * t198) * qJD(1), t180 * t263, -pkin(3) * t115 - t284 + t228 * t180 - t136 * t162 + (pkin(8) * t280 + t287) * qJD(4) + (-t80 * t198 + (-pkin(8) * t259 - t120 * t201) * t197) * qJD(1), -pkin(3) * t114 + t285 - t269 * t180 - t136 * t164 + (-pkin(8) * t180 * t197 + t286) * qJD(4) + (-t120 * t273 + (-pkin(8) * t253 + t81) * t198) * qJD(1), t295 * t105 - t109 * t70 + t110 * t69 - t6 * t156 - t7 * t218 + t294 * t219 - t309 * t22 + t310 * t23, t6 * t109 + t7 * t110 + t295 * t22 + t294 * t23 + t82 * t245 + t302 * t99, t10 * t104 + t297 * t306, t10 * t220 - t104 * t11 - t296 * t306 + t297 * t53, -t297 * t176 + (qJD(3) * t104 - t306) * t263, t296 * t176 + (qJD(3) * t220 - t53) * t263, t176 * t263, -t37 * t220 + t127 * t11 + t296 * t56 - t289 * t53 + (t196 * t224 + t199 * t225) * t176 + ((-t196 * t95 + t199 * t94) * qJD(3) - t1) * t263, t127 * t10 + t37 * t104 + t297 * t56 + t289 * t306 + (-t196 * t225 + t199 * t224) * t176 + (-(t196 * t94 + t199 * t95) * qJD(3) + t2) * t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164 * t162, -t162 ^ 2 + t164 ^ 2, t114 - t283, -t115 - t281, t185, -t120 * t164 - t81 * t180 + t204, t120 * t162 - t180 * t80 - t209 (t192 * t69 - t194 * t70) * pkin(4) + (-t29 + t22) * t219 + (-t23 - t28) * t105, -t22 * t28 - t23 * t29 + (-t164 * t99 + t192 * t7 + t194 * t6) * pkin(4), -t319, t318, t320, t315, t185, t212 * t185 + (t17 * t199 - t18 * t196) * t176 + t75 * t53 + (t176 * t213 - t2) * qJD(6) + t312, -t213 * t185 - t199 * t5 - (t17 * t196 + t18 * t199) * t176 - t75 * t306 + (t176 * t212 - t291) * qJD(6) + t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105 ^ 2 - t219 ^ 2, -t105 * t22 - t219 * t23 + t82, 0, 0, 0, 0, 0, t11 - t290, t10 - t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t319, t318, t320, t315, t185 (-qJD(6) - t176) * t2 + t312, -t1 * t176 - t199 * t235 + t317;];
tauc_reg  = t3;

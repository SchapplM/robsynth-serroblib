% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRR10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:21:56
% EndTime: 2019-03-09 19:22:15
% DurationCPUTime: 7.32s
% Computational Cost: add. (7552->473), mult. (16903->796), div. (0->0), fcn. (14669->8), ass. (0->220)
t147 = sin(qJ(3));
t149 = cos(qJ(3));
t284 = sin(qJ(5));
t227 = qJD(5) * t284;
t286 = cos(qJ(5));
t228 = qJD(5) * t286;
t299 = t147 * t227 + t149 * t228;
t150 = cos(qJ(2));
t278 = t149 * pkin(2);
t280 = t147 * pkin(7);
t197 = pkin(3) + t278 + t280;
t148 = sin(qJ(2));
t287 = pkin(8) - pkin(9);
t292 = t287 * t148 + pkin(1);
t170 = t292 * t149 + (pkin(4) + t197) * t150;
t271 = t149 * t150;
t130 = pkin(7) * t271;
t279 = t148 * pkin(8);
t282 = pkin(2) * t150;
t209 = -t279 - t282;
t198 = pkin(1) - t209;
t275 = -t147 * t198 + t130;
t68 = -qJ(4) * t150 + t275;
t190 = t147 * t148 * pkin(9) + t68;
t28 = t286 * t170 - t284 * t190;
t261 = qJD(3) * t148;
t238 = qJ(4) * t261;
t263 = qJD(2) * t150;
t244 = pkin(7) * t263;
t231 = t149 * t263;
t259 = qJD(4) * t149;
t270 = qJ(4) * t231 + t148 * t259;
t137 = qJD(3) * t149;
t90 = t148 * t137 + t147 * t263;
t40 = pkin(3) * t90 + t147 * t238 + t244 - t270;
t132 = t284 * t147;
t191 = t286 * t149 + t132;
t294 = t191 * qJD(3);
t298 = -t294 + t299;
t193 = -t282 - t292;
t273 = t147 * t150;
t251 = pkin(7) * t273;
t257 = t150 * qJD(4);
t277 = t149 * pkin(7);
t283 = pkin(2) * t148;
t161 = -t257 + (t193 * t149 - t251) * qJD(3) + ((qJ(4) - t277) * t148 + (-t287 * t150 + t283) * t147) * qJD(2);
t281 = pkin(8) * t150;
t208 = -t281 + t283;
t138 = qJD(2) * t148;
t235 = t147 * t138;
t264 = qJD(2) * t149;
t276 = pkin(7) * t235 + t208 * t264;
t288 = pkin(3) + pkin(4);
t162 = (t193 * t147 + t130) * qJD(3) + (-pkin(9) * t271 - t288 * t148) * qJD(2) - t276;
t29 = t284 * t170 + t286 * t190;
t152 = t29 * qJD(5) + t284 * t161 - t286 * t162;
t229 = qJD(2) * t284;
t230 = qJD(2) * t286;
t187 = t299 * t148 + t229 * t271 - t230 * t273;
t164 = t148 * t294 - t187;
t139 = t147 * qJ(4);
t296 = -t149 * pkin(3) - t139;
t142 = t147 ^ 2;
t144 = t149 ^ 2;
t266 = t142 - t144;
t213 = t149 * t227;
t240 = t284 * t149;
t269 = -qJD(3) * t240 - t147 * t228;
t295 = t213 + t269;
t180 = t149 * t198 + t251;
t172 = t180 * qJD(3);
t106 = t266 * qJD(3);
t293 = qJD(5) + qJD(6);
t219 = t287 * t286;
t102 = t149 * t219;
t291 = (-t287 * t132 - t102) * qJD(3);
t105 = -pkin(2) + t296;
t272 = t148 * t149;
t128 = qJ(4) * t272;
t81 = -t128 + (pkin(3) * t147 + pkin(7)) * t148;
t262 = qJD(3) * t147;
t267 = qJ(4) * t137 + t147 * qJD(4);
t85 = pkin(3) * t262 - t267;
t290 = qJD(2) * (-t105 * t150 + t279) - qJD(3) * t81 - t148 * t85;
t42 = -qJD(3) * t275 + t276;
t104 = t286 * qJ(4) - t284 * t288;
t169 = t284 * qJD(4) + t104 * qJD(5);
t289 = 0.2e1 * qJD(4);
t285 = cos(qJ(6));
t146 = sin(qJ(6));
t274 = t146 * t104;
t217 = t284 * t287;
t65 = t147 * t217 + t102;
t143 = t148 ^ 2;
t265 = -t150 ^ 2 + t143;
t260 = qJD(3) * t150;
t258 = qJD(6) * t146;
t256 = -0.2e1 * pkin(1) * qJD(2);
t255 = -0.2e1 * pkin(2) * qJD(3);
t218 = t286 * t288;
t188 = t284 * qJ(4) + t218;
t184 = -pkin(5) - t188;
t175 = t285 * t184;
t82 = qJ(4) * t227 - t286 * qJD(4) + qJD(5) * t218;
t254 = -qJD(6) * t175 + t146 * t169 + t285 * t82;
t200 = qJD(5) * t219;
t201 = t149 * t217;
t242 = t286 * t147;
t212 = qJD(3) * t242;
t253 = -qJD(3) * t201 - t147 * t200 + t287 * t212;
t252 = pkin(8) * t273;
t250 = pkin(7) * t272;
t249 = pkin(3) * t138;
t248 = pkin(5) * t138;
t247 = pkin(8) * t262;
t246 = pkin(8) * t137;
t245 = pkin(5) * t258;
t243 = t288 * t147;
t241 = t285 * t104;
t239 = t147 * t260;
t236 = t149 * t260;
t233 = t147 * t137;
t232 = t148 * t263;
t226 = qJD(6) * t285;
t224 = -t146 * t82 + t285 * t169;
t223 = t146 * t191;
t222 = t265 * qJD(2);
t221 = 0.2e1 * t232;
t95 = t149 * pkin(4) - t105;
t220 = pkin(5) * t226;
t216 = t147 * t231;
t215 = t143 * t233;
t214 = t285 * t286;
t211 = -pkin(7) - t243;
t210 = t147 * pkin(2) - t277;
t207 = t285 * t191;
t70 = -t149 * (-pkin(1) - t279) + t197 * t150;
t205 = -t147 * t68 + t149 * t70;
t204 = -t147 * t275 + t149 * t180;
t199 = qJD(5) * t217;
t114 = t148 * t240;
t196 = -t148 * t242 + t114;
t34 = t191 * t263 + (-t212 - t295) * t148;
t151 = -t34 * pkin(10) - t152 - t248;
t8 = -t28 * qJD(5) - t286 * t161 - t284 * t162;
t157 = pkin(10) * t164 - t8;
t83 = t191 * t148;
t160 = t150 * pkin(5) - t83 * pkin(10) + t28;
t158 = t285 * t160;
t18 = -t196 * pkin(10) + t29;
t1 = -qJD(6) * t158 - t146 * t151 - t285 * t157 + t18 * t258;
t195 = t147 * t208;
t194 = (-t286 * pkin(5) - t288) * t147;
t192 = t212 + t269;
t189 = t285 * t196;
t99 = t146 * t286 + t285 * t284;
t183 = t296 * qJD(3) + t259;
t46 = -t146 * t196 + t285 * t83;
t64 = t147 * t219 - t201;
t176 = -t213 - t192;
t37 = -t257 - t172 + (-t252 + (qJ(4) + t210) * t148) * qJD(2);
t39 = -t249 - t42;
t174 = t205 * qJD(3) + t147 * t39 + t149 * t37;
t41 = t172 + (-t195 + t250) * qJD(2);
t173 = t204 * qJD(3) - t147 * t42 - t149 * t41;
t100 = t242 - t240;
t171 = -t100 * pkin(10) + t64;
t57 = t146 * t184 + t241;
t168 = t146 * t171;
t167 = t285 * t171;
t163 = t192 * pkin(10) + (t284 * pkin(10) - t217) * qJD(5) * t149 - t253;
t159 = t146 * t160;
t156 = -t146 * t157 + t285 * t151;
t11 = t285 * t18 + t159;
t155 = -pkin(10) * t298 + t147 * t199 + t149 * t200 + t291;
t153 = qJD(6) * t159 + t18 * t226 - t156;
t119 = -0.2e1 * t232;
t118 = -0.2e1 * t233;
t117 = 0.2e1 * t233;
t116 = pkin(8) * t236;
t98 = -t146 * t284 + t214;
t91 = t235 - t236;
t89 = t148 * t264 + t239;
t88 = -t147 * t261 + t231;
t73 = -qJD(3) * t243 + t267;
t72 = 0.2e1 * t144 * t232 - 0.2e1 * t215;
t71 = 0.2e1 * t142 * t232 + 0.2e1 * t215;
t69 = t266 * t261 - t216;
t67 = -t147 * t222 + t148 * t236;
t66 = 0.4e1 * t148 * t233 + t266 * t263;
t63 = t211 * t148 + t128;
t62 = 0.2e1 * t148 * t239 + 0.2e1 * t265 * t264;
t61 = t143 * t106 - 0.2e1 * t148 * t216;
t60 = t191 * pkin(5) + t95;
t59 = t293 * t99;
t58 = (t284 * qJD(6) + t227) * t146 - t293 * t214;
t56 = t175 - t274;
t55 = t285 * t100 - t223;
t54 = t100 * t146 + t207;
t48 = -t191 * pkin(10) + t65;
t45 = t146 * t83 + t189;
t43 = t114 * pkin(5) + t128 + (-pkin(7) + t194) * t148;
t38 = -t295 * pkin(5) + qJD(3) * t194 + t267;
t36 = t65 * qJD(5) + t291;
t35 = t149 * t199 + t253;
t33 = (-t288 * t149 - t139) * t261 + t211 * t263 + t270;
t27 = qJD(6) * t57 + t224;
t26 = t104 * t258 + t254;
t23 = t285 * t48 + t168;
t22 = -t146 * t48 + t167;
t17 = -qJD(6) * t223 + t100 * t226 - t146 * t298 + t285 * t176;
t16 = qJD(6) * t207 + t100 * t258 + t146 * t176 + t285 * t298;
t14 = -t90 * pkin(4) - pkin(5) * t164 - t40;
t13 = t46 * qJD(6) + t146 * t34 - t285 * t164;
t12 = qJD(6) * t189 - t146 * t164 + t83 * t258 - t285 * t34;
t10 = -t146 * t18 + t158;
t7 = qJD(6) * t168 + t146 * t163 + t285 * t155 + t48 * t226;
t6 = -qJD(6) * t167 + t146 * t155 - t285 * t163 + t48 * t258;
t2 = -qJD(6) * t11 + t156;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, -0.2e1 * t222, 0, t119, 0, 0, t148 * t256, t150 * t256, 0, 0, t72, 0.2e1 * t61, t62, t71, 0.2e1 * t67, t119, -0.2e1 * t180 * t138 - 0.2e1 * t150 * t42 + 0.2e1 * (t137 * t143 + t147 * t221) * pkin(7), -0.2e1 * t275 * t138 - 0.2e1 * t150 * t41 + 0.2e1 * (-t143 * t262 + t149 * t221) * pkin(7), 0.2e1 * t204 * t263 + 0.2e1 * (t147 * t41 - t149 * t42 + (-t147 * t180 - t149 * t275) * qJD(3)) * t148, 0.2e1 * pkin(7) ^ 2 * t232 - 0.2e1 * t180 * t42 - 0.2e1 * t275 * t41, t72, t62, -0.2e1 * t61, t119, -0.2e1 * t67, t71, 0.2e1 * (qJD(2) * t147 * t81 + t39) * t150 + 0.2e1 * (-qJD(2) * t70 + t137 * t81 + t40 * t147) * t148, 0.2e1 * t205 * t263 + 0.2e1 * (-t147 * t37 + t149 * t39 + (-t147 * t70 - t149 * t68) * qJD(3)) * t148, 0.2e1 * (-t264 * t81 - t37) * t150 + 0.2e1 * (qJD(2) * t68 - t40 * t149 + t262 * t81) * t148, 0.2e1 * t37 * t68 + 0.2e1 * t39 * t70 + 0.2e1 * t40 * t81, 0.2e1 * t83 * t34, 0.2e1 * t164 * t83 - 0.2e1 * t196 * t34, -0.2e1 * t138 * t83 + 0.2e1 * t150 * t34, -0.2e1 * t196 * t164, 0.2e1 * t138 * t196 + 0.2e1 * t150 * t164, t119, -0.2e1 * t152 * t150 + 0.2e1 * t33 * t114 + 0.2e1 * t63 * t187 + 0.2e1 * (-t28 * qJD(2) - t242 * t33 - t294 * t63) * t148, 0.2e1 * t138 * t29 + 0.2e1 * t150 * t8 + 0.2e1 * t33 * t83 + 0.2e1 * t63 * t34, 0.2e1 * t152 * t83 + 0.2e1 * t164 * t29 + 0.2e1 * t196 * t8 - 0.2e1 * t28 * t34, -0.2e1 * t152 * t28 - 0.2e1 * t29 * t8 + 0.2e1 * t33 * t63, -0.2e1 * t46 * t12, 0.2e1 * t12 * t45 - 0.2e1 * t13 * t46, -0.2e1 * t12 * t150 - 0.2e1 * t138 * t46, 0.2e1 * t45 * t13, -0.2e1 * t13 * t150 + 0.2e1 * t138 * t45, t119, -0.2e1 * t10 * t138 + 0.2e1 * t43 * t13 + 0.2e1 * t14 * t45 + 0.2e1 * t150 * t2, 0.2e1 * t1 * t150 + 0.2e1 * t11 * t138 - 0.2e1 * t43 * t12 + 0.2e1 * t14 * t46, 0.2e1 * t1 * t45 + 0.2e1 * t10 * t12 - 0.2e1 * t11 * t13 - 0.2e1 * t2 * t46, -0.2e1 * t1 * t11 + 0.2e1 * t10 * t2 + 0.2e1 * t14 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, 0, -t138, 0, -t244, pkin(7) * t138, 0, 0, -t69, -t66, t91, t69, t89, 0, t116 + (-t278 + t280) * t261 + (t147 * t209 - t130) * qJD(2) (t195 + t250) * qJD(3) + (t149 * t209 + t251) * qJD(2), t173, -pkin(2) * t244 + pkin(8) * t173, -t69, t91, t66, 0, -t89, t69, t116 + (t105 * t261 - t40) * t149 - t290 * t147, t174 (-t40 + (t105 * t148 + t281) * qJD(3)) * t147 + t290 * t149, pkin(8) * t174 + t105 * t40 + t81 * t85, t34 * t100 - t298 * t83, t100 * t164 - t176 * t83 - t191 * t34 + t196 * t298, -t100 * t138 - t150 * t298, -t164 * t191 + t176 * t196, t138 * t191 - t150 * t176, 0, -t138 * t64 - t36 * t150 - t164 * t95 + t176 * t63 + t191 * t33 + t196 * t73, t33 * t100 + t138 * t65 + t35 * t150 - t298 * t63 + t95 * t34 + t73 * t83, t100 * t152 + t164 * t65 - t176 * t29 + t191 * t8 + t196 * t35 + t28 * t298 - t64 * t34 + t36 * t83, -t152 * t64 - t28 * t36 - t29 * t35 + t33 * t95 + t63 * t73 - t65 * t8, -t12 * t55 - t16 * t46, t12 * t54 - t13 * t55 + t16 * t45 - t17 * t46, -t138 * t55 - t150 * t16, t13 * t54 + t17 * t45, t138 * t54 - t150 * t17, 0, t60 * t13 - t138 * t22 + t14 * t54 - t150 * t7 + t43 * t17 + t38 * t45, -t60 * t12 + t138 * t23 + t14 * t55 + t150 * t6 - t43 * t16 + t38 * t46, t1 * t54 + t10 * t16 - t11 * t17 + t12 * t22 - t13 * t23 - t2 * t55 + t45 * t6 + t46 * t7, -t1 * t23 - t10 * t7 - t11 * t6 + t14 * t60 + t2 * t22 + t38 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, -0.2e1 * t106, 0, t118, 0, 0, t147 * t255, t149 * t255, 0, 0, t117, 0, 0.2e1 * t106, 0, 0, t118, 0.2e1 * t105 * t262 - 0.2e1 * t149 * t85, 0, -0.2e1 * t105 * t137 - 0.2e1 * t147 * t85, 0.2e1 * t105 * t85, -0.2e1 * t100 * t298, -0.2e1 * t100 * t176 + 0.2e1 * t191 * t298, 0, 0.2e1 * t191 * t176, 0, 0, 0.2e1 * t176 * t95 + 0.2e1 * t191 * t73, 0.2e1 * t73 * t100 - 0.2e1 * t298 * t95, 0.2e1 * t36 * t100 - 0.2e1 * t176 * t65 + 0.2e1 * t191 * t35 + 0.2e1 * t298 * t64, -0.2e1 * t35 * t65 - 0.2e1 * t36 * t64 + 0.2e1 * t73 * t95, -0.2e1 * t55 * t16, 0.2e1 * t16 * t54 - 0.2e1 * t17 * t55, 0, 0.2e1 * t54 * t17, 0, 0, 0.2e1 * t17 * t60 + 0.2e1 * t38 * t54, -0.2e1 * t16 * t60 + 0.2e1 * t38 * t55, 0.2e1 * t16 * t22 - 0.2e1 * t17 * t23 + 0.2e1 * t54 * t6 + 0.2e1 * t55 * t7, -0.2e1 * t22 * t7 - 0.2e1 * t23 * t6 + 0.2e1 * t38 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, -t90, t138, t42, t41, 0, 0, 0, t88, 0, t138, t90, 0, t42 + 0.2e1 * t249 (-pkin(3) * t263 - t238) * t149 + (-qJ(4) * t263 + (pkin(3) * qJD(3) - qJD(4)) * t148) * t147, -0.2e1 * t257 - t172 + (-t252 + (0.2e1 * qJ(4) + t210) * t148) * qJD(2), -pkin(3) * t39 + qJ(4) * t37 + qJD(4) * t68, 0, 0, -t34, 0, -t164, t138, t138 * t188 - t150 * t169 + t152, t104 * t138 + t82 * t150 - t8, t104 * t164 + t169 * t83 + t188 * t34 + t196 * t82, -t8 * t104 + t152 * t188 - t169 * t28 - t29 * t82, 0, 0, t12, 0, t13, t138, -t138 * t56 - t27 * t150 + t153, t138 * t57 + t150 * t26 - t1, t12 * t56 - t13 * t57 + t26 * t45 + t27 * t46, -t1 * t57 - t10 * t27 - t11 * t26 + t2 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, -t262, 0, -t246, t247, 0, 0, 0, t137, 0, 0, t262, 0, -t246, t183, -t247, t183 * pkin(8), 0, 0, t298, 0, t176, 0, t36, -t35, t100 * t169 - t104 * t176 - t188 * t298 + t191 * t82, -t35 * t104 - t169 * t64 + t188 * t36 - t65 * t82, 0, 0, t16, 0, t17, 0, t7, -t6, t16 * t56 - t17 * t57 + t26 * t54 + t27 * t55, -t22 * t27 - t23 * t26 - t56 * t7 - t57 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, qJ(4) * t289, 0, 0, 0, 0, 0, 0, 0.2e1 * t169, -0.2e1 * t82, 0, -0.2e1 * t104 * t82 + 0.2e1 * t169 * t188, 0, 0, 0, 0, 0, 0, 0.2e1 * t27, -0.2e1 * t26, 0, -0.2e1 * t26 * t57 - 0.2e1 * t27 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t88, 0, t39, 0, 0, 0, 0, 0, 0, -t148 * t230 - t150 * t227, t148 * t229 - t150 * t228, t284 * t164 - t286 * t34 + (-t286 * t196 + t83 * t284) * qJD(5), -t152 * t286 - t8 * t284 + (-t284 * t28 + t286 * t29) * qJD(5), 0, 0, 0, 0, 0, 0, -t138 * t98 - t150 * t59, t138 * t99 + t150 * t58, t12 * t98 - t13 * t99 + t45 * t58 + t46 * t59, -t1 * t99 - t10 * t59 - t11 * t58 + t2 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, t246, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t227 - t284 * t176 - t191 * t228 + t286 * t298, -t36 * t286 - t35 * t284 + (-t284 * t64 + t286 * t65) * qJD(5), 0, 0, 0, 0, 0, 0, 0, 0, t16 * t98 - t17 * t99 + t54 * t58 + t55 * t59, -t22 * t59 - t23 * t58 - t6 * t99 - t7 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t228, 0, t104 * t228 - t169 * t286 + t188 * t227 - t82 * t284, 0, 0, 0, 0, 0, 0, t59, -t58, 0, -t26 * t99 - t27 * t98 - t56 * t59 - t57 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t58 * t99 - 0.2e1 * t59 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t164, -t138, -t152, t8, 0, 0, 0, 0, -t12, 0, -t13, -t138, -t150 * t245 - t285 * t248 - t153 (t138 * t146 - t150 * t226) * pkin(5) + t1 (t285 * t12 - t13 * t146 + (t146 * t46 - t285 * t45) * qJD(6)) * pkin(5) (t285 * t2 - t1 * t146 + (-t10 * t146 + t285 * t11) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t298, 0, -t176, 0, -t36, t35, 0, 0, 0, 0, -t16, 0, -t17, 0, -t7, t6 (t285 * t16 - t146 * t17 + (t146 * t55 - t285 * t54) * qJD(6)) * pkin(5) (-t285 * t7 - t146 * t6 + (-t146 * t22 + t285 * t23) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, t82, 0, 0, 0, 0, 0, 0, 0, 0 (-t241 + (0.2e1 * pkin(5) + t188) * t146) * qJD(6) - t224 (t285 * pkin(5) + t274) * qJD(6) + t254, 0 (-t285 * t27 - t146 * t26 + (-t146 * t56 + t285 * t57) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, -t228, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0 (-t285 * t59 - t146 * t58 + (-t146 * t98 + t285 * t99) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t245, -0.2e1 * t220, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t13, -t138, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, -t17, 0, -t7, t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, -t220, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;

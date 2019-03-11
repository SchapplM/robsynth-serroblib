% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:57
% EndTime: 2019-03-08 20:39:10
% DurationCPUTime: 5.40s
% Computational Cost: add. (10001->445), mult. (25333->616), div. (0->0), fcn. (20460->12), ass. (0->216)
t188 = cos(pkin(12));
t296 = cos(qJ(4));
t232 = t296 * t188;
t186 = sin(pkin(12));
t192 = sin(qJ(4));
t253 = t192 * t186;
t205 = t232 - t253;
t187 = sin(pkin(6));
t196 = cos(qJ(2));
t259 = t187 * t196;
t199 = t205 * t259;
t293 = pkin(8) + qJ(3);
t164 = t293 * t186;
t165 = t293 * t188;
t303 = -t296 * t164 - t192 * t165;
t272 = qJD(1) * t199 - t205 * qJD(3) - qJD(4) * t303;
t151 = t205 * qJD(4);
t157 = t296 * t186 + t192 * t188;
t152 = t157 * qJD(4);
t193 = sin(qJ(2));
t250 = qJD(1) * t187;
t228 = t193 * t250;
t322 = t152 * pkin(4) - t151 * pkin(9) - t228;
t177 = -t188 * pkin(3) - pkin(2);
t107 = -pkin(4) * t205 - t157 * pkin(9) + t177;
t121 = -t192 * t164 + t296 * t165;
t191 = sin(qJ(5));
t195 = cos(qJ(5));
t243 = qJD(5) * t195;
t244 = qJD(5) * t191;
t290 = t107 * t243 - t121 * t244 + t322 * t191 - t272 * t195;
t321 = t272 * t191 + t322 * t195;
t111 = t195 * t121;
t252 = t195 * t151;
t320 = -pkin(10) * t252 + t152 * pkin(5) + (-t111 + (pkin(10) * t157 - t107) * t191) * qJD(5) + t321;
t255 = t191 * t151;
t204 = t157 * t243 + t255;
t319 = -pkin(10) * t204 + t290;
t247 = qJD(2) * t157;
t126 = -t195 * qJD(4) + t191 * t247;
t128 = t191 * qJD(4) + t195 * t247;
t190 = sin(qJ(6));
t194 = cos(qJ(6));
t211 = t190 * t126 - t194 * t128;
t69 = t194 * t126 + t190 * t128;
t294 = t69 * t211;
t172 = qJD(2) * t232;
t226 = qJD(2) * t253;
t147 = -t172 + t226;
t257 = t190 * t191;
t159 = -t194 * t195 + t257;
t238 = qJD(5) + qJD(6);
t241 = qJD(6) * t194;
t275 = t159 * t147 - t194 * t243 - t195 * t241 + t238 * t257;
t160 = t190 * t195 + t194 * t191;
t113 = t238 * t160;
t274 = t160 * t147 + t113;
t299 = -pkin(10) - pkin(9);
t234 = qJD(5) * t299;
t266 = t147 * t191;
t105 = pkin(4) * t247 + t147 * pkin(9);
t163 = qJD(2) * qJ(3) + t228;
t189 = cos(pkin(6));
t249 = qJD(1) * t189;
t171 = t188 * t249;
t288 = pkin(8) * qJD(2);
t114 = t171 + (-t163 - t288) * t186;
t130 = t188 * t163 + t186 * t249;
t115 = t188 * t288 + t130;
t304 = t296 * t114 - t192 * t115;
t44 = t191 * t105 + t195 * t304;
t318 = -pkin(10) * t266 + t191 * t234 - t44;
t43 = t195 * t105 - t191 * t304;
t317 = pkin(5) * t247 + t43 + (pkin(10) * t147 - t234) * t195;
t169 = qJD(4) * t172;
t140 = qJD(4) * t226 - t169;
t316 = -qJD(4) * qJD(5) + t140;
t315 = t244 + t266;
t314 = t211 ^ 2 - t69 ^ 2;
t143 = qJD(5) + t147;
t63 = t192 * t114 + t296 * t115;
t59 = qJD(4) * pkin(9) + t63;
t248 = qJD(1) * t196;
t227 = t187 * t248;
t214 = qJD(3) - t227;
t246 = qJD(2) * t177;
t142 = t214 + t246;
t73 = t147 * pkin(4) - pkin(9) * t247 + t142;
t35 = -t191 * t59 + t195 * t73;
t30 = -t128 * pkin(10) + t35;
t23 = t143 * pkin(5) + t30;
t36 = t191 * t73 + t195 * t59;
t31 = -t126 * pkin(10) + t36;
t283 = t194 * t31;
t11 = t190 * t23 + t283;
t141 = qJD(2) * t152;
t155 = (qJD(3) + t227) * qJD(2);
t306 = t205 * t155;
t45 = qJD(4) * t304 + t306;
t168 = qJD(2) * t228;
t82 = t141 * pkin(4) + t140 * pkin(9) + t168;
t15 = -qJD(5) * t36 - t191 * t45 + t195 * t82;
t83 = t316 * t195 + t247 * t244;
t8 = t141 * pkin(5) + t83 * pkin(10) + t15;
t14 = t191 * t82 + t195 * t45 + t73 * t243 - t244 * t59;
t235 = -t316 * t191 + t247 * t243;
t9 = -pkin(10) * t235 + t14;
t2 = -qJD(6) * t11 - t190 * t9 + t194 * t8;
t58 = -qJD(4) * pkin(4) - t304;
t48 = t126 * pkin(5) + t58;
t313 = t48 * t211 + t2;
t139 = qJD(6) + t143;
t242 = qJD(6) * t190;
t24 = t126 * t241 + t128 * t242 + t190 * t235 + t194 * t83;
t312 = t69 * t139 - t24;
t1 = (qJD(6) * t23 + t9) * t194 + t190 * t8 - t31 * t242;
t311 = t48 * t69 - t1;
t198 = qJD(6) * t211 + t190 * t83 - t194 * t235;
t310 = -t139 * t211 + t198;
t212 = -t35 * t143 + t14;
t309 = t36 * t143 + t15;
t200 = t157 * t259;
t270 = -qJD(1) * t200 + qJD(3) * t157 + qJD(4) * t121;
t220 = t143 * t191;
t308 = t128 * t220;
t307 = t155 * t157;
t210 = (-t186 * t163 + t171) * t186 - t130 * t188;
t305 = t210 * t196;
t56 = t191 * t107 + t111;
t203 = -t157 * t244 + t252;
t302 = -t159 * t24 - t211 * t274;
t301 = -t139 * t275 + t160 * t141;
t300 = t247 ^ 2;
t262 = t157 * t195;
t55 = t195 * t107 - t191 * t121;
t47 = -pkin(5) * t205 - pkin(10) * t262 + t55;
t263 = t157 * t191;
t51 = -pkin(10) * t263 + t56;
t20 = -t190 * t51 + t194 * t47;
t298 = qJD(6) * t20 + t320 * t190 + t319 * t194;
t21 = t190 * t47 + t194 * t51;
t297 = -qJD(6) * t21 - t319 * t190 + t320 * t194;
t260 = t187 * t193;
t145 = -t186 * t260 + t189 * t188;
t146 = t189 * t186 + t188 * t260;
t206 = t296 * t145 - t192 * t146;
t46 = t63 * qJD(4) + t307;
t295 = t46 * t206;
t166 = t299 * t191;
t167 = t299 * t195;
t125 = t190 * t166 - t194 * t167;
t292 = qJD(6) * t125 + t318 * t190 + t317 * t194;
t124 = t194 * t166 + t190 * t167;
t291 = -qJD(6) * t124 + t317 * t190 - t318 * t194;
t289 = -t56 * qJD(5) + t321;
t287 = qJD(2) * pkin(2);
t286 = t247 * t69;
t284 = t190 * t31;
t280 = t46 * t303;
t279 = t46 * t191;
t278 = t46 * t195;
t277 = t211 * t247;
t276 = t83 * t191;
t273 = pkin(5) * t204 + t270;
t240 = t126 * qJD(5);
t79 = t191 * t235;
t271 = -t195 * t240 - t79;
t269 = t126 * t147;
t268 = t128 * t126;
t267 = t128 * t247;
t109 = t141 * t205;
t265 = t247 * t126;
t264 = t247 * t147;
t197 = qJD(2) ^ 2;
t258 = t187 * t197;
t256 = t191 * t141;
t133 = t195 * t141;
t251 = t186 ^ 2 + t188 ^ 2;
t245 = qJD(2) * t193;
t237 = t193 * t258;
t236 = t196 * t258;
t231 = t187 ^ 2 * t248;
t230 = t187 * t245;
t221 = t251 * t155;
t219 = t143 * t195;
t218 = t160 * t198 + t275 * t69;
t217 = t315 * pkin(5) - t63;
t216 = -t274 * t139 - t159 * t141;
t215 = t235 * t195;
t94 = t192 * t145 + t296 * t146;
t208 = t191 * t259 - t195 * t94;
t76 = -t191 * t94 - t195 * t259;
t40 = t190 * t208 + t194 * t76;
t41 = t190 * t76 - t194 * t208;
t213 = t36 * t191 + t35 * t195;
t209 = -t315 * t143 + t133;
t202 = -pkin(9) * t141 + t143 * t58;
t181 = -t195 * pkin(5) - pkin(4);
t158 = t214 - t287;
t144 = t147 ^ 2;
t101 = t159 * t157;
t100 = t160 * t157;
t88 = pkin(5) * t263 - t303;
t61 = qJD(2) * t200 + qJD(4) * t94;
t60 = qJD(2) * t199 + qJD(4) * t206;
t39 = -t242 * t263 + (t238 * t262 + t255) * t194 + t203 * t190;
t38 = t113 * t157 + t190 * t255 - t194 * t252;
t34 = qJD(5) * t208 - t191 * t60 + t195 * t230;
t33 = qJD(5) * t76 + t191 * t230 + t195 * t60;
t29 = pkin(5) * t235 + t46;
t13 = t194 * t30 - t284;
t12 = -t190 * t30 - t283;
t10 = t194 * t23 - t284;
t6 = -qJD(6) * t41 - t190 * t33 + t194 * t34;
t5 = qJD(6) * t40 + t190 * t34 + t194 * t33;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, -t236, 0, 0, 0, 0, 0, 0, 0, 0, -t188 * t237, t186 * t237, t251 * t236 (-t145 * t186 + t146 * t188) * t155 + (-t193 * t231 + (t158 * t193 - t305) * t187) * qJD(2), 0, 0, 0, 0, 0, 0, -t61 * qJD(4) + (-t141 * t196 + t147 * t245) * t187, -t60 * qJD(4) + (t140 * t196 + t245 * t247) * t187, t140 * t206 - t94 * t141 - t60 * t147 + t247 * t61, t45 * t94 - t295 + t63 * t60 - t304 * t61 + (t142 * t187 - t231) * t245, 0, 0, 0, 0, 0, 0, t61 * t126 + t76 * t141 + t34 * t143 - t206 * t235, t128 * t61 + t141 * t208 - t143 * t33 + t206 * t83, -t33 * t126 - t34 * t128 + t208 * t235 + t76 * t83, -t14 * t208 + t15 * t76 + t33 * t36 + t34 * t35 + t58 * t61 - t295, 0, 0, 0, 0, 0, 0, t139 * t6 + t141 * t40 + t198 * t206 + t61 * t69, -t139 * t5 - t141 * t41 + t206 * t24 - t211 * t61, t198 * t41 + t211 * t6 + t24 * t40 - t5 * t69, t1 * t41 + t10 * t6 + t11 * t5 + t2 * t40 - t206 * t29 + t48 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t214 * t251 + t221, -t210 * qJD(3) + qJ(3) * t221 + (t305 + (-t158 - t287) * t193) * t250, -t140 * t157 + t151 * t247, -t140 * t205 - t157 * t141 - t151 * t147 - t152 * t247, t151 * qJD(4), t147 * t152 - t109, -t152 * qJD(4), 0, t177 * t141 + t142 * t152 - t270 * qJD(4) + (-qJD(2) * t205 - t147) * t228, t272 * qJD(4) - t177 * t140 + t142 * t151, -t121 * t141 + t140 * t303 + t147 * t272 - t151 * t304 - t63 * t152 + t46 * t157 + t205 * t45 + t247 * t270, -t280 + t45 * t121 - t272 * t63 - t270 * t304 + (-t142 + t246) * t228, t128 * t203 - t262 * t83 -(t195 * t126 + t128 * t191) * t151 + (-t215 + t276 + (t126 * t191 - t128 * t195) * qJD(5)) * t157, t128 * t152 + t133 * t157 + t143 * t203 + t205 * t83, t126 * t204 + t157 * t79, -t126 * t152 - t143 * t204 - t157 * t256 + t205 * t235, t143 * t152 - t109, t55 * t141 - t15 * t205 + t35 * t152 - t303 * t235 + t58 * t255 + (t243 * t58 + t279) * t157 + t289 * t143 + t270 * t126, t58 * t252 + t303 * t83 + t14 * t205 - t56 * t141 - t36 * t152 + (-t244 * t58 + t278) * t157 - t290 * t143 + t270 * t128, -t56 * t235 + t55 * t83 - t213 * t151 - t289 * t128 - t290 * t126 + (-t14 * t191 - t15 * t195 + (t191 * t35 - t195 * t36) * qJD(5)) * t157, t14 * t56 + t15 * t55 + t270 * t58 + t289 * t35 + t290 * t36 - t280, t101 * t24 + t211 * t38, t100 * t24 - t101 * t198 + t211 * t39 + t38 * t69, -t101 * t141 - t38 * t139 - t152 * t211 + t205 * t24, -t100 * t198 + t39 * t69, -t100 * t141 - t39 * t139 - t69 * t152 - t198 * t205, t139 * t152 - t109, t10 * t152 + t29 * t100 + t297 * t139 + t20 * t141 - t198 * t88 - t2 * t205 + t273 * t69 + t48 * t39, t1 * t205 - t29 * t101 - t11 * t152 - t298 * t139 - t21 * t141 - t211 * t273 - t88 * t24 - t48 * t38, -t1 * t100 + t10 * t38 + t101 * t2 - t11 * t39 + t198 * t21 + t20 * t24 + t211 * t297 - t298 * t69, t1 * t21 + t297 * t10 + t298 * t11 + t2 * t20 + t273 * t48 + t29 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251 * t197, qJD(2) * t210 + t168, 0, 0, 0, 0, 0, 0, 0.2e1 * t247 * qJD(4), t169 + (-t147 - t226) * qJD(4), -t144 - t300, t63 * t147 + t247 * t304 + t168, 0, 0, 0, 0, 0, 0, t209 - t265, -t143 ^ 2 * t195 - t256 - t267 (t83 - t269) * t195 + t308 + t271, t212 * t191 + t309 * t195 - t247 * t58, 0, 0, 0, 0, 0, 0, t216 - t286, t277 - t301, t218 + t302, t1 * t160 - t10 * t274 - t11 * t275 - t2 * t159 - t247 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, -t144 + t300, t169 + (t147 - t226) * qJD(4), -t264, 0, 0, -t142 * t247 - t307, t142 * t147 - t306, 0, 0, t128 * t219 - t276 (-t83 - t269) * t195 - t308 + t271, t143 * t219 + t256 - t267, t126 * t220 - t215, t209 + t265, -t143 * t247, -pkin(4) * t235 - t278 - t35 * t247 - t63 * t126 + (-pkin(9) * t243 - t43) * t143 + t202 * t191, pkin(4) * t83 - t63 * t128 + t36 * t247 + t279 + (pkin(9) * t244 + t44) * t143 + t202 * t195, t44 * t126 + t43 * t128 + ((t128 * qJD(5) - t235) * pkin(9) + t212) * t195 + ((-t83 + t240) * pkin(9) - t309) * t191, -t46 * pkin(4) - t35 * t43 - t36 * t44 - t58 * t63 + (-qJD(5) * t213 + t14 * t195 - t15 * t191) * pkin(9), -t24 * t160 + t211 * t275, t218 - t302, t277 + t301, -t159 * t198 + t274 * t69, t216 + t286, -t139 * t247, -t10 * t247 + t124 * t141 - t139 * t292 + t29 * t159 - t181 * t198 + t217 * t69 + t274 * t48, t11 * t247 - t125 * t141 + t139 * t291 + t29 * t160 - t181 * t24 - t211 * t217 - t275 * t48, -t1 * t159 + t275 * t10 - t274 * t11 + t124 * t24 + t125 * t198 - t2 * t160 - t211 * t292 + t291 * t69, t1 * t125 - t10 * t292 - t291 * t11 + t2 * t124 + t29 * t181 + t217 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, -t126 ^ 2 + t128 ^ 2, t126 * t143 - t83, -t268, t128 * t143 - t235, t141, -t58 * t128 + t309, t58 * t126 - t212, 0, 0, -t294, t314, t312, t294, t310, t141, -t12 * t139 + (-t128 * t69 - t139 * t242 + t141 * t194) * pkin(5) + t313, t13 * t139 + (t128 * t211 - t139 * t241 - t141 * t190) * pkin(5) + t311, -t10 * t69 - t11 * t211 - t12 * t211 + t13 * t69 + (t190 * t198 + t194 * t24 + (-t190 * t211 - t194 * t69) * qJD(6)) * pkin(5), -t10 * t12 - t11 * t13 + (t1 * t190 - t128 * t48 + t194 * t2 + (-t10 * t190 + t11 * t194) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t294, t314, t312, t294, t310, t141, t11 * t139 + t313, t10 * t139 + t311, 0, 0;];
tauc_reg  = t3;

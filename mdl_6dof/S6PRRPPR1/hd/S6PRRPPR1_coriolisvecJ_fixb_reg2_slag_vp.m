% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:22
% EndTime: 2019-03-08 21:02:32
% DurationCPUTime: 6.16s
% Computational Cost: add. (7465->440), mult. (19555->621), div. (0->0), fcn. (15379->12), ass. (0->231)
t190 = sin(pkin(11));
t195 = sin(qJ(3));
t198 = cos(qJ(3));
t279 = cos(pkin(11));
t164 = t190 * t198 + t195 * t279;
t153 = t164 * qJD(2);
t189 = sin(pkin(12));
t192 = cos(pkin(12));
t131 = qJD(3) * t189 + t153 * t192;
t194 = sin(qJ(6));
t197 = cos(qJ(6));
t265 = t189 * t153;
t317 = qJD(3) * t192 - t265;
t212 = t197 * t317;
t73 = -t131 * t194 + t212;
t323 = t73 ^ 2;
t226 = t279 * t198;
t176 = qJD(2) * t226;
t248 = qJD(2) * t195;
t150 = t190 * t248 - t176;
t146 = qJD(6) + t150;
t322 = t73 * t146;
t72 = t197 * t131 + t194 * t317;
t321 = t72 ^ 2;
t196 = sin(qJ(2));
t191 = sin(pkin(6));
t251 = qJD(1) * t191;
t237 = t196 * t251;
t245 = qJD(3) * t195;
t320 = pkin(3) * t245 - t237;
t293 = -qJ(4) - pkin(8);
t228 = qJD(3) * t293;
t243 = t198 * qJD(4);
t147 = t195 * t228 + t243;
t208 = -t195 * qJD(4) + t198 * t228;
t210 = -t190 * t195 + t226;
t199 = cos(qJ(2));
t232 = t199 * t251;
t283 = -t147 * t279 - t190 * t208 + t210 * t232;
t152 = t164 * qJD(3);
t155 = t210 * qJD(3);
t319 = pkin(4) * t152 - qJ(5) * t155 - qJD(5) * t164 + t320;
t242 = qJD(2) * qJD(3);
t229 = t195 * t242;
t173 = t190 * t229;
t144 = qJD(3) * t176 - t173;
t261 = t192 * t144;
t318 = t150 * t317 + t261;
t292 = t189 * t283 + t192 * t319;
t291 = t189 * t319 - t192 * t283;
t258 = t197 * t192;
t316 = -t189 * t194 + t258;
t269 = t155 * t192;
t315 = pkin(5) * t152 - pkin(9) * t269 + t292;
t270 = t155 * t189;
t314 = -pkin(9) * t270 + t291;
t307 = t316 * qJD(6);
t282 = -t150 * t316 - t307;
t165 = t189 * t197 + t192 * t194;
t157 = t165 * qJD(6);
t281 = t150 * t165 + t157;
t274 = t144 * t189;
t311 = t131 * t150 + t274;
t280 = t147 * t190 - t164 * t232 - t208 * t279;
t310 = t144 * t165;
t309 = t153 * t317;
t143 = qJD(2) * t152;
t112 = t143 * t210;
t308 = t150 * t152 - t112;
t306 = -qJD(6) + t146;
t34 = t194 * (qJD(6) * t131 + t274) - qJD(6) * t212 - t144 * t258;
t305 = -t281 * t72 - t316 * t34;
t304 = t165 * t143 - t146 * t282;
t149 = t150 ^ 2;
t303 = -t143 * t189 - t149 * t192;
t302 = t153 ^ 2;
t296 = pkin(9) * t192;
t184 = -pkin(3) * t198 - pkin(2);
t108 = -pkin(4) * t210 - qJ(5) * t164 + t184;
t171 = t293 * t195;
t172 = t293 * t198;
t129 = t171 * t190 - t172 * t279;
t57 = t108 * t192 - t129 * t189;
t40 = -pkin(5) * t210 - t164 * t296 + t57;
t268 = t164 * t189;
t58 = t108 * t189 + t129 * t192;
t44 = -pkin(9) * t268 + t58;
t14 = -t194 * t44 + t197 * t40;
t301 = qJD(6) * t14 + t194 * t315 + t197 * t314;
t15 = t194 * t40 + t197 * t44;
t300 = -qJD(6) * t15 - t194 * t314 + t197 * t315;
t180 = pkin(3) * t190 + qJ(5);
t294 = pkin(9) + t180;
t160 = t294 * t189;
t161 = t294 * t192;
t105 = -t160 * t194 + t161 * t197;
t168 = qJD(2) * pkin(8) + t237;
t225 = qJ(4) * qJD(2) + t168;
t193 = cos(pkin(6));
t250 = qJD(1) * t195;
t231 = t193 * t250;
t124 = t198 * t225 + t231;
t114 = t190 * t124;
t260 = t193 * t198;
t177 = qJD(1) * t260;
t123 = -t195 * t225 + t177;
t65 = t123 * t279 - t114;
t241 = pkin(3) * t248;
t95 = pkin(4) * t153 + qJ(5) * t150 + t241;
t29 = -t189 * t65 + t192 * t95;
t20 = pkin(5) * t153 + t150 * t296 + t29;
t271 = t150 * t189;
t30 = t189 * t95 + t192 * t65;
t22 = pkin(9) * t271 + t30;
t299 = qJD(5) * t165 + qJD(6) * t105 - t194 * t22 + t197 * t20;
t104 = -t160 * t197 - t161 * t194;
t298 = qJD(5) * t316 + qJD(6) * t104 - t194 * t20 - t197 * t22;
t297 = pkin(3) * t195;
t295 = t72 * t73;
t249 = qJD(2) * t191;
t230 = qJD(1) * t249;
t220 = t199 * t230;
t253 = qJD(3) * t177 + t198 * t220;
t106 = -t168 * t245 + t253;
t85 = (-qJ(4) * t245 + t243) * qJD(2) + t106;
t136 = t168 * t198 + t231;
t244 = qJD(3) * t198;
t86 = -t136 * qJD(3) + (-qJ(4) * t244 + (-qJD(4) - t232) * t195) * qJD(2);
t37 = t190 * t86 + t279 * t85;
t33 = qJD(3) * qJD(5) + t37;
t148 = pkin(3) * t229 + t196 * t230;
t52 = pkin(4) * t143 - qJ(5) * t144 - qJD(5) * t153 + t148;
t17 = t189 * t52 + t192 * t33;
t118 = qJD(3) * pkin(3) + t123;
t227 = t279 * t124;
t61 = t118 * t190 + t227;
t56 = qJD(3) * qJ(5) + t61;
t145 = qJD(2) * t184 + qJD(4) - t232;
t82 = t150 * pkin(4) - t153 * qJ(5) + t145;
t24 = t189 * t82 + t192 * t56;
t290 = qJD(2) * pkin(2);
t289 = t153 * t73;
t264 = t191 * t196;
t158 = -t195 * t264 + t260;
t159 = t193 * t195 + t198 * t264;
t100 = -t158 * t279 + t159 * t190;
t36 = t190 * t85 - t279 * t86;
t287 = t36 * t100;
t128 = -t171 * t279 - t172 * t190;
t286 = t36 * t128;
t285 = t72 * t153;
t284 = pkin(5) * t270 + t280;
t278 = t131 * t153;
t277 = t131 * t189;
t275 = t144 * t164;
t272 = t150 * t153;
t266 = t168 * t195;
t263 = t191 * t199;
t201 = qJD(2) ^ 2;
t262 = t191 * t201;
t200 = qJD(3) ^ 2;
t257 = t200 * t195;
t256 = t200 * t198;
t254 = t143 * t192 - t149 * t189;
t187 = t195 ^ 2;
t188 = t198 ^ 2;
t252 = t187 - t188;
t247 = qJD(2) * t196;
t239 = t196 * t262;
t238 = t195 * t201 * t198;
t236 = t191 * t247;
t235 = t199 * t249;
t16 = -t189 * t33 + t192 * t52;
t23 = -t189 * t56 + t192 * t82;
t63 = t123 * t190 + t227;
t35 = qJD(6) * t72 + t310;
t223 = -t165 * t35 - t282 * t73;
t222 = t195 * t235;
t221 = t198 * t235;
t219 = t198 * t229;
t25 = pkin(5) * t274 + t36;
t183 = -pkin(3) * t279 - pkin(4);
t218 = t143 * t316 - t146 * t281;
t169 = -t232 - t290;
t217 = -t169 - t232;
t12 = pkin(5) * t143 - pkin(9) * t261 + t16;
t13 = -pkin(9) * t274 + t17;
t216 = t194 * t12 + t197 * t13;
t18 = pkin(5) * t150 - pkin(9) * t131 + t23;
t19 = pkin(9) * t317 + t24;
t5 = t18 * t197 - t19 * t194;
t6 = t18 * t194 + t19 * t197;
t215 = -t189 * t23 + t192 * t24;
t101 = t158 * t190 + t159 * t279;
t83 = -t101 * t189 - t192 * t263;
t84 = t101 * t192 - t189 * t263;
t31 = -t194 * t84 + t197 * t83;
t32 = t194 * t83 + t197 * t84;
t60 = t118 * t279 - t114;
t214 = t128 * t144 + t36 * t164;
t213 = t192 * t317;
t55 = -qJD(3) * pkin(4) + qJD(5) - t60;
t209 = t155 * t55 + t214;
t207 = t143 * t164 - t144 * t210 + t150 * t155;
t206 = t213 - t277;
t205 = qJD(3) * (-t217 - t290);
t204 = -t180 * t143 + t183 * t144 + (-qJD(5) + t55) * t150;
t2 = -qJD(6) * t6 + t12 * t197 - t194 * t13;
t107 = -t168 * t244 + (-qJD(3) * t193 - t235) * t250;
t135 = t177 - t266;
t202 = t106 * t198 - t107 * t195 + (-t135 * t198 - t136 * t195) * qJD(3);
t186 = t192 ^ 2;
t185 = t189 ^ 2;
t170 = -pkin(5) * t192 + t183;
t122 = -qJD(3) * t159 - t222;
t121 = qJD(3) * t158 + t221;
t99 = t316 * t164;
t98 = t165 * t164;
t92 = pkin(5) * t268 + t128;
t64 = t121 * t279 + t122 * t190;
t62 = t121 * t190 - t122 * t279;
t49 = t189 * t236 + t192 * t64;
t48 = -t189 * t64 + t192 * t236;
t47 = -pkin(5) * t271 + t63;
t46 = t155 * t165 + t164 * t307;
t45 = -t155 * t316 + t157 * t164;
t41 = -pkin(5) * t317 + t55;
t10 = -qJD(6) * t32 - t194 * t49 + t197 * t48;
t9 = qJD(6) * t31 + t194 * t48 + t197 * t49;
t1 = qJD(6) * t5 + t216;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239, -t199 * t262, 0, 0, 0, 0, 0, 0, 0, 0, -t198 * t239 + (t122 - t222) * qJD(3), t195 * t239 + (-t121 - t221) * qJD(3) (t121 * t198 - t122 * t195 + (-t158 * t198 - t159 * t195) * qJD(3)) * qJD(2), t106 * t159 + t107 * t158 + t136 * t121 + t135 * t122 + (t169 - t232) * t236, 0, 0, 0, 0, 0, 0, -t62 * qJD(3) + (-t143 * t199 + t150 * t247) * t191, -t64 * qJD(3) + (-t144 * t199 + t153 * t247) * t191, t100 * t144 - t101 * t143 - t150 * t64 + t153 * t62, t287 + t37 * t101 - t60 * t62 + t61 * t64 + (t145 * t247 - t148 * t199) * t191, 0, 0, 0, 0, 0, 0, t100 * t274 + t143 * t83 + t150 * t48 - t317 * t62, t100 * t261 + t131 * t62 - t143 * t84 - t150 * t49, t49 * t317 - t48 * t131 + (-t189 * t84 - t192 * t83) * t144, t16 * t83 + t17 * t84 + t23 * t48 + t24 * t49 + t55 * t62 + t287, 0, 0, 0, 0, 0, 0, t10 * t146 + t100 * t35 + t143 * t31 - t62 * t73, -t100 * t34 - t143 * t32 - t146 * t9 + t62 * t72, -t10 * t72 + t31 * t34 - t32 * t35 + t73 * t9, t1 * t32 + t10 * t5 + t100 * t25 + t2 * t31 + t41 * t62 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t219, -0.2e1 * t252 * t242, t256, -0.2e1 * t219, -t257, 0, -pkin(8) * t256 + t195 * t205, pkin(8) * t257 + t198 * t205 (-t187 - t188) * t220 + t202 ((t135 * t195 - t136 * t198) * t199 + (-t169 - t290) * t196) * t251 + t202 * pkin(8), t153 * t155 + t275, -t152 * t153 - t207, t155 * qJD(3), t308, -t152 * qJD(3), 0, -t150 * t237 + t184 * t143 + t145 * t152 - t148 * t210 + (t150 * t297 - t280) * qJD(3), -t153 * t237 + t184 * t144 + t145 * t155 + t148 * t164 + (t153 * t297 + t283) * qJD(3), -t129 * t143 + t150 * t283 - t61 * t152 + t153 * t280 - t60 * t155 + t210 * t37 + t214, t37 * t129 + t145 * t320 + t148 * t184 - t280 * t60 - t283 * t61 + t286, t131 * t269 + t186 * t275, t155 * t206 - 0.2e1 * t261 * t268, t131 * t152 + t192 * t207, t185 * t275 - t270 * t317, t152 * t317 - t189 * t207, t308, t57 * t143 + t292 * t150 + t23 * t152 - t16 * t210 + t209 * t189 - t280 * t317, t131 * t280 - t58 * t143 - t150 * t291 - t24 * t152 + t17 * t210 + t192 * t209, -t291 * t265 - t292 * t131 + (qJD(3) * t291 - t57 * t144 - t23 * t155 - t16 * t164) * t192 + (-t144 * t58 - t155 * t24 - t164 * t17) * t189, t16 * t57 + t17 * t58 + t23 * t292 + t24 * t291 + t280 * t55 + t286, -t34 * t99 - t45 * t72, t34 * t98 - t35 * t99 - t45 * t73 - t46 * t72, t143 * t99 - t146 * t45 + t152 * t72 + t210 * t34, t35 * t98 - t46 * t73, -t143 * t98 - t146 * t46 + t152 * t73 + t210 * t35, t146 * t152 - t112, t14 * t143 + t146 * t300 + t5 * t152 - t2 * t210 + t25 * t98 - t284 * t73 + t92 * t35 + t41 * t46, t1 * t210 - t15 * t143 - t146 * t301 - t6 * t152 + t25 * t99 + t284 * t72 - t92 * t34 - t41 * t45, -t1 * t98 + t14 * t34 - t15 * t35 - t2 * t99 - t300 * t72 + t301 * t73 + t5 * t45 - t6 * t46, t1 * t15 + t2 * t14 + t25 * t92 + t284 * t41 + t300 * t5 + t301 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t252 * t201, 0, t238, 0, 0, t217 * t248, -t169 * t198 * qJD(2) + (t135 + t266) * qJD(3) - t253, 0, 0, t272, -t149 + t302, -t173 + (t176 + t150) * qJD(3), -t272, 0, 0, qJD(3) * t63 - t145 * t153 - t150 * t241 - t36, qJD(3) * t65 + t145 * t150 - t153 * t241 - t37 (t61 - t63) * t153 + (-t60 + t65) * t150 + (-t143 * t190 - t144 * t279) * pkin(3), t60 * t63 - t61 * t65 + (-t145 * t248 + t190 * t37 - t279 * t36) * pkin(3), t311 * t192, t206 * t150 + (-t185 + t186) * t144, -t278 - t303, -t318 * t189, t254 - t309, -t272, -t29 * t150 - t23 * t153 + t189 * t204 - t36 * t192 + t317 * t63, -t63 * t131 + t30 * t150 + t24 * t153 + t36 * t189 + t192 * t204, t29 * t131 + t30 * t265 + (-qJD(5) * t265 - t23 * t150 + t17 + (qJD(5) * t192 - t30) * qJD(3)) * t192 + (qJD(5) * t131 - t150 * t24 - t16) * t189, t36 * t183 - t23 * t29 - t24 * t30 - t55 * t63 + (-t16 * t189 + t17 * t192) * t180 + t215 * qJD(5), -t165 * t34 - t282 * t72, t223 + t305, -t285 + t304, -t281 * t73 - t316 * t35, t218 - t289, -t146 * t153, t104 * t143 - t146 * t299 - t5 * t153 + t170 * t35 - t25 * t316 + t281 * t41 + t47 * t73, -t105 * t143 - t146 * t298 + t6 * t153 + t25 * t165 - t170 * t34 - t282 * t41 - t47 * t72, t1 * t316 + t104 * t34 - t105 * t35 - t2 * t165 - t281 * t6 + t282 * t5 + t298 * t73 + t299 * t72, t1 * t105 + t2 * t104 + t25 * t170 + t298 * t6 - t299 * t5 - t41 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t153 * qJD(3), -t173 + (t176 - t150) * qJD(3), -t149 - t302, t150 * t61 + t153 * t60 + t148, 0, 0, 0, 0, 0, 0, t254 + t309, -t278 + t303 (t213 + t277) * t150 + (-t185 - t186) * t144, t150 * t215 - t55 * t153 + t16 * t192 + t17 * t189, 0, 0, 0, 0, 0, 0, t218 + t289, -t285 - t304, t223 - t305, t1 * t165 - t41 * t153 + t2 * t316 - t281 * t5 - t282 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, t318, -t131 ^ 2 - t317 ^ 2, t131 * t23 - t24 * t317 + t36, 0, 0, 0, 0, 0, 0, t72 * t146 + t35, -t34 + t322, -t321 - t323, t5 * t72 - t6 * t73 + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295, t321 - t323, -t34 - t322, t295, t306 * t72 - t310, t143, t6 * t146 - t41 * t72 + t2, t306 * t5 - t41 * t73 - t216, 0, 0;];
tauc_reg  = t3;
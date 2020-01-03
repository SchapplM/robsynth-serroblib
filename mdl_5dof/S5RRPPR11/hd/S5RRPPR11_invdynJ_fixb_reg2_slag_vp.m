% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:48:04
% EndTime: 2019-12-31 19:48:11
% DurationCPUTime: 3.40s
% Computational Cost: add. (3853->449), mult. (8229->566), div. (0->0), fcn. (5181->10), ass. (0->242)
t169 = cos(pkin(8));
t168 = sin(pkin(8));
t250 = t168 * qJD(2);
t176 = cos(qJ(2));
t255 = qJD(1) * t176;
t100 = -t169 * t255 - t250;
t233 = t168 * t255;
t249 = t169 * qJD(2);
t101 = -t233 + t249;
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t45 = -t100 * t175 + t101 * t172;
t306 = t45 ^ 2;
t48 = t100 * t172 + t101 * t175;
t305 = t48 ^ 2;
t298 = pkin(3) + pkin(6);
t104 = t168 * t175 + t169 * t172;
t173 = sin(qJ(2));
t194 = t104 * t173;
t90 = t104 * qJD(5);
t283 = qJD(1) * t194 + t90;
t248 = t173 * qJD(1);
t136 = qJD(5) + t248;
t304 = t136 * t45;
t174 = sin(qJ(1));
t177 = cos(qJ(1));
t214 = g(1) * t177 + g(2) * t174;
t146 = pkin(6) * t248;
t303 = qJD(3) + t146;
t251 = qJD(5) * t175;
t252 = qJD(5) * t172;
t245 = qJD(1) * qJD(2);
t232 = t173 * t245;
t242 = t176 * qJDD(1);
t244 = t168 * qJDD(2);
t65 = t244 + (-t232 + t242) * t169;
t125 = t168 * t242;
t222 = t169 * qJDD(2) - t125;
t66 = t168 * t232 + t222;
t16 = -t100 * t251 + t101 * t252 + t172 * t65 - t175 * t66;
t205 = t168 * t172 - t169 * t175;
t302 = t16 * t205 - t283 * t48;
t226 = t172 * t66 + t175 * t65;
t17 = qJD(5) * t48 + t226;
t235 = t169 * t248;
t236 = t168 * t248;
t284 = t168 * t252 - t169 * t251 + t172 * t236 - t175 * t235;
t301 = t104 * t17 - t284 * t45;
t253 = qJD(2) * t176;
t289 = pkin(2) + qJ(4);
t147 = pkin(6) * t255;
t111 = pkin(3) * t255 + t147;
t165 = qJD(2) * qJ(3);
t88 = qJD(4) + t165 + t111;
t300 = t289 * t253 + (qJD(4) - t88) * t173;
t154 = t173 * qJ(3);
t227 = -pkin(1) - t154;
t191 = -t176 * t289 + t227;
t70 = t191 * qJD(1);
t246 = pkin(3) * t248 + t303;
t79 = -qJD(2) * t289 + t246;
t28 = -t168 * t70 + t169 * t79;
t18 = pkin(4) * t248 - pkin(7) * t101 + t28;
t29 = t168 * t79 + t169 * t70;
t19 = pkin(7) * t100 + t29;
t209 = t172 * t19 - t175 * t18;
t135 = pkin(2) * t232;
t277 = qJ(3) * t176;
t208 = qJ(4) * t173 - t277;
t247 = t173 * qJD(3);
t182 = qJD(2) * t208 - t176 * qJD(4) - t247;
t24 = qJD(1) * t182 + qJDD(1) * t191 + t135;
t229 = t176 * t245;
t243 = t173 * qJDD(1);
t193 = t229 + t243;
t134 = pkin(6) * t229;
t143 = pkin(6) * t243;
t228 = qJDD(3) + t134 + t143;
t41 = pkin(3) * t193 - qJD(2) * qJD(4) - qJDD(2) * t289 + t228;
t12 = -t168 * t24 + t169 * t41;
t8 = pkin(4) * t193 - t66 * pkin(7) + t12;
t13 = t168 * t41 + t169 * t24;
t9 = -pkin(7) * t65 + t13;
t1 = -qJD(5) * t209 + t172 * t8 + t175 * t9;
t7 = t172 * t18 + t175 * t19;
t2 = -qJD(5) * t7 - t172 * t9 + t175 * t8;
t161 = g(3) * t176;
t271 = t173 * t177;
t272 = t173 * t174;
t238 = -g(1) * t271 - g(2) * t272 + t161;
t299 = t1 * t104 - t2 * t205 + t209 * t283 - t284 * t7 + t238;
t297 = t65 * pkin(4);
t296 = g(1) * t174;
t293 = g(2) * t177;
t292 = g(3) * t173;
t291 = t168 * pkin(4);
t158 = t176 * pkin(2);
t290 = t48 * t45;
t287 = -pkin(7) - t289;
t274 = t168 * t173;
t203 = pkin(4) * t176 - pkin(7) * t274;
t150 = pkin(2) * t248;
t84 = qJD(1) * t208 + t150;
t39 = t111 * t169 - t168 * t84;
t25 = qJD(1) * t203 + t39;
t40 = t111 * t168 + t169 * t84;
t32 = pkin(7) * t235 + t40;
t113 = t287 * t168;
t114 = t287 * t169;
t53 = -t113 * t172 + t114 * t175;
t286 = -qJD(4) * t104 + qJD(5) * t53 - t172 * t25 - t175 * t32;
t54 = t113 * t175 + t114 * t172;
t285 = qJD(4) * t205 - qJD(5) * t54 + t172 * t32 - t175 * t25;
t112 = t298 * t253;
t254 = qJD(2) * t173;
t149 = pkin(2) * t254;
t55 = t149 + t182;
t31 = t112 * t168 + t169 * t55;
t280 = t65 * t168;
t279 = t66 * t169;
t122 = t298 * t173;
t259 = t158 + t154;
t266 = t176 * qJ(4);
t211 = t259 + t266;
t97 = -pkin(1) - t211;
t50 = t122 * t168 + t169 * t97;
t278 = pkin(6) * qJDD(2);
t276 = qJDD(2) * pkin(2);
t166 = t173 ^ 2;
t179 = qJD(1) ^ 2;
t275 = t166 * t179;
t273 = t169 * t176;
t270 = t173 * t179;
t162 = pkin(8) + qJ(5);
t152 = sin(t162);
t269 = t174 * t152;
t153 = cos(t162);
t268 = t174 * t153;
t267 = t174 * t176;
t170 = -pkin(7) - qJ(4);
t265 = t176 * t170;
t264 = t176 * t177;
t263 = t177 * t152;
t262 = t177 * t153;
t142 = pkin(4) * t169 + pkin(3);
t261 = t142 * t248 + t303;
t123 = t298 * t176;
t258 = pkin(1) * t177 + pkin(6) * t174;
t167 = t176 ^ 2;
t257 = t166 - t167;
t256 = t166 + t167;
t241 = t176 * t291;
t240 = pkin(4) * t274;
t144 = pkin(6) * t242;
t163 = qJDD(2) * qJ(3);
t164 = qJD(2) * qJD(3);
t239 = t144 + t163 + t164;
t237 = t298 * qJD(2);
t231 = t168 * t243;
t126 = t169 * t243;
t230 = t169 * t242;
t30 = t112 * t169 - t168 * t55;
t225 = -qJD(2) * pkin(2) + qJD(3);
t224 = qJD(1) * t257;
t223 = t100 + t250;
t221 = pkin(2) * t264 + qJ(3) * t271 + t258;
t220 = -t143 - t238;
t103 = qJDD(5) + t193;
t219 = -t205 * t103 - t136 * t283;
t218 = qJD(1) * t237;
t217 = t173 * t229;
t216 = t256 * qJDD(1) * pkin(6);
t178 = qJD(2) ^ 2;
t215 = pkin(6) * t178 + t293;
t213 = -t293 + t296;
t210 = t12 * t169 + t13 * t168;
t107 = t169 * t122;
t33 = t173 * pkin(4) + t107 + (pkin(7) * t176 - t97) * t168;
t38 = -pkin(7) * t273 + t50;
t14 = -t172 * t38 + t175 * t33;
t15 = t172 * t33 + t175 * t38;
t207 = t100 * t168 + t101 * t169;
t115 = t146 + t225;
t121 = -t147 - t165;
t206 = t115 * t176 + t121 * t173;
t204 = pkin(3) * t242 + qJDD(4) + t239;
t202 = t227 - t158;
t201 = -t168 * t275 + t169 * t229 + t126;
t200 = t214 * t176;
t96 = t202 * qJD(1);
t199 = t248 * t96 + qJDD(3) - t220;
t198 = -0.2e1 * pkin(1) * t245 - t278;
t197 = (-t168 * t28 + t169 * t29) * t173;
t196 = -t103 * t104 + t136 * t284;
t195 = -qJ(3) * t253 - t247;
t80 = t205 * t176;
t43 = -t173 * t218 + t204;
t192 = -t176 * t43 + t214;
t190 = 0.2e1 * qJDD(1) * pkin(1) - t215;
t189 = -t169 * t275 - t231;
t116 = -pkin(1) - t259;
t187 = t278 + (-qJD(1) * t116 - t96) * qJD(2);
t186 = -t200 - t292;
t185 = t43 + t186;
t42 = qJD(1) * t195 + qJDD(1) * t202 + t135;
t86 = t149 + t195;
t184 = qJD(1) * t86 + qJDD(1) * t116 + t215 + t42;
t71 = pkin(6) * t232 - t239;
t83 = t228 - t276;
t183 = qJD(2) * t206 + t83 * t173 - t71 * t176;
t181 = (-g(3) - t218) * t173 - t200 + t204;
t159 = t177 * pkin(6);
t141 = qJ(3) + t291;
t139 = g(1) * t267;
t133 = qJ(3) * t264;
t131 = qJ(3) * t267;
t130 = t176 * t270;
t120 = t257 * t179;
t118 = qJDD(2) * t176 - t173 * t178;
t117 = qJDD(2) * t173 + t176 * t178;
t110 = t173 * t237;
t108 = -qJ(3) * t255 + t150;
t99 = qJDD(1) * t167 - 0.2e1 * t217;
t98 = qJDD(1) * t166 + 0.2e1 * t217;
t89 = pkin(4) * t273 + t123;
t81 = t104 * t176;
t78 = -t173 * t269 + t262;
t77 = t173 * t268 + t263;
t76 = t173 * t263 + t268;
t75 = t173 * t262 - t269;
t73 = (-pkin(6) - t142) * t254;
t64 = 0.2e1 * t173 * t242 - 0.2e1 * t245 * t257;
t51 = -pkin(4) * t100 + t88;
t49 = -t168 * t97 + t107;
t35 = t176 * t90 - t205 * t254;
t34 = qJD(2) * t194 + qJD(5) * t80;
t23 = pkin(7) * t173 * t249 + t31;
t22 = t43 + t297;
t20 = qJD(2) * t203 + t30;
t4 = -qJD(5) * t15 - t172 * t23 + t175 * t20;
t3 = qJD(5) * t14 + t172 * t20 + t175 * t23;
t5 = [0, 0, 0, 0, 0, qJDD(1), t213, t214, 0, 0, t98, t64, t117, t99, t118, 0, t173 * t198 + t176 * t190 + t139, t198 * t176 + (-t190 - t296) * t173, -t214 + 0.2e1 * t216, -g(1) * (-t174 * pkin(1) + t159) - g(2) * t258 + (pkin(6) ^ 2 * t256 + pkin(1) ^ 2) * qJDD(1), 0, -t117, -t118, t98, t64, t99, t216 + t183 - t214, t173 * t187 + t176 * t184 - t139, t187 * t176 + (-t184 + t296) * t173, pkin(6) * t183 - g(1) * t159 - g(2) * t221 + t42 * t116 - t202 * t296 + t96 * t86, (t101 * t254 - t176 * t66) * t168, (-t279 + t280) * t176 + t207 * t254, (t66 - t125) * t173 + (t101 * t176 + t168 * t224) * qJD(2), (t100 * t254 + t176 * t65) * t169, (-t65 - t230) * t173 + (t100 * t176 + t169 * t224) * qJD(2), t98, t110 * t100 + t123 * t65 + (qJD(1) * t49 + t28) * t253 - t192 * t169 + (qJD(1) * t30 + qJDD(1) * t49 + t168 * t213 - t249 * t88 + t12) * t173, -t110 * t101 + t123 * t66 + (-qJD(1) * t50 - t29) * t253 + t192 * t168 + (-t31 * qJD(1) - t50 * qJDD(1) + t169 * t213 + t250 * t88 - t13) * t173, t31 * t100 - t30 * t101 - t49 * t66 - t50 * t65 + t139 + qJD(2) * t197 + (t12 * t168 - t13 * t169 - t293) * t176, t13 * t50 + t29 * t31 + t12 * t49 + t28 * t30 + t43 * t123 - t88 * t110 - g(1) * (t177 * pkin(3) + t159) - g(2) * (qJ(4) * t264 + t221) + (-g(1) * (t202 - t266) - g(2) * pkin(3)) * t174, t16 * t81 + t34 * t48, -t16 * t80 + t17 * t81 - t34 * t45 + t35 * t48, -t103 * t81 + t136 * t34 - t16 * t173 + t253 * t48, -t17 * t80 - t35 * t45, t103 * t80 + t136 * t35 - t17 * t173 - t253 * t45, t103 * t173 + t136 * t253, -g(1) * t78 - g(2) * t76 + t103 * t14 + t136 * t4 + t17 * t89 + t173 * t2 - t209 * t253 - t22 * t80 - t35 * t51 + t45 * t73, g(1) * t77 - g(2) * t75 - t1 * t173 - t103 * t15 - t136 * t3 - t16 * t89 - t22 * t81 - t253 * t7 + t34 * t51 + t48 * t73, -g(2) * t264 + t1 * t80 + t14 * t16 - t15 * t17 + t2 * t81 + t209 * t34 - t3 * t45 + t35 * t7 - t4 * t48 + t139, t1 * t15 + t7 * t3 + t2 * t14 - t209 * t4 + t22 * t89 + t51 * t73 - g(1) * (t177 * t142 + t159) - g(2) * (-t170 * t264 + t177 * t240 + t221) + (-g(1) * (t202 - t240 + t265) - g(2) * t142) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t120, t243, t130, t242, qJDD(2), pkin(1) * t270 + t220, t292 - t144 + (pkin(1) * t179 + t214) * t176, 0, 0, qJDD(2), -t243, -t242, -t130, t120, t130, (-pkin(2) * t173 + t277) * qJDD(1) + ((-t121 - t165) * t173 + (-t115 + t225) * t176) * qJD(1), -t108 * t255 + t199 - 0.2e1 * t276, t144 + 0.2e1 * t163 + 0.2e1 * t164 + (qJD(1) * t108 - g(3)) * t173 + (qJD(1) * t96 - t214) * t176, -t71 * qJ(3) - t121 * qJD(3) - t83 * pkin(2) - t96 * t108 - g(1) * (-pkin(2) * t271 + t133) - g(2) * (-pkin(2) * t272 + t131) - g(3) * t259 - t206 * qJD(1) * pkin(6), -t101 * t236 + t279, -t66 * t168 - t169 * t65 - t207 * t248, -t101 * t255 + t201, -t100 * t235 + t280, -t223 * t255 + t189, -t130, -t289 * t126 + qJ(3) * t65 - t246 * t100 + t185 * t168 + (-t169 * t300 - t173 * t39 - t28 * t176) * qJD(1), t289 * t231 + qJ(3) * t66 + t246 * t101 + t185 * t169 + (t168 * t300 + t173 * t40 + t29 * t176) * qJD(1), -t40 * t100 + t39 * t101 + (qJD(4) * t101 - t248 * t29 + t289 * t66 - t12) * t169 + (-qJD(4) * t100 + t248 * t28 + t289 * t65 - t13) * t168 - t238, t43 * qJ(3) - t29 * t40 - t28 * t39 - g(1) * t133 - g(2) * t131 - g(3) * t211 + t246 * t88 + (-t168 * t29 - t169 * t28) * qJD(4) + (t173 * t214 - t210) * t289, t302, t16 * t104 + t17 * t205 + t283 * t45 + t284 * t48, -t255 * t48 + t219, t301, t255 * t45 + t196, -t136 * t255, t53 * t103 + t22 * t104 + t136 * t285 + t141 * t17 + t152 * t186 + t209 * t255 + t261 * t45 - t284 * t51, -t54 * t103 - t136 * t286 - t141 * t16 + t153 * t186 - t205 * t22 + t255 * t7 + t261 * t48 - t283 * t51, t53 * t16 - t54 * t17 - t285 * t48 - t286 * t45 - t299, t1 * t54 + t2 * t53 + t22 * t141 - g(1) * (t177 * t241 + t133) - g(2) * (t174 * t241 + t131) - g(3) * (t259 - t265) + t286 * t7 - t285 * t209 + t261 * t51 + (-g(3) * t291 + t214 * (pkin(2) - t170)) * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, qJDD(2) + t130, -t178 - t275, qJD(2) * t121 + t134 + t199 - t276, 0, 0, 0, 0, 0, 0, qJD(2) * t100 + t201, (-t101 - t233) * qJD(2) + t189, -t280 - t279 + (t100 * t169 + t101 * t168) * t248, qJD(1) * t197 - t88 * qJD(2) + t210 + t238, 0, 0, 0, 0, 0, 0, -qJD(2) * t45 + t219, -qJD(2) * t48 + t196, -t301 - t302, -t51 * qJD(2) + t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230 + t244 + (t101 - t249) * t248, t223 * t248 + t222, -t100 ^ 2 - t101 ^ 2, -t29 * t100 + t28 * t101 + t181, 0, 0, 0, 0, 0, 0, t48 * t136 + t17, -t16 - t304, -t305 - t306, -t209 * t48 + t45 * t7 + t181 + t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t305 - t306, -t16 + t304, -t290, -t226 + (-qJD(5) + t136) * t48, t103, -g(1) * t75 - g(2) * t77 + t7 * t136 + t153 * t161 - t51 * t48 + t2, g(1) * t76 - g(2) * t78 - t136 * t209 - t152 * t161 + t51 * t45 - t1, 0, 0;];
tau_reg = t5;

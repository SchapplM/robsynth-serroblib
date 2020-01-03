% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR14_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:43
% EndTime: 2019-12-31 19:20:02
% DurationCPUTime: 6.51s
% Computational Cost: add. (15032->470), mult. (50046->696), div. (0->0), fcn. (42094->12), ass. (0->208)
t170 = sin(pkin(5));
t169 = sin(pkin(11));
t171 = cos(pkin(11));
t174 = sin(qJ(3));
t268 = cos(pkin(6));
t236 = t174 * t268;
t286 = cos(qJ(3));
t193 = t169 * t236 - t286 * t171;
t187 = t170 * t193;
t135 = qJD(1) * t187;
t267 = sin(pkin(6));
t221 = t267 * t286;
t302 = qJD(3) * t221 + t135;
t269 = cos(pkin(5));
t214 = t269 * t267;
t222 = t268 * t286;
t261 = t170 * t171;
t289 = -t286 * t214 - t222 * t261;
t253 = qJD(1) * t170;
t241 = t169 * t253;
t299 = t289 * qJD(1);
t115 = t174 * t241 + t299;
t196 = qJD(4) + t115;
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t173 = sin(qJ(4));
t176 = cos(qJ(4));
t124 = t170 * (t286 * t169 + t171 * t236) + t174 * t214;
t119 = t124 * qJD(1);
t231 = qJD(1) * t269;
t224 = pkin(1) * t231;
t160 = t171 * t224;
t262 = t169 * t170;
t183 = t269 * pkin(2) + (-t268 * pkin(8) - qJ(2)) * t262;
t114 = qJD(1) * t183 + t160;
t137 = (-t267 * pkin(8) * t169 - pkin(2) * t171 - pkin(1)) * t170;
t129 = qJD(1) * t137 + qJD(2);
t88 = -t267 * t114 + t268 * t129;
t47 = t115 * pkin(3) - t119 * pkin(9) + t88;
t215 = t269 * t268;
t230 = qJD(1) * t267;
t220 = t170 * t230;
t153 = t171 * t220;
t245 = qJD(3) - t153;
t190 = -qJD(1) * t215 - t245;
t161 = qJ(2) * t261;
t142 = qJD(1) * t161 + t169 * t224;
t182 = (t268 * t261 + t214) * pkin(8);
t108 = qJD(1) * t182 + t142;
t62 = (t268 * t114 + t267 * t129) * t174 + t286 * t108;
t49 = -pkin(9) * t190 + t62;
t27 = t173 * t47 + t176 * t49;
t24 = pkin(10) * t196 + t27;
t290 = -t174 * t108 + t114 * t222 + t129 * t221;
t48 = pkin(3) * t190 - t290;
t136 = t176 * t190;
t90 = t119 * t173 + t136;
t92 = t176 * t119 - t173 * t190;
t30 = t90 * pkin(4) - t92 * pkin(10) + t48;
t213 = t172 * t24 - t175 * t30;
t191 = t169 * t222 + t171 * t174;
t186 = t170 * t191;
t184 = qJD(2) * t186;
t45 = qJD(1) * t184 + qJD(3) * t62;
t244 = t174 * t262;
t225 = qJD(3) * t244;
t106 = qJD(1) * t225 + t299 * qJD(3);
t251 = qJD(4) * t173;
t64 = qJD(4) * t136 + t176 * t106 + t119 * t251;
t259 = t173 * t106;
t266 = qJD(4) * t92;
t65 = -t259 + t266;
t22 = t65 * pkin(4) + t64 * pkin(10) + t45;
t118 = t124 * qJD(3);
t107 = qJD(1) * t118;
t249 = qJD(4) * t176;
t185 = qJD(2) * t187;
t44 = -qJD(1) * t185 + qJD(3) * t290;
t229 = qJD(2) * t267;
t206 = t229 * t262;
t75 = t107 * pkin(3) + t106 * pkin(9) + qJD(1) * t206;
t201 = -t173 * t75 - t176 * t44 - t47 * t249 + t49 * t251;
t7 = pkin(10) * t107 - t201;
t1 = -t213 * qJD(5) + t172 * t22 + t175 * t7;
t89 = qJD(5) + t90;
t301 = t213 * t89 + t1;
t188 = qJD(4) * t196;
t300 = pkin(9) * t188 + t45;
t235 = t174 * t267;
t144 = t173 * t235 - t176 * t268;
t207 = t169 * t220;
t257 = -qJD(4) * t144 - t173 * t207 + t302 * t176;
t26 = -t173 * t49 + t176 * t47;
t298 = -t196 * t26 - t201;
t134 = qJD(1) * t186;
t202 = qJD(3) * t235 - t134;
t297 = pkin(9) * qJD(5) * t176 + t62 - t196 * (pkin(4) * t173 - pkin(10) * t176);
t243 = pkin(1) * t269;
t255 = t169 * t243 + t161;
t121 = t182 + t255;
t163 = t171 * t243;
t125 = t163 + t183;
t70 = t286 * t121 + (t268 * t125 + t267 * t137) * t174;
t6 = t172 * t30 + t175 * t24;
t2 = -qJD(5) * t6 - t172 * t7 + t175 * t22;
t295 = -t6 * t89 - t2;
t293 = t196 * t90;
t292 = t196 * t92;
t291 = t170 ^ 2 * (t169 ^ 2 + t171 ^ 2);
t69 = -t174 * t121 + t125 * t222 + t137 * t221;
t145 = t173 * t268 + t176 * t235;
t256 = qJD(4) * t145 + t302 * t173 + t176 * t207;
t68 = t172 * t196 + t175 * t92;
t265 = qJD(5) * t68;
t32 = -t175 * t107 - t172 * t64 + t265;
t237 = t173 * t44 - t176 * t75;
t10 = -qJD(4) * t27 - t237;
t123 = t244 + t289;
t93 = -t267 * t125 + t268 * t137;
t56 = t123 * pkin(3) - t124 * pkin(9) + t93;
t143 = t267 * t261 - t215;
t60 = -pkin(9) * t143 + t70;
t278 = t173 * t56 + t176 * t60;
t51 = t69 * qJD(3) - t185;
t117 = qJD(3) * t289 + t225;
t79 = t118 * pkin(3) + t117 * pkin(9) + t206;
t16 = -qJD(4) * t278 - t173 * t51 + t176 * t79;
t94 = t124 * t173 + t143 * t176;
t285 = t65 * t94;
t189 = t175 * t196;
t66 = t172 * t92 - t189;
t284 = t66 * t89;
t283 = t68 * t66;
t282 = t68 * t89;
t281 = t92 * t90;
t156 = -pkin(4) * t176 - pkin(10) * t173 - pkin(3);
t248 = qJD(5) * t172;
t83 = pkin(3) * t119 + pkin(9) * t115;
t41 = t173 * t83 + t176 * t290;
t36 = pkin(10) * t119 + t41;
t280 = t156 * t248 + (-t251 * pkin(9) - t36) * t172 + t297 * t175;
t247 = qJD(5) * t175;
t250 = qJD(4) * t175;
t279 = t173 * t250 * pkin(9) - t156 * t247 + t297 * t172 + t175 * t36;
t194 = -t175 * t145 + t172 * t221;
t277 = -t194 * qJD(5) + t257 * t172 - t202 * t175;
t130 = -t172 * t145 - t175 * t221;
t276 = -t130 * qJD(5) - t202 * t172 - t257 * t175;
t31 = -qJD(5) * t189 - t172 * t107 + t175 * t64 + t92 * t248;
t275 = t172 * t31;
t274 = t172 * t65;
t273 = t172 * t89;
t272 = t175 * t32;
t271 = t175 * t65;
t270 = t176 * t65;
t264 = t107 * t123;
t263 = t119 * t115;
t260 = t172 * t176;
t258 = t175 * t176;
t252 = qJD(2) * t170;
t240 = t169 * t252;
t239 = t267 * t88;
t233 = t267 * t119;
t228 = t175 * t89;
t227 = t173 * t196;
t8 = -pkin(4) * t107 - t10;
t226 = pkin(10) * qJD(5) * t89 + t8;
t177 = qJD(1) ^ 2;
t223 = t170 * t177 * t269;
t217 = -t172 * t6 + t175 * t213;
t29 = pkin(10) * t123 + t278;
t59 = t143 * pkin(3) - t69;
t95 = t124 * t176 - t143 * t173;
t37 = t94 * pkin(4) - t95 * pkin(10) + t59;
t12 = t172 * t37 + t175 * t29;
t11 = -t172 * t29 + t175 * t37;
t33 = -t173 * t60 + t176 * t56;
t40 = -t173 * t290 + t176 * t83;
t210 = t123 * t175 - t172 * t95;
t77 = t123 * t172 + t175 * t95;
t209 = (-qJ(2) * t241 + t160) * t169 - t142 * t171;
t208 = t45 * t221;
t23 = -pkin(4) * t196 - t26;
t203 = -pkin(10) * t65 + t89 * t23;
t15 = t173 * t79 + t176 * t51 + t56 * t249 - t60 * t251;
t200 = -0.2e1 * t231 * t252;
t199 = -pkin(9) * t107 + t196 * t48;
t178 = t115 * t196 + t188;
t52 = t70 * qJD(3) + t184;
t139 = pkin(9) * t258 + t156 * t172;
t138 = -pkin(9) * t260 + t156 * t175;
t81 = -t115 * t258 + t119 * t172;
t80 = -t115 * t260 - t175 * t119;
t74 = -qJD(4) * t94 - t117 * t176;
t73 = qJD(4) * t95 - t117 * t173;
t54 = pkin(4) * t92 + pkin(10) * t90;
t39 = qJD(5) * t210 + t118 * t172 + t175 * t74;
t38 = qJD(5) * t77 - t118 * t175 + t172 * t74;
t35 = -pkin(4) * t119 - t40;
t28 = -pkin(4) * t123 - t33;
t25 = t73 * pkin(4) - t74 * pkin(10) + t52;
t20 = t172 * t54 + t175 * t26;
t19 = -t172 * t26 + t175 * t54;
t14 = -pkin(4) * t118 - t16;
t13 = pkin(10) * t118 + t15;
t4 = -t12 * qJD(5) - t13 * t172 + t175 * t25;
t3 = t11 * qJD(5) + t13 * t175 + t172 * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169 * t200, t171 * t200, 0.2e1 * qJD(2) * qJD(1) * t291, ((t171 * t255 + (qJ(2) * t262 - t163) * t169) * qJD(1) - t209) * t252, -t106 * t124 - t117 * t119, t106 * t123 - t107 * t124 + t115 * t117 - t118 * t119, t106 * t143 + t117 * t190, t115 * t118 + t264, t107 * t143 + t118 * t190, 0, t52 * t190 + t45 * t143 + t93 * t107 + t88 * t118 + (t267 * t115 + t123 * t230) * t240, t51 * t190 + t44 * t143 - t93 * t106 - t88 * t117 + (t124 * t230 + t233) * t240, t106 * t69 - t107 * t70 - t115 * t51 + t117 * t290 - t118 * t62 + t119 * t52 - t123 * t44 + t124 * t45, t44 * t70 - t45 * t69 + t62 * t51 - t290 * t52 + (t93 * t230 + t239) * t240, -t64 * t95 + t74 * t92, t64 * t94 - t65 * t95 - t73 * t92 - t74 * t90, t95 * t107 + t92 * t118 - t64 * t123 + t196 * t74, t73 * t90 + t285, -t94 * t107 - t90 * t118 - t65 * t123 - t196 * t73, t118 * t196 + t264, t10 * t123 + t33 * t107 + t26 * t118 + t16 * t196 + t45 * t94 + t48 * t73 + t52 * t90 + t59 * t65, -t107 * t278 - t27 * t118 + t123 * t201 - t15 * t196 + t45 * t95 + t48 * t74 + t52 * t92 - t59 * t64, -t10 * t95 - t15 * t90 - t16 * t92 + t201 * t94 - t26 * t74 - t27 * t73 - t278 * t65 + t33 * t64, t10 * t33 + t15 * t27 + t16 * t26 - t201 * t278 + t45 * t59 + t48 * t52, -t31 * t77 + t39 * t68, -t210 * t31 - t32 * t77 - t38 * t68 - t39 * t66, -t31 * t94 + t39 * t89 + t65 * t77 + t68 * t73, -t210 * t32 + t38 * t66, t210 * t65 - t32 * t94 - t38 * t89 - t66 * t73, t73 * t89 + t285, t11 * t65 + t14 * t66 + t2 * t94 - t210 * t8 - t213 * t73 + t23 * t38 + t28 * t32 + t4 * t89, -t1 * t94 - t12 * t65 + t14 * t68 + t23 * t39 - t28 * t31 - t3 * t89 - t6 * t73 + t77 * t8, t1 * t210 + t11 * t31 - t12 * t32 - t2 * t77 + t213 * t39 - t3 * t66 - t38 * t6 - t4 * t68, t1 * t12 + t11 * t2 + t14 * t23 - t213 * t4 + t28 * t8 + t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169 * t223, t171 * t223, -t177 * t291, t209 * t253, 0, 0, 0, 0, 0, 0, t268 * t107 - t115 * t207 + t190 * t202, -t268 * t106 - t119 * t207 + t302 * t190, t106 * t221 - t107 * t235 - t135 * t115 - t134 * t119 + (-t115 * t221 + t174 * t233) * qJD(3), -t208 + t44 * t235 + t290 * t134 + t62 * t135 + (t62 * t221 - t235 * t290) * qJD(3) + (t268 * t229 - t239) * t241, 0, 0, 0, 0, 0, 0, -t144 * t107 - t196 * t256 + t202 * t90 - t221 * t65, -t145 * t107 - t257 * t196 + t202 * t92 + t221 * t64, -t144 * t64 - t145 * t65 + t256 * t92 - t257 * t90, -t10 * t144 - t145 * t201 + t202 * t48 - t256 * t26 + t257 * t27 - t208, 0, 0, 0, 0, 0, 0, t130 * t65 + t144 * t32 + t256 * t66 - t277 * t89, -t144 * t31 + t194 * t65 + t256 * t68 + t276 * t89, t130 * t31 + t194 * t32 + t276 * t66 + t277 * t68, -t1 * t194 + t130 * t2 + t144 * t8 + t213 * t277 + t23 * t256 - t276 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, -t115 ^ 2 + t119 ^ 2, -t115 * t190 - t106, -t263, t119 * t245 + (t119 * t215 - t118) * qJD(1), 0, -t88 * t119 - t62 * t153 + (-t191 * t252 + t62 * t215) * qJD(1), t88 * t115 - t290 * t153 + (t193 * t252 + t215 * t290) * qJD(1), 0, 0, -t173 * t64 + t176 * t292, (-t64 - t293) * t176 + (-t65 - t292) * t173, t173 * t107 - t92 * t119 + t176 * t178, t227 * t90 - t270, t176 * t107 + t90 * t119 - t173 * t178, -t196 * t119, -pkin(3) * t65 - t26 * t119 + t199 * t173 - t300 * t176 - t40 * t196 - t62 * t90, pkin(3) * t64 + t27 * t119 + t300 * t173 + t199 * t176 + t41 * t196 - t62 * t92, t40 * t92 + t41 * t90 + ((-t65 + t266) * pkin(9) + t298) * t176 + (-t10 - t196 * t27 + (qJD(4) * t90 - t64) * pkin(9)) * t173, -pkin(3) * t45 - t26 * t40 - t27 * t41 - t48 * t62 + (-t10 * t173 - t176 * t201 + (-t173 * t27 - t176 * t26) * qJD(4)) * pkin(9), -t173 * t175 * t31 + (-t173 * t248 + t175 * t249 - t81) * t68, t66 * t81 + t68 * t80 + (-t172 * t68 - t175 * t66) * t249 + (t275 - t272 + (t172 * t66 - t175 * t68) * qJD(5)) * t173, -t81 * t89 + (t250 * t89 + t31) * t176 + (t196 * t68 - t248 * t89 + t271) * t173, t172 * t173 * t32 + (t172 * t249 + t173 * t247 - t80) * t66, t80 * t89 + (-qJD(4) * t273 + t32) * t176 + (-t196 * t66 - t247 * t89 - t274) * t173, t227 * t89 - t270, t138 * t65 - t23 * t80 - t35 * t66 - t280 * t89 + (-t2 + (pkin(9) * t66 + t172 * t23) * qJD(4)) * t176 + (pkin(9) * t32 + t172 * t8 - t196 * t213 + t23 * t247) * t173, -t139 * t65 - t23 * t81 - t35 * t68 + t279 * t89 + (t1 + (pkin(9) * t68 + t175 * t23) * qJD(4)) * t176 + (-pkin(9) * t31 + t175 * t8 - t196 * t6 - t23 * t248) * t173, t138 * t31 - t139 * t32 - t213 * t81 + t6 * t80 + t280 * t68 + t279 * t66 + t217 * t249 + (-t1 * t172 - t175 * t2 + (-t172 * t213 - t175 * t6) * qJD(5)) * t173, t1 * t139 + t138 * t2 - t23 * t35 - t279 * t6 + t280 * t213 + (t173 * t8 + t23 * t249) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, -t90 ^ 2 + t92 ^ 2, -t64 + t293, -t281, t115 * t92 + t259, t107, t115 * t27 - t48 * t92 - t237, t48 * t90 - t298, 0, 0, t228 * t68 - t275, (-t31 - t284) * t175 + (-t32 - t282) * t172, t228 * t89 - t68 * t92 + t274, t273 * t66 - t272, -t172 * t89 ^ 2 + t66 * t92 + t271, -t89 * t92, -pkin(4) * t32 + t172 * t203 - t175 * t226 - t19 * t89 + t213 * t92 - t27 * t66, pkin(4) * t31 + t172 * t226 + t175 * t203 + t20 * t89 - t27 * t68 + t6 * t92, t19 * t68 + t20 * t66 + ((-t32 + t265) * pkin(10) + t301) * t175 + ((qJD(5) * t66 - t31) * pkin(10) + t295) * t172, -pkin(4) * t8 + t19 * t213 - t20 * t6 - t23 * t27 + (qJD(5) * t217 + t1 * t175 - t172 * t2) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, -t66 ^ 2 + t68 ^ 2, -t31 + t284, -t283, t282 - t32, t65, -t23 * t68 - t295, t23 * t66 - t301, 0, 0;];
tauc_reg = t5;

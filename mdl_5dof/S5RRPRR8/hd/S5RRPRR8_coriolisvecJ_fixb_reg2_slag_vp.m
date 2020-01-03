% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:24
% DurationCPUTime: 4.04s
% Computational Cost: add. (7506->358), mult. (19511->481), div. (0->0), fcn. (14596->8), ass. (0->196)
t156 = cos(qJ(5));
t153 = sin(qJ(5));
t199 = qJD(5) * t156;
t152 = sin(pkin(9));
t155 = sin(qJ(2));
t157 = cos(qJ(2));
t212 = cos(pkin(9));
t130 = t152 * t157 + t212 * t155;
t242 = t130 * qJD(2);
t110 = qJD(1) * t242;
t198 = qJD(1) * qJD(2);
t190 = t155 * t198;
t140 = t152 * t190;
t187 = t212 * t157;
t179 = qJD(1) * t187;
t111 = qJD(2) * t179 - t140;
t202 = qJD(1) * t155;
t123 = -t152 * t202 + t179;
t124 = t130 * qJD(1);
t154 = sin(qJ(4));
t236 = cos(qJ(4));
t191 = qJD(4) * t236;
t201 = qJD(4) * t154;
t172 = -t154 * t110 + t236 * t111 + t123 * t191 - t124 * t201;
t173 = -t154 * t123 - t236 * t124;
t197 = qJD(2) + qJD(4);
t65 = t153 * t197 - t156 * t173;
t211 = qJD(5) * t65;
t27 = t153 * t172 + t211;
t182 = t156 * t197;
t62 = -t153 * t173 - t182;
t222 = -t153 * t27 - t62 * t199;
t77 = t236 * t123 - t154 * t124;
t256 = qJD(5) - t77;
t266 = t153 * t256;
t246 = t65 * t266;
t258 = t156 * t77;
t200 = qJD(5) * t153;
t26 = -qJD(5) * t182 - t156 * t172 - t173 * t200;
t268 = -t156 * t26 + t258 * t62 + t222 - t246;
t146 = t212 * pkin(2) + pkin(3);
t235 = pkin(2) * t152;
t118 = t154 * t146 + t236 * t235;
t233 = t123 * pkin(7);
t223 = -qJ(3) - pkin(6);
t138 = t223 * t157;
t135 = qJD(1) * t138;
t128 = t212 * t135;
t137 = t223 * t155;
t134 = qJD(1) * t137;
t83 = -t152 * t134 + t128;
t164 = t83 - t233;
t232 = t124 * pkin(7);
t126 = t152 * t135;
t84 = t212 * t134 + t126;
t66 = t84 - t232;
t213 = t118 * qJD(4) - t154 * t66 + t236 * t164;
t219 = qJD(2) * pkin(2);
t129 = t134 + t219;
t79 = t212 * t129 + t126;
t56 = qJD(2) * pkin(3) - t232 + t79;
t80 = t152 * t129 - t128;
t60 = t80 + t233;
t33 = t154 * t56 + t236 * t60;
t267 = -t33 - t213;
t31 = t197 * pkin(8) + t33;
t147 = -t157 * pkin(2) - pkin(1);
t203 = qJD(1) * t147;
t136 = qJD(3) + t203;
t85 = -t123 * pkin(3) + t136;
t34 = -pkin(4) * t77 + pkin(8) * t173 + t85;
t175 = t153 * t31 - t156 * t34;
t265 = t256 * t175;
t245 = t77 * t197;
t264 = t172 - t245;
t226 = t173 ^ 2;
t227 = t77 ^ 2;
t263 = t226 - t227;
t10 = t153 * t34 + t156 * t31;
t188 = qJD(2) * t223;
t119 = t157 * qJD(3) + t155 * t188;
t97 = t119 * qJD(1);
t120 = -t155 * qJD(3) + t157 * t188;
t98 = t120 * qJD(1);
t61 = -t152 * t97 + t212 * t98;
t51 = -t111 * pkin(7) + t61;
t64 = t152 * t98 + t212 * t97;
t52 = -t110 * pkin(7) + t64;
t161 = -t154 * t51 - t56 * t191 + t60 * t201 - t236 * t52;
t186 = t236 * t110 + t154 * t111;
t254 = qJD(4) * t173;
t45 = t186 - t254;
t144 = pkin(2) * t190;
t86 = t110 * pkin(3) + t144;
t18 = t45 * pkin(4) - pkin(8) * t172 + t86;
t3 = -qJD(5) * t10 + t153 * t161 + t156 * t18;
t247 = t256 * t10 + t3;
t23 = t26 * t153;
t261 = -t23 + (t199 - t258) * t65;
t41 = t153 * t45;
t221 = t199 * t256 + t41;
t229 = t65 * t173;
t260 = -t256 * t258 + t221 + t229;
t32 = -t154 * t60 + t236 * t56;
t30 = -t197 * pkin(4) - t32;
t259 = t30 * t77;
t228 = t173 * t62;
t257 = t256 * t173;
t225 = t77 * t173;
t206 = t173 * qJD(2);
t255 = -t206 - t186;
t253 = -t85 * t77 + t161;
t29 = t30 * t199;
t189 = t154 * t52 - t236 * t51;
t8 = qJD(4) * t33 + t189;
t252 = -t10 * t173 + t8 * t153 + t29;
t28 = t30 * t200;
t251 = -t175 * t173 + t28;
t250 = t85 * t173 - t189;
t47 = -pkin(4) * t173 - t77 * pkin(8);
t249 = -0.2e1 * t198;
t2 = -qJD(5) * t175 + t153 * t18 - t156 * t161;
t248 = t2 + t265;
t117 = t236 * t146 - t154 * t235;
t103 = t117 * qJD(4);
t36 = t154 * t164 + t236 * t66;
t214 = t103 - t36;
t43 = t156 * t45;
t244 = t200 * t256 - t43;
t171 = t152 * t155 - t187;
t243 = qJD(2) * t171;
t241 = t124 ^ 2;
t158 = qJD(2) ^ 2;
t87 = t212 * t137 + t152 * t138;
t71 = -t130 * pkin(7) + t87;
t88 = t152 * t137 - t212 * t138;
t72 = -pkin(7) * t171 + t88;
t39 = t154 * t72 - t236 * t71;
t240 = t8 * t39;
t169 = t154 * t171;
t82 = t236 * t130 - t169;
t239 = t8 * t82;
t1 = t2 * t156;
t163 = t236 * t171;
t81 = t154 * t130 + t163;
t231 = t45 * t81;
t230 = t65 * t62;
t224 = t82 * t45;
t217 = t153 * t62;
t25 = t27 * t156;
t210 = t124 * t123;
t159 = qJD(1) ^ 2;
t209 = t157 * t159;
t208 = t158 * t155;
t207 = t158 * t157;
t70 = t212 * t119 + t152 * t120;
t204 = t155 ^ 2 - t157 ^ 2;
t149 = t155 * t219;
t148 = pkin(2) * t202;
t194 = t82 * t200;
t193 = t82 * t199;
t192 = t155 * t209;
t92 = t124 * pkin(3) + t148;
t183 = pkin(1) * t249;
t181 = t157 * t190;
t48 = t130 * t201 + t154 * t242 + t197 * t163;
t180 = -t30 * t48 + t239;
t178 = -t256 * t48 + t224;
t177 = t10 * t153 - t156 * t175;
t100 = pkin(3) * t171 + t147;
t38 = t81 * pkin(4) - t82 * pkin(8) + t100;
t40 = t154 * t71 + t236 * t72;
t19 = -t153 * t40 + t156 * t38;
t20 = t153 * t38 + t156 * t40;
t174 = t266 * t77 - t244;
t69 = -t152 * t119 + t212 * t120;
t115 = pkin(8) + t118;
t170 = -t103 * t256 - t115 * t45 - t259;
t93 = pkin(3) * t242 + t149;
t162 = -qJD(5) * t177 - t3 * t153 + t1;
t160 = pkin(7) * t243 + t69;
t121 = t123 ^ 2;
t114 = -pkin(4) - t117;
t54 = -pkin(7) * t242 + t70;
t49 = -qJD(4) * t169 + t130 * t191 - t154 * t243 + t236 * t242;
t37 = t47 + t92;
t21 = t49 * pkin(4) + t48 * pkin(8) + t93;
t17 = t153 * t47 + t156 * t32;
t16 = -t153 * t32 + t156 * t47;
t14 = t40 * qJD(4) + t154 * t54 - t236 * t160;
t13 = t154 * t160 + t71 * t191 - t72 * t201 + t236 * t54;
t12 = t153 * t37 + t156 * t36;
t11 = -t153 * t36 + t156 * t37;
t5 = -qJD(5) * t20 - t153 * t13 + t156 * t21;
t4 = qJD(5) * t19 + t156 * t13 + t153 * t21;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t181, t204 * t249, t207, -0.2e1 * t181, -t208, 0, -pkin(6) * t207 + t155 * t183, pkin(6) * t208 + t157 * t183, 0, 0, t111 * t130 - t124 * t243, -t130 * t110 - t111 * t171 - t123 * t243 - t124 * t242, -t171 * t158, t110 * t171 - t123 * t242, -t130 * t158, 0, t69 * qJD(2) + t147 * t110 - t123 * t149 + t136 * t242 + t148 * t243, -t70 * qJD(2) + t147 * t111 + t124 * t149 + t130 * t144 - t136 * t243, -t88 * t110 - t87 * t111 + t70 * t123 - t69 * t124 - t61 * t130 - t171 * t64 - t242 * t80 + t243 * t79, t61 * t87 + t64 * t88 + t79 * t69 + t80 * t70 + (t136 + t203) * t149, t172 * t82 + t173 * t48, -t172 * t81 + t173 * t49 - t48 * t77 - t224, -t48 * t197, -t49 * t77 + t231, -t49 * t197, 0, t100 * t45 - t14 * t197 + t85 * t49 - t77 * t93 + t86 * t81, t100 * t172 - t13 * t197 - t173 * t93 - t85 * t48 + t86 * t82, t13 * t77 - t14 * t173 + t161 * t81 + t172 * t39 + t32 * t48 - t33 * t49 - t40 * t45 + t239, t86 * t100 + t33 * t13 - t32 * t14 - t161 * t40 + t85 * t93 + t240, -t65 * t194 + (-t26 * t82 - t48 * t65) * t156, (t153 * t65 + t156 * t62) * t48 + (t23 - t25 + (-t156 * t65 + t217) * qJD(5)) * t82, t156 * t178 - t194 * t256 - t26 * t81 + t65 * t49, t62 * t193 + (t27 * t82 - t48 * t62) * t153, -t153 * t178 - t193 * t256 - t27 * t81 - t62 * t49, t256 * t49 + t231, t14 * t62 + t153 * t180 - t175 * t49 + t19 * t45 + t256 * t5 + t39 * t27 + t29 * t82 + t3 * t81, -t10 * t49 + t14 * t65 + t156 * t180 - t2 * t81 - t20 * t45 - t256 * t4 - t39 * t26 - t28 * t82, t19 * t26 - t20 * t27 - t4 * t62 - t5 * t65 + t177 * t48 + (-t2 * t153 - t3 * t156 + (-t10 * t156 - t153 * t175) * qJD(5)) * t82, t10 * t4 + t30 * t14 - t175 * t5 + t3 * t19 + t2 * t20 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, t204 * t159, 0, t192, 0, 0, t159 * pkin(1) * t155, pkin(1) * t209, 0, 0, -t210, -t121 + t241, -t140 + (t179 - t123) * qJD(2), t210, 0, 0, -t83 * qJD(2) + t123 * t148 - t136 * t124 + t61, t84 * qJD(2) - t136 * t123 - t124 * t148 - t64, (t80 + t83) * t124 + (t79 - t84) * t123 + (-t110 * t152 - t111 * t212) * pkin(2), -t79 * t83 - t80 * t84 + (-t136 * t202 + t152 * t64 + t212 * t61) * pkin(2), t225, t263, t264, -t225, t255, 0, -t213 * qJD(2) + t267 * qJD(4) + t92 * t77 + t250, t173 * t92 - t214 * t197 + t253, -t117 * t172 - t118 * t45 + (t214 + t32) * t77 + t267 * t173, -t8 * t117 - t118 * t161 - t213 * t32 + t214 * t33 - t85 * t92, t261, t268, t260, t266 * t62 - t25, t174 - t228, t257, -t11 * t256 + t114 * t27 + t213 * t62 + (-qJD(5) * t115 * t256 - t8) * t156 + t170 * t153 + t251, -t114 * t26 + (t115 * t200 + t12) * t256 + t213 * t65 + t170 * t156 + t252, t11 * t65 + t12 * t62 + t1 + (-t103 * t62 - t115 * t27 - t77 * t175 + (t115 * t65 + t175) * qJD(5)) * t156 + (t10 * t77 + t103 * t65 - t115 * t26 - t3 + (t115 * t62 - t10) * qJD(5)) * t153, t8 * t114 - (-t103 * t153 - t11) * t175 + t213 * t30 + (t103 * t156 - t12) * t10 + t162 * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t124 * qJD(2), -t140 + (t179 + t123) * qJD(2), -t121 - t241, -t80 * t123 + t79 * t124 + t144, 0, 0, 0, 0, 0, 0, t186 - t206 - 0.2e1 * t254, t172 + t245, -t226 - t227, -t173 * t32 - t33 * t77 + t86, 0, 0, 0, 0, 0, 0, t174 + t228, -t156 * t256 ^ 2 + t229 - t41, (t62 * t77 + t26) * t156 + t246 + t222, t248 * t153 + t247 * t156 + t173 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, t263, t264, -t225, t255, 0, t33 * qJD(2) + t250, t32 * t197 + t253, 0, 0, t261, t268, t260, t217 * t256 - t25, -t256 * t266 - t228 + t43, t257, -pkin(4) * t27 - pkin(8) * t221 - t153 * t259 - t8 * t156 - t16 * t256 - t33 * t62 + t251, pkin(4) * t26 + t244 * pkin(8) + t17 * t256 - t258 * t30 - t33 * t65 + t252, t16 * t65 + t17 * t62 + t1 + (t265 + (-t27 + t211) * pkin(8)) * t156 + ((qJD(5) * t62 - t26) * pkin(8) - t247) * t153, -t8 * pkin(4) + pkin(8) * t162 - t10 * t17 + t16 * t175 - t30 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, -t62 ^ 2 + t65 ^ 2, t256 * t62 - t26, -t230, t256 * t65 - t27, t45, -t30 * t65 + t247, t30 * t62 - t248, 0, 0;];
tauc_reg = t6;

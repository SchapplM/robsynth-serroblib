% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:48:02
% EndTime: 2019-12-31 19:48:08
% DurationCPUTime: 2.35s
% Computational Cost: add. (1956->364), mult. (4209->473), div. (0->0), fcn. (2689->10), ass. (0->205)
t256 = pkin(3) + pkin(6);
t148 = sin(qJ(2));
t144 = sin(pkin(8));
t145 = cos(pkin(8));
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t88 = t150 * t144 + t147 * t145;
t165 = t88 * t148;
t75 = t88 * qJD(5);
t247 = -qJD(1) * t165 - t75;
t216 = t148 * qJD(1);
t114 = qJD(5) + t216;
t151 = cos(qJ(2));
t223 = qJD(1) * t151;
t202 = t145 * t223;
t218 = t144 * qJD(2);
t83 = t202 + t218;
t203 = t144 * t223;
t217 = t145 * qJD(2);
t85 = -t203 + t217;
t180 = t147 * t83 - t150 * t85;
t260 = t114 * t180;
t149 = sin(qJ(1));
t152 = cos(qJ(1));
t187 = g(1) * t152 + g(2) * t149;
t123 = pkin(6) * t216;
t259 = qJD(3) + t123;
t258 = -qJD(5) + t114;
t221 = qJD(2) * t151;
t249 = pkin(2) + qJ(4);
t141 = qJD(2) * qJ(3);
t124 = pkin(6) * t223;
t95 = pkin(3) * t223 + t124;
t73 = qJD(4) + t141 + t95;
t257 = t249 * t221 + (qJD(4) - t73) * t148;
t214 = qJD(1) * qJD(2);
t201 = t148 * t214;
t212 = t151 * qJDD(1);
t50 = t144 * qJDD(2) + (-t201 + t212) * t145;
t193 = t145 * qJDD(2) - t144 * t212;
t51 = t144 * t201 + t193;
t9 = -qJD(5) * t180 + t147 * t51 + t150 * t50;
t255 = pkin(7) * t148;
t254 = g(1) * t149;
t251 = g(2) * t152;
t250 = g(3) * t148;
t137 = g(3) * t151;
t134 = t151 * pkin(2);
t248 = -pkin(7) - t249;
t113 = pkin(2) * t201;
t243 = qJ(3) * t151;
t178 = qJ(4) * t148 - t243;
t215 = t148 * qJD(3);
t155 = qJD(2) * t178 - t151 * qJD(4) - t215;
t130 = t148 * qJ(3);
t196 = -pkin(1) - t130;
t162 = -t151 * t249 + t196;
t16 = qJD(1) * t155 + qJDD(1) * t162 + t113;
t198 = t151 * t214;
t213 = t148 * qJDD(1);
t164 = t198 + t213;
t112 = pkin(6) * t198;
t120 = pkin(6) * t213;
t197 = qJDD(3) + t112 + t120;
t31 = pkin(3) * t164 - qJD(2) * qJD(4) - qJDD(2) * t249 + t197;
t7 = t144 * t31 + t145 * t16;
t222 = qJD(2) * t148;
t126 = pkin(2) * t222;
t41 = t126 + t155;
t96 = t256 * t221;
t21 = t144 * t96 + t145 * t41;
t55 = t162 * qJD(1);
t227 = pkin(3) * t216 + t259;
t64 = -qJD(2) * t249 + t227;
t19 = t144 * t64 + t145 * t55;
t127 = pkin(2) * t216;
t69 = qJD(1) * t178 + t127;
t30 = t144 * t95 + t145 * t69;
t103 = t256 * t148;
t225 = t134 + t130;
t232 = t151 * qJ(4);
t183 = t225 + t232;
t82 = -pkin(1) - t183;
t38 = t103 * t144 + t145 * t82;
t205 = t145 * t216;
t219 = qJD(5) * t150;
t220 = qJD(5) * t147;
t239 = t147 * t144;
t246 = -t144 * t220 + t145 * t219 + t150 * t205 - t216 * t239;
t34 = t147 * t85 + t150 * t83;
t245 = t34 * t114;
t244 = pkin(6) * qJDD(2);
t242 = qJDD(2) * pkin(2);
t142 = t148 ^ 2;
t154 = qJD(1) ^ 2;
t241 = t142 * t154;
t240 = t145 * t151;
t238 = t148 * t149;
t237 = t148 * t152;
t236 = t148 * t154;
t138 = pkin(8) + qJ(5);
t128 = sin(t138);
t235 = t149 * t128;
t129 = cos(t138);
t234 = t149 * t129;
t233 = t149 * t151;
t231 = t151 * t152;
t230 = t152 * t128;
t229 = t152 * t129;
t208 = -pkin(4) * t145 - pkin(3);
t228 = -t208 * t216 + t259;
t104 = t256 * t151;
t143 = t151 ^ 2;
t224 = t142 - t143;
t211 = t151 * t236;
t121 = pkin(6) * t212;
t139 = qJDD(2) * qJ(3);
t140 = qJD(2) * qJD(3);
t210 = t121 + t139 + t140;
t209 = -g(1) * t237 - g(2) * t238 + t137;
t6 = -t144 * t16 + t145 * t31;
t4 = pkin(4) * t164 - t51 * pkin(7) + t6;
t5 = -pkin(7) * t50 + t7;
t207 = -t147 * t5 + t150 * t4;
t206 = t256 * qJD(2);
t200 = t144 * t213;
t199 = t145 * t213;
t20 = -t144 * t41 + t145 * t96;
t18 = -t144 * t55 + t145 * t64;
t29 = -t144 * t69 + t145 * t95;
t194 = -qJD(2) * pkin(2) + qJD(3);
t192 = pkin(1) * t152 + pkin(2) * t231 + pkin(6) * t149 + qJ(3) * t237;
t191 = -t120 - t209;
t177 = -t145 * t150 + t239;
t87 = qJDD(5) + t164;
t190 = t114 * t247 - t177 * t87;
t189 = qJD(1) * t206;
t153 = qJD(2) ^ 2;
t188 = pkin(6) * t153 + t251;
t186 = -t251 + t254;
t185 = t7 * t144 + t6 * t145;
t184 = t147 * t4 + t150 * t5;
t10 = pkin(4) * t216 - pkin(7) * t85 + t18;
t11 = -pkin(7) * t83 + t19;
t2 = t10 * t150 - t11 * t147;
t3 = t10 * t147 + t11 * t150;
t91 = t145 * t103;
t23 = t148 * pkin(4) + t91 + (pkin(7) * t151 - t82) * t144;
t28 = -pkin(7) * t240 + t38;
t182 = -t147 * t28 + t150 * t23;
t181 = t147 * t23 + t150 * t28;
t102 = -t124 - t141;
t99 = t123 + t194;
t179 = t102 * t148 + t151 * t99;
t176 = pkin(3) * t212 + qJDD(4) + t210;
t175 = pkin(4) * t151 - t144 * t255;
t174 = t196 - t134;
t173 = t187 * t151;
t172 = -t114 * t246 - t88 * t87;
t81 = t174 * qJD(1);
t171 = t216 * t81 + qJDD(3) - t191;
t170 = -0.2e1 * pkin(1) * t214 - t244;
t98 = t248 * t145;
t169 = pkin(7) * t205 + qJD(4) * t144 - qJD(5) * t98 + t30;
t97 = t248 * t144;
t168 = qJD(1) * t175 + qJD(4) * t145 + qJD(5) * t97 + t29;
t8 = -t147 * t50 + t150 * t51 - t219 * t83 - t220 * t85;
t167 = (-t144 * t18 + t145 * t19) * t148;
t166 = -qJ(3) * t221 - t215;
t65 = t177 * t151;
t33 = -t148 * t189 + t176;
t163 = -t151 * t33 + t187;
t161 = 0.2e1 * qJDD(1) * pkin(1) - t188;
t100 = -pkin(1) - t225;
t160 = t244 + (-qJD(1) * t100 - t81) * qJD(2);
t159 = -t173 - t250;
t158 = t33 + t159;
t32 = qJD(1) * t166 + qJDD(1) * t174 + t113;
t71 = t126 + t166;
t157 = qJD(1) * t71 + qJDD(1) * t100 + t188 + t32;
t56 = pkin(6) * t201 - t210;
t68 = t197 - t242;
t156 = qJD(2) * t179 + t68 * t148 - t56 * t151;
t135 = t152 * pkin(6);
t119 = pkin(4) * t144 + qJ(3);
t117 = g(1) * t233;
t111 = qJ(3) * t231;
t109 = qJ(3) * t233;
t94 = t148 * t206;
t92 = -qJ(3) * t223 + t127;
t74 = pkin(4) * t240 + t104;
t66 = t88 * t151;
t63 = -t148 * t235 + t229;
t62 = t148 * t234 + t230;
t61 = t148 * t230 + t234;
t60 = t148 * t229 - t235;
t58 = (-pkin(6) + t208) * t222;
t39 = pkin(4) * t83 + t73;
t37 = -t144 * t82 + t91;
t25 = t151 * t75 - t177 * t222;
t24 = qJD(2) * t165 + qJD(5) * t65;
t15 = t217 * t255 + t21;
t14 = t50 * pkin(4) + t33;
t12 = qJD(2) * t175 + t20;
t1 = [qJDD(1), t186, t187, qJDD(1) * t142 + 0.2e1 * t148 * t198, 0.2e1 * t148 * t212 - 0.2e1 * t214 * t224, qJDD(2) * t148 + t151 * t153, qJDD(2) * t151 - t148 * t153, 0, t148 * t170 + t151 * t161 + t117, t170 * t151 + (-t161 - t254) * t148, (t142 + t143) * qJDD(1) * pkin(6) + t156 - t187, t148 * t160 + t151 * t157 - t117, t160 * t151 + (-t157 + t254) * t148, pkin(6) * t156 - g(1) * t135 - g(2) * t192 + t32 * t100 - t174 * t254 + t81 * t71, t104 * t50 - t94 * t83 + (qJD(1) * t37 + t18) * t221 - t163 * t145 + (t20 * qJD(1) + t37 * qJDD(1) + t144 * t186 - t217 * t73 + t6) * t148, t104 * t51 - t94 * t85 + (-qJD(1) * t38 - t19) * t221 + t163 * t144 + (-qJD(1) * t21 - qJDD(1) * t38 + t145 * t186 + t218 * t73 - t7) * t148, -t20 * t85 - t21 * t83 - t37 * t51 - t38 * t50 + t117 + qJD(2) * t167 + (t144 * t6 - t145 * t7 - t251) * t151, t7 * t38 + t19 * t21 + t6 * t37 + t18 * t20 + t33 * t104 - t73 * t94 - g(1) * (t152 * pkin(3) + t135) - g(2) * (qJ(4) * t231 + t192) + (-g(1) * (t174 - t232) - g(2) * pkin(3)) * t149, -t180 * t24 - t66 * t8, -t180 * t25 - t24 * t34 + t65 * t8 + t66 * t9, t114 * t24 + t148 * t8 - t180 * t221 - t66 * t87, t114 * t25 - t148 * t9 - t221 * t34 + t65 * t87, t114 * t221 + t148 * t87, (t150 * t12 - t147 * t15) * t114 + t182 * t87 + t207 * t148 + t2 * t221 + t58 * t34 + t74 * t9 - t14 * t65 - t39 * t25 - g(1) * t63 - g(2) * t61 + (-t114 * t181 - t148 * t3) * qJD(5), -(t147 * t12 + t150 * t15) * t114 - t181 * t87 - t184 * t148 - t3 * t221 - t58 * t180 + t74 * t8 - t14 * t66 + t39 * t24 + g(1) * t62 - g(2) * t60 + (-t114 * t182 - t148 * t2) * qJD(5); 0, 0, 0, -t211, t224 * t154, t213, t212, qJDD(2), pkin(1) * t236 + t191, t250 - t121 + (pkin(1) * t154 + t187) * t151, (-pkin(2) * t148 + t243) * qJDD(1) + ((-t102 - t141) * t148 + (t194 - t99) * t151) * qJD(1), -t223 * t92 + t171 - 0.2e1 * t242, t121 + 0.2e1 * t139 + 0.2e1 * t140 + (qJD(1) * t92 - g(3)) * t148 + (qJD(1) * t81 - t187) * t151, -t56 * qJ(3) - t102 * qJD(3) - t68 * pkin(2) - t81 * t92 - g(1) * (-pkin(2) * t237 + t111) - g(2) * (-pkin(2) * t238 + t109) - g(3) * t225 - t179 * qJD(1) * pkin(6), -t249 * t199 + qJ(3) * t50 + t227 * t83 + t158 * t144 + (-t145 * t257 - t148 * t29 - t18 * t151) * qJD(1), t249 * t200 + qJ(3) * t51 + t227 * t85 + t158 * t145 + (t144 * t257 + t148 * t30 + t19 * t151) * qJD(1), t29 * t85 + t30 * t83 + (qJD(4) * t85 - t19 * t216 + t249 * t51 - t6) * t145 + (qJD(4) * t83 + t18 * t216 + t249 * t50 - t7) * t144 - t209, t33 * qJ(3) - t19 * t30 - t18 * t29 - g(1) * t111 - g(2) * t109 - g(3) * t183 + t227 * t73 + (-t144 * t19 - t145 * t18) * qJD(4) + (t148 * t187 - t185) * t249, -t177 * t8 - t180 * t247, t177 * t9 + t180 * t246 - t247 * t34 - t8 * t88, t180 * t223 + t190, t223 * t34 + t172, -t114 * t223, (-t147 * t97 + t150 * t98) * t87 + t119 * t9 + t14 * t88 - t2 * t223 + t246 * t39 + t228 * t34 + (t147 * t169 - t150 * t168) * t114 + t159 * t128, -(t147 * t98 + t150 * t97) * t87 + t119 * t8 - t14 * t177 + t3 * t223 + t247 * t39 - t228 * t180 + (t147 * t168 + t150 * t169) * t114 + t159 * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, qJDD(2) + t211, -t153 - t241, qJD(2) * t102 + t112 + t171 - t242, t199 - t144 * t241 + (-t83 + t202) * qJD(2), -t200 - t145 * t241 + (-t85 - t203) * qJD(2), -t144 * t50 - t145 * t51 + (t144 * t85 - t145 * t83) * t216, qJD(1) * t167 - t73 * qJD(2) + t185 + t209, 0, 0, 0, 0, 0, -qJD(2) * t34 + t190, qJD(2) * t180 + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216 * t85 + t50, (-t83 + t218) * t216 + t193, -t83 ^ 2 - t85 ^ 2, t18 * t85 + t19 * t83 - t173 + (-g(3) - t189) * t148 + t176, 0, 0, 0, 0, 0, t9 - t260, t8 - t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180 * t34, t180 ^ 2 - t34 ^ 2, t8 + t245, -t9 - t260, t87, -g(1) * t60 - g(2) * t62 + t129 * t137 + t180 * t39 + t258 * t3 + t207, g(1) * t61 - g(2) * t63 - t128 * t137 + t2 * t258 + t39 * t34 - t184;];
tau_reg = t1;

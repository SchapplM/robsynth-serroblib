% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:30:01
% EndTime: 2019-03-09 04:30:06
% DurationCPUTime: 2.36s
% Computational Cost: add. (4400->350), mult. (10303->467), div. (0->0), fcn. (6629->8), ass. (0->175)
t163 = cos(qJ(3));
t207 = t163 * qJD(1);
t142 = -qJD(4) + t207;
t160 = sin(qJ(4));
t161 = sin(qJ(3));
t215 = qJD(1) * t161;
t191 = t160 * t215;
t162 = cos(qJ(4));
t209 = t162 * qJD(3);
t119 = t191 - t209;
t210 = t160 * qJD(3);
t121 = t162 * t215 + t210;
t156 = sin(pkin(10));
t158 = cos(pkin(10));
t70 = t158 * t119 + t156 * t121;
t256 = t142 * t70;
t115 = t156 * t162 + t158 * t160;
t104 = t115 * qJD(4);
t235 = -t115 * t207 + t104;
t211 = qJD(4) * t162;
t212 = qJD(4) * t160;
t105 = -t156 * t212 + t158 * t211;
t197 = t160 * t207;
t226 = t158 * t162;
t234 = t156 * t197 - t207 * t226 + t105;
t192 = t163 * t210;
t195 = t161 * t211;
t255 = t192 + t195;
t176 = -t156 * t119 + t158 * t121;
t254 = t176 ^ 2;
t248 = -qJ(5) - pkin(8);
t185 = qJD(4) * t248;
t208 = t162 * qJD(5);
t101 = t160 * t185 + t208;
t167 = -t160 * qJD(5) + t162 * t185;
t222 = t162 * t163;
t172 = pkin(4) * t161 - qJ(5) * t222;
t178 = pkin(3) * t161 - pkin(8) * t163;
t123 = t178 * qJD(1);
t147 = sin(pkin(9)) * pkin(1) + pkin(7);
t131 = t147 * qJD(1);
t252 = t163 * qJD(2) - t161 * t131;
t183 = t162 * t123 - t160 * t252;
t53 = t172 * qJD(1) + t183;
t231 = t160 * t123 + t162 * t252;
t58 = -qJ(5) * t197 + t231;
t244 = (t167 - t53) * t158 + (-t101 + t58) * t156;
t83 = -qJD(3) * pkin(3) - t252;
t67 = t119 * pkin(4) + qJD(5) + t83;
t25 = t70 * pkin(5) - qJ(6) * t176 + t67;
t253 = t25 * t176;
t152 = t161 * qJD(2);
t95 = t163 * t131 + t152;
t251 = -t95 + (-t197 + t212) * pkin(4);
t154 = t161 ^ 2;
t175 = qJD(1) * t154 - t142 * t163;
t196 = t161 * t212;
t250 = -t142 * t196 - t175 * t209;
t204 = qJD(3) * qJD(4);
t80 = t255 * qJD(1) + t160 * t204;
t149 = -cos(pkin(9)) * pkin(1) - pkin(2);
t111 = -t163 * pkin(3) - t161 * pkin(8) + t149;
t89 = t111 * qJD(1);
t241 = t160 * t89;
t84 = qJD(3) * pkin(8) + t95;
t56 = t162 * t84 + t241;
t47 = -t119 * qJ(5) + t56;
t44 = t158 * t47;
t55 = -t160 * t84 + t162 * t89;
t46 = -t121 * qJ(5) + t55;
t19 = t156 * t46 + t44;
t249 = t19 * t176;
t124 = t178 * qJD(3);
t110 = qJD(1) * t124;
t87 = t252 * qJD(3);
t186 = -t162 * t110 + t160 * t87;
t166 = -t56 * qJD(4) - t186;
t205 = qJD(1) * qJD(3);
t187 = t161 * t205;
t190 = t163 * t209;
t79 = qJD(1) * t190 - qJD(4) * t191 + t162 * t204;
t16 = pkin(4) * t187 - t79 * qJ(5) - t121 * qJD(5) + t166;
t203 = t160 * t110 + t162 * t87 + t89 * t211;
t171 = -t84 * t212 + t203;
t22 = -t80 * qJ(5) - t119 * qJD(5) + t171;
t3 = -t156 * t22 + t158 * t16;
t4 = t156 * t16 + t158 * t22;
t27 = t156 * t53 + t158 * t58;
t23 = qJ(6) * t215 + t27;
t64 = t158 * t101 + t156 * t167;
t247 = t23 - t64;
t246 = pkin(5) * t215 - t244;
t122 = t147 * t222;
t134 = t161 * t147;
t218 = t162 * t124 + t210 * t134;
t30 = -t161 * t208 + t172 * qJD(3) + (-t122 + (qJ(5) * t161 - t111) * t160) * qJD(4) + t218;
t223 = t161 * t162;
t233 = t111 * t211 + t160 * t124;
t35 = (-qJ(5) * qJD(4) - qJD(3) * t147) * t223 + (-qJD(5) * t161 + (-qJ(5) * qJD(3) - qJD(4) * t147) * t163) * t160 + t233;
t9 = t156 * t30 + t158 * t35;
t245 = -t235 * pkin(5) + t234 * qJ(6) + t115 * qJD(6) - t251;
t43 = -t142 * pkin(4) + t46;
t15 = t156 * t43 + t44;
t99 = t162 * t111;
t62 = -qJ(5) * t223 + t99 + (-t147 * t160 - pkin(4)) * t163;
t225 = t160 * t161;
t232 = t160 * t111 + t122;
t66 = -qJ(5) * t225 + t232;
t34 = t156 * t62 + t158 * t66;
t243 = t156 * t47;
t242 = t160 * t83;
t240 = t162 * t83;
t239 = t163 * t80;
t238 = t79 * t160;
t213 = qJD(3) * t163;
t88 = qJD(3) * t152 + t131 * t213;
t237 = t88 * t160;
t236 = t88 * t162;
t230 = t119 * t142;
t229 = t119 * t161;
t228 = t121 * t142;
t227 = t142 * t162;
t224 = t160 * t163;
t164 = qJD(3) ^ 2;
t221 = t164 * t161;
t220 = t164 * t163;
t20 = t158 * t46 - t243;
t219 = qJD(6) - t20;
t217 = pkin(4) * t225 + t134;
t216 = -t163 ^ 2 + t154;
t132 = qJD(1) * t149;
t214 = qJD(3) * t161;
t202 = qJ(6) * t187 + t4;
t199 = t255 * pkin(4) + t147 * t213;
t198 = -t162 * pkin(4) - pkin(3);
t189 = t248 * t160;
t50 = t156 * t79 + t158 * t80;
t184 = t121 * t214 - t79 * t163;
t182 = t142 * t147 + t84;
t51 = -t156 * t80 + t158 * t79;
t133 = t248 * t162;
t75 = -t156 * t133 - t158 * t189;
t76 = -t158 * t133 + t156 * t189;
t181 = -t76 * t50 + t75 * t51 - t64 * t70;
t179 = t142 * t195;
t57 = t80 * pkin(4) + t88;
t177 = -t70 ^ 2 - t254;
t8 = -t156 * t35 + t158 * t30;
t14 = t158 * t43 - t243;
t33 = -t156 * t66 + t158 * t62;
t173 = 0.2e1 * qJD(3) * t132;
t2 = -pkin(5) * t187 - t3;
t60 = -t105 * t161 - t156 * t190 - t158 * t192;
t61 = t161 * t104 + t156 * t192 - t158 * t190;
t92 = t115 * t161;
t93 = -t156 * t225 + t158 * t223;
t170 = -t176 * t60 - t93 * t50 + t92 * t51 + t61 * t70;
t169 = t175 * t160;
t5 = t50 * pkin(5) - t51 * qJ(6) - qJD(6) * t176 + t57;
t165 = qJD(1) ^ 2;
t148 = -t158 * pkin(4) - pkin(5);
t145 = t156 * pkin(4) + qJ(6);
t114 = t156 * t160 - t226;
t68 = t114 * pkin(5) - t115 * qJ(6) + t198;
t48 = t92 * pkin(5) - t93 * qJ(6) + t217;
t36 = t121 * pkin(4) + pkin(5) * t176 + qJ(6) * t70;
t31 = t163 * pkin(5) - t33;
t29 = -t163 * qJ(6) + t34;
t21 = -t60 * pkin(5) + t61 * qJ(6) - t93 * qJD(6) + t199;
t11 = -t142 * qJ(6) + t15;
t10 = t142 * pkin(5) + qJD(6) - t14;
t7 = -pkin(5) * t214 - t8;
t6 = qJ(6) * t214 - t163 * qJD(6) + t9;
t1 = -t142 * qJD(6) + t202;
t12 = [0, 0, 0, 0, 0.2e1 * t163 * t187, -0.2e1 * t216 * t205, t220, -t221, 0, -t147 * t220 + t161 * t173, t147 * t221 + t163 * t173, t79 * t223 + (t190 - t196) * t121 (-t119 * t162 - t121 * t160) * t213 + (-t238 - t162 * t80 + (t119 * t160 - t121 * t162) * qJD(4)) * t161, t184 - t250, t179 + t239 + (-t169 - t229) * qJD(3) (-t142 - t207) * t214 -(-t111 * t212 + t218) * t142 + ((t119 * t147 + t242) * qJD(3) + (t182 * t162 + t241) * qJD(4) + t186) * t163 + (t83 * t211 + t147 * t80 + t237 + ((-t147 * t224 + t99) * qJD(1) + t55) * qJD(3)) * t161, t233 * t142 + (-t182 * t212 + (t121 * t147 + t240) * qJD(3) + t203) * t163 + (-t83 * t212 + t147 * t79 + t236 + (-t232 * qJD(1) - t147 * t227 - t56) * qJD(3)) * t161, t14 * t61 + t15 * t60 - t176 * t8 - t3 * t93 - t33 * t51 - t34 * t50 - t4 * t92 - t9 * t70, t14 * t8 + t15 * t9 + t67 * t199 + t57 * t217 + t3 * t33 + t4 * t34, t7 * t142 + t2 * t163 + t21 * t70 - t25 * t60 + t48 * t50 + t5 * t92 + (-qJD(1) * t31 - t10) * t214, -t1 * t92 - t10 * t61 + t11 * t60 + t176 * t7 + t2 * t93 - t29 * t50 + t31 * t51 - t6 * t70, -t1 * t163 - t6 * t142 - t21 * t176 + t25 * t61 - t48 * t51 - t5 * t93 + (qJD(1) * t29 + t11) * t214, t1 * t29 + t10 * t7 + t11 * t6 + t2 * t31 + t25 * t21 + t5 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, -t220, 0, 0, 0, 0, 0, t179 - t239 + (-t169 + t229) * qJD(3), t184 + t250, t170, t14 * t60 - t15 * t61 - t57 * t163 + t67 * t214 - t3 * t92 + t4 * t93, -t60 * t142 - t163 * t50 + (-qJD(1) * t92 + t70) * t214, t170, t61 * t142 + t163 * t51 + (qJD(1) * t93 - t176) * t214, t1 * t93 - t10 * t60 - t11 * t61 - t5 * t163 + t2 * t92 + t25 * t214; 0, 0, 0, 0, -t163 * t165 * t161, t216 * t165, 0, 0, 0, t95 * qJD(3) - t132 * t215 - t88, -t132 * t207, -t121 * t227 + t238 (t79 + t230) * t162 + (-t80 + t228) * t160, -t142 * t211 + (t142 * t222 + (-t121 + t210) * t161) * qJD(1), t142 * t212 + (-t142 * t224 + (t119 + t209) * t161) * qJD(1), t142 * t215, -pkin(3) * t80 - t236 + t183 * t142 - t95 * t119 + (pkin(8) * t227 + t242) * qJD(4) + (-t55 * t161 + (-pkin(8) * t214 - t163 * t83) * t160) * qJD(1), -pkin(3) * t79 + t237 - t231 * t142 - t95 * t121 + (-t160 * pkin(8) * t142 + t240) * qJD(4) + (-t83 * t222 + (-pkin(8) * t209 + t56) * t161) * qJD(1), -t4 * t114 - t3 * t115 - t234 * t14 - t235 * t15 - t176 * t244 + t27 * t70 + t181, t4 * t76 - t3 * t75 + t57 * t198 + t251 * t67 + (t64 - t27) * t15 + t244 * t14, t5 * t114 + t68 * t50 - t245 * t70 + t235 * t25 + t246 * t142 + (-qJD(3) * t75 + t10) * t215, -t1 * t114 + t234 * t10 - t235 * t11 + t2 * t115 + t176 * t246 + t23 * t70 + t181, -t5 * t115 - t68 * t51 + t245 * t176 - t234 * t25 + t247 * t142 + (qJD(3) * t76 - t11) * t215, t1 * t76 + t246 * t10 - t247 * t11 + t2 * t75 - t245 * t25 + t5 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t119, -t119 ^ 2 + t121 ^ 2, t79 - t230, -t228 - t80, t187, -t83 * t121 - t56 * t142 + t166, t83 * t119 - t55 * t142 - t171, t15 * t176 - t249 + (-t156 * t50 - t158 * t51) * pkin(4) + (-t14 + t20) * t70, t14 * t19 - t15 * t20 + (-t121 * t67 + t156 * t4 + t158 * t3) * pkin(4), -t19 * t142 - t253 - t36 * t70 + (pkin(5) - t148) * t187 + t3, t11 * t176 - t145 * t50 + t148 * t51 - t249 + (t10 - t219) * t70, t145 * t187 - t25 * t70 + t36 * t176 + (-0.2e1 * qJD(6) + t20) * t142 + t202, t1 * t145 - t10 * t19 + t219 * t11 + t2 * t148 - t25 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t14 * t176 + t15 * t70 + t57, -t142 * t176 + t50, t177, -t51 - t256, -t10 * t176 + t11 * t70 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 * t70 - t187, t51 - t256, -t142 ^ 2 - t254, t11 * t142 + t2 + t253;];
tauc_reg  = t12;

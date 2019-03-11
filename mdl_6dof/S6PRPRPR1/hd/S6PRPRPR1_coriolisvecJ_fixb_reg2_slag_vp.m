% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:57
% EndTime: 2019-03-08 19:28:04
% DurationCPUTime: 2.92s
% Computational Cost: add. (4812->346), mult. (12196->488), div. (0->0), fcn. (9645->12), ass. (0->191)
t134 = sin(qJ(4));
t189 = qJD(4) * t134;
t129 = sin(pkin(11));
t130 = sin(pkin(6));
t131 = cos(pkin(11));
t135 = sin(qJ(2));
t138 = cos(qJ(2));
t91 = (t129 * t138 + t131 * t135) * t130;
t86 = qJD(1) * t91;
t242 = pkin(4) * t189 - t86;
t128 = sin(pkin(12));
t137 = cos(qJ(4));
t121 = pkin(2) * t129 + pkin(8);
t194 = qJ(5) + t121;
t165 = qJD(4) * t194;
t145 = -t134 * qJD(5) - t137 * t165;
t202 = cos(pkin(12));
t170 = t202 * t137;
t146 = -t128 * t134 + t170;
t183 = t137 * qJD(5);
t83 = -t134 * t165 + t183;
t192 = qJD(1) * t130;
t174 = t135 * t192;
t112 = t129 * t174;
t173 = t138 * t192;
t89 = t131 * t173 - t112;
t219 = -t128 * t145 + t146 * t89 - t202 * t83;
t107 = t128 * t137 + t202 * t134;
t100 = t107 * qJD(4);
t103 = t146 * qJD(4);
t241 = -pkin(5) * t100 + pkin(9) * t103 - t242;
t132 = cos(pkin(6));
t117 = qJD(1) * t132 + qJD(3);
t188 = qJD(4) * t137;
t90 = (t129 * t135 - t131 * t138) * t130;
t88 = qJD(2) * t90;
t79 = qJD(1) * t88;
t203 = t117 * t188 - t137 * t79;
t111 = qJD(2) * pkin(2) + t173;
t163 = t131 * t174;
t76 = t129 * t111 + t163;
t73 = qJD(2) * pkin(8) + t76;
t30 = -t73 * t189 + t203;
t24 = (-qJ(5) * t189 + t183) * qJD(2) + t30;
t172 = t202 * t24;
t164 = -qJD(2) * qJD(5) + t79;
t166 = qJ(5) * qJD(2) + t73;
t197 = t134 * t117;
t53 = t166 * t137 + t197;
t236 = -t53 * qJD(4) + t164 * t134;
t10 = t236 * t128 + t172;
t133 = sin(qJ(6));
t136 = cos(qJ(6));
t44 = t202 * t53;
t113 = t137 * t117;
t52 = -t166 * t134 + t113;
t47 = qJD(4) * pkin(4) + t52;
t18 = t128 * t47 + t44;
t16 = qJD(4) * pkin(9) + t18;
t101 = t107 * qJD(2);
t175 = -pkin(4) * t137 - pkin(3);
t75 = t131 * t111 - t112;
t64 = t175 * qJD(2) + qJD(5) - t75;
t116 = qJD(2) * t170;
t191 = qJD(2) * t134;
t98 = t128 * t191 - t116;
t33 = t98 * pkin(5) - t101 * pkin(9) + t64;
t156 = t133 * t16 - t136 * t33;
t182 = qJD(2) * qJD(4);
t171 = t134 * t182;
t78 = (t129 * t173 + t163) * qJD(2);
t68 = pkin(4) * t171 + t78;
t92 = qJD(2) * t100;
t115 = t128 * t171;
t93 = qJD(4) * t116 - t115;
t32 = pkin(5) * t92 - pkin(9) * t93 + t68;
t1 = -qJD(6) * t156 + t136 * t10 + t133 * t32;
t96 = qJD(6) + t98;
t240 = t156 * t96 + t1;
t6 = t133 * t33 + t136 * t16;
t2 = -qJD(6) * t6 - t133 * t10 + t136 * t32;
t239 = t6 * t96 + t2;
t169 = t133 * t96;
t82 = qJD(4) * t133 + t101 * t136;
t238 = t82 * t169;
t187 = qJD(6) * t133;
t150 = -t136 * t92 + t96 * t187;
t237 = qJD(2) * t86 - t78;
t199 = t103 * t136;
t235 = -t107 * t150 + t96 * t199;
t234 = t101 ^ 2;
t69 = t132 * t137 - t134 * t91;
t70 = t132 * t134 + t137 * t91;
t34 = t128 * t70 - t202 * t69;
t9 = t128 * t24 - t202 * t236;
t233 = t34 * t9;
t105 = t194 * t137;
t167 = t194 * t134;
t61 = t128 * t105 + t202 * t167;
t230 = t9 * t61;
t227 = pkin(2) * t131;
t114 = t175 - t227;
t60 = -pkin(5) * t146 - pkin(9) * t107 + t114;
t62 = t202 * t105 - t128 * t167;
t28 = t133 * t60 + t136 * t62;
t229 = qJD(6) * t28 - t219 * t133 + t241 * t136;
t27 = -t133 * t62 + t136 * t60;
t228 = -qJD(6) * t27 + t241 * t133 + t219 * t136;
t226 = pkin(4) * t134;
t225 = t78 * t90;
t184 = t136 * qJD(4);
t80 = t101 * t133 - t184;
t224 = t80 * t98;
t223 = t82 * t80;
t222 = t9 * t146;
t209 = t133 * t93;
t51 = qJD(6) * t82 + t209;
t205 = t51 * t136;
t221 = -t107 * t205 - t80 * t199;
t220 = t107 * t89 - t128 * t83 + t202 * t145;
t186 = qJD(6) * t136;
t218 = -t133 * t51 - t80 * t186;
t50 = -qJD(6) * t184 + t101 * t187 - t136 * t93;
t217 = t82 * t100 + t146 * t50;
t216 = -t103 * t98 - t107 * t92;
t215 = t101 * t80;
t214 = t101 * t98;
t213 = t146 * t92;
t212 = t128 * t53;
t211 = t133 * t80;
t210 = t133 * t92;
t208 = t134 * t73;
t207 = t136 * t82;
t206 = t50 * t133;
t204 = t82 * t101;
t200 = t103 * t133;
t140 = qJD(2) ^ 2;
t198 = t130 * t140;
t139 = qJD(4) ^ 2;
t196 = t139 * t134;
t195 = t139 * t137;
t126 = t134 ^ 2;
t127 = t137 ^ 2;
t193 = t126 - t127;
t190 = qJD(2) * t137;
t185 = t103 * qJD(4);
t181 = pkin(4) * t191;
t178 = t82 * t200;
t176 = t134 * t140 * t137;
t168 = t136 * t96;
t162 = t137 * t171;
t120 = pkin(4) * t128 + pkin(9);
t161 = qJD(6) * t120 * t96 + t9;
t160 = -t133 * t6 + t136 * t156;
t159 = -t133 * t156 - t136 * t6;
t157 = -t100 * t80 + t146 * t51;
t35 = t128 * t69 + t202 * t70;
t22 = t133 * t90 + t136 * t35;
t21 = -t133 * t35 + t136 * t90;
t55 = t113 - t208;
t56 = t137 * t73 + t197;
t155 = t134 * t55 - t137 * t56;
t154 = t100 * t101 - t146 * t93;
t151 = -t98 * t169 - t150;
t149 = -t121 * t139 + t237;
t123 = -pkin(3) - t227;
t72 = -qJD(2) * pkin(3) - t75;
t148 = qJD(4) * (qJD(2) * t123 + t72 + t89);
t17 = t202 * t47 - t212;
t15 = -qJD(4) * pkin(5) - t17;
t147 = -t120 * t92 + t96 * t15;
t143 = t160 * qJD(6) + t1 * t136 - t2 * t133;
t142 = (-t96 * t186 - t210) * t107 - t96 * t200;
t31 = -t56 * qJD(4) + t134 * t79;
t141 = -t31 * t134 + t30 * t137 + (-t134 * t56 - t137 * t55) * qJD(4);
t122 = -t202 * pkin(4) - pkin(5);
t97 = t98 ^ 2;
t95 = t100 * qJD(4);
t87 = qJD(2) * t91;
t57 = pkin(5) * t101 + pkin(9) * t98 + t181;
t40 = qJD(4) * t69 - t88 * t137;
t39 = -qJD(4) * t70 + t88 * t134;
t20 = t202 * t52 - t212;
t19 = t128 * t52 + t44;
t14 = t128 * t39 + t202 * t40;
t13 = t128 * t40 - t202 * t39;
t12 = t133 * t57 + t136 * t20;
t11 = -t133 * t20 + t136 * t57;
t4 = -qJD(6) * t22 - t133 * t14 + t87 * t136;
t3 = qJD(6) * t21 + t87 * t133 + t136 * t14;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135 * t198, -t138 * t198, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * qJD(2), t88 * qJD(2), 0, -t75 * t87 - t76 * t88 - t79 * t91 + t225, 0, 0, 0, 0, 0, 0, t39 * qJD(4) + (-t137 * t87 + t189 * t90) * qJD(2), -t40 * qJD(4) + (t134 * t87 + t188 * t90) * qJD(2) (-t134 * t39 + t137 * t40 + (-t134 * t70 - t137 * t69) * qJD(4)) * qJD(2), t30 * t70 + t31 * t69 + t39 * t55 + t40 * t56 + t72 * t87 + t225, 0, 0, 0, 0, 0, 0, -qJD(4) * t13 + t87 * t98 + t90 * t92, -qJD(4) * t14 + t101 * t87 + t90 * t93, t101 * t13 - t14 * t98 + t34 * t93 - t35 * t92, t10 * t35 - t13 * t17 + t14 * t18 + t64 * t87 + t68 * t90 + t233, 0, 0, 0, 0, 0, 0, t13 * t80 + t21 * t92 + t34 * t51 + t4 * t96, t13 * t82 - t22 * t92 - t3 * t96 - t34 * t50, t21 * t50 - t22 * t51 - t3 * t80 - t4 * t82, t1 * t22 + t13 * t15 - t156 * t4 + t2 * t21 + t3 * t6 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237 (qJD(1) * t90 + t89) * qJD(2), 0, t75 * t86 - t76 * t89 + (-t129 * t79 - t131 * t78) * pkin(2), 0.2e1 * t162, -0.2e1 * t193 * t182, t195, -0.2e1 * t162, -t196, 0, t134 * t148 + t137 * t149, -t134 * t149 + t137 * t148 (-t126 - t127) * t89 * qJD(2) + t141, t121 * t141 + t78 * t123 + t155 * t89 - t72 * t86, t101 * t103 + t107 * t93, -t154 + t216, t185, t100 * t98 - t213, -t95, 0, t64 * t100 - t68 * t146 + t114 * t92 - t86 * t98 + (t98 * t226 + t220) * qJD(4), -t86 * t101 + t64 * t103 + t68 * t107 + t114 * t93 + (t101 * t226 + t219) * qJD(4), t10 * t146 - t18 * t100 - t220 * t101 - t17 * t103 + t9 * t107 + t219 * t98 + t61 * t93 - t62 * t92, t10 * t62 + t68 * t114 + t220 * t17 - t219 * t18 + t242 * t64 + t230, t82 * t199 + (-t136 * t50 - t187 * t82) * t107, -t178 + (t206 + (-t207 + t211) * qJD(6)) * t107 + t221, t217 + t235, -t218 * t107 + t80 * t200, t142 + t157, t100 * t96 - t213, t15 * t200 - t156 * t100 - t2 * t146 + t27 * t92 + t61 * t51 - t229 * t96 - t220 * t80 + (t9 * t133 + t15 * t186) * t107, t15 * t199 + t1 * t146 - t6 * t100 - t28 * t92 - t61 * t50 + t228 * t96 - t220 * t82 + (t9 * t136 - t15 * t187) * t107, t27 * t50 - t28 * t51 + t229 * t82 + t228 * t80 + t160 * t103 + (qJD(6) * t159 - t1 * t133 - t2 * t136) * t107, t1 * t28 - t220 * t15 + t156 * t229 + t2 * t27 - t228 * t6 + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, -t195, 0, -qJD(4) * t155 + t30 * t134 + t31 * t137, 0, 0, 0, 0, 0, 0, -t95, -t185, t154 + t216, t10 * t107 - t100 * t17 + t103 * t18 - t222, 0, 0, 0, 0, 0, 0, t142 - t157, t217 - t235, t178 + (-t206 + (t207 + t211) * qJD(6)) * t107 + t221, t15 * t100 - t159 * t103 + t107 * t143 - t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t193 * t140, 0, t176, 0, 0 (-qJD(2) * t72 + t79) * t134, -t72 * t190 + (t55 + t208) * qJD(4) - t203, 0, 0, t214, -t97 + t234, -t115 + (t116 + t98) * qJD(4), -t214, 0, 0, qJD(4) * t19 - t101 * t64 - t181 * t98 - t9, -t172 + t64 * t98 + (-qJD(2) * pkin(4) * t101 - t128 * t164) * t134 + (-t128 * (-qJ(5) * t190 - t56) + t20) * qJD(4) (-t17 + t20) * t98 + (t18 - t19) * t101 + (-t128 * t92 - t202 * t93) * pkin(4), t17 * t19 - t18 * t20 + (t10 * t128 - t64 * t191 - t202 * t9) * pkin(4), t168 * t82 - t206 (-t50 - t224) * t136 - t238 + t218, t168 * t96 - t204 + t210, t169 * t80 - t205, t151 + t215, -t96 * t101, t101 * t156 - t11 * t96 + t122 * t51 + t133 * t147 - t136 * t161 - t19 * t80, t6 * t101 + t12 * t96 - t122 * t50 + t133 * t161 + t136 * t147 - t19 * t82, t11 * t82 + t12 * t80 + (-t120 * t51 + t156 * t98 + t1 + (t120 * t82 + t156) * qJD(6)) * t136 + (-t120 * t50 - t6 * t98 - t2 + (t120 * t80 - t6) * qJD(6)) * t133, t11 * t156 - t6 * t12 + t120 * t143 + t9 * t122 - t15 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101 * qJD(4), -t115 + (t116 - t98) * qJD(4), -t97 - t234, t101 * t17 + t18 * t98 + t68, 0, 0, 0, 0, 0, 0, t151 - t215, -t136 * t96 ^ 2 - t204 - t210 (t50 - t224) * t136 + t238 + t218, -t15 * t101 + t240 * t133 + t239 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -t80 ^ 2 + t82 ^ 2, t80 * t96 - t50, -t223, -t209 + (-qJD(6) + t96) * t82, t92, -t15 * t82 + t239, t15 * t80 - t240, 0, 0;];
tauc_reg  = t5;

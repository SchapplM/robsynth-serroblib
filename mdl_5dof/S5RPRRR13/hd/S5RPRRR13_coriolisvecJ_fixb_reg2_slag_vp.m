% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR13_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:27
% EndTime: 2019-12-31 19:15:36
% DurationCPUTime: 2.91s
% Computational Cost: add. (3927->352), mult. (8443->514), div. (0->0), fcn. (5229->6), ass. (0->177)
t131 = cos(qJ(4));
t132 = cos(qJ(3));
t228 = -pkin(8) - pkin(7);
t170 = qJD(4) * t228;
t129 = sin(qJ(3));
t201 = t131 * t129;
t174 = pkin(8) * t201;
t128 = sin(qJ(4));
t133 = -pkin(1) - pkin(6);
t114 = t133 * qJD(1) + qJD(2);
t199 = t132 * t114;
t150 = pkin(3) * t132 + pkin(7) * t129;
t96 = t150 * qJD(1);
t47 = -t128 * t199 + t131 * t96;
t237 = (pkin(4) * t132 + t174) * qJD(1) + t47 - t131 * t170;
t180 = t129 * qJD(1);
t169 = t128 * t180;
t48 = t128 * t96 + t131 * t199;
t236 = pkin(8) * t169 - t128 * t170 + t48;
t127 = sin(qJ(5));
t130 = cos(qJ(5));
t179 = t131 * qJD(3);
t190 = qJD(1) * t132;
t91 = t128 * t190 - t179;
t162 = t131 * t190;
t181 = t128 * qJD(3);
t93 = t162 + t181;
t148 = t127 * t91 - t130 * t93;
t41 = t127 * t93 + t130 * t91;
t227 = t41 * t148;
t235 = t148 ^ 2 - t41 ^ 2;
t149 = t129 * pkin(3) - t132 * pkin(7);
t100 = qJ(2) + t149;
t80 = t100 * qJD(1);
t99 = t129 * t114;
t83 = qJD(3) * pkin(7) + t99;
t38 = t128 * t80 + t131 * t83;
t89 = qJD(3) * t150 + qJD(2);
t67 = t89 * qJD(1);
t140 = -qJD(4) * t38 + t131 * t67;
t188 = qJD(3) * t132;
t184 = qJD(4) * t132;
t164 = t128 * t184;
t166 = t129 * t179;
t142 = t164 + t166;
t177 = qJD(3) * qJD(4);
t56 = qJD(1) * t142 - t131 * t177;
t10 = t56 * pkin(8) + (pkin(4) * qJD(1) - t114 * t128) * t188 + t140;
t167 = t129 * t181;
t193 = qJD(4) * t162 + t128 * t177;
t143 = qJD(1) * t167 - t193;
t168 = t114 * t188;
t185 = qJD(4) * t131;
t186 = qJD(4) * t128;
t16 = t128 * t67 + t131 * t168 + t80 * t185 - t186 * t83;
t11 = pkin(8) * t143 + t16;
t33 = -t91 * pkin(8) + t38;
t215 = t130 * t33;
t117 = qJD(4) + t180;
t37 = -t128 * t83 + t131 * t80;
t32 = -t93 * pkin(8) + t37;
t26 = t117 * pkin(4) + t32;
t6 = t127 * t26 + t215;
t2 = -qJD(5) * t6 + t130 * t10 - t127 * t11;
t219 = qJD(3) * pkin(3);
t84 = -t199 - t219;
t50 = t91 * pkin(4) + t84;
t234 = t50 * t148 + t2;
t113 = qJD(5) + t117;
t182 = qJD(5) * t130;
t183 = qJD(5) * t127;
t12 = -t127 * t143 + t130 * t56 + t91 * t182 + t183 * t93;
t233 = t41 * t113 - t12;
t1 = (qJD(5) * t26 + t11) * t130 + t127 * t10 - t33 * t183;
t232 = t50 * t41 - t1;
t137 = qJD(5) * t148 + t127 * t56 + t130 * t143;
t231 = -t113 * t148 + t137;
t175 = 0.2e1 * qJD(1);
t203 = t130 * t131;
t206 = t127 * t128;
t94 = -t203 + t206;
t70 = t94 * t129;
t95 = t127 * t131 + t130 * t128;
t68 = t95 * t129;
t230 = -t37 * t117 + t16;
t17 = -t128 * t168 + t140;
t229 = -t38 * t117 - t17;
t126 = t132 ^ 2;
t192 = t129 ^ 2 - t126;
t176 = qJD(4) + qJD(5);
t226 = t93 * t91;
t108 = t228 * t128;
t109 = t228 * t131;
t55 = t127 * t108 - t130 * t109;
t225 = qJD(5) * t55 - t236 * t127 + t237 * t130;
t54 = t130 * t108 + t127 * t109;
t224 = -qJD(5) * t54 + t237 * t127 + t236 * t130;
t71 = t94 * t132;
t78 = t95 * qJD(1);
t223 = -qJD(3) * t71 - t176 * t68 - t78;
t222 = t94 * qJD(1) + t176 * t70 - t95 * t188;
t221 = t127 * t169 - t130 * t185 - t131 * t182 + t176 * t206 - t180 * t203;
t46 = t176 * t95;
t220 = t129 * t78 + t46;
t218 = t127 * t33;
t217 = t128 * t84;
t216 = t129 * t84;
t214 = t131 * t37;
t213 = t132 * t56;
t210 = t56 * t128;
t209 = t91 * t117;
t112 = t133 * t201;
t59 = t128 * t100 + t112;
t208 = t117 * t128;
t207 = t117 * t129;
t205 = t128 * t132;
t204 = t129 * t133;
t202 = t131 * t117;
t200 = t131 * t132;
t134 = qJD(3) ^ 2;
t198 = t134 * t129;
t197 = t134 * t132;
t135 = qJD(1) ^ 2;
t196 = t135 * qJ(2);
t195 = t91 * qJD(4);
t194 = t93 * qJD(4);
t191 = -t134 - t135;
t189 = qJD(3) * t129;
t187 = qJD(3) * t133;
t178 = qJD(1) * qJD(3);
t173 = qJD(2) * t175;
t172 = t128 * t204;
t171 = t132 * t135 * t129;
t165 = t132 * t187;
t163 = t131 * t184;
t118 = t132 * t178;
t161 = -t128 * t133 + pkin(4);
t158 = qJD(1) * t59 + t38;
t157 = -t56 + t195;
t155 = -t84 + t199;
t154 = -t93 + t181;
t153 = t91 + t179;
t152 = qJD(4) * t129 + qJD(1);
t151 = -t99 + (t169 + t186) * pkin(4);
t88 = t131 * t100;
t39 = -pkin(8) * t200 + t129 * t161 + t88;
t49 = -pkin(8) * t205 + t59;
t19 = -t127 * t49 + t130 * t39;
t20 = t127 * t39 + t130 * t49;
t147 = -t128 * t38 - t214;
t146 = t128 * t37 - t131 * t38;
t145 = t129 * t154;
t144 = t155 * qJD(3);
t141 = -t163 + t167;
t30 = -qJD(4) * t172 + t100 * t185 + t128 * t89 + t131 * t165;
t138 = t143 + t194;
t136 = qJD(4) * t147 - t17 * t128 + t16 * t131;
t123 = qJ(2) * t173;
t122 = -t131 * pkin(4) - pkin(3);
t110 = t129 * t118;
t90 = (pkin(4) * t128 - t133) * t132;
t74 = t131 * t89;
t69 = t95 * t132;
t58 = t88 - t172;
t57 = -pkin(4) * t141 + t129 * t187;
t35 = -pkin(4) * t143 + t114 * t189;
t31 = -t59 * qJD(4) - t128 * t165 + t74;
t25 = -t183 * t205 + (t176 * t200 - t167) * t130 - t142 * t127;
t23 = -t127 * t167 + t130 * t166 + t132 * t46;
t21 = pkin(8) * t141 + t30;
t18 = t74 + (-t112 + (pkin(8) * t132 - t100) * t128) * qJD(4) + (t132 * t161 + t174) * qJD(3);
t9 = t130 * t32 - t218;
t8 = -t127 * t32 - t215;
t5 = t130 * t26 - t218;
t4 = -qJD(5) * t20 - t127 * t21 + t130 * t18;
t3 = qJD(5) * t19 + t127 * t18 + t130 * t21;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t123, -0.2e1 * t110, 0.2e1 * t192 * t178, -t198, 0.2e1 * t110, -t197, 0, -t133 * t198 + (qJ(2) * t188 + qJD(2) * t129) * t175, -t133 * t197 + (-qJ(2) * t189 + qJD(2) * t132) * t175, 0, t123, -t93 * t164 + (-t189 * t93 - t213) * t131, (t93 * t128 + t131 * t91) * t189 + (t143 * t131 + t210 + (t91 * t128 - t93 * t131) * qJD(4)) * t132, -t117 * t164 - t56 * t129 + (t132 * t93 + (qJD(1) * t126 - t207) * t131) * qJD(3), -t91 * t167 + (-t128 * t143 + t185 * t91) * t132, -t117 * t163 - t193 * t129 + (-t91 * t132 + (t192 * qJD(1) + t207) * t128) * qJD(3), t117 * t188 + t110, t31 * t117 + (t17 + (t133 * t91 - t217) * qJD(3)) * t129 + (-t133 * t193 + t84 * t185 + (t128 * t99 + t37 + (t58 + t172) * qJD(1)) * qJD(3)) * t132, -t30 * t117 - t16 * t129 + (t133 * t56 - t186 * t84) * t132 + (-t158 * t132 + (t131 * t155 + t133 * t93) * t129) * qJD(3), -t30 * t91 - t59 * t193 - t31 * t93 + t58 * t56 + (t128 * t158 + t214) * t189 + (qJD(4) * t146 - t16 * t128 - t17 * t131) * t132, -t144 * t204 + t16 * t59 + t17 * t58 + t38 * t30 + t37 * t31, t12 * t71 + t148 * t23, t12 * t69 - t137 * t71 + t148 * t25 + t23 * t41, -t23 * t113 - t12 * t129 + (-qJD(1) * t71 - t148) * t188, -t137 * t69 + t41 * t25, -t25 * t113 + t137 * t129 + (-qJD(1) * t69 - t41) * t188, t113 * t188 + t110, t4 * t113 + t2 * t129 - t90 * t137 + t50 * t25 + t35 * t69 + t57 * t41 + (qJD(1) * t19 + t5) * t188, -t1 * t129 - t3 * t113 - t90 * t12 - t50 * t23 - t35 * t71 - t57 * t148 + (-qJD(1) * t20 - t6) * t188, -t1 * t69 + t19 * t12 + t137 * t20 + t148 * t4 + t2 * t71 + t5 * t23 - t6 * t25 - t3 * t41, t1 * t20 + t2 * t19 + t6 * t3 + t35 * t90 + t5 * t4 + t50 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t196, 0, 0, 0, 0, 0, 0, t191 * t129, t191 * t132, 0, -t196, 0, 0, 0, 0, 0, 0, -t132 * t193 - t152 * t202 + (-t117 * t205 + t129 * t91) * qJD(3), t213 + t152 * t208 + (-t117 * t200 + (t93 - t162) * t129) * qJD(3), (qJD(1) * t91 + t129 * t157 + t188 * t93) * t128 + (qJD(1) * t93 + t129 * t138 - t188 * t91) * t131, -t146 * t188 + t147 * qJD(1) + (-t144 + t136) * t129, 0, 0, 0, 0, 0, 0, t132 * t137 + t222 * t113 + (t129 * t41 - t190 * t68) * qJD(3), t132 * t12 - t223 * t113 + (-t129 * t148 + t190 * t70) * qJD(3), -t68 * t12 - t137 * t70 + t148 * t222 - t223 * t41, -t1 * t70 - t35 * t132 + t50 * t189 - t2 * t68 + t222 * t5 + t223 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t192 * t135, 0, -t171, 0, 0, -t132 * t196, t129 * t196, 0, 0, t93 * t202 - t210, (-t56 - t209) * t131 + (qJD(1) * t145 - t193 - t194) * t128, t117 * t185 + (t117 * t201 + t132 * t154) * qJD(1), -t193 * t131 + (t153 * t180 + t195) * t128, -t117 * t186 + (-t128 * t207 + t132 * t153) * qJD(1), -t117 * t190, -pkin(3) * t193 - t47 * t117 - t153 * t99 + (-pkin(7) * t202 + t217) * qJD(4) + (-t37 * t132 + (qJD(3) * t149 + t216) * t128) * qJD(1), pkin(3) * t56 + t48 * t117 + t114 * t145 + (pkin(7) * t208 + t131 * t84) * qJD(4) + (t132 * t38 + (-pkin(7) * t188 + t216) * t131) * qJD(1), t47 * t93 + t48 * t91 + (pkin(7) * t157 + t229) * t128 + (pkin(7) * t138 + t230) * t131, -t37 * t47 - t38 * t48 + (-t84 - t219) * t99 + t136 * pkin(7), -t12 * t95 + t148 * t221, t12 * t94 + t137 * t95 + t148 * t220 + t221 * t41, -t221 * t113 + (qJD(3) * t95 + t148) * t190, -t137 * t94 + t220 * t41, -t220 * t113 + (-qJD(3) * t94 + t41) * t190, -t113 * t190, -t122 * t137 + t35 * t94 + t220 * t50 + t151 * t41 - t225 * t113 + (qJD(3) * t54 - t5) * t190, -t122 * t12 + t35 * t95 - t221 * t50 - t151 * t148 + t224 * t113 + (-qJD(3) * t55 + t6) * t190, -t1 * t94 + t54 * t12 + t137 * t55 - t148 * t225 - t2 * t95 - t220 * t6 + t221 * t5 + t224 * t41, t1 * t55 + t35 * t122 + t151 * t50 + t2 * t54 - t224 * t6 - t225 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t91 ^ 2 + t93 ^ 2, t209 - t56, -t226, t93 * t117 + t143, t118, -t84 * t93 - t229, t84 * t91 - t230, 0, 0, -t227, t235, t233, t227, t231, t118, -t8 * t113 + (-t113 * t183 + t118 * t130 - t41 * t93) * pkin(4) + t234, t9 * t113 + (-t113 * t182 - t118 * t127 + t148 * t93) * pkin(4) + t232, -t6 * t148 + t9 * t41 - t5 * t41 - t8 * t148 + (t12 * t130 + t127 * t137 + (-t127 * t148 - t130 * t41) * qJD(5)) * pkin(4), -t5 * t8 - t6 * t9 + (t1 * t127 + t130 * t2 - t50 * t93 + (-t127 * t5 + t130 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, t235, t233, t227, t231, t118, t6 * t113 + t234, t5 * t113 + t232, 0, 0;];
tauc_reg = t7;

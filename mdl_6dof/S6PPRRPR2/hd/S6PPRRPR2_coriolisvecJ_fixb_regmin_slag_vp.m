% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:51:04
% EndTime: 2019-03-08 18:51:07
% DurationCPUTime: 1.69s
% Computational Cost: add. (1563->255), mult. (4282->378), div. (0->0), fcn. (3601->12), ass. (0->157)
t89 = cos(pkin(12));
t90 = cos(pkin(7));
t182 = t89 * t90;
t86 = sin(pkin(12));
t88 = sin(pkin(6));
t94 = sin(qJ(3));
t97 = cos(qJ(3));
t104 = (t94 * t182 + t86 * t97) * t88;
t87 = sin(pkin(7));
t184 = t87 * t94;
t91 = cos(pkin(6));
t78 = qJD(1) * t91 + qJD(2);
t36 = qJD(1) * t104 + t78 * t184;
t34 = qJD(3) * pkin(9) + t36;
t172 = qJD(1) * t88;
t146 = t89 * t172;
t51 = -t87 * t146 + t78 * t90;
t93 = sin(qJ(4));
t96 = cos(qJ(4));
t178 = -t93 * t34 + t96 * t51;
t204 = -qJD(5) + t178;
t14 = -qJD(4) * pkin(4) - t204;
t183 = t87 * t97;
t201 = (t97 * t182 - t86 * t94) * t88;
t205 = t91 * t183 + t201;
t203 = qJD(1) * t201 + t78 * t183;
t164 = qJD(4) * t93;
t202 = -pkin(4) * t164 + t36;
t157 = qJD(4) * qJ(5);
t17 = t96 * t34 + t93 * t51;
t15 = -t157 - t17;
t130 = t90 * t146;
t168 = qJD(3) * t94;
t145 = t87 * t168;
t147 = t86 * t172;
t166 = qJD(3) * t97;
t32 = t130 * t168 + t78 * t145 + t147 * t166;
t200 = qJD(3) * t36 - t32;
t169 = qJD(3) * t93;
t80 = qJD(6) + t169;
t199 = -qJD(6) + t80;
t149 = pkin(5) * t169;
t158 = t149 - t204;
t35 = (t78 * t87 + t130) * t97 - t94 * t147;
t41 = t91 * t184 + t104;
t39 = t41 * qJD(3);
t54 = -t87 * t88 * t89 + t90 * t91;
t23 = t41 * t96 + t54 * t93;
t38 = t205 * qJD(3);
t8 = t23 * qJD(4) + t38 * t93;
t198 = (-t164 * t205 - t39 * t96) * qJD(3) - t8 * qJD(4);
t162 = qJD(4) * t96;
t9 = -t41 * t164 + (qJD(4) * t54 + t38) * t96;
t197 = (-t162 * t205 + t39 * t93) * qJD(3) - t9 * qJD(4);
t159 = t93 * qJD(5);
t108 = -t96 * t157 - t159;
t179 = -t108 + t202;
t156 = qJD(3) * qJD(4);
t137 = t93 * t156;
t134 = pkin(4) * t137 + t32;
t19 = t108 * qJD(3) + t134;
t99 = qJD(4) ^ 2;
t192 = pkin(9) * t99;
t196 = t179 * qJD(3) - t19 - t192;
t144 = t87 * t166;
t128 = t93 * t144;
t100 = qJD(3) ^ 2;
t175 = t100 * t87;
t150 = t94 * t175;
t56 = t96 * t184 + t90 * t93;
t43 = t56 * qJD(4) + t128;
t195 = (t43 + t128) * qJD(4) + t96 * t150;
t129 = t96 * t144;
t155 = t93 * t184;
t42 = qJD(4) * t155 - t90 * t162 - t129;
t194 = (-t42 + t129) * qJD(4) - t93 * t150;
t98 = -pkin(4) - pkin(10);
t193 = pkin(5) + pkin(9);
t31 = t203 * qJD(3);
t151 = t51 * t162 - t34 * t164 + t96 * t31;
t3 = (qJD(5) - t149) * qJD(4) + t151;
t92 = sin(qJ(6));
t191 = t3 * t92;
t95 = cos(qJ(6));
t190 = t3 * t95;
t165 = qJD(4) * t92;
t167 = qJD(3) * t96;
t60 = t95 * t167 + t165;
t48 = -t60 * qJD(6) + t92 * t137;
t189 = t48 * t95;
t188 = t60 * t80;
t140 = t92 * t167;
t163 = qJD(4) * t95;
t62 = -t140 + t163;
t187 = t62 * t80;
t186 = t80 * t93;
t185 = t80 * t98;
t181 = t96 * t62;
t127 = pkin(10) * t93 - qJ(5) * t96;
t102 = t127 * qJD(4) - t159;
t180 = -t102 + t202;
t84 = t93 ^ 2;
t85 = t96 ^ 2;
t177 = t84 - t85;
t176 = qJD(3) * pkin(3);
t135 = -t93 * qJ(5) - pkin(3);
t68 = -pkin(4) * t96 + t135;
t170 = qJD(3) * t68;
t161 = qJD(6) * t92;
t160 = qJD(6) * t95;
t152 = t95 * t186;
t7 = t34 * t162 + t51 * t164 + t93 * t31;
t148 = t93 * t100 * t96;
t71 = t193 * t96;
t143 = t80 * t161;
t142 = t80 * t160;
t141 = t96 * t160;
t18 = t102 * qJD(3) + t134;
t136 = t96 * t156;
t5 = pkin(5) * t136 + t7;
t138 = -t92 * t18 + t95 * t5;
t10 = t98 * qJD(4) + t158;
t59 = t98 * t96 + t135;
t20 = t59 * qJD(3) - t35;
t1 = t10 * t95 - t20 * t92;
t2 = t10 * t92 + t20 * t95;
t22 = t41 * t93 - t54 * t96;
t126 = t205 * t92 + t22 * t95;
t125 = -t205 * t95 + t22 * t92;
t70 = t193 * t93;
t124 = t59 * t95 + t70 * t92;
t123 = -qJD(3) * t85 + t186;
t121 = t80 * t92;
t118 = qJD(4) * t17 - t7;
t55 = -t90 * t96 + t155;
t117 = t95 * t183 - t55 * t92;
t116 = t92 * t183 + t55 * t95;
t49 = qJD(4) * t160 - qJD(6) * t140 - t95 * t137;
t114 = t192 - t200;
t81 = pkin(5) * t167;
t11 = -t15 + t81;
t113 = t11 * t93 + t98 * t162;
t33 = -t35 - t176;
t110 = qJD(4) * (t33 + t35 - t176);
t24 = -t35 + t170;
t109 = qJD(4) * (-t24 - t35 - t170);
t6 = -qJD(4) * qJD(5) - t151;
t101 = -t6 * t96 + t7 * t93 + (t14 * t96 + t15 * t93) * qJD(4);
t83 = pkin(4) * t169;
t74 = t95 * t136;
t66 = qJD(4) * t71;
t65 = t193 * t164;
t63 = -qJ(5) * t167 + t83;
t52 = t127 * qJD(3) + t83;
t21 = t24 * t169;
t13 = t81 + t17;
t4 = [0, 0, 0, -t39 * qJD(3), -t38 * qJD(3), 0, 0, 0, 0, 0, t198, t197 (t8 * t93 + t9 * t96 + (t22 * t96 - t23 * t93) * qJD(4)) * qJD(3), -t198, -t197, t14 * t8 - t15 * t9 - t19 * t205 + t22 * t7 - t23 * t6 + t24 * t39, 0, 0, 0, 0, 0 (-t125 * qJD(6) - t39 * t92 + t8 * t95) * t80 + t126 * t136 + t9 * t60 + t23 * t49 -(t126 * qJD(6) + t39 * t95 + t8 * t92) * t80 - t125 * t136 + t9 * t62 + t23 * t48; 0, 0, 0, -t150, -t97 * t175, 0, 0, 0, 0, 0, -t195, -t194 (-t42 * t96 + t43 * t93 + (t55 * t96 - t56 * t93) * qJD(4)) * qJD(3), t195, t194, t14 * t43 + t15 * t42 + t7 * t55 - t6 * t56 + (t24 * t168 - t19 * t97) * t87, 0, 0, 0, 0, 0 (t117 * qJD(6) - t92 * t145 + t95 * t43) * t80 + t116 * t136 - t42 * t60 + t56 * t49 -(t116 * qJD(6) + t95 * t145 + t92 * t43) * t80 + t117 * t136 - t42 * t62 + t56 * t48; 0, 0, 0, t200 (-t203 + t35) * qJD(3), 0.2e1 * t93 * t136, -0.2e1 * t177 * t156, t99 * t96, -t99 * t93, 0, t93 * t110 - t114 * t96, t96 * t110 + t114 * t93 (-t84 - t85) * t35 * qJD(3) + t101, t93 * t109 - t196 * t96, t96 * t109 + t196 * t93, t19 * t68 + (-t14 * t93 + t15 * t96) * t35 - t179 * t24 + t101 * pkin(9), -t48 * t92 * t96 + (t92 * t164 - t141) * t62 (-t60 * t92 + t62 * t95) * t164 + (-t189 + t49 * t92 + (t60 * t95 + t62 * t92) * qJD(6)) * t96, -t80 * t141 + t48 * t93 + (t123 * t92 + t181) * qJD(4), t96 * t143 - t49 * t93 + (t123 * t95 - t96 * t60) * qJD(4) (t80 + t169) * t162, t71 * t49 - t65 * t60 + (-t11 * t163 + t138) * t93 + ((-t35 * t93 + t66) * t95 + t180 * t92) * t80 + (-t124 * t80 - t2 * t93) * qJD(6) + (-t11 * t161 + t190 - t35 * t60 + ((-t59 * t92 + t70 * t95) * qJD(3) + t1) * qJD(4)) * t96, t71 * t48 - t65 * t62 + (-(qJD(6) * t10 + t18) * t93 + (-qJD(6) * t70 + t180) * t80) * t95 + (-(-qJD(6) * t59 + t66) * t80 + (qJD(4) * t11 + qJD(6) * t20 + t35 * t80 - t5) * t93) * t92 + (-t11 * t160 - t191 - t35 * t62 + (-qJD(3) * t124 - t2) * qJD(4)) * t96; 0, 0, 0, 0, 0, -t148, t177 * t100, 0, 0, 0, -t33 * t169 + t118, qJD(4) * t178 - t33 * t167 - t151, 0, -t63 * t167 - t118 + t21 (0.2e1 * qJD(5) - t178) * qJD(4) + (t24 * t96 + t63 * t93) * qJD(3) + t151, -t7 * pkin(4) - t6 * qJ(5) - t14 * t17 + t204 * t15 - t24 * t63, -t62 * t121 + t189 (-t49 - t187) * t95 + (-t48 + t188) * t92, -t143 + t74 + (-t92 * t186 - t181) * qJD(3), -t142 + (-t152 + (t60 - t165) * t96) * qJD(3), -t80 * t167, qJ(5) * t49 + t191 - (t13 * t95 - t52 * t92) * t80 + t158 * t60 + (t11 * t95 - t92 * t185) * qJD(6) + (-t1 * t96 + t113 * t95) * qJD(3), qJ(5) * t48 + t190 + (t13 * t92 + t52 * t95) * t80 + t158 * t62 + (-t11 * t92 - t95 * t185) * qJD(6) + (-t113 * t92 + t2 * t96) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t100 * t84 - t99, qJD(4) * t15 + t21 + t7, 0, 0, 0, 0, 0, -qJD(4) * t60 - t80 * t121 + t74, -t142 - qJD(4) * t62 + (-t92 * t162 - t152) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, t48 + t188, -t49 + t187, t136, -t11 * t62 + t199 * t2 + t138, t1 * t199 + t11 * t60 - t95 * t18 - t92 * t5;];
tauc_reg  = t4;

% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:51
% EndTime: 2019-12-05 18:38:58
% DurationCPUTime: 2.15s
% Computational Cost: add. (3272->245), mult. (8891->342), div. (0->0), fcn. (6610->8), ass. (0->163)
t152 = sin(qJ(5));
t155 = cos(qJ(5));
t190 = qJD(5) * t152;
t156 = cos(qJ(3));
t157 = cos(qJ(2));
t192 = qJD(1) * t157;
t185 = t156 * t192;
t153 = sin(qJ(3));
t154 = sin(qJ(2));
t193 = qJD(1) * t154;
t186 = t153 * t193;
t111 = -t185 + t186;
t113 = -t153 * t192 - t156 * t193;
t150 = sin(pkin(9));
t151 = cos(pkin(9));
t169 = -t151 * t111 + t150 * t113;
t211 = t155 * t169;
t147 = qJD(2) + qJD(3);
t189 = qJD(1) * qJD(2);
t182 = t157 * t189;
t95 = qJD(3) * t185 - t147 * t186 + t156 * t182;
t125 = t153 * t157 + t156 * t154;
t100 = t147 * t125;
t96 = t100 * qJD(1);
t54 = -t150 * t95 - t151 * t96;
t55 = -t150 * t96 + t151 * t95;
t91 = t150 * t111 + t151 * t113;
t11 = qJD(5) * t211 + t152 * t54 + t155 * t55 + t190 * t91;
t146 = qJD(5) + t147;
t49 = t152 * t91 + t211;
t209 = t49 * t146;
t4 = t11 - t209;
t225 = t152 * t169 - t155 * t91;
t215 = t225 * t49;
t6 = t225 ^ 2 - t49 ^ 2;
t218 = t169 * pkin(8);
t219 = pkin(6) + pkin(7);
t133 = t219 * t157;
t128 = qJD(1) * t133;
t118 = t156 * t128;
t132 = t219 * t154;
t126 = qJD(1) * t132;
t213 = qJD(2) * pkin(2);
t120 = -t126 + t213;
t168 = -t153 * t120 - t118;
t204 = t111 * qJ(4);
t78 = -t168 - t204;
t212 = t151 * t78;
t109 = t113 * qJ(4);
t114 = t153 * t128;
t177 = t156 * t120 - t114;
t77 = t109 + t177;
t66 = t147 * pkin(3) + t77;
t33 = t150 * t66 + t212;
t19 = t33 + t218;
t143 = -t157 * pkin(2) - pkin(1);
t131 = t143 * qJD(1);
t101 = t111 * pkin(3) + qJD(4) + t131;
t61 = -pkin(4) * t169 + t101;
t183 = t19 * t190 - t61 * t49;
t12 = qJD(5) * t225 + t152 * t55 - t155 * t54;
t210 = t225 * t146;
t5 = -t12 + t210;
t187 = qJD(2) * t219;
t172 = qJD(1) * t187;
t122 = t157 * t172;
t191 = qJD(3) * t153;
t175 = -t153 * t122 - t128 * t191;
t121 = t154 * t172;
t222 = (qJD(3) * t120 - t121) * t156;
t27 = -t96 * qJ(4) - t111 * qJD(4) + t175 + t222;
t176 = t153 * t121 - t156 * t122;
t161 = qJD(3) * t168 + t176;
t28 = -t95 * qJ(4) + t113 * qJD(4) + t161;
t9 = -t150 * t27 + t151 * t28;
t2 = -t55 * pkin(8) + t9;
t10 = t150 * t28 + t151 * t27;
t3 = t54 * pkin(8) + t10;
t165 = -t152 * t3 + t155 * t2 - t225 * t61;
t88 = t91 * pkin(8);
t224 = -0.2e1 * t189;
t200 = t151 * t153;
t214 = pkin(2) * qJD(3);
t174 = t153 * t126 - t118;
t79 = t174 + t204;
t195 = -t156 * t126 - t114;
t80 = t109 + t195;
t207 = t150 * t80 - t151 * t79 + (-t150 * t156 - t200) * t214;
t201 = t150 * t153;
t205 = -t150 * t79 - t151 * t80 + (t151 * t156 - t201) * t214;
t221 = qJD(1) * t125;
t220 = qJD(5) - t146;
t217 = pkin(3) * t150;
t216 = t113 * pkin(3);
t124 = t153 * t154 - t156 * t157;
t127 = t154 * t187;
t129 = t157 * t187;
t199 = t156 * t132;
t163 = -qJD(3) * t199 - t156 * t127 - t153 * t129 - t133 * t191;
t40 = -t100 * qJ(4) - t124 * qJD(4) + t163;
t167 = t153 * t132 - t156 * t133;
t160 = qJD(3) * t167 + t153 * t127 - t156 * t129;
t99 = t147 * t124;
t41 = t99 * qJ(4) - t125 * qJD(4) + t160;
t16 = t150 * t41 + t151 * t40;
t69 = t150 * t78;
t35 = t151 * t77 - t69;
t89 = -t125 * qJ(4) - t153 * t133 - t199;
t90 = -t124 * qJ(4) - t167;
t45 = t150 * t89 + t151 * t90;
t208 = t218 + t207;
t206 = t88 - t205;
t203 = t113 * t111;
t202 = t131 * t113;
t159 = qJD(1) ^ 2;
t198 = t157 * t159;
t158 = qJD(2) ^ 2;
t197 = t158 * t154;
t196 = t158 * t157;
t194 = t154 ^ 2 - t157 ^ 2;
t145 = t154 * t213;
t144 = pkin(2) * t193;
t184 = t96 * pkin(3) + qJD(2) * t144;
t181 = -pkin(2) * t147 - t120;
t180 = t100 * pkin(3) + t145;
t15 = -t150 * t40 + t151 * t41;
t32 = t151 * t66 - t69;
t34 = -t150 * t77 - t212;
t44 = -t150 * t90 + t151 * t89;
t178 = pkin(1) * t224;
t142 = t156 * pkin(2) + pkin(3);
t107 = -pkin(2) * t201 + t151 * t142;
t65 = -t91 * pkin(4) - t216;
t171 = t169 * t32 - t33 * t91;
t17 = t147 * pkin(4) + t32 + t88;
t170 = -t152 * t17 - t155 * t19;
t97 = -t151 * t124 - t150 * t125;
t98 = -t150 * t124 + t151 * t125;
t56 = t152 * t98 - t155 * t97;
t57 = t152 * t97 + t155 * t98;
t166 = t124 * pkin(3) + t143;
t164 = t131 * t111 - t175;
t140 = t151 * pkin(3) + pkin(4);
t108 = pkin(2) * t200 + t150 * t142;
t102 = pkin(4) + t107;
t81 = -t111 ^ 2 + t113 ^ 2;
t73 = -t97 * pkin(4) + t166;
t68 = (-t113 - t221) * t147;
t67 = t111 * t147 + t95;
t62 = t144 + t65;
t60 = -t150 * t100 - t151 * t99;
t59 = -t151 * t100 + t150 * t99;
t39 = -t59 * pkin(4) + t180;
t31 = t97 * pkin(8) + t45;
t30 = -t98 * pkin(8) + t44;
t29 = -t54 * pkin(4) + t184;
t21 = t88 + t35;
t20 = t34 - t218;
t14 = qJD(5) * t57 + t152 * t60 - t155 * t59;
t13 = -qJD(5) * t56 + t152 * t59 + t155 * t60;
t8 = t59 * pkin(8) + t16;
t7 = -t60 * pkin(8) + t15;
t1 = [0, 0, 0, 0.2e1 * t154 * t182, t194 * t224, t196, -t197, 0, -pkin(6) * t196 + t154 * t178, pkin(6) * t197 + t157 * t178, t113 * t99 + t95 * t125, t113 * t100 + t99 * t111 - t95 * t124 - t125 * t96, -t99 * t147, -t100 * t147, 0, t143 * t96 + t131 * t100 + t160 * t147 + (qJD(1) * t124 + t111) * t145, t143 * t95 - t131 * t99 - t163 * t147 + (-t113 + t221) * t145, t10 * t97 + t15 * t91 + t16 * t169 - t32 * t60 + t33 * t59 - t44 * t55 + t45 * t54 - t9 * t98, t10 * t45 + t101 * t180 + t32 * t15 + t33 * t16 + t166 * t184 + t9 * t44, t11 * t57 + t13 * t225, -t11 * t56 - t57 * t12 + t13 * t49 - t14 * t225, t13 * t146, -t14 * t146, 0, -t39 * t49 + t73 * t12 + t29 * t56 + t61 * t14 + (-t152 * t8 + t155 * t7 + (-t152 * t30 - t155 * t31) * qJD(5)) * t146, t39 * t225 + t73 * t11 + t29 * t57 + t61 * t13 - (t152 * t7 + t155 * t8 + (-t152 * t31 + t155 * t30) * qJD(5)) * t146; 0, 0, 0, -t154 * t198, t194 * t159, 0, 0, 0, t159 * pkin(1) * t154, pkin(1) * t198, -t203, t81, t67, t68, 0, -t111 * t144 + t202 - t174 * t147 + (t153 * t181 - t118) * qJD(3) + t176, t113 * t144 + t195 * t147 + (qJD(3) * t181 + t121) * t156 + t164, -t107 * t55 + t108 * t54 + t169 * t205 + t207 * t91 + t171, t10 * t108 + t9 * t107 - t101 * (t144 - t216) + t205 * t33 + t207 * t32, -t215, t6, t4, t5, 0, t62 * t49 + (t152 * t206 + t155 * t208) * t146 + ((-t102 * t152 - t108 * t155) * t146 + t170) * qJD(5) + t165, -t62 * t225 + (-t2 + (qJD(5) * t108 - t208) * t146) * t152 + (-qJD(5) * t17 - t3 + (-qJD(5) * t102 + t206) * t146) * t155 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203, t81, t67, t68, 0, -t147 * t168 + t161 + t202, t147 * t177 + t164 - t222, -t34 * t91 - t35 * t169 + (t150 * t54 - t151 * t55) * pkin(3) + t171, -t32 * t34 - t33 * t35 + (t10 * t150 + t101 * t113 + t151 * t9) * pkin(3), -t215, t6, t4, t5, 0, t65 * t49 - (-t152 * t21 + t155 * t20) * t146 + ((-t140 * t152 - t155 * t217) * t146 + t170) * qJD(5) + t165, -t155 * t3 - t152 * t2 - t65 * t225 + (t152 * t20 + t155 * t21) * t146 + (-(t140 * t155 - t152 * t217) * t146 - t155 * t17) * qJD(5) + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169 ^ 2 - t91 ^ 2, -t169 * t33 - t32 * t91 + t184, 0, 0, 0, 0, 0, t12 + t210, t11 + t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t6, t4, t5, 0, t220 * t170 + t165, (-t19 * t146 - t2) * t152 + (-t220 * t17 - t3) * t155 + t183;];
tauc_reg = t1;

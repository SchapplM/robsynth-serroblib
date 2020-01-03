% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:58
% EndTime: 2019-12-31 19:30:03
% DurationCPUTime: 1.57s
% Computational Cost: add. (2523->249), mult. (6512->322), div. (0->0), fcn. (4525->6), ass. (0->156)
t131 = sin(qJ(5));
t133 = cos(qJ(5));
t130 = sin(pkin(8));
t132 = sin(qJ(2));
t134 = cos(qJ(2));
t180 = cos(pkin(8));
t107 = t130 * t134 + t180 * t132;
t142 = qJD(1) * t107;
t158 = t180 * t134;
t152 = qJD(1) * t158;
t172 = qJD(1) * t132;
t94 = t130 * t172 - t152;
t51 = t131 * t94 + t133 * t142;
t96 = t107 * qJD(2);
t84 = qJD(1) * t96;
t168 = qJD(1) * qJD(2);
t163 = t132 * t168;
t116 = t130 * t163;
t85 = qJD(2) * t152 - t116;
t13 = qJD(5) * t51 + t131 * t85 - t133 * t84;
t125 = qJD(2) - qJD(5);
t214 = t51 * t125;
t217 = t13 + t214;
t186 = -qJ(3) - pkin(6);
t114 = t186 * t132;
t111 = qJD(1) * t114;
t115 = t186 * t134;
t112 = qJD(1) * t115;
t178 = t130 * t112;
t61 = t180 * t111 + t178;
t174 = qJD(4) - t61;
t150 = t131 * t142 - t133 * t94;
t189 = t150 ^ 2;
t190 = t51 ^ 2;
t216 = t189 - t190;
t201 = pkin(3) + pkin(4);
t164 = t134 * pkin(2) + pkin(1);
t210 = t164 * qJD(1);
t113 = qJD(3) - t210;
t212 = -t142 * qJ(4) + t113;
t25 = -t201 * t94 - t212;
t215 = t25 * t51;
t188 = t51 * t150;
t198 = t142 * pkin(7);
t213 = -t198 + t174;
t169 = qJD(5) * t133;
t170 = qJD(5) * t131;
t144 = t131 * t84 + t133 * t85 - t142 * t170 + t94 * t169;
t204 = t125 * t150;
t211 = t144 - t204;
t209 = -0.2e1 * t168;
t126 = qJD(2) * qJD(4);
t160 = qJD(2) * t186;
t89 = t134 * qJD(3) + t132 * t160;
t72 = t89 * qJD(1);
t143 = -t132 * qJD(3) + t134 * t160;
t73 = t143 * qJD(1);
t36 = t130 * t73 + t180 * t72;
t33 = t126 + t36;
t18 = pkin(7) * t84 + t33;
t35 = t130 * t72 - t180 * t73;
t21 = -pkin(7) * t85 + t35;
t105 = qJD(2) * pkin(2) + t111;
t56 = t180 * t105 + t178;
t145 = qJD(4) - t56;
t26 = -t201 * qJD(2) + t145 - t198;
t199 = pkin(7) * t94;
t159 = t180 * t112;
t57 = t130 * t105 - t159;
t52 = qJD(2) * qJ(4) + t57;
t32 = t52 + t199;
t6 = t131 * t26 + t133 * t32;
t2 = -qJD(5) * t6 - t131 * t18 + t133 * t21;
t208 = -t125 * t6 + t2;
t202 = t142 ^ 2;
t90 = t94 ^ 2;
t207 = -t90 - t202;
t206 = -t90 + t202;
t205 = 0.2e1 * t142;
t1 = t131 * t21 + t133 * t18 + t26 * t169 - t32 * t170;
t203 = t150 * t25 - t1;
t120 = pkin(2) * t163;
t157 = t85 * qJ(4) - t120;
t146 = qJD(4) * t142 + t157;
t11 = -t201 * t84 + t146;
t200 = pkin(3) * t84;
t60 = t111 * t130 - t159;
t37 = t60 + t199;
t123 = -t180 * pkin(2) - pkin(3);
t119 = -pkin(4) + t123;
t121 = pkin(2) * t130 + qJ(4);
t74 = t119 * t133 - t121 * t131;
t197 = t74 * qJD(5) - t131 * t37 + t213 * t133;
t75 = t119 * t131 + t121 * t133;
t196 = -t75 * qJD(5) - t213 * t131 - t133 * t37;
t195 = pkin(2) * t132;
t5 = -t131 * t32 + t133 * t26;
t194 = t125 * t5;
t64 = -t180 * t114 - t115 * t130;
t192 = t35 * t64;
t40 = t94 * pkin(3) + t212;
t191 = t40 * t142;
t187 = t142 * t94;
t43 = t130 * t143 + t180 * t89;
t179 = qJD(2) * t96;
t136 = qJD(1) ^ 2;
t177 = t134 * t136;
t135 = qJD(2) ^ 2;
t176 = t135 * t132;
t175 = t135 * t134;
t65 = t130 * t114 - t180 * t115;
t173 = t132 ^ 2 - t134 ^ 2;
t171 = qJD(2) * t132;
t167 = pkin(2) * t171;
t166 = pkin(2) * t172;
t165 = t132 * t177;
t42 = t130 * t89 - t180 * t143;
t156 = pkin(1) * t209;
t155 = t125 ^ 2;
t154 = t134 * t163;
t153 = qJD(2) * t61 - t36;
t106 = t130 * t132 - t158;
t151 = t106 * t84 + t94 * t96;
t44 = -pkin(7) * t107 + t64;
t45 = pkin(7) * t106 + t65;
t9 = -t131 * t45 + t133 * t44;
t10 = t131 * t44 + t133 * t45;
t59 = t106 * t131 + t107 * t133;
t149 = t107 * qJ(4) + t164;
t148 = -t94 * qJ(4) - t166;
t147 = qJD(2) * t60 - t35;
t99 = qJD(2) * t158 - t130 * t171;
t141 = t99 * qJ(4) + t107 * qJD(4) - t167;
t140 = t107 * t35 + t142 * t42 - t43 * t94 + t64 * t85 - t65 * t84;
t139 = t106 * t85 + t107 * t84 + t142 * t96 + t94 * t99;
t137 = t205 * qJD(2);
t87 = t99 * qJD(2);
t58 = -t133 * t106 + t107 * t131;
t55 = pkin(3) * t106 - t149;
t54 = -t116 + (t152 + t94) * qJD(2);
t53 = -t116 + (t152 - t94) * qJD(2);
t46 = -qJD(2) * pkin(3) + t145;
t41 = pkin(3) * t142 - t148;
t34 = -t201 * t106 + t149;
t31 = pkin(3) * t96 - t141;
t29 = pkin(7) * t96 + t43;
t28 = -pkin(7) * t99 + t42;
t27 = -t142 * t201 + t148;
t23 = t107 * t85 + t142 * t99;
t22 = -t146 + t200;
t17 = -t201 * t96 + t141;
t16 = t59 * qJD(5) + t131 * t99 - t133 * t96;
t15 = -t106 * t169 + t107 * t170 - t131 * t96 - t133 * t99;
t4 = -t10 * qJD(5) - t131 * t29 + t133 * t28;
t3 = t9 * qJD(5) + t131 * t28 + t133 * t29;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t154, t173 * t209, t175, -0.2e1 * t154, -t176, 0, -pkin(6) * t175 + t132 * t156, pkin(6) * t176 + t134 * t156, 0, 0, t23, -t139, t87, t151, -t179, 0, t113 * t96 - t164 * t84 + (-t42 + (qJD(1) * t106 + t94) * t195) * qJD(2), t113 * t99 - t164 * t85 + (t205 * t195 - t43) * qJD(2), -t106 * t36 - t56 * t99 - t57 * t96 + t140, t192 + t36 * t65 - t42 * t56 + t43 * t57 + (t113 - t210) * t167, t23, t87, t139, 0, t179, t151, -t42 * qJD(2) + t106 * t22 + t31 * t94 + t40 * t96 + t55 * t84, -t106 * t33 + t46 * t99 - t52 * t96 + t140, t43 * qJD(2) - t107 * t22 - t142 * t31 - t40 * t99 - t55 * t85, t22 * t55 + t31 * t40 + t33 * t65 + t42 * t46 + t43 * t52 + t192, t144 * t59 - t15 * t51, -t13 * t59 - t144 * t58 + t15 * t150 - t16 * t51, t15 * t125, t13 * t58 + t150 * t16, t16 * t125, 0, t11 * t58 - t125 * t4 + t13 * t34 + t150 * t17 + t16 * t25, t11 * t59 + t125 * t3 + t144 * t34 - t15 * t25 + t17 * t51, -t1 * t58 - t10 * t13 - t144 * t9 + t15 * t5 - t150 * t3 - t16 * t6 - t2 * t59 - t4 * t51, t1 * t10 + t11 * t34 + t17 * t25 + t2 * t9 + t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t173 * t136, 0, t165, 0, 0, t136 * pkin(1) * t132, pkin(1) * t177, 0, 0, t187, t206, t54, -t187, 0, 0, -t113 * t142 - t94 * t166 + t147, t113 * t94 - t142 * t166 + t153, (t57 - t60) * t142 + (-t56 + t61) * t94 + (-t130 * t84 - t180 * t85) * pkin(2), t56 * t60 - t57 * t61 + (-t113 * t172 + t130 * t36 - t180 * t35) * pkin(2), t187, t54, -t206, 0, 0, -t187, -t41 * t94 + t147 - t191, -t121 * t84 + t123 * t85 + (t52 - t60) * t142 + (t46 - t174) * t94, t142 * t41 - t40 * t94 + 0.2e1 * t126 - t153, t121 * t33 + t123 * t35 + t174 * t52 - t40 * t41 - t46 * t60, -t188, t216, -t211, t188, t217, 0, -t196 * t125 - t150 * t27 - t2 + t215, t197 * t125 - t27 * t51 - t203, -t144 * t74 - t13 * t75 + (-t197 + t5) * t150 + (-t196 - t6) * t51, t1 * t75 + t196 * t5 + t197 * t6 + t2 * t74 - t25 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t53, t207, t142 * t56 + t57 * t94 + t120, 0, 0, 0, 0, 0, 0, t137, t207, -t53, t200 + t52 * t94 + (-qJD(4) - t46) * t142 - t157, 0, 0, 0, 0, 0, 0, -t13 + t214, -t144 - t204, t189 + t190, -t150 * t6 - t5 * t51 - t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t54, -t202 - t135, -qJD(2) * t52 + t191 + t35, 0, 0, 0, 0, 0, 0, -t131 * t155 - t142 * t150, -t133 * t155 - t142 * t51, -t217 * t131 - t211 * t133, -t25 * t142 + t208 * t133 + (t1 + t194) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, -t216, t211, -t188, -t217, 0, t208 - t215, -t194 + t203, 0, 0;];
tauc_reg = t7;

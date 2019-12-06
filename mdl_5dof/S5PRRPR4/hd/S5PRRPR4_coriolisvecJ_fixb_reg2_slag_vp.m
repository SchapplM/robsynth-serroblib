% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:22
% EndTime: 2019-12-05 16:23:30
% DurationCPUTime: 1.71s
% Computational Cost: add. (2482->242), mult. (6414->347), div. (0->0), fcn. (4743->8), ass. (0->145)
t125 = sin(pkin(9));
t126 = cos(pkin(9));
t128 = sin(qJ(3));
t131 = cos(qJ(3));
t107 = t125 * t131 + t126 * t128;
t132 = cos(qJ(2));
t159 = t132 * qJD(1);
t183 = -qJ(4) - pkin(6);
t149 = qJD(3) * t183;
t94 = t131 * qJD(4) + t128 * t149;
t95 = -t128 * qJD(4) + t131 * t149;
t181 = t107 * t159 - t125 * t94 + t126 * t95;
t173 = t126 * t131;
t106 = t125 * t128 - t173;
t138 = t106 * t132;
t180 = qJD(1) * t138 + t125 * t95 + t126 * t94;
t127 = sin(qJ(5));
t130 = cos(qJ(5));
t101 = t107 * qJD(2);
t154 = qJD(2) * t173;
t165 = qJD(2) * t128;
t98 = t125 * t165 - t154;
t141 = -t130 * t101 + t127 * t98;
t100 = t107 * qJD(3);
t85 = qJD(2) * t100;
t158 = qJD(2) * qJD(3);
t153 = t128 * t158;
t113 = t125 * t153;
t152 = t131 * t158;
t86 = t126 * t152 - t113;
t135 = t141 * qJD(5) - t127 * t86 - t130 * t85;
t122 = qJD(3) + qJD(5);
t176 = t122 * t141;
t204 = t135 - t176;
t160 = qJD(5) * t130;
t161 = qJD(5) * t127;
t139 = -t101 * t161 - t127 * t85 + t130 * t86 - t98 * t160;
t49 = -t101 * t127 - t130 * t98;
t175 = t122 * t49;
t203 = t139 - t175;
t185 = t49 ^ 2;
t186 = t141 ^ 2;
t202 = -t185 + t186;
t103 = t106 * qJD(3);
t201 = pkin(7) * t103 + t181;
t200 = -pkin(7) * t100 + t180;
t184 = t49 * t141;
t129 = sin(qJ(2));
t166 = qJD(1) * t129;
t178 = qJD(2) * pkin(6);
t114 = t166 + t178;
t143 = qJD(4) + t159;
t163 = qJD(3) * t128;
t51 = -t114 * t163 + (-qJ(4) * t163 + t143 * t131) * qJD(2);
t162 = qJD(3) * t131;
t52 = -t114 * t162 + (-qJ(4) * t162 - t143 * t128) * qJD(2);
t18 = -t125 * t51 + t126 * t52;
t12 = -pkin(7) * t86 + t18;
t19 = t125 * t52 + t126 * t51;
t13 = -pkin(7) * t85 + t19;
t187 = pkin(7) * t101;
t146 = qJ(4) * qJD(2) + t114;
t90 = t146 * t131;
t66 = t125 * t90;
t89 = t146 * t128;
t70 = qJD(3) * pkin(3) - t89;
t31 = t126 * t70 - t66;
t23 = qJD(3) * pkin(4) - t187 + t31;
t193 = pkin(7) * t98;
t174 = t126 * t90;
t32 = t125 * t70 + t174;
t24 = t32 - t193;
t6 = t127 * t23 + t130 * t24;
t2 = -t6 * qJD(5) + t130 * t12 - t127 * t13;
t121 = -pkin(3) * t131 - pkin(2);
t104 = t121 * qJD(2) + qJD(4) - t159;
t56 = pkin(4) * t98 + t104;
t199 = t141 * t56 + t2;
t1 = (qJD(5) * t23 + t13) * t130 + t127 * t12 - t24 * t161;
t198 = -t49 * t56 - t1;
t197 = -0.2e1 * t158;
t133 = qJD(3) ^ 2;
t134 = qJD(2) ^ 2;
t196 = (t133 + t134) * t129;
t195 = pkin(3) * t163 - t166;
t194 = t101 ^ 2;
t111 = t183 * t128;
t112 = t183 * t131;
t57 = t126 * t111 + t112 * t125;
t41 = -pkin(7) * t107 + t57;
t58 = t125 * t111 - t126 * t112;
t42 = -pkin(7) * t106 + t58;
t15 = t127 * t41 + t130 * t42;
t192 = t15 * qJD(5) + t200 * t127 - t201 * t130;
t14 = -t127 * t42 + t130 * t41;
t191 = t14 * qJD(5) + t201 * t127 + t200 * t130;
t33 = t125 * t89 - t174;
t25 = t33 + t193;
t34 = -t126 * t89 - t66;
t26 = t34 - t187;
t119 = pkin(3) * t126 + pkin(4);
t189 = pkin(3) * t125;
t93 = t119 * t127 + t130 * t189;
t190 = -t93 * qJD(5) + t127 * t26 - t130 * t25;
t188 = pkin(3) * t128;
t92 = t119 * t130 - t127 * t189;
t182 = -t92 * qJD(5) + t127 * t25 + t130 * t26;
t179 = qJD(2) * pkin(2);
t177 = t101 * t98;
t172 = t133 * t128;
t171 = t133 * t131;
t170 = t134 * t132;
t164 = qJD(2) * t129;
t105 = pkin(3) * t153 + qJD(1) * t164;
t123 = t128 ^ 2;
t124 = t131 ^ 2;
t169 = t123 - t124;
t168 = t123 + t124;
t157 = pkin(3) * t165;
t155 = t128 * t134 * t131;
t147 = pkin(4) * t100 + t195;
t115 = -t159 - t179;
t145 = -t115 - t159;
t144 = t132 * t197;
t53 = pkin(4) * t85 + t105;
t142 = t128 * t152;
t87 = t107 * t129;
t88 = t106 * t129;
t35 = t127 * t88 - t130 * t87;
t36 = -t127 * t87 - t130 * t88;
t55 = -t106 * t127 + t107 * t130;
t140 = qJD(2) * t145;
t137 = qJD(3) * (-t145 - t179);
t96 = t98 ^ 2;
t74 = pkin(4) * t106 + t121;
t61 = pkin(4) * t101 + t157;
t54 = t130 * t106 + t107 * t127;
t38 = -qJD(2) * t138 - qJD(3) * t87;
t37 = -t132 * t101 + t129 * t103;
t21 = t55 * qJD(5) + t130 * t100 - t127 * t103;
t20 = t127 * t100 + t130 * t103 + t106 * t160 + t107 * t161;
t8 = -t36 * qJD(5) - t127 * t38 + t130 * t37;
t7 = t35 * qJD(5) + t127 * t37 + t130 * t38;
t5 = -t127 * t24 + t130 * t23;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t129, -t170, 0, 0, 0, 0, 0, 0, 0, 0, t128 * t144 - t131 * t196, t128 * t196 + t131 * t144, t168 * t170, (t115 * t129 + (-t166 + (t114 + t166) * t168) * t132) * qJD(2), 0, 0, 0, 0, 0, 0, qJD(3) * t37 - t132 * t85 + t164 * t98, -qJD(3) * t38 + t101 * t164 - t132 * t86, -t101 * t37 - t38 * t98 + t85 * t88 + t86 * t87, t104 * t164 - t105 * t132 - t18 * t87 - t19 * t88 + t31 * t37 + t32 * t38, 0, 0, 0, 0, 0, 0, t122 * t8 + t132 * t135 - t164 * t49, -t122 * t7 - t132 * t139 - t141 * t164, t135 * t36 - t139 * t35 + t141 * t8 + t49 * t7, t1 * t36 - t132 * t53 + t164 * t56 + t2 * t35 + t5 * t8 + t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t142, t169 * t197, t171, -0.2e1 * t142, -t172, 0, -pkin(6) * t171 + t128 * t137, pkin(6) * t172 + t131 * t137, 0, ((-t115 - t179) * t129 + (-t114 + t178) * t132 * t168) * qJD(1), -t101 * t103 + t107 * t86, -t100 * t101 + t103 * t98 - t106 * t86 - t107 * t85, -t103 * qJD(3), t100 * t98 + t106 * t85, -t100 * qJD(3), 0, -t98 * t166 + t100 * t104 + t105 * t106 + t121 * t85 + (t98 * t188 + t181) * qJD(3), -t101 * t166 - t103 * t104 + t105 * t107 + t121 * t86 + (t101 * t188 - t180) * qJD(3), -t100 * t32 - t181 * t101 + t103 * t31 - t106 * t19 - t107 * t18 - t180 * t98 - t57 * t86 - t58 * t85, t195 * t104 + t105 * t121 + t18 * t57 + t180 * t32 + t181 * t31 + t19 * t58, t139 * t55 + t141 * t20, t135 * t55 - t139 * t54 + t141 * t21 - t20 * t49, -t20 * t122, -t135 * t54 - t21 * t49, -t21 * t122, 0, -t192 * t122 - t135 * t74 - t147 * t49 + t21 * t56 + t53 * t54, -t191 * t122 + t139 * t74 - t141 * t147 - t20 * t56 + t53 * t55, -t1 * t54 + t135 * t15 - t139 * t14 - t141 * t192 + t191 * t49 - t2 * t55 + t20 * t5 - t21 * t6, t1 * t15 + t14 * t2 + t147 * t56 + t191 * t6 - t192 * t5 + t53 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t169 * t134, 0, t155, 0, 0, t128 * t140, t131 * t140, 0, 0, t177, -t96 + t194, -t113 + (t98 + t154) * qJD(3), -t177, 0, 0, -qJD(3) * t33 - t101 * t104 - t157 * t98 + t18, qJD(3) * t34 - t101 * t157 + t104 * t98 - t19, (-t31 + t34) * t98 + (t32 + t33) * t101 + (-t125 * t85 - t126 * t86) * pkin(3), -t31 * t33 - t32 * t34 + (-t104 * t165 + t125 * t19 + t126 * t18) * pkin(3), t184, t202, t203, -t184, t204, 0, t190 * t122 + t49 * t61 + t199, t182 * t122 + t141 * t61 + t198, -t139 * t92 + t135 * t93 + (-t182 + t5) * t49 + (t190 - t6) * t141, t1 * t93 - t182 * t6 + t190 * t5 + t2 * t92 - t56 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101 * qJD(3), -t113 + (-t98 + t154) * qJD(3), -t96 - t194, t101 * t31 + t32 * t98 + t105, 0, 0, 0, 0, 0, 0, -t135 - t176, t139 + t175, -t185 - t186, -t141 * t5 - t49 * t6 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t202, t203, -t184, t204, 0, t122 * t6 + t199, t122 * t5 + t198, 0, 0;];
tauc_reg = t3;

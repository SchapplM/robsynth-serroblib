% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:32
% EndTime: 2019-12-05 16:13:41
% DurationCPUTime: 1.84s
% Computational Cost: add. (1183->254), mult. (3122->361), div. (0->0), fcn. (1995->6), ass. (0->157)
t90 = sin(qJ(3));
t86 = t90 ^ 2;
t92 = cos(qJ(3));
t87 = t92 ^ 2;
t124 = qJD(2) * (t86 - 0.2e1 * t87);
t89 = cos(pkin(8));
t147 = t89 * qJD(3);
t154 = qJD(2) * t90;
t88 = sin(pkin(8));
t60 = t88 * t154 - t147;
t201 = (t124 * t88 + t60 * t90) * qJD(3);
t151 = qJD(3) * t92;
t91 = sin(qJ(2));
t131 = t91 * t151;
t171 = t91 * t92;
t139 = t89 * t171;
t93 = cos(qJ(2));
t170 = t93 * t89;
t178 = t88 * t91;
t102 = t92 * t170 + t178;
t152 = qJD(3) * t90;
t132 = t91 * t152;
t30 = t102 * qJD(2) - t89 * t132;
t53 = -t93 * t88 + t139;
t130 = t89 * t154;
t148 = t88 * qJD(3);
t62 = t130 + t148;
t200 = ((qJD(3) * (-t53 + t139) + t62 * t93) * t90 + t30 * t92) * qJD(2) + t62 * t131;
t153 = qJD(2) * t92;
t182 = t60 * t89;
t84 = t88 ^ 2;
t85 = t89 ^ 2;
t199 = ((t84 - t85) * qJD(3) - t62 * t88 - t182) * t153;
t198 = ((t62 + 0.2e1 * t130) * t88 + t182) * t151;
t100 = t93 * t154 + t131;
t141 = qJD(2) * qJD(3);
t195 = -0.2e1 * t141;
t109 = pkin(3) * t90 - qJ(4) * t92;
t50 = t109 * qJD(3) - t90 * qJD(4);
t194 = -t102 * qJD(1) + t88 * t50;
t58 = t62 ^ 2;
t193 = -t60 ^ 2 - t58;
t144 = t93 * qJD(1);
t146 = t91 * qJD(1);
t177 = t88 * t92;
t191 = t144 * t177 + (-t146 + t50) * t89;
t126 = t92 * t141;
t115 = t89 * t126;
t127 = qJD(2) * t144;
t76 = t90 * t127;
t156 = qJD(2) * pkin(6);
t77 = t146 + t156;
t37 = t77 * t151 + t76;
t75 = t88 * t126;
t98 = pkin(4) * t75 - qJ(5) * t115 + t37;
t5 = -t62 * qJD(5) + t98;
t187 = t5 * t88;
t186 = t5 * t89;
t184 = t37 * t88;
t183 = t37 * t89;
t95 = qJD(2) ^ 2;
t179 = t87 * t95;
t65 = t109 * qJD(2);
t175 = t89 * t65;
t155 = t90 * qJ(4);
t110 = -t92 * pkin(3) - t155;
t70 = -pkin(2) + t110;
t174 = t89 * t70;
t173 = t89 * t92;
t67 = t90 * t77;
t172 = t90 * t91;
t68 = t92 * t77;
t94 = qJD(3) ^ 2;
t169 = t94 * t90;
t168 = t94 * t92;
t167 = t95 * t93;
t145 = t92 * qJD(5);
t166 = -t145 + (-pkin(6) * t89 + qJ(5)) * t152 + t194;
t129 = pkin(6) * t88 + pkin(4);
t165 = -t129 * t152 - t191;
t138 = pkin(6) * t152;
t120 = t88 * t138;
t164 = t120 + t191;
t119 = t89 * t138;
t163 = -t119 + t194;
t32 = (t50 + t146) * qJD(2);
t35 = t92 * t127 + (qJD(4) - t67) * qJD(3);
t4 = t88 * t32 + t89 * t35;
t43 = t70 * qJD(2) - t144;
t142 = qJD(3) * qJ(4);
t55 = t68 + t142;
t10 = t88 * t43 + t89 * t55;
t39 = pkin(6) * t173 + t88 * t70;
t160 = t86 - t87;
t159 = t86 + t87;
t158 = t94 + t95;
t157 = qJD(2) * pkin(2);
t150 = qJD(4) * t62;
t149 = qJD(4) * t89;
t143 = qJ(5) * qJD(2);
t140 = t88 * t179;
t137 = t60 * t144;
t136 = t62 * t144;
t135 = t90 * t144;
t134 = t88 * t153;
t128 = t90 * t142;
t125 = t90 * t141;
t3 = t89 * t32 - t88 * t35;
t78 = -t144 - t157;
t123 = -t78 - t144;
t122 = t60 + t147;
t108 = pkin(4) * t88 - qJ(5) * t89;
t24 = t108 * t153 + t68;
t121 = qJD(5) * t88 + t24;
t118 = -qJD(3) * pkin(3) + qJD(4);
t116 = t93 * t195;
t114 = t92 * t125;
t112 = -t62 * t153 + t75;
t51 = t118 + t67;
t9 = t89 * t43 - t88 * t55;
t106 = qJD(2) * t122;
t105 = qJD(2) * (-t62 + t148);
t104 = t60 * t92 * t148 + t84 * t114;
t103 = pkin(6) + t108;
t101 = t106 * t177;
t99 = qJD(3) * (-t123 - t157);
t29 = -qJD(2) * t91 * t89 - t88 * t132 + t93 * t134;
t52 = t88 * t171 + t170;
t97 = t114 * t178 + (-t52 * t152 + t29 * t92) * qJD(2) + t100 * t60;
t96 = t29 * t62 - t30 * t60 + (t52 * t89 - t53 * t88) * t126;
t81 = t90 * t95 * t92;
t74 = qJD(4) * t134;
t73 = -0.2e1 * t114;
t69 = -t89 * pkin(4) - t88 * qJ(5) - pkin(3);
t54 = t88 * t65;
t49 = t60 * t153;
t48 = t60 * t149;
t44 = t103 * t90;
t38 = -pkin(6) * t177 + t174;
t34 = t129 * t92 - t174;
t33 = -t92 * qJ(5) + t39;
t31 = t49 + t115;
t22 = -t89 * t67 + t54;
t21 = t88 * t67 + t175;
t20 = -t89 * t90 * qJD(5) + t103 * t151;
t19 = (t85 * t154 + t62 * t89) * t151;
t18 = t105 * t173;
t16 = t90 * t105 + t89 * t179;
t15 = -t175 + (-pkin(4) * qJD(2) - t77 * t88) * t90;
t14 = (t89 * t124 + t62 * t90) * qJD(3);
t13 = t54 + (-t77 * t89 + t143) * t90;
t8 = t60 * pkin(4) - t62 * qJ(5) + t51;
t7 = -t92 * t143 + t10;
t6 = pkin(4) * t153 + qJD(5) - t9;
t2 = -pkin(4) * t125 - t3;
t1 = (qJ(5) * t152 - t145) * qJD(2) + t4;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95 * t91, -t167, 0, 0, 0, 0, 0, 0, 0, 0, t90 * t116 - t158 * t171, t92 * t116 + t158 * t172, t159 * t167, (t78 * t91 + (-t146 + (t77 + t146) * t159) * t93) * qJD(2), 0, 0, 0, 0, 0, 0, t97, t200, t96, t10 * t30 + t100 * t51 + t37 * t172 - t9 * t29 - t3 * t52 + t4 * t53, 0, 0, 0, 0, 0, 0, t97, t96, -t200, t1 * t53 + t100 * t8 + t5 * t172 + t2 * t52 + t6 * t29 + t7 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t114, t160 * t195, t168, t73, -t169, 0, -pkin(6) * t168 + t90 * t99, pkin(6) * t169 + t92 * t99, 0, ((-t78 - t157) * t91 + (t156 - t77) * t93 * t159) * qJD(1), t19, -t198, t14, t104, -t201, t73, (-t137 + t184 + (qJD(2) * t38 + t9) * qJD(3)) * t90 + (-t3 + (pkin(6) * t60 + t51 * t88) * qJD(3) + (t120 - t164) * qJD(2)) * t92, (-t136 + t183 + (-qJD(2) * t39 - t10) * qJD(3)) * t90 + (t4 + (pkin(6) * t62 + t51 * t89) * qJD(3) + (t119 + t163) * qJD(2)) * t92, (-t3 * t89 - t4 * t88) * t90 - t164 * t62 - t163 * t60 + (-t10 * t88 - t89 * t9 + (-t38 * t89 - t39 * t88) * qJD(2)) * t151, -t51 * t135 + t3 * t38 + t4 * t39 + t164 * t9 + t163 * t10 + (t151 * t51 + t37 * t90) * pkin(6), t19, t14, t198, t73, t201, t104, t20 * t60 + (-t137 + t187 + (-qJD(2) * t34 - t6) * qJD(3)) * t90 + (t8 * t148 + t2 + (t148 * t44 + t165) * qJD(2)) * t92, (-t1 * t88 + t2 * t89) * t90 + t165 * t62 - t166 * t60 + (t6 * t89 - t7 * t88 + (-t33 * t88 + t34 * t89) * qJD(2)) * t151, -t20 * t62 + (t136 - t186 + (qJD(2) * t33 + t7) * qJD(3)) * t90 + (-t8 * t147 - t1 + (-t147 * t44 - t166) * qJD(2)) * t92, t1 * t33 + t2 * t34 + t5 * t44 + (t20 - t135) * t8 + t166 * t7 + t165 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t160 * t95, 0, t81, 0, 0, -t78 * t154 - t76, t123 * t153, 0, 0, t18, -t199, t16, -t101, t122 * t154 - t140, t81, -t60 * t68 - t183 + t74 + (t21 * t92 - t9 * t90 + (t110 * qJD(3) - t51 * t92) * t88) * qJD(2), -t62 * t68 + t184 + (t10 * t90 - t22 * t92 + (-t128 + (t118 - t51) * t92) * t89) * qJD(2), t21 * t62 + t22 * t60 - t48 + (t153 * t9 + t4) * t89 + (t10 * t153 + t150 - t3) * t88, -t51 * t68 - t37 * pkin(3) - t10 * t22 - t9 * t21 + (t10 * t89 - t88 * t9) * qJD(4) + (-t3 * t88 + t4 * t89) * qJ(4), t18, t16, t199, t81, -t106 * t90 + t140, -t101, -t186 + t74 - t121 * t60 + (-t15 * t92 + t6 * t90 + (-t8 * t92 + (t69 * t92 - t155) * qJD(3)) * t88) * qJD(2), t13 * t60 - t15 * t62 - t48 + (-t153 * t6 + t1) * t89 + (t153 * t7 + t150 + t2) * t88, -t187 + t121 * t62 + (t13 * t92 - t7 * t90 + (t128 + (-qJD(3) * t69 - qJD(4) + t8) * t92) * t89) * qJD(2), t1 * t89 * qJ(4) - t6 * t15 - t8 * t24 + t5 * t69 + (-t13 + t149) * t7 + (qJ(4) * t2 + qJD(4) * t6 - qJD(5) * t8) * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t31, t193, t10 * t60 + t9 * t62 + t37, 0, 0, 0, 0, 0, 0, t112, t193, -t31, t7 * t60 + (-qJD(5) - t6) * t62 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60 - t125, -t49 + t115, -t58 - t179, t8 * t62 + (-pkin(4) * t152 + t7 * t92) * qJD(2) - t3;];
tauc_reg = t11;

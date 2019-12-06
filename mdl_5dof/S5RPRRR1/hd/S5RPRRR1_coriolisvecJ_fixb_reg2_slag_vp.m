% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:51
% EndTime: 2019-12-05 18:10:02
% DurationCPUTime: 2.46s
% Computational Cost: add. (1403->284), mult. (4282->439), div. (0->0), fcn. (3099->6), ass. (0->169)
t60 = cos(qJ(4));
t152 = qJD(4) * t60;
t61 = cos(qJ(3));
t166 = t60 * t61;
t56 = sin(qJ(5));
t58 = sin(qJ(3));
t175 = t56 * t58;
t59 = cos(qJ(5));
t82 = t59 * t166 + t175;
t26 = t82 * qJD(1);
t210 = t59 * t152 - t26;
t140 = qJD(2) * qJD(4);
t118 = t60 * t140;
t146 = t60 * qJD(2);
t126 = t61 * t146;
t151 = qJD(4) * t61;
t57 = sin(qJ(4));
t132 = t57 * t151;
t51 = t60 * qJD(3);
t10 = t118 + (t126 + (-t58 * t51 - t132) * qJ(2)) * qJD(1);
t150 = qJD(5) * t56;
t144 = qJ(2) * qJD(1);
t123 = t61 * t144;
t147 = t57 * qJD(2);
t35 = t60 * t123 + t147;
t149 = qJD(5) * t59;
t154 = qJD(3) * t61;
t70 = t58 * t149 + t56 * t154;
t1 = (t70 * qJ(2) + qJD(2) * t175) * qJD(1) + t59 * t10 - t35 * t150;
t53 = t61 * qJD(1);
t201 = t53 - qJD(4);
t113 = t59 * t201;
t156 = qJD(3) * t57;
t157 = qJD(1) * t58;
t38 = t60 * t157 + t156;
t16 = t38 * t56 + t113;
t124 = t58 * t144;
t19 = t59 * t124 - t35 * t56;
t36 = t57 * t157 - t51;
t32 = qJD(5) + t36;
t34 = t57 * t123 - t146;
t209 = -t34 * t16 - t19 * t32 + t1;
t86 = t201 * t34;
t208 = t35 * t201;
t203 = t201 * t61;
t207 = t60 * t203;
t206 = t57 * t150 - t210;
t54 = t58 ^ 2;
t148 = t54 * qJD(1);
t75 = t203 + t148;
t205 = t206 * t32;
t204 = t201 * t16;
t18 = -t201 * t56 + t59 * t38;
t142 = qJD(1) * qJD(3);
t115 = t61 * t142;
t141 = qJD(2) * qJD(1);
t120 = t58 * t141;
t20 = t56 * t124 + t35 * t59;
t2 = -t20 * qJD(5) - t56 * t10 + (qJ(2) * t115 + t120) * t59;
t200 = t34 * t18 - t20 * t32 - t2;
t176 = t56 * t57;
t116 = t58 * t142;
t139 = qJD(3) * qJD(4);
t153 = qJD(4) * t57;
t133 = t58 * t153;
t72 = t61 * t51 - t133;
t21 = -t72 * qJD(1) - t60 * t139;
t6 = qJD(5) * t113 - t56 * t116 + t38 * t150 + t59 * t21;
t136 = t56 * t166;
t25 = qJD(1) * t136 - t59 * t157;
t97 = t56 * t152 - t25;
t69 = t57 * t149 + t97;
t199 = t6 * t176 - t18 * t69;
t101 = -qJD(5) + t51;
t197 = t101 * t58 + t132;
t71 = t58 * t152 + t57 * t154;
t22 = t71 * qJD(1) + t57 * t139;
t7 = t18 * qJD(5) - t59 * t116 - t56 * t21;
t196 = t60 * t7;
t100 = t57 * t124;
t93 = -qJD(3) * t100 + t57 * t140;
t11 = (qJ(2) * t152 + t147) * t53 + t93;
t195 = t11 * t56;
t194 = t11 * t59;
t193 = t16 * t32;
t192 = t18 * t16;
t191 = t18 * t32;
t188 = t21 * t57;
t187 = t22 * t60;
t184 = t36 * t57;
t183 = t36 * t58;
t182 = t36 * t60;
t181 = t38 * t36;
t180 = t38 * t57;
t179 = t38 * t58;
t178 = t38 * t60;
t177 = t56 * t22;
t174 = t57 * t58;
t173 = t57 * t59;
t172 = t57 * t61;
t171 = t58 * t22;
t170 = t58 * t60;
t169 = t59 * t22;
t168 = t59 * t61;
t167 = t60 * t21;
t63 = qJD(1) ^ 2;
t165 = t61 * t63;
t164 = -t36 * t152 - t57 * t22;
t55 = t61 ^ 2;
t163 = t54 - t55;
t162 = t54 + t55;
t161 = qJ(2) * t58;
t160 = qJ(2) * t61;
t62 = qJD(3) ^ 2;
t159 = qJ(2) * t62;
t158 = t63 * qJ(2);
t155 = qJD(3) * t58;
t143 = qJ(2) * qJD(3);
t137 = t57 * t171;
t135 = t58 * t165;
t134 = t18 * t53;
t130 = t60 * t151;
t128 = t61 * t147;
t127 = t36 * t157;
t125 = 0.2e1 * t141;
t122 = t18 * t153 + t6 * t60;
t121 = t162 * t63;
t111 = t32 * t59;
t109 = t36 + t51;
t108 = -t38 + t156;
t107 = qJ(2) * t135;
t106 = qJD(1) * t201;
t105 = qJD(4) * t201;
t104 = 0.2e1 * t116;
t102 = -qJD(5) * t60 + qJD(3);
t99 = t58 * t115;
t95 = t57 * t105;
t94 = t60 * t105;
t92 = -t19 * t59 - t20 * t56;
t91 = t19 * t56 - t20 * t59;
t90 = t34 * t60 - t35 * t57;
t89 = t34 * t57 + t35 * t60;
t88 = -t180 + t182;
t87 = t180 + t182;
t85 = t201 * t57;
t84 = t102 * t61;
t83 = t162 * t125;
t31 = t59 * t170 - t56 * t61;
t81 = t58 * t59 - t136;
t30 = t56 * t170 + t168;
t80 = t97 * t32;
t79 = -t38 * t153 - t167;
t77 = -t32 * t149 - t177;
t76 = (-t201 + t53) * t58;
t74 = -t203 + t148;
t73 = -t106 * t172 + t60 * t116 + t95;
t67 = t206 * t16 - t7 * t173;
t66 = t75 * t60;
t65 = t90 * qJD(4) + t10 * t60 + t11 * t57;
t64 = qJ(2) ^ 2;
t28 = t82 * qJ(2);
t27 = t81 * qJ(2);
t24 = t31 * t144;
t23 = t30 * t144;
t9 = t101 * t168 + (t102 * t56 - t59 * t153) * t58;
t8 = -t56 * t133 - t61 * t150 - t59 * t155 + t70 * t60;
t5 = t81 * qJD(2) + (t197 * t56 + t59 * t84) * qJ(2);
t4 = t82 * qJD(2) + (-t197 * t59 + t56 * t84) * qJ(2);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, qJ(2) * t125, 0.2e1 * t99, -0.2e1 * t163 * t142, t62 * t61, -0.2e1 * t99, -t62 * t58, 0, -t61 * t159, t58 * t159, t83, qJ(2) * t83, -t58 * t167 + t38 * t72, -t87 * t154 + (t188 - t187 + (-t178 + t184) * qJD(4)) * t58, t58 * t95 + t21 * t61 + (t60 * t74 + t179) * qJD(3), t36 * t71 + t137, t58 * t94 + t22 * t61 + (-t57 * t74 - t183) * qJD(3), (-t201 - t53) * t155, -t34 * t155 + t11 * t61 + (t57 * t75 + t183) * qJD(2) + (t171 + qJD(4) * t66 + (t61 * t36 + t57 * t76) * qJD(3)) * qJ(2), -t35 * t155 + t10 * t61 + (t66 + t179) * qJD(2) + (-t58 * t21 - t75 * t153 + (t61 * t38 + t60 * t76) * qJD(3)) * qJ(2), ((-qJD(4) * t35 + t143 * t36 + t11) * t60 + (-qJD(4) * t34 - t143 * t38 - t10) * t57) * t58 + (t90 * qJD(3) - t88 * qJD(2) + (-t188 - t187 + (t178 + t184) * qJD(4)) * qJ(2)) * t61, (qJD(2) * t89 + t104 * t64) * t61 + (t125 * t54 - t155 * t89 + t61 * t65) * qJ(2), t18 * t9 - t31 * t6, -t16 * t9 - t18 * t8 + t30 * t6 - t31 * t7, -t174 * t6 + t18 * t71 + t31 * t22 + t9 * t32, t16 * t8 + t30 * t7, -t16 * t71 - t174 * t7 - t30 * t22 - t8 * t32, t32 * t71 + t137, t11 * t30 + t22 * t27 + t32 * t5 + t34 * t8 + (t16 * t160 + t19 * t58) * t152 + ((-t143 * t16 + t2) * t58 + (qJ(2) * t7 + qJD(2) * t16 + qJD(3) * t19) * t61) * t57, t11 * t31 - t22 * t28 - t32 * t4 + t34 * t9 + (t160 * t18 - t20 * t58) * t152 + ((-t143 * t18 - t1) * t58 + (-qJ(2) * t6 + qJD(2) * t18 - qJD(3) * t20) * t61) * t57, -t1 * t30 - t16 * t4 - t18 * t5 - t19 * t9 - t2 * t31 - t20 * t8 + t27 * t6 - t28 * t7, t34 * t128 + t1 * t28 + t19 * t5 + t2 * t27 + t20 * t4 + (t11 * t172 + (-t155 * t57 + t130) * t34) * qJ(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t158, 0, 0, 0, 0, 0, 0, t104, 0.2e1 * t115, -t121, -qJ(2) * t121, 0, 0, 0, 0, 0, 0, t73 - t127, t201 * t152 + (-t207 + (-t38 - t156) * t58) * qJD(1), t53 * t88 + t164 - t79, -t54 * t158 + (-t11 - t208) * t60 + (t10 - t86) * t57, 0, 0, 0, 0, 0, 0, -t196 - t80 + (t77 - t204) * t57, (-t134 - t169) * t57 + t205 + t122, -t199 + t67, t19 * t25 - t20 * t26 + (-qJD(4) * t91 - t11) * t60 + (qJD(5) * t92 + t1 * t59 - t2 * t56 - t86) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, t163 * t63, 0, t135, 0, 0, -0.2e1 * t120, -0.2e1 * t61 * t141, 0, 0, -t178 * t201 - t188, t87 * t53 + t164 + t79, -t94 + (t108 * t58 + t207) * qJD(1), -t36 * t85 - t187, t73 + t127, t58 * t106, -t57 * t107 + ((t34 - t146) * t58 + (-t109 * t61 + t174 * t53) * qJ(2)) * qJD(1), -t60 * t107 + ((t35 + t147) * t58 + (t108 * t61 + t170 * t53) * qJ(2)) * qJD(1), (-t161 * t88 - t61 * t90) * qJD(1) + t65, (t144 * t89 - t165 * t64) * t58, -t173 * t6 - t18 * t206, t199 + t67, (-t134 + t169) * t57 - t205 + t122, t16 * t69 + t176 * t7, t196 - t80 + (t77 + t204) * t57, -t32 * t85 - t187, -t2 * t60 - t23 * t32 + t97 * t34 + (t34 * t149 + qJD(4) * t19 + t195 + (t16 * t161 - t19 * t61) * qJD(1)) * t57, t1 * t60 - t24 * t32 + t210 * t34 + (-t34 * t150 - qJD(4) * t20 + t194 + (t161 * t18 + t20 * t61) * qJD(1)) * t57, -t16 * t24 + t18 * t23 + t19 * t26 + t20 * t25 + t92 * t152 + (qJD(5) * t91 - t1 * t56 - t2 * t59) * t57, t100 * t34 - t19 * t23 + t20 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, -t36 ^ 2 + t38 ^ 2, -t201 * t36 - t21, -t181, -t201 * t38 - t22, t116, -t208 + (-t128 + (-t130 - t179) * qJ(2)) * qJD(1) - t93, -t118 + t86 + (-t126 + (t109 * t58 + t132) * qJ(2)) * qJD(1), 0, 0, t111 * t18 - t6 * t56, (-t6 - t193) * t59 + (-t7 - t191) * t56, t111 * t32 - t18 * t38 + t177, t56 * t193 - t7 * t59, -t32 ^ 2 * t56 + t16 * t38 + t169, -t32 * t38, -t35 * t16 - t19 * t38 - t194, -t35 * t18 + t20 * t38 + t195, t200 * t56 + t209 * t59, (-t35 - t91) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, -t16 ^ 2 + t18 ^ 2, -t6 + t193, -t192, t191 - t7, t22, -t200, -t209, 0, 0;];
tauc_reg = t3;

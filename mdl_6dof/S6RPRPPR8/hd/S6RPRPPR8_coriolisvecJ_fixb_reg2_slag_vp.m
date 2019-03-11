% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:00:00
% EndTime: 2019-03-09 03:00:05
% DurationCPUTime: 1.94s
% Computational Cost: add. (2230->286), mult. (4209->375), div. (0->0), fcn. (2031->4), ass. (0->169)
t83 = cos(qJ(3));
t150 = qJD(1) * t83;
t136 = qJ(5) * qJD(1);
t85 = -pkin(1) - pkin(7);
t60 = t85 * qJD(1) + qJD(2);
t34 = (t60 + t136) * t83;
t190 = qJD(4) - t34;
t84 = -pkin(3) - pkin(4);
t21 = qJD(3) * t84 + t190;
t82 = cos(qJ(6));
t148 = qJD(3) * t82;
t81 = sin(qJ(3));
t151 = qJD(1) * t81;
t80 = sin(qJ(6));
t43 = t80 * t151 + t148;
t64 = qJD(6) + t150;
t101 = t43 * t64;
t135 = qJD(1) * qJD(3);
t126 = t83 * t135;
t24 = qJD(6) * t43 - t82 * t126;
t191 = t24 - t101;
t160 = t83 * t60;
t48 = t81 * t60;
t76 = qJD(3) * qJ(4);
t40 = t48 + t76;
t172 = t40 * t83;
t36 = (qJD(4) + t160) * qJD(3);
t154 = qJD(3) * pkin(3);
t119 = -qJD(4) + t154;
t37 = -t119 - t160;
t189 = ((-t37 + t160) * t81 - t172) * qJD(3) - t36 * t81;
t131 = 0.2e1 * qJD(1);
t72 = t83 * qJD(4);
t65 = qJD(1) * t72;
t75 = -pkin(8) + t84;
t79 = qJ(4) + pkin(5);
t93 = t75 * t83 - t79 * t81;
t89 = t93 * qJD(3) - qJD(2);
t10 = t89 * qJD(1) + t65;
t134 = qJD(1) * qJD(5);
t127 = t81 * t135;
t149 = qJD(3) * t81;
t46 = t60 * t149;
t157 = qJ(5) * t127 + t46;
t20 = -t83 * t134 + t157;
t68 = qJ(4) * t150;
t138 = qJD(5) + t68;
t94 = t83 * pkin(5) + t75 * t81 - qJ(2);
t14 = t94 * qJD(1) + t138;
t19 = qJD(3) * t75 + t190;
t5 = t14 * t82 - t19 * t80;
t1 = qJD(6) * t5 + t80 * t10 + t82 * t20;
t188 = -t5 * t64 + t1;
t6 = t14 * t80 + t19 * t82;
t2 = -qJD(6) * t6 + t82 * t10 - t80 * t20;
t187 = -t6 * t64 - t2;
t141 = t80 * qJD(3);
t45 = t82 * t151 - t141;
t100 = t45 * t64;
t145 = qJD(6) * t82;
t95 = t83 * t141 + t81 * t145;
t25 = t95 * qJD(1) - qJD(6) * t141;
t97 = -t25 + t100;
t110 = t5 * t80 - t6 * t82;
t183 = t110 * qJD(6) - t1 * t80 - t2 * t82;
t111 = t5 * t82 + t6 * t80;
t182 = -qJD(6) * t111 + t1 * t82 - t2 * t80;
t16 = t81 * t134 + (qJD(4) + t34) * qJD(3);
t152 = qJ(5) + t85;
t49 = t152 * t81;
t179 = t16 * t49;
t178 = t16 * t80;
t177 = t16 * t82;
t176 = t24 * t83;
t175 = t25 * t82;
t174 = t25 * t83;
t171 = t43 * t81;
t170 = t45 * t43;
t169 = t45 * t81;
t168 = t64 * t80;
t167 = t64 * t82;
t166 = t64 * t83;
t77 = t81 ^ 2;
t87 = qJD(1) ^ 2;
t165 = t77 * t87;
t164 = t80 * t24;
t163 = t81 * t24;
t162 = t81 * t25;
t161 = t82 * t83;
t86 = qJD(3) ^ 2;
t159 = t86 * t81;
t158 = t86 * t83;
t156 = t86 + t87;
t155 = qJ(4) * t81;
t153 = t87 * qJ(2);
t147 = qJD(3) * t83;
t146 = qJD(6) * t80;
t67 = t81 * t136;
t29 = -t67 - t40;
t23 = qJD(3) * pkin(5) - t29;
t144 = t23 * qJD(3);
t143 = t29 * qJD(3);
t125 = pkin(3) * t81 + qJ(2);
t38 = t125 * qJD(1) - t68;
t142 = t38 * qJD(1);
t109 = t84 * t81 - qJ(2);
t28 = t109 * qJD(1) + t138;
t139 = qJD(5) + t28;
t137 = qJ(2) * qJD(3);
t133 = t80 * t166;
t132 = t64 * t161;
t130 = t81 * t146;
t129 = t64 * t145;
t128 = qJD(2) * t131;
t123 = t139 * t83;
t98 = t84 * t83 - t155;
t90 = t98 * qJD(3) - qJD(2);
t15 = t90 * qJD(1) + t65;
t26 = t72 + t90;
t122 = qJD(1) * t26 + t15;
t73 = t83 * qJ(4);
t41 = t73 + t109;
t121 = qJD(1) * t41 + t28;
t120 = qJD(3) * t152;
t117 = qJD(6) * t83 + qJD(1);
t116 = t81 * t126;
t57 = t80 * t127;
t115 = t57 - t129;
t58 = t82 * t127;
t114 = -t64 * t146 - t58;
t108 = pkin(3) * t83 + t155;
t27 = t73 + t94;
t50 = t152 * t83;
t11 = t27 * t82 + t50 * t80;
t12 = t27 * t80 - t50 * t82;
t106 = (-t64 - t150) * t81;
t105 = qJD(1) * t77 - t166;
t51 = t125 - t73;
t104 = qJD(3) * (qJD(1) * t51 + t38);
t92 = t108 * qJD(3) + qJD(2);
t22 = t92 * qJD(1) - t65;
t32 = -t72 + t92;
t99 = qJD(1) * t32 - t85 * t86 + t22;
t91 = t16 * t81 + t21 * t149 + (-t20 - t143) * t83;
t78 = t83 ^ 2;
t74 = 0.2e1 * qJD(3) * qJD(4);
t71 = t78 * t87;
t69 = qJ(2) * t128;
t63 = t83 * t87 * t81;
t61 = -t71 - t86;
t56 = -0.2e1 * t116;
t55 = 0.2e1 * t116;
t54 = t156 * t83;
t53 = t156 * t81;
t52 = -t71 + t165;
t47 = t108 * qJD(1);
t42 = (t77 - t78) * t135;
t39 = 0.2e1 * t42;
t35 = t98 * qJD(1);
t33 = t48 + t67;
t31 = t81 * qJD(5) + t83 * t120;
t30 = -t83 * qJD(5) + t81 * t120;
t18 = t93 * qJD(1);
t13 = t72 + t89;
t8 = t18 * t80 + t33 * t82;
t7 = t18 * t82 - t33 * t80;
t4 = -t12 * qJD(6) + t82 * t13 - t80 * t30;
t3 = t11 * qJD(6) + t80 * t13 + t82 * t30;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t69, t56, t39, -t159, t55, -t158, 0, -t85 * t159 + (qJD(2) * t81 + t83 * t137) * t131, -t85 * t158 + (qJD(2) * t83 - t81 * t137) * t131, 0, t69, t56, -t159, -0.2e1 * t42, 0, t158, t55, t83 * t104 + t99 * t81, t189, t81 * t104 - t99 * t83, -t189 * t85 + t22 * t51 + t38 * t32, t55, t39, -t158, t56, -t159, 0, t122 * t83 + (-t121 * t81 + t31) * qJD(3), t122 * t81 + (t121 * t83 + t30) * qJD(3) (-t30 * t83 + t31 * t81 + (t49 * t83 - t50 * t81) * qJD(3)) * qJD(1) + t91, t15 * t41 - t20 * t50 + t21 * t30 + t26 * t28 - t29 * t31 + t179, -t82 * t163 + (t82 * t147 - t130) * t45 (-t43 * t82 - t45 * t80) * t147 + (t164 - t175 + (t43 * t80 - t45 * t82) * qJD(6)) * t81, -t64 * t130 - t176 + (-t105 * t82 - t169) * qJD(3), t80 * t162 + t95 * t43, -t81 * t129 - t174 + (t105 * t80 + t171) * qJD(3), qJD(3) * t106, t25 * t49 + t31 * t43 + t4 * t64 + (t23 * t141 + t2) * t83 + (t23 * t145 + t178 + (-qJD(1) * t11 - t5) * qJD(3)) * t81, -t24 * t49 - t3 * t64 + t31 * t45 + (t82 * t144 - t1) * t83 + (-t23 * t146 + t177 + (qJD(1) * t12 + t6) * qJD(3)) * t81, t11 * t24 - t111 * t147 - t12 * t25 + t183 * t81 - t3 * t43 - t4 * t45, t1 * t12 + t11 * t2 + t23 * t31 + t3 * t6 + t4 * t5 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t153, 0, 0, 0, 0, 0, 0, -t53, -t54, 0, -t153, 0, 0, 0, 0, 0, 0, -t53, 0, t54, -t189 - t142, 0, 0, 0, 0, 0, 0, t54, t53, 0, t28 * qJD(1) + t91, 0, 0, 0, 0, 0, 0, t162 + t117 * t167 + (t80 * t106 + t43 * t83) * qJD(3), -t163 - t117 * t168 + (t82 * t106 + t45 * t83) * qJD(3) (-t117 * t45 - t43 * t149 + t174) * t82 + (-t117 * t43 + t45 * t149 + t176) * t80, t111 * qJD(1) + (-qJD(3) * t110 + t16) * t81 + (t144 - t182) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t52, 0, -t63, 0, 0, -t83 * t153, t81 * t153, 0, 0, t63, 0, t52, 0, 0, -t63 (-t38 * t83 - t47 * t81) * qJD(1) ((t40 - t76) * t83 + (t119 + t37) * t81) * qJD(1), t74 + (-t38 * t81 + t47 * t83) * qJD(1), t36 * qJ(4) + t40 * qJD(4) - t38 * t47 + (-t172 + (-t37 - t154) * t81) * t60, -t63, -t52, 0, t63, 0, 0, t74 + (-t34 + t160) * qJD(3) + ((qJ(5) * qJD(3) - t35) * t83 + t139 * t81) * qJD(1), -t33 * qJD(3) + (-t35 * t81 - t123) * qJD(1) + t157 (t29 + t33 + t76) * t150, t16 * qJ(4) - t190 * t29 + t20 * t84 - t21 * t33 - t28 * t35, -t167 * t45 + t164 (t24 + t101) * t82 + (t25 + t100) * t80 (-t132 + t169) * qJD(1) + t115, -t101 * t80 + t175 (t133 - t171) * qJD(1) - t114, t64 * t151, t177 + t25 * t79 - t64 * t7 + t190 * t43 + (-t75 * t167 - t23 * t80) * qJD(6) + (t5 * t81 + (t75 * t149 - t23 * t83) * t80) * qJD(1), -t178 - t24 * t79 + t64 * t8 + t190 * t45 + (t75 * t168 - t23 * t82) * qJD(6) + (-t23 * t161 + (t75 * t148 - t6) * t81) * qJD(1), t43 * t8 + t45 * t7 + (t5 * t150 - t25 * t75 - t1 + (t45 * t75 + t5) * qJD(6)) * t82 + (t6 * t150 - t24 * t75 + t2 + (t43 * t75 + t6) * qJD(6)) * t80, t16 * t79 + t182 * t75 + t190 * t23 - t5 * t7 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t61, -qJD(3) * t40 + t83 * t142 + t46, 0, 0, 0, 0, 0, 0, t61, -t63, 0, -qJD(1) * t123 + t143 + t157, 0, 0, 0, 0, 0, 0, -qJD(3) * t43 - t167 * t64 + t57, t64 ^ 2 * t80 - qJD(3) * t45 + t58, -t191 * t80 + t97 * t82, t187 * t80 + t188 * t82 - t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t127, 0.2e1 * t126, -t71 - t165, t65 + (t21 * t83 + t29 * t81 + t90) * qJD(1), 0, 0, 0, 0, 0, 0 (-t133 - t171) * qJD(1) + t114 (-t132 - t169) * qJD(1) + t115, t191 * t82 + t97 * t80 (-t110 * t83 - t23 * t81) * qJD(1) - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, -t43 ^ 2 + t45 ^ 2, -t191, -t170, t97, -t127, -t23 * t45 - t187, t23 * t43 - t188, 0, 0;];
tauc_reg  = t9;

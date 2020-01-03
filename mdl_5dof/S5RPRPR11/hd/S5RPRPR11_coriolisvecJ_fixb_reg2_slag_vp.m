% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR11_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:27:59
% DurationCPUTime: 1.38s
% Computational Cost: add. (2342->224), mult. (6206->276), div. (0->0), fcn. (4530->6), ass. (0->138)
t116 = sin(qJ(5));
t118 = cos(qJ(5));
t114 = sin(pkin(8));
t115 = cos(pkin(8));
t117 = sin(qJ(3));
t177 = cos(qJ(3));
t94 = t177 * t114 + t117 * t115;
t185 = t94 * qJD(1);
t148 = t177 * t115;
t137 = qJD(1) * t148;
t155 = t114 * t117;
t147 = qJD(1) * t155;
t83 = -t137 + t147;
t49 = t116 * t83 + t118 * t185;
t103 = qJD(3) * t137;
t146 = qJD(3) * t155;
t70 = qJD(1) * t146 - t103;
t88 = t94 * qJD(3);
t71 = qJD(1) * t88;
t12 = t49 * qJD(5) - t116 * t70 - t118 * t71;
t111 = qJD(3) - qJD(5);
t193 = t49 * t111;
t196 = t12 + t193;
t168 = pkin(6) + qJ(2);
t101 = t168 * t115;
t96 = qJD(1) * t101;
t81 = t117 * t96;
t100 = t168 * t114;
t95 = qJD(1) * t100;
t58 = -t177 * t95 - t81;
t154 = qJD(4) - t58;
t135 = t116 * t185 - t118 * t83;
t172 = t135 ^ 2;
t173 = t49 ^ 2;
t195 = t172 - t173;
t181 = pkin(3) + pkin(4);
t149 = t115 * pkin(2) + pkin(1);
t97 = -t149 * qJD(1) + qJD(2);
t192 = -t185 * qJ(4) + t97;
t23 = -t181 * t83 - t192;
t194 = t23 * t49;
t171 = t49 * t135;
t151 = qJD(5) * t118;
t152 = qJD(5) * t116;
t130 = t116 * t71 - t118 * t70 + t83 * t151 - t152 * t185;
t186 = t111 * t135;
t191 = t130 - t186;
t190 = -pkin(7) * t185 + t154;
t112 = qJD(3) * qJD(4);
t150 = qJD(1) * qJD(2);
t144 = t114 * t150;
t104 = qJD(2) * t148;
t145 = qJD(3) * t177;
t159 = qJD(1) * t104 - t95 * t145;
t31 = (-qJD(3) * t96 - t144) * t117 + t159;
t30 = t112 + t31;
t15 = t71 * pkin(7) + t30;
t127 = t94 * qJD(2);
t125 = qJD(1) * t127;
t59 = -t117 * t95 + t177 * t96;
t32 = t59 * qJD(3) + t125;
t20 = t70 * pkin(7) + t32;
t29 = -t181 * qJD(3) + t190;
t113 = qJD(3) * qJ(4);
t38 = t83 * pkin(7) + t59;
t33 = t113 + t38;
t6 = t116 * t29 + t118 * t33;
t2 = -qJD(5) * t6 - t116 * t15 + t118 * t20;
t189 = -t6 * t111 + t2;
t182 = t185 ^ 2;
t77 = t83 ^ 2;
t188 = -t77 - t182;
t187 = -t77 + t182;
t1 = t116 * t20 + t118 * t15 + t29 * t151 - t33 * t152;
t184 = t135 * t23 - t1;
t162 = t70 * qJ(4);
t134 = -qJD(4) * t185 + t162;
t14 = -t181 * t71 - t134;
t183 = (t58 + t81) * qJD(3) + t117 * t144 - t159;
t180 = t71 * pkin(3);
t99 = qJ(4) * t118 - t116 * t181;
t179 = t99 * qJD(5) + t190 * t116 + t118 * t38;
t98 = -qJ(4) * t116 - t118 * t181;
t178 = t98 * qJD(5) - t116 * t38 + t190 * t118;
t5 = -t116 * t33 + t118 * t29;
t176 = t111 * t5;
t60 = t177 * t100 + t101 * t117;
t175 = t32 * t60;
t40 = t83 * pkin(3) + t192;
t174 = t40 * t185;
t169 = t185 * t83;
t53 = t113 + t59;
t167 = t53 - t59;
t61 = -t117 * t100 + t177 * t101;
t161 = t83 * qJ(4);
t41 = t104 - t100 * t145 + (-qJD(2) * t114 - qJD(3) * t101) * t117;
t158 = qJD(3) * t41;
t42 = t61 * qJD(3) + t127;
t157 = qJD(3) * t42;
t156 = qJD(3) * t88;
t153 = t114 ^ 2 + t115 ^ 2;
t140 = t153 * qJD(1) ^ 2;
t139 = t111 ^ 2;
t93 = -t148 + t155;
t136 = t71 * t93 + t83 * t88;
t43 = -pkin(7) * t94 + t60;
t44 = pkin(7) * t93 + t61;
t9 = -t116 * t44 + t118 * t43;
t10 = t116 * t43 + t118 * t44;
t57 = t116 * t93 + t118 * t94;
t87 = -t115 * t145 + t146;
t133 = -t87 * qJ(4) + t94 * qJD(4);
t132 = 0.2e1 * t153 * t150;
t131 = t94 * qJ(4) + t149;
t126 = t185 * t42 + t32 * t94 - t41 * t83 - t60 * t70 - t61 * t71;
t124 = -t185 * t88 + t70 * t93 - t71 * t94 + t83 * t87;
t122 = 0.2e1 * t185 * qJD(3);
t121 = qJD(2) * t185;
t74 = t87 * qJD(3);
t56 = t116 * t94 - t118 * t93;
t55 = pkin(3) * t93 - t131;
t54 = pkin(3) * t185 + t161;
t52 = t103 + (t83 - t147) * qJD(3);
t51 = -t103 + (t83 + t147) * qJD(3);
t50 = -qJD(3) * pkin(3) + t154;
t36 = -t181 * t93 + t131;
t35 = pkin(3) * t88 - t133;
t34 = -t181 * t185 - t161;
t28 = t134 + t180;
t25 = t87 * pkin(7) + t42;
t24 = t88 * pkin(7) + t41;
t22 = -t185 * t87 - t70 * t94;
t21 = -t181 * t88 + t133;
t19 = t57 * qJD(5) - t116 * t87 - t118 * t88;
t18 = -t116 * t88 + t118 * t87 - t93 * t151 + t94 * t152;
t4 = -t10 * qJD(5) - t116 * t24 + t118 * t25;
t3 = t9 * qJD(5) + t116 * t25 + t118 * t24;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, qJ(2) * t132, t22, t124, -t74, t136, -t156, 0, -t149 * t71 + t88 * t97 - t157, t149 * t70 - t87 * t97 - t158, -t31 * t93 + t58 * t87 - t59 * t88 + t126, t31 * t61 + t41 * t59 - t42 * t58 + t175, t22, -t74, -t124, 0, t156, t136, t28 * t93 + t35 * t83 + t40 * t88 + t55 * t71 - t157, -t30 * t93 - t50 * t87 - t53 * t88 + t126, -t185 * t35 - t28 * t94 + t40 * t87 + t55 * t70 + t158, t28 * t55 + t30 * t61 + t35 * t40 + t41 * t53 + t42 * t50 + t175, t130 * t57 - t18 * t49, -t12 * t57 - t130 * t56 + t135 * t18 - t19 * t49, t18 * t111, t12 * t56 + t135 * t19, t19 * t111, 0, -t111 * t4 + t12 * t36 + t135 * t21 + t14 * t56 + t19 * t23, t111 * t3 + t130 * t36 + t14 * t57 - t18 * t23 + t21 * t49, -t1 * t56 - t10 * t12 - t130 * t9 - t135 * t3 + t18 * t5 - t19 * t6 - t2 * t57 - t4 * t49, t1 * t10 + t14 * t36 + t2 * t9 + t21 * t23 + t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -qJ(2) * t140, 0, 0, 0, 0, 0, 0, t122, -t51, t188, t185 * t58 + t59 * t83, 0, 0, 0, 0, 0, 0, t122, t188, t51, t180 + t162 + t53 * t83 + (-qJD(4) - t50) * t185, 0, 0, 0, 0, 0, 0, -t12 + t193, -t130 - t186, t172 + t173, -t135 * t6 - t49 * t5 - t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t187, t52, -t169, 0, 0, -t185 * t97 - t121, t97 * t83 + t183, 0, 0, t169, t52, -t187, 0, 0, -t169, -t54 * t83 - t121 - t174, pkin(3) * t70 - qJ(4) * t71 + t167 * t185 + (t50 - t154) * t83, t185 * t54 - t40 * t83 + 0.2e1 * t112 - t183, -t32 * pkin(3) + t30 * qJ(4) + t154 * t53 - t40 * t54 - t50 * t59, -t171, t195, -t191, t171, t196, 0, t111 * t179 - t135 * t34 + t194 - t2, t111 * t178 - t34 * t49 - t184, -t130 * t98 - t12 * t99 + (-t178 + t5) * t135 + (t179 - t6) * t49, t1 * t99 + t178 * t6 - t179 * t5 + t2 * t98 - t23 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t52, -qJD(3) ^ 2 - t182, -t167 * qJD(3) + t125 + t174, 0, 0, 0, 0, 0, 0, -t116 * t139 - t135 * t185, -t118 * t139 - t185 * t49, -t196 * t116 - t191 * t118, -t23 * t185 + t189 * t118 + (t1 + t176) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t195, t191, -t171, -t196, 0, t189 - t194, -t176 + t184, 0, 0;];
tauc_reg = t7;

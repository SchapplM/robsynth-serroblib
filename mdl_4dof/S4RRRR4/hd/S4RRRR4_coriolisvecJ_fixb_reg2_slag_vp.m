% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:19
% DurationCPUTime: 1.43s
% Computational Cost: add. (2726->248), mult. (7005->357), div. (0->0), fcn. (4686->6), ass. (0->141)
t160 = cos(qJ(3));
t91 = cos(qJ(2));
t118 = t160 * t91;
t88 = sin(qJ(3));
t89 = sin(qJ(2));
t137 = t88 * t89;
t97 = t118 - t137;
t166 = qJD(1) * t97;
t130 = qJD(2) * pkin(2);
t120 = t89 * t130;
t165 = 0.2e1 * t120;
t125 = qJD(1) * qJD(2);
t164 = -0.2e1 * t125;
t129 = qJD(1) * t89;
t59 = -qJD(1) * t118 + t88 * t129;
t136 = t88 * t91;
t68 = t160 * t89 + t136;
t61 = qJD(1) * t68;
t84 = -pkin(2) * t91 - pkin(1);
t73 = qJD(1) * t84;
t36 = pkin(3) * t59 - pkin(7) * t61 + t73;
t124 = qJD(2) + qJD(3);
t161 = -pkin(6) - pkin(5);
t75 = t161 * t91;
t72 = qJD(1) * t75;
t119 = t160 * t72;
t74 = t161 * t89;
t70 = qJD(1) * t74;
t64 = t70 + t130;
t46 = t88 * t64 - t119;
t39 = t124 * pkin(7) + t46;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t13 = t36 * t90 - t39 * t87;
t14 = t36 * t87 + t39 * t90;
t163 = -t13 * t87 + t14 * t90;
t162 = t124 * qJD(1);
t41 = t97 * t162;
t53 = t87 * t124 + t90 * t61;
t26 = qJD(4) * t53 + t87 * t41;
t116 = t89 * t125;
t50 = t124 * t68;
t42 = t50 * qJD(1);
t16 = pkin(2) * t116 + pkin(3) * t42 - pkin(7) * t41;
t114 = qJD(3) * t160;
t128 = qJD(3) * t88;
t117 = qJD(2) * t161;
t105 = qJD(1) * t117;
t65 = t89 * t105;
t66 = t91 * t105;
t22 = t64 * t114 + t72 * t128 + t160 * t65 + t88 * t66;
t103 = t16 * t87 + t22 * t90;
t3 = qJD(4) * t13 + t103;
t2 = t3 * t90;
t56 = t160 * t66;
t115 = t88 * t65 - t56;
t23 = t46 * qJD(3) + t115;
t98 = t160 * t74 + t88 * t75;
t157 = t23 * t98;
t111 = t90 * t124;
t127 = qJD(4) * t87;
t25 = -qJD(4) * t111 + t61 * t127 - t90 * t41;
t156 = t25 * t87;
t155 = t26 * t90;
t154 = t42 * t97;
t153 = t42 * t87;
t152 = t42 * t90;
t51 = t61 * t87 - t111;
t57 = qJD(4) + t59;
t151 = t51 * t57;
t150 = t51 * t87;
t149 = t53 * t51;
t148 = t53 * t57;
t147 = t53 * t90;
t146 = t57 * t61;
t145 = t59 * t87;
t144 = t59 * t90;
t143 = t61 * t59;
t142 = t68 * t87;
t141 = t68 * t90;
t140 = t73 * t61;
t138 = t88 * t72;
t93 = qJD(1) ^ 2;
t135 = t91 * t93;
t92 = qJD(2) ^ 2;
t134 = t92 * t89;
t133 = t92 * t91;
t132 = t89 ^ 2 - t91 ^ 2;
t131 = pkin(2) * qJD(3);
t126 = qJD(4) * t90;
t123 = -t13 * t144 - t14 * t145 + t2;
t122 = t89 * t135;
t121 = pkin(2) * t129;
t113 = t57 * t90;
t112 = pkin(1) * t164;
t110 = pkin(2) * t114;
t109 = t91 * t116;
t45 = t160 * t64 + t138;
t38 = -t124 * pkin(3) - t45;
t108 = t38 * t126 + t14 * t61 + t23 * t87;
t107 = qJD(2) * t118;
t47 = t70 * t88 - t119;
t106 = pkin(2) * t128 - t47;
t43 = pkin(3) * t61 + pkin(7) * t59;
t104 = t13 * t90 + t14 * t87;
t82 = pkin(2) * t88 + pkin(7);
t102 = t38 * t59 - t42 * t82;
t44 = -pkin(3) * t97 - pkin(7) * t68 + t84;
t55 = -t160 * t75 + t88 * t74;
t27 = t44 * t90 - t55 * t87;
t28 = t44 * t87 + t55 * t90;
t101 = t38 * t127 - t13 * t61 - t23 * t90;
t49 = -t91 * t114 + t124 * t137 - t107;
t100 = t68 * t126 - t49 * t87;
t99 = -t68 * t127 - t49 * t90;
t4 = -qJD(4) * t14 + t90 * t16 - t22 * t87;
t96 = -t104 * qJD(4) - t4 * t87;
t95 = t96 + t2;
t94 = t73 * t59 - t22;
t83 = -t160 * pkin(2) - pkin(3);
t71 = t89 * t117;
t48 = t160 * t70 + t138;
t37 = t43 + t121;
t33 = -t59 ^ 2 + t61 ^ 2;
t32 = t61 * t124 - t68 * t162;
t31 = (t59 + t166) * t124;
t30 = t55 * qJD(3) - t161 * t107 + t88 * t71;
t29 = t98 * qJD(3) + t117 * t136 + t160 * t71;
t24 = pkin(3) * t50 + pkin(7) * t49 + t120;
t21 = t43 * t87 + t45 * t90;
t20 = t43 * t90 - t45 * t87;
t18 = t37 * t87 + t48 * t90;
t17 = t37 * t90 - t48 * t87;
t10 = t57 * t113 - t53 * t61 + t153;
t9 = -t57 ^ 2 * t87 + t51 * t61 + t152;
t8 = t57 * t150 - t155;
t7 = t53 * t113 - t156;
t6 = -qJD(4) * t28 + t24 * t90 - t29 * t87;
t5 = qJD(4) * t27 + t24 * t87 + t29 * t90;
t1 = (-t25 - t151) * t90 + (-t26 - t148) * t87;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t109, t132 * t164, t133, -0.2e1 * t109, -t134, 0, -pkin(5) * t133 + t89 * t112, pkin(5) * t134 + t91 * t112, 0, 0, t41 * t68 - t49 * t61, t41 * t97 - t42 * t68 + t49 * t59 - t50 * t61, -t49 * t124, t50 * t59 - t154, -t50 * t124, 0, t84 * t42 + t73 * t50 - t30 * t124 + (t59 - t166) * t120, -t29 * t124 + t61 * t165 + t84 * t41 - t73 * t49, t22 * t97 + t23 * t68 - t29 * t59 + t30 * t61 - t41 * t98 - t42 * t55 + t45 * t49 - t46 * t50, t73 * t165 + t22 * t55 + t29 * t46 - t30 * t45 - t157, -t25 * t141 + t53 * t99, (t51 * t90 + t53 * t87) * t49 + (t156 - t155 + (-t147 + t150) * qJD(4)) * t68, t42 * t141 + t25 * t97 + t50 * t53 + t57 * t99, t100 * t51 + t26 * t142, -t100 * t57 - t42 * t142 + t26 * t97 - t50 * t51, t50 * t57 - t154, t100 * t38 + t13 * t50 + t23 * t142 - t26 * t98 + t27 * t42 + t30 * t51 - t4 * t97 + t57 * t6, -t14 * t50 + t23 * t141 + t25 * t98 - t28 * t42 + t3 * t97 + t30 * t53 + t38 * t99 - t5 * t57, t25 * t27 - t26 * t28 - t5 * t51 - t53 * t6 + t104 * t49 + (-t163 * qJD(4) - t3 * t87 - t4 * t90) * t68, t13 * t6 + t14 * t5 + t27 * t4 + t28 * t3 + t30 * t38 - t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t132 * t93, 0, t122, 0, 0, t93 * pkin(1) * t89, pkin(1) * t135, 0, 0, t143, t33, t31, -t143, t32, 0, t72 * t114 + t56 - t59 * t121 - t140 + t47 * t124 + (-qJD(3) * t64 - t124 * t131 - t65) * t88, t48 * t124 + (-t124 * t114 - t61 * t129) * pkin(2) + t94, (t46 - t47) * t61 + (-t45 + t48) * t59 + (-t160 * t41 - t42 * t88 + (-t160 * t59 + t61 * t88) * qJD(3)) * pkin(2), t45 * t47 - t46 * t48 + (-t73 * t129 - t160 * t23 + t22 * t88 + (t160 * t46 - t45 * t88) * qJD(3)) * pkin(2), t7, t1, t10, t8, t9, -t146, t83 * t26 + t102 * t87 + t106 * t51 + (-t87 * t110 - t82 * t126 - t17) * t57 + t101, -t83 * t25 + t102 * t90 + t106 * t53 + (-t90 * t110 + t82 * t127 + t18) * t57 + t108, t17 * t53 + t18 * t51 + (-t51 * t110 - t26 * t82 + (t53 * t82 - t13) * qJD(4)) * t90 + (t53 * t110 - t25 * t82 - t4 + (t51 * t82 - t14) * qJD(4)) * t87 + t123, -t13 * t17 - t14 * t18 + t23 * t83 - t38 * t47 + (t163 * t160 + t38 * t88) * t131 + t95 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t33, t31, -t143, t32, 0, t46 * qJD(2) - t115 - t140, t45 * t124 + t94, 0, 0, t7, t1, t10, t8, t9, -t146, t38 * t145 - pkin(3) * t26 - t20 * t57 - t46 * t51 + (-t57 * t126 - t153) * pkin(7) + t101, t38 * t144 + pkin(3) * t25 + t21 * t57 - t46 * t53 + (t57 * t127 - t152) * pkin(7) + t108, t20 * t53 + t21 * t51 + (-t156 - t155 + (t147 + t150) * qJD(4)) * pkin(7) + t96 + t123, -pkin(3) * t23 + pkin(7) * t95 - t13 * t20 - t14 * t21 - t38 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t51 ^ 2 + t53 ^ 2, -t25 + t151, -t149, t148 - t26, t42, t14 * t57 - t38 * t53 + t4, t38 * t51 - t103 + (-qJD(4) + t57) * t13, 0, 0;];
tauc_reg = t11;

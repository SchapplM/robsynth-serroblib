% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:32
% EndTime: 2019-12-05 17:31:38
% DurationCPUTime: 1.33s
% Computational Cost: add. (1938->200), mult. (5250->336), div. (0->0), fcn. (4017->8), ass. (0->124)
t94 = cos(pkin(9));
t96 = cos(pkin(7));
t143 = t96 * t94;
t93 = sin(pkin(7));
t95 = cos(pkin(8));
t144 = t93 * t95;
t91 = sin(pkin(9));
t63 = t91 * t144 + t143;
t51 = t63 * qJD(1);
t154 = qJD(5) + t51;
t103 = (-qJD(4) * t95 + qJD(2)) * t93;
t100 = qJD(1) * t103;
t132 = qJD(3) * t93;
t92 = sin(pkin(8));
t118 = t92 * t132;
t101 = -t96 * qJD(4) - t118;
t133 = qJD(2) * t96;
t81 = t95 * t133;
t74 = qJD(1) * t81;
t41 = t101 * qJD(1) + t74;
t18 = t91 * t100 + t94 * t41;
t65 = t95 * t132 + t92 * t133;
t57 = t65 * qJD(1);
t135 = qJD(1) * t93;
t123 = t92 * t135;
t134 = qJD(1) * t96;
t129 = qJ(2) * qJD(1);
t116 = t96 * t129;
t71 = -t96 * pkin(2) - t93 * qJ(3) - pkin(1);
t62 = t71 * qJD(1) + qJD(2);
t35 = t95 * t116 + t92 * t62;
t23 = -qJ(4) * t134 + t35;
t107 = pkin(3) * t92 - qJ(4) * t95;
t78 = t93 * t129 + qJD(3);
t40 = t107 * t135 + t78;
t13 = t94 * t23 + t91 * t40;
t11 = pkin(6) * t123 + t13;
t34 = -t92 * t116 + t95 * t62;
t22 = pkin(3) * t134 + qJD(4) - t34;
t117 = t94 * t135;
t119 = t91 * t134;
t54 = t95 * t117 - t119;
t9 = t51 * pkin(4) - t54 * pkin(6) + t22;
t97 = sin(qJ(5));
t98 = cos(qJ(5));
t6 = t98 * t11 + t97 * t9;
t2 = -t6 * qJD(5) - t97 * t18 + t98 * t57;
t157 = t154 * t6 + t2;
t108 = t97 * t11 - t98 * t9;
t1 = -t108 * qJD(5) + t98 * t18 + t97 * t57;
t156 = t108 * t154 + t1;
t136 = qJD(1) * t92;
t121 = t98 * t136;
t110 = t93 * t121;
t131 = qJD(5) * t97;
t24 = -qJD(5) * t110 + t54 * t131;
t26 = t97 * t54 - t110;
t155 = t154 * t26 - t24;
t17 = -t94 * t100 + t91 * t41;
t153 = t17 * t94;
t122 = t97 * t136;
t29 = t93 * t122 + t98 * t54;
t152 = t29 * t26;
t151 = t57 * t95;
t150 = t78 * t93;
t89 = t93 ^ 2;
t99 = qJD(1) ^ 2;
t149 = t89 * t99;
t148 = t91 * t92;
t147 = t92 * t93;
t146 = t92 * t97;
t145 = t92 * t98;
t56 = (t95 * t143 + t91 * t93) * qJD(1);
t68 = t94 * t145 - t97 * t95;
t142 = t68 * qJD(5) + t96 * t121 - t97 * t56;
t102 = t94 * t146 + t98 * t95;
t141 = t102 * qJD(5) + t96 * t122 + t98 * t56;
t137 = qJ(2) * t96;
t139 = t95 * t137 + t92 * t71;
t39 = -t96 * qJ(4) + t139;
t48 = (qJ(2) + t107) * t93;
t140 = t94 * t39 + t91 * t48;
t90 = t96 ^ 2;
t138 = t89 + t90;
t128 = qJD(1) * qJD(2);
t88 = t92 ^ 2;
t127 = t88 * t149;
t126 = t93 * t145;
t125 = t93 * t96 * t99;
t124 = 0.2e1 * qJD(2) * t89;
t120 = t95 * t135;
t115 = t138 * t99;
t114 = -t92 * t137 + t95 * t71;
t113 = qJ(2) * t128;
t112 = t154 ^ 2;
t111 = t91 * t123;
t109 = t96 * pkin(3) - t114;
t64 = t94 * t144 - t96 * t91;
t14 = t63 * pkin(4) - t64 * pkin(6) + t109;
t16 = pkin(6) * t147 + t140;
t7 = t98 * t14 - t97 * t16;
t8 = t97 * t14 + t98 * t16;
t12 = -t91 * t23 + t94 * t40;
t105 = -t91 * t39 + t94 * t48;
t58 = -qJD(1) * t118 + t74;
t104 = t58 * t92 - t151;
t37 = t93 * t146 + t64 * t98;
t84 = t89 * t113;
t66 = t81 - t118;
t53 = t95 * t119 - t117;
t46 = t101 + t81;
t43 = t68 * t135;
t42 = t102 * t135;
t36 = t64 * t97 - t126;
t32 = t37 * qJD(5);
t31 = -qJD(5) * t126 + t64 * t131;
t25 = t29 * qJD(5);
t20 = t91 * t103 + t94 * t46;
t19 = -t94 * t103 + t91 * t46;
t15 = -pkin(4) * t147 - t105;
t10 = -pkin(4) * t123 - t12;
t4 = -t8 * qJD(5) - t97 * t20 + t98 * t65;
t3 = t7 * qJD(5) + t98 * t20 + t97 * t65;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t138 * t128, 0.2e1 * t90 * t113 + 0.2e1 * t84, 0, 0, 0, 0, 0, 0, t57 * t96 + (t92 * t124 + t65 * t96) * qJD(1), t58 * t96 + (t95 * t124 + t66 * t96) * qJD(1), ((t65 * t95 - t66 * t92) * qJD(1) - t104) * t93, qJD(2) * t150 - t57 * t114 + t58 * t139 - t34 * t65 + t35 * t66 + t84, 0, 0, 0, 0, 0, 0, t65 * t51 + t57 * t63 + (-qJD(1) * t19 - t17) * t147, t65 * t54 + t57 * t64 + (-qJD(1) * t20 - t18) * t147, t17 * t64 - t18 * t63 + t19 * t54 - t20 * t51, -t17 * t105 + t57 * t109 - t12 * t19 + t13 * t20 + t18 * t140 + t22 * t65, -t24 * t37 - t29 * t31, t24 * t36 - t37 * t25 + t31 * t26 - t29 * t32, -t154 * t31 - t24 * t63, t25 * t36 + t26 * t32, -t154 * t32 - t25 * t63, 0, t10 * t32 + t15 * t25 + t154 * t4 + t17 * t36 + t19 * t26 + t2 * t63, -t1 * t63 - t10 * t31 - t15 * t24 - t154 * t3 + t17 * t37 + t19 * t29, -t1 * t36 - t108 * t31 - t2 * t37 + t7 * t24 - t8 * t25 - t3 * t26 - t4 * t29 - t6 * t32, t1 * t8 + t10 * t19 - t108 * t4 + t17 * t15 + t2 * t7 + t6 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -qJ(2) * t115, 0, 0, 0, 0, 0, 0, -t92 * t115, -t95 * t115, 0, (-t150 + (t34 * t92 - t35 * t95) * t96) * qJD(1) + t104, 0, 0, 0, 0, 0, 0, (-t51 * t96 + t53 * t93) * t136, (-t54 * t96 + t56 * t93) * t136, t56 * t51 - t53 * t54, t12 * t53 - t13 * t56 - t151 + (-t22 * t134 + t17 * t91 + t18 * t94) * t92, 0, 0, 0, 0, 0, 0, -t142 * t154 + t25 * t148 - t53 * t26, t141 * t154 - t24 * t148 - t53 * t29, -t102 * t24 + t141 * t26 + t142 * t29 - t68 * t25, t1 * t68 - t10 * t53 - t102 * t2 + t108 * t142 - t141 * t6 + t17 * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95 * t125, t92 * t125, (-t95 ^ 2 - t88) * t149, (t34 * t95 + t35 * t92 + qJD(2)) * t135, 0, 0, 0, 0, 0, 0, -t51 * t120 - t91 * t127, -t54 * t120 - t94 * t127, (-t51 * t94 + t54 * t91) * t123, -t153 + t18 * t91 + (-t22 * t95 + (-t12 * t91 + t13 * t94) * t92) * t135, 0, 0, 0, 0, 0, 0, t26 * t111 - t94 * t25 + (-qJD(5) * t91 * t98 - t42) * t154, t29 * t111 + t94 * t24 + (t131 * t91 - t43) * t154, -t43 * t26 + t42 * t29 + (-t24 * t97 - t25 * t98 + (t26 * t97 + t29 * t98) * qJD(5)) * t91, -t153 + t108 * t42 + t6 * t43 + (t10 * t123 + t1 * t98 - t2 * t97 + (t108 * t98 - t6 * t97) * qJD(5)) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t123, -t51 * t123, -t51 ^ 2 - t54 ^ 2, t12 * t54 + t13 * t51 + t57, 0, 0, 0, 0, 0, 0, -t112 * t97 - t54 * t26, -t112 * t98 - t54 * t29, -t155 * t98 + (t154 * t29 - t25) * t97, -t10 * t54 + t156 * t97 + t157 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t26 ^ 2 + t29 ^ 2, t155, -t152, (-qJD(5) + t154) * t29, 0, -t10 * t29 + t157, t10 * t26 - t156, 0, 0;];
tauc_reg = t5;

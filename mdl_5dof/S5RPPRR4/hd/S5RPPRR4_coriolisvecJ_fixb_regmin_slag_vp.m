% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:23
% EndTime: 2020-01-03 11:31:30
% DurationCPUTime: 1.24s
% Computational Cost: add. (1328->177), mult. (3933->281), div. (0->0), fcn. (3046->8), ass. (0->119)
t96 = sin(qJ(5));
t134 = qJD(5) * t96;
t92 = sin(pkin(9));
t94 = cos(pkin(9));
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t73 = t99 * t92 + t97 * t94;
t104 = qJD(1) * t73;
t93 = sin(pkin(8));
t48 = t93 * t104;
t98 = cos(qJ(5));
t146 = t98 * t48;
t69 = t73 * qJD(4);
t53 = t93 * t69;
t43 = qJD(1) * t53;
t139 = qJD(1) * t93;
t126 = t92 * t139;
t145 = t99 * t94;
t129 = t93 * t145;
t51 = qJD(1) * t129 - t97 * t126;
t44 = t51 * qJD(4);
t10 = -qJD(5) * t146 - t51 * t134 - t98 * t43 - t96 * t44;
t17 = t96 * t51 + t146;
t95 = cos(pkin(8));
t159 = -t95 * qJD(1) + qJD(4);
t81 = -qJD(5) - t159;
t150 = t17 * t81;
t162 = t10 - t150;
t109 = -t96 * t48 + t98 * t51;
t161 = t109 * t17;
t11 = t109 * qJD(5) - t96 * t43 + t98 * t44;
t151 = t109 * t81;
t160 = -t11 - t151;
t158 = t109 ^ 2 - t17 ^ 2;
t131 = qJ(2) * qJD(1);
t82 = t93 * t131 + qJD(3);
t65 = pkin(3) * t126 + t82;
t26 = t48 * pkin(4) + t65;
t103 = -t94 * t93 * pkin(6) + (-qJ(2) * t92 - pkin(3)) * t95;
t76 = -t95 * pkin(2) - t93 * qJ(3) - pkin(1);
t64 = t76 * qJD(1) + qJD(2);
t56 = t94 * t64;
t23 = t103 * qJD(1) + t56;
t123 = t95 * t131;
t33 = t94 * t123 + t92 * t64;
t28 = -pkin(6) * t126 + t33;
t111 = -t97 * t23 - t99 * t28;
t137 = qJD(3) * t93;
t138 = qJD(2) * t95;
t66 = -t94 * t137 - t92 * t138;
t60 = t66 * qJD(1);
t67 = -t92 * t137 + t94 * t138;
t61 = t67 * qJD(1);
t120 = t99 * t60 - t97 * t61;
t102 = t111 * qJD(4) + t120;
t3 = t43 * pkin(7) + t102;
t9 = -t48 * pkin(7) - t111;
t7 = t9 * t134;
t157 = t26 * t17 + t7 + (t9 * t81 - t3) * t96;
t155 = qJD(5) + t81;
t135 = qJD(4) * t99;
t136 = qJD(4) * t97;
t106 = -t23 * t135 + t28 * t136 - t97 * t60 - t99 * t61;
t2 = -t44 * pkin(7) - t106;
t124 = -t96 * t2 + t98 * t3;
t154 = -t26 * t109 + t124;
t153 = pkin(4) * t51;
t152 = t98 * t9;
t149 = t48 * t159;
t148 = t51 * t159;
t147 = t92 * t93;
t144 = -t95 * t104 + t69;
t72 = -t97 * t92 + t145;
t143 = t159 * t72;
t140 = qJ(2) * t95;
t142 = t94 * t140 + t92 * t76;
t74 = pkin(3) * t147 + t93 * qJ(2);
t90 = t93 ^ 2;
t91 = t95 ^ 2;
t141 = t90 + t91;
t133 = t93 * qJD(2);
t130 = qJD(1) * qJD(2);
t128 = 0.2e1 * qJD(2) * t90;
t100 = qJD(1) ^ 2;
t127 = t100 * t93 * t95;
t122 = t99 * t23 - t97 * t28;
t8 = -t51 * pkin(7) + t122;
t6 = pkin(4) * t159 + t8;
t125 = pkin(4) * t81 - t6;
t119 = t99 * t66 - t97 * t67;
t118 = t141 * t100;
t117 = qJ(2) * t130;
t114 = -t96 * t6 - t152;
t113 = qJD(5) * t73 + t144;
t112 = qJD(5) * t72 + t143;
t71 = t94 * t76;
t29 = t71 + t103;
t34 = -pkin(6) * t147 + t142;
t110 = -t97 * t29 - t99 * t34;
t108 = -t60 * t94 - t61 * t92;
t62 = t73 * t93;
t63 = t72 * t93;
t24 = t98 * t62 + t96 * t63;
t25 = -t96 * t62 + t98 * t63;
t105 = t29 * t135 - t34 * t136 + t97 * t66 + t99 * t67;
t87 = t93 * t130;
t84 = t90 * t117;
t54 = qJD(4) * t129 - t136 * t147;
t35 = t54 * pkin(4) + t133;
t32 = -t92 * t123 + t56;
t31 = t44 * pkin(4) + t87;
t30 = t62 * pkin(4) + t74;
t15 = -t62 * pkin(7) - t110;
t14 = -t95 * pkin(4) - t63 * pkin(7) + t99 * t29 - t97 * t34;
t13 = t25 * qJD(5) - t96 * t53 + t98 * t54;
t12 = -t24 * qJD(5) - t98 * t53 - t96 * t54;
t5 = t53 * pkin(7) + t110 * qJD(4) + t119;
t4 = -t54 * pkin(7) + t105;
t1 = [0, 0, 0, 0, 0, 0.2e1 * t141 * t130, 0.2e1 * t91 * t117 + 0.2e1 * t84, -t60 * t95 + (t92 * t128 - t66 * t95) * qJD(1), t61 * t95 + (t94 * t128 + t67 * t95) * qJD(1), ((-t66 * t94 - t67 * t92) * qJD(1) + t108) * t93, t61 * t142 + t33 * t67 + t60 * (-t92 * t140 + t71) + t32 * t66 + t84 + t82 * t133, -t43 * t63 - t51 * t53, t43 * t62 - t63 * t44 + t53 * t48 - t51 * t54, -t159 * t53 + t43 * t95, -t159 * t54 + t44 * t95, 0, t119 * t159 - t120 * t95 + t74 * t44 + t65 * t54 + (t110 * t159 - t111 * t95) * qJD(4) + (qJD(1) * t62 + t48) * t133, -t105 * t159 - t106 * t95 - t74 * t43 - t65 * t53 + (qJD(1) * t63 + t51) * t133, t10 * t25 + t109 * t12, -t10 * t24 - t109 * t13 - t25 * t11 - t12 * t17, -t10 * t95 - t12 * t81, t11 * t95 + t13 * t81, 0, -(-t96 * t4 + t98 * t5) * t81 - t124 * t95 + t35 * t17 + t30 * t11 + t31 * t24 + t26 * t13 + (-(-t14 * t96 - t15 * t98) * t81 - t114 * t95) * qJD(5), t30 * t10 + t26 * t12 + t35 * t109 + t31 * t25 - t7 * t95 + ((-qJD(5) * t15 + t5) * t81 + t3 * t95) * t96 + ((qJD(5) * t14 + t4) * t81 + (qJD(5) * t6 + t2) * t95) * t98; 0, 0, 0, 0, 0, -t118, -qJ(2) * t118, -t92 * t118, -t94 * t118, 0, (-t82 * t93 + (t32 * t92 - t33 * t94) * t95) * qJD(1) - t108, 0, 0, 0, 0, 0, -t48 * t139 - t144 * t159, -t51 * t139 - t143 * t159, 0, 0, 0, 0, 0, -t17 * t139 + (t112 * t96 + t113 * t98) * t81, -t109 * t139 + (t112 * t98 - t113 * t96) * t81; 0, 0, 0, 0, 0, 0, 0, -t94 * t127, t92 * t127, (-t92 ^ 2 - t94 ^ 2) * t90 * t100, t87 + (t32 * t94 + t33 * t92) * t139, 0, 0, 0, 0, 0, t44 + t148, -t43 - t149, 0, 0, 0, 0, 0, t11 - t151, t10 + t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t48, -t48 ^ 2 + t51 ^ 2, -t43 + t149, -t44 + t148, 0, -t111 * t159 - t65 * t51 + t102, t122 * t159 + t65 * t48 + t106, t161, t158, t162, t160, 0, (-t96 * t8 - t152) * t81 - t17 * t153 + (t125 * t96 - t152) * qJD(5) + t154, -t109 * t153 + (t125 * qJD(5) - t8 * t81 - t2) * t98 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t158, t162, t160, 0, t155 * t114 + t154, (-t155 * t6 - t2) * t98 + t157;];
tauc_reg = t1;

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
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:16:58
% EndTime: 2022-01-23 09:17:02
% DurationCPUTime: 1.07s
% Computational Cost: add. (1320->174), mult. (3893->273), div. (0->0), fcn. (3018->8), ass. (0->118)
t96 = sin(qJ(5));
t133 = qJD(5) * t96;
t92 = sin(pkin(9));
t94 = cos(pkin(9));
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t73 = t99 * t92 + t97 * t94;
t104 = qJD(1) * t73;
t93 = sin(pkin(8));
t48 = t93 * t104;
t98 = cos(qJ(5));
t145 = t98 * t48;
t69 = t73 * qJD(4);
t53 = t93 * t69;
t43 = qJD(1) * t53;
t138 = qJD(1) * t93;
t125 = t92 * t138;
t144 = t99 * t94;
t128 = t93 * t144;
t51 = qJD(1) * t128 - t97 * t125;
t44 = t51 * qJD(4);
t10 = -qJD(5) * t145 - t51 * t133 - t98 * t43 - t96 * t44;
t17 = t96 * t51 + t145;
t95 = cos(pkin(8));
t158 = -t95 * qJD(1) + qJD(4);
t81 = -qJD(5) - t158;
t149 = t17 * t81;
t161 = t10 - t149;
t108 = -t96 * t48 + t98 * t51;
t160 = t108 * t17;
t11 = t108 * qJD(5) - t96 * t43 + t98 * t44;
t150 = t108 * t81;
t159 = -t11 - t150;
t157 = t108 ^ 2 - t17 ^ 2;
t130 = qJ(2) * qJD(1);
t82 = t93 * t130 + qJD(3);
t65 = pkin(3) * t125 + t82;
t26 = t48 * pkin(4) + t65;
t103 = -t94 * t93 * pkin(6) + (-qJ(2) * t92 - pkin(3)) * t95;
t76 = -pkin(2) * t95 - t93 * qJ(3) - pkin(1);
t64 = t76 * qJD(1) + qJD(2);
t56 = t94 * t64;
t23 = t103 * qJD(1) + t56;
t122 = t95 * t130;
t33 = t94 * t122 + t92 * t64;
t28 = -pkin(6) * t125 + t33;
t110 = -t97 * t23 - t99 * t28;
t136 = qJD(3) * t93;
t137 = qJD(2) * t95;
t66 = -t94 * t136 - t92 * t137;
t60 = t66 * qJD(1);
t67 = -t92 * t136 + t94 * t137;
t61 = t67 * qJD(1);
t119 = t99 * t60 - t97 * t61;
t102 = t110 * qJD(4) + t119;
t3 = t43 * pkin(7) + t102;
t9 = -t48 * pkin(7) - t110;
t7 = t9 * t133;
t156 = t26 * t17 + t7 + (t9 * t81 - t3) * t96;
t154 = qJD(5) + t81;
t134 = qJD(4) * t99;
t135 = qJD(4) * t97;
t106 = -t23 * t134 + t28 * t135 - t97 * t60 - t99 * t61;
t2 = -t44 * pkin(7) - t106;
t123 = -t96 * t2 + t98 * t3;
t153 = -t26 * t108 + t123;
t152 = pkin(4) * t51;
t151 = t98 * t9;
t148 = t48 * t158;
t147 = t51 * t158;
t146 = t92 * t93;
t143 = -t95 * t104 + t69;
t72 = -t97 * t92 + t144;
t142 = t158 * t72;
t139 = qJ(2) * t95;
t141 = t94 * t139 + t92 * t76;
t74 = pkin(3) * t146 + t93 * qJ(2);
t90 = t93 ^ 2;
t91 = t95 ^ 2;
t140 = t90 + t91;
t132 = t93 * qJD(2);
t129 = qJD(1) * qJD(2);
t127 = 0.2e1 * qJD(2) * t90;
t100 = qJD(1) ^ 2;
t126 = t100 * t93 * t95;
t121 = t99 * t23 - t97 * t28;
t8 = -t51 * pkin(7) + t121;
t6 = pkin(4) * t158 + t8;
t124 = pkin(4) * t81 - t6;
t118 = t99 * t66 - t97 * t67;
t117 = t140 * t100;
t116 = qJ(2) * t129;
t113 = -t96 * t6 - t151;
t112 = qJD(5) * t73 + t143;
t111 = qJD(5) * t72 + t142;
t71 = t94 * t76;
t29 = t71 + t103;
t34 = -pkin(6) * t146 + t141;
t109 = -t97 * t29 - t99 * t34;
t62 = t73 * t93;
t63 = t72 * t93;
t24 = t98 * t62 + t96 * t63;
t25 = -t96 * t62 + t98 * t63;
t105 = t29 * t134 - t34 * t135 + t97 * t66 + t99 * t67;
t87 = t93 * t129;
t84 = t90 * t116;
t54 = qJD(4) * t128 - t135 * t146;
t35 = t54 * pkin(4) + t132;
t32 = -t92 * t122 + t56;
t31 = t44 * pkin(4) + t87;
t30 = t62 * pkin(4) + t74;
t15 = -t62 * pkin(7) - t109;
t14 = -t95 * pkin(4) - t63 * pkin(7) + t99 * t29 - t97 * t34;
t13 = t25 * qJD(5) - t96 * t53 + t98 * t54;
t12 = -t24 * qJD(5) - t98 * t53 - t96 * t54;
t5 = t53 * pkin(7) + t109 * qJD(4) + t118;
t4 = -t54 * pkin(7) + t105;
t1 = [0, 0, 0, 0, 0.2e1 * t140 * t129, 0.2e1 * t91 * t116 + 0.2e1 * t84, -t60 * t95 + (t92 * t127 - t66 * t95) * qJD(1), t61 * t95 + (t94 * t127 + t67 * t95) * qJD(1), t61 * t141 + t33 * t67 + t60 * (-t92 * t139 + t71) + t32 * t66 + t84 + t82 * t132, -t43 * t63 - t51 * t53, t43 * t62 - t63 * t44 + t53 * t48 - t51 * t54, -t158 * t53 + t43 * t95, -t158 * t54 + t44 * t95, 0, t118 * t158 - t119 * t95 + t74 * t44 + t65 * t54 + (t109 * t158 - t110 * t95) * qJD(4) + (qJD(1) * t62 + t48) * t132, -t105 * t158 - t106 * t95 - t74 * t43 - t65 * t53 + (qJD(1) * t63 + t51) * t132, t10 * t25 + t108 * t12, -t10 * t24 - t108 * t13 - t25 * t11 - t12 * t17, -t10 * t95 - t12 * t81, t11 * t95 + t13 * t81, 0, -(-t96 * t4 + t98 * t5) * t81 - t123 * t95 + t35 * t17 + t30 * t11 + t31 * t24 + t26 * t13 + (-(-t14 * t96 - t15 * t98) * t81 - t113 * t95) * qJD(5), t30 * t10 + t26 * t12 + t35 * t108 + t31 * t25 - t7 * t95 + ((-qJD(5) * t15 + t5) * t81 + t3 * t95) * t96 + ((qJD(5) * t14 + t4) * t81 + (qJD(5) * t6 + t2) * t95) * t98; 0, 0, 0, 0, -t117, -qJ(2) * t117, -t92 * t117, -t94 * t117, t60 * t94 + t61 * t92 + (-t82 * t93 + (t32 * t92 - t33 * t94) * t95) * qJD(1), 0, 0, 0, 0, 0, -t48 * t138 - t143 * t158, -t51 * t138 - t142 * t158, 0, 0, 0, 0, 0, -t17 * t138 + (t111 * t96 + t112 * t98) * t81, -t108 * t138 + (t111 * t98 - t112 * t96) * t81; 0, 0, 0, 0, 0, 0, -t94 * t126, t92 * t126, t87 + (t32 * t94 + t33 * t92) * t138, 0, 0, 0, 0, 0, t44 + t147, -t43 - t148, 0, 0, 0, 0, 0, t11 - t150, t10 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t48, -t48 ^ 2 + t51 ^ 2, -t43 + t148, -t44 + t147, 0, -t110 * t158 - t65 * t51 + t102, t121 * t158 + t65 * t48 + t106, t160, t157, t161, t159, 0, (-t96 * t8 - t151) * t81 - t17 * t152 + (t124 * t96 - t151) * qJD(5) + t153, -t108 * t152 + (t124 * qJD(5) - t8 * t81 - t2) * t98 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t157, t161, t159, 0, t154 * t113 + t153, (-t154 * t6 - t2) * t98 + t156;];
tauc_reg = t1;

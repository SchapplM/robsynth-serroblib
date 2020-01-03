% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR6
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:08
% EndTime: 2019-12-31 17:58:11
% DurationCPUTime: 0.98s
% Computational Cost: add. (2498->200), mult. (6020->274), div. (0->0), fcn. (4312->8), ass. (0->116)
t142 = cos(qJ(4));
t118 = pkin(6) * qJD(1);
t70 = sin(pkin(8)) * pkin(1) + qJ(3);
t66 = t70 * qJD(1);
t80 = cos(pkin(9));
t74 = t80 * qJD(2);
t78 = sin(pkin(9));
t38 = t74 + (-t66 - t118) * t78;
t46 = qJD(2) * t78 + t66 * t80;
t39 = t118 * t80 + t46;
t83 = sin(qJ(4));
t18 = t142 * t39 + t38 * t83;
t16 = qJD(4) * pkin(7) + t18;
t65 = -cos(pkin(8)) * pkin(1) - pkin(3) * t80 - pkin(2);
t52 = qJD(1) * t65 + qJD(3);
t124 = t83 * t78;
t110 = qJD(1) * t124;
t108 = t142 * t80;
t68 = qJD(1) * t108;
t54 = -t68 + t110;
t64 = t142 * t78 + t80 * t83;
t56 = t64 * qJD(1);
t22 = t54 * pkin(4) - t56 * pkin(7) + t52;
t82 = sin(qJ(5));
t84 = cos(qJ(5));
t100 = t16 * t82 - t22 * t84;
t91 = t108 - t124;
t88 = t91 * qJD(3);
t148 = qJD(1) * t88;
t93 = -t142 * t38 + t39 * t83;
t11 = -qJD(4) * t93 + t148;
t109 = qJD(4) * t124;
t67 = qJD(4) * t68;
t47 = qJD(1) * t109 - t67;
t59 = t64 * qJD(4);
t48 = qJD(1) * t59;
t27 = pkin(4) * t48 + pkin(7) * t47;
t1 = -qJD(5) * t100 + t84 * t11 + t82 * t27;
t51 = qJD(5) + t54;
t102 = t100 * t51 + t1;
t6 = t16 * t84 + t22 * t82;
t2 = -qJD(5) * t6 - t82 * t11 + t27 * t84;
t150 = t6 * t51 + t2;
t106 = t51 * t82;
t42 = qJD(4) * t82 + t56 * t84;
t149 = t42 * t106;
t117 = qJD(5) * t42;
t24 = -t82 * t47 + t117;
t44 = t84 * t48;
t116 = qJD(5) * t82;
t58 = -qJD(4) * t108 + t109;
t128 = t58 * t84;
t94 = t116 * t64 + t128;
t147 = t64 * t44 - t51 * t94;
t146 = t56 ^ 2;
t143 = pkin(6) + t70;
t89 = t64 * qJD(3);
t12 = qJD(1) * t89 + t18 * qJD(4);
t60 = t143 * t78;
t61 = t143 * t80;
t92 = -t142 * t60 - t61 * t83;
t141 = t12 * t92;
t140 = t12 * t91;
t139 = t12 * t82;
t113 = t84 * qJD(4);
t23 = -qJD(5) * t113 + t116 * t56 + t47 * t84;
t138 = t23 * t82;
t40 = t56 * t82 - t113;
t137 = t40 * t54;
t136 = t40 * t56;
t135 = t40 * t82;
t134 = t42 * t40;
t133 = t42 * t56;
t132 = t42 * t84;
t131 = t48 * t91;
t130 = t56 * t54;
t129 = t58 * t82;
t127 = t64 * t84;
t21 = t82 * t24;
t125 = t82 * t48;
t115 = qJD(5) * t84;
t123 = -t115 * t40 - t21;
t122 = -t127 * t24 + t128 * t40;
t121 = t23 * t91 + t42 * t59;
t120 = -t48 * t64 + t54 * t58;
t119 = t78 ^ 2 + t80 ^ 2;
t114 = t58 * qJD(4);
t112 = t42 * t129;
t107 = qJD(1) * t119;
t105 = t51 * t84;
t104 = -t100 * t84 + t6 * t82;
t103 = -t100 * t82 - t6 * t84;
t99 = t24 * t91 - t59 * t40;
t29 = -pkin(4) * t91 - pkin(7) * t64 + t65;
t31 = t142 * t61 - t60 * t83;
t9 = t29 * t84 - t31 * t82;
t10 = t29 * t82 + t31 * t84;
t98 = (-t66 * t78 + t74) * t78 - t46 * t80;
t97 = t47 * t91 + t56 * t59;
t96 = -t106 * t54 - t116 * t51 + t44;
t95 = t115 * t64 - t129;
t15 = -qJD(4) * pkin(4) + t93;
t90 = -pkin(7) * t48 + t15 * t51;
t87 = -t125 * t64 - t51 * t95;
t86 = -qJD(5) * t104 + t1 * t84 - t2 * t82;
t53 = t54 ^ 2;
t50 = t59 * qJD(4);
t34 = pkin(4) * t59 + pkin(7) * t58;
t33 = pkin(4) * t56 + pkin(7) * t54;
t20 = qJD(4) * t31 + t89;
t19 = qJD(4) * t92 + t88;
t8 = t33 * t82 - t84 * t93;
t7 = t33 * t84 + t82 * t93;
t4 = -qJD(5) * t10 - t82 * t19 + t84 * t34;
t3 = qJD(5) * t9 + t84 * t19 + t82 * t34;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t107, (t107 * t70 - t98) * qJD(3), -t47 * t64 - t56 * t58, -t97 + t120, -t114, t54 * t59 - t131, -t50, 0, -qJD(4) * t20 + t48 * t65 + t52 * t59, -qJD(4) * t19 - t47 * t65 - t52 * t58, t11 * t91 + t12 * t64 - t18 * t59 - t19 * t54 + t20 * t56 - t31 * t48 + t47 * t92 - t58 * t93, t11 * t31 + t18 * t19 + t20 * t93 - t141, -t127 * t23 - t42 * t94, t112 + (t138 + (-t132 + t135) * qJD(5)) * t64 + t122, t121 + t147, t21 * t64 + t40 * t95, t87 + t99, t51 * t59 - t131, -t100 * t59 + t139 * t64 + t15 * t95 - t2 * t91 + t20 * t40 - t24 * t92 + t4 * t51 + t9 * t48, t1 * t91 - t10 * t48 + t12 * t127 - t15 * t94 + t20 * t42 + t23 * t92 - t3 * t51 - t6 * t59, -t10 * t24 + t23 * t9 - t3 * t40 - t4 * t42 + t104 * t58 + (qJD(5) * t103 - t1 * t82 - t2 * t84) * t64, t1 * t10 - t100 * t4 + t15 * t20 + t2 * t9 + t3 * t6 - t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t114, t97 + t120, t11 * t64 - t18 * t58 + t59 * t93 - t140, 0, 0, 0, 0, 0, 0, t87 - t99, t121 - t147, -t112 + (-t138 + (t132 + t135) * qJD(5)) * t64 + t122, t103 * t58 + t15 * t59 + t64 * t86 - t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 * qJD(1) ^ 2, t98 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t56 * qJD(4), t67 + (-t54 - t110) * qJD(4), -t53 - t146, t18 * t54 - t56 * t93, 0, 0, 0, 0, 0, 0, t96 - t136, -t51 ^ 2 * t84 - t125 - t133, (t23 - t137) * t84 + t149 + t123, t102 * t82 - t15 * t56 + t150 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t53 + t146, t67 + (t54 - t110) * qJD(4), -t130, 0, 0, -(qJD(3) + t52) * t56, t52 * t54 - t148, 0, 0, t105 * t42 - t138, (-t23 - t137) * t84 - t149 + t123, t105 * t51 + t125 - t133, t106 * t40 - t24 * t84, t96 + t136, -t51 * t56, -pkin(4) * t24 - t12 * t84 - t18 * t40 + t100 * t56 + (-pkin(7) * t115 - t7) * t51 + t90 * t82, pkin(4) * t23 + t139 - t18 * t42 + t6 * t56 + (pkin(7) * t116 + t8) * t51 + t90 * t84, t40 * t8 + t42 * t7 + ((-t24 + t117) * pkin(7) + t102) * t84 + ((qJD(5) * t40 - t23) * pkin(7) - t150) * t82, -t12 * pkin(4) + pkin(7) * t86 + t100 * t7 - t15 * t18 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t40 ^ 2 + t42 ^ 2, t40 * t51 - t23, -t134, t42 * t51 - t24, t48, -t15 * t42 + t150, t15 * t40 - t102, 0, 0;];
tauc_reg = t5;

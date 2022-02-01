% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:32
% EndTime: 2022-01-20 10:48:35
% DurationCPUTime: 0.59s
% Computational Cost: add. (889->130), mult. (1719->191), div. (0->0), fcn. (1106->8), ass. (0->111)
t83 = sin(qJ(5));
t84 = sin(qJ(4));
t86 = cos(qJ(5));
t87 = cos(qJ(4));
t56 = t83 * t87 + t86 * t84;
t78 = qJD(1) + qJD(2);
t43 = t56 * t78;
t77 = qJD(4) + qJD(5);
t139 = qJD(5) - t77;
t119 = pkin(1) * qJD(1);
t88 = cos(qJ(2));
t109 = t88 * t119;
t58 = t78 * pkin(2) + t109;
t85 = sin(qJ(2));
t111 = t85 * t119;
t82 = cos(pkin(9));
t67 = t82 * t111;
t81 = sin(pkin(9));
t34 = t81 * t58 + t67;
t105 = t34 + (pkin(7) + pkin(8)) * t78;
t15 = t87 * qJD(3) - t105 * t84;
t16 = t84 * qJD(3) + t105 * t87;
t138 = t82 * pkin(2);
t137 = t87 * pkin(4);
t129 = t82 * t85;
t72 = t88 * pkin(1) + pkin(2);
t121 = pkin(1) * t129 + t81 * t72;
t45 = pkin(7) + t121;
t136 = -pkin(8) - t45;
t70 = t81 * pkin(2) + pkin(7);
t135 = -pkin(8) - t70;
t107 = -pkin(3) - t137;
t66 = t81 * t111;
t33 = t82 * t58 - t66;
t21 = t107 * t78 - t33;
t27 = t77 * t56;
t116 = t84 * qJD(4);
t110 = pkin(4) * t116;
t118 = pkin(1) * qJD(2);
t47 = (t81 * t88 + t129) * t118;
t38 = qJD(1) * t47;
t28 = t78 * t110 + t38;
t124 = t86 * t87;
t128 = t83 * t84;
t55 = -t124 + t128;
t134 = t21 * t27 + t28 * t55;
t26 = t77 * t55;
t133 = -t21 * t26 + t28 * t56;
t22 = t26 * t77;
t112 = t78 * t124;
t113 = t78 * t128;
t41 = -t112 + t113;
t132 = t43 * t41;
t48 = t82 * t109 - t66;
t131 = t48 * t77;
t130 = t81 * t85;
t126 = t86 * t16;
t89 = qJD(4) ^ 2;
t123 = t89 * t84;
t115 = t87 * qJD(4);
t29 = -t78 * pkin(3) - t33;
t122 = t29 * t115 + t38 * t84;
t120 = t84 ^ 2 - t87 ^ 2;
t114 = pkin(4) * t78 * t84;
t108 = t78 * t115;
t14 = qJD(4) * pkin(4) + t15;
t106 = -pkin(4) * t77 - t14;
t49 = (t82 * t88 - t130) * t118;
t39 = qJD(1) * t49;
t104 = -t29 * t78 - t39;
t103 = -pkin(1) * t130 + t82 * t72;
t102 = qJD(4) * t136;
t101 = qJD(4) * t135;
t44 = -pkin(3) - t103;
t46 = t81 * t109 + t67;
t100 = -t46 + t110;
t98 = (-qJD(2) + t78) * t119;
t97 = (-qJD(1) - t78) * t118;
t95 = t45 * t89 + t47 * t78;
t94 = -t46 * t78 + t70 * t89;
t93 = qJD(4) * (t44 * t78 - t49);
t92 = qJD(4) * ((-pkin(3) - t138) * t78 + t48);
t4 = t15 * qJD(4) + t87 * t39;
t5 = -t16 * qJD(4) - t84 * t39;
t91 = -t21 * t43 - t83 * t4 + t86 * t5;
t17 = qJD(5) * t112 + t86 * t108 - t77 * t113;
t90 = t21 * t41 + (t139 * t16 - t5) * t83;
t76 = t78 ^ 2;
t75 = t87 * pkin(8);
t74 = t89 * t87;
t61 = t107 - t138;
t57 = 0.2e1 * t84 * t108;
t54 = t87 * t70 + t75;
t53 = t135 * t84;
t51 = t87 * t101;
t50 = t84 * t101;
t40 = -0.2e1 * t120 * t78 * qJD(4);
t37 = t44 - t137;
t36 = t47 + t110;
t32 = t87 * t45 + t75;
t31 = t136 * t84;
t24 = t29 * t116;
t23 = t27 * t77;
t18 = t27 * t78;
t13 = t87 * t102 - t84 * t49;
t12 = t84 * t102 + t87 * t49;
t10 = -t41 ^ 2 + t43 ^ 2;
t6 = t41 * t77 + t17;
t2 = t17 * t56 - t43 * t26;
t1 = -t17 * t55 - t56 * t18 + t26 * t41 - t43 * t27;
t3 = [0, 0, 0, 0, t85 * t97, t88 * t97, -t38 * t103 + t39 * t121 - t33 * t47 + t34 * t49, t57, t40, t74, -t123, 0, t24 + t84 * t93 + (-t38 - t95) * t87, t95 * t84 + t87 * t93 + t122, t2, t1, -t22, -t23, 0, t36 * t41 + t37 * t18 + (-t83 * t12 + t86 * t13 + (-t31 * t83 - t32 * t86) * qJD(5)) * t77 + t134, t36 * t43 + t37 * t17 - (t86 * t12 + t83 * t13 + (t31 * t86 - t32 * t83) * qJD(5)) * t77 + t133; 0, 0, 0, 0, t85 * t98, t88 * t98, t33 * t46 - t34 * t48 + (-t38 * t82 + t39 * t81) * pkin(2), t57, t40, t74, -t123, 0, t24 + t84 * t92 + (-t38 - t94) * t87, t94 * t84 + t87 * t92 + t122, t2, t1, -t22, -t23, 0, t61 * t18 + (-t83 * t50 + t86 * t51 + (-t53 * t83 - t54 * t86) * qJD(5)) * t77 + t56 * t131 + t100 * t41 + t134, t61 * t17 - (t86 * t50 + t83 * t51 + (t53 * t86 - t54 * t83) * qJD(5)) * t77 - t55 * t131 + t100 * t43 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t74, 0, 0, 0, 0, 0, -t23, t22; 0, 0, 0, 0, 0, 0, 0, -t84 * t76 * t87, t120 * t76, 0, 0, 0, t104 * t84, t104 * t87, t132, t10, t6, 0, 0, -t41 * t114 - (-t83 * t15 - t126) * t77 + (t106 * t83 - t126) * qJD(5) + t91, -t43 * t114 + (t106 * qJD(5) + t15 * t77 - t4) * t86 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t10, t6, 0, 0, t91 + t139 * (-t83 * t14 - t126), (-t139 * t14 - t4) * t86 + t90;];
tauc_reg = t3;

% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:44
% DurationCPUTime: 0.83s
% Computational Cost: add. (1536->163), mult. (4297->226), div. (0->0), fcn. (3225->6), ass. (0->97)
t84 = qJD(3) + qJD(4);
t86 = cos(pkin(7));
t90 = cos(qJ(3));
t114 = t90 * t86;
t105 = qJD(1) * t114;
t85 = sin(pkin(7));
t88 = sin(qJ(3));
t116 = t88 * t85;
t106 = qJD(1) * t116;
t59 = -t105 + t106;
t70 = t90 * t85 + t88 * t86;
t61 = t70 * qJD(1);
t87 = sin(qJ(4));
t89 = cos(qJ(4));
t97 = t87 * t59 - t89 * t61;
t124 = t97 * t84;
t77 = qJD(3) * t105;
t52 = qJD(3) * t106 - t77;
t64 = t70 * qJD(3);
t53 = qJD(1) * t64;
t92 = t97 * qJD(4) + t87 * t52 - t89 * t53;
t131 = t92 - t124;
t33 = -t89 * t59 - t87 * t61;
t122 = t33 * t84;
t108 = qJD(4) * t89;
t109 = qJD(4) * t87;
t95 = -t59 * t108 - t61 * t109 - t89 * t52 - t87 * t53;
t130 = t95 - t122;
t121 = t97 ^ 2;
t123 = t33 ^ 2;
t129 = t121 - t123;
t120 = t33 * t97;
t81 = -t86 * pkin(2) - pkin(1);
t73 = t81 * qJD(1) + qJD(2);
t39 = t59 * pkin(3) + t73;
t128 = t39 * t33;
t127 = t39 * t97;
t113 = pkin(5) + qJ(2);
t74 = t113 * t85;
t75 = t113 * t86;
t41 = -t88 * t74 + t90 * t75;
t71 = qJD(1) * t74;
t72 = qJD(1) * t75;
t38 = -t88 * t71 + t90 * t72;
t94 = t70 * qJD(2);
t93 = qJD(1) * t94;
t22 = -t38 * qJD(3) - t93;
t17 = t52 * pkin(6) + t22;
t25 = -t59 * pkin(6) + t38;
t100 = -t25 * t109 + t87 * t17;
t107 = qJD(1) * qJD(2);
t102 = t85 * t107;
t110 = qJD(3) * t90;
t78 = qJD(2) * t114;
t112 = qJD(1) * t78 - t71 * t110;
t21 = (-qJD(3) * t72 - t102) * t88 + t112;
t16 = -t53 * pkin(6) + t21;
t117 = t88 * t72;
t37 = -t90 * t71 - t117;
t24 = -t61 * pkin(6) + t37;
t23 = qJD(3) * pkin(3) + t24;
t1 = (qJD(4) * t23 + t16) * t89 + t100;
t126 = t61 ^ 2;
t125 = pkin(3) * t61;
t119 = t61 * t59;
t118 = t87 * t25;
t115 = t89 * t25;
t111 = t85 ^ 2 + t86 ^ 2;
t104 = -pkin(3) * t84 - t23;
t103 = t111 * qJD(1) ^ 2;
t101 = -t87 * t16 + t89 * t17;
t40 = -t90 * t74 - t88 * t75;
t6 = t87 * t23 + t115;
t28 = -t70 * pkin(6) + t40;
t69 = t114 - t116;
t29 = t69 * pkin(6) + t41;
t9 = t89 * t28 - t87 * t29;
t10 = t87 * t28 + t89 * t29;
t36 = t87 * t69 + t89 * t70;
t96 = 0.2e1 * t111 * t107;
t26 = -t74 * t110 + t78 + (-qJD(2) * t85 - qJD(3) * t75) * t88;
t2 = -t6 * qJD(4) + t101;
t27 = -t41 * qJD(3) - t94;
t63 = qJD(3) * t116 - t86 * t110;
t56 = t59 ^ 2;
t45 = -t69 * pkin(3) + t81;
t35 = -t89 * t69 + t87 * t70;
t19 = t63 * pkin(6) + t27;
t18 = -t64 * pkin(6) + t26;
t15 = t36 * qJD(4) - t87 * t63 + t89 * t64;
t14 = -t69 * t108 + t70 * t109 + t89 * t63 + t87 * t64;
t8 = t89 * t24 - t118;
t7 = -t87 * t24 - t115;
t5 = t89 * t23 - t118;
t4 = -t10 * qJD(4) - t87 * t18 + t89 * t19;
t3 = t9 * qJD(4) + t89 * t18 + t87 * t19;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, qJ(2) * t96, -t52 * t70 - t61 * t63, -t52 * t69 - t70 * t53 + t63 * t59 - t61 * t64, -t63 * qJD(3), -t53 * t69 + t59 * t64, -t64 * qJD(3), 0, t27 * qJD(3) + t81 * t53 + t73 * t64, -t26 * qJD(3) - t81 * t52 - t73 * t63, t21 * t69 - t22 * t70 - t26 * t59 - t27 * t61 + t37 * t63 - t38 * t64 + t40 * t52 - t41 * t53, t21 * t41 + t22 * t40 + t38 * t26 + t37 * t27, t14 * t97 + t36 * t95, -t14 * t33 + t15 * t97 - t35 * t95 + t36 * t92, -t14 * t84, -t15 * t33 - t35 * t92, -t15 * t84, 0, -t45 * t92 + t39 * t15 + t4 * t84 + (-t33 * t64 + t35 * t53) * pkin(3), t45 * t95 - t39 * t14 - t3 * t84 + (t36 * t53 - t64 * t97) * pkin(3), -t1 * t35 + t10 * t92 + t5 * t14 - t6 * t15 - t2 * t36 + t3 * t33 + t4 * t97 - t9 * t95, t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4 + (t39 * t64 + t45 * t53) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -qJ(2) * t103, 0, 0, 0, 0, 0, 0, 0.2e1 * t61 * qJD(3), t77 + (-t59 - t106) * qJD(3), -t56 - t126, t37 * t61 + t38 * t59, 0, 0, 0, 0, 0, 0, -t92 - t124, t95 + t122, -t121 - t123, t53 * pkin(3) - t6 * t33 - t5 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t56 + t126, t77 + (t59 - t106) * qJD(3), -t119, 0, 0, -t73 * t61 - t93, t88 * t102 + t73 * t59 + (t37 + t117) * qJD(3) - t112, 0, 0, t120, t129, t130, -t120, t131, 0, t33 * t125 + t127 - t7 * t84 + (t104 * t87 - t115) * qJD(4) + t101, t97 * t125 - t128 + t8 * t84 + (t104 * qJD(4) - t16) * t89 - t100, -t6 * t97 - t8 * t33 + t5 * t33 - t7 * t97 + (-t95 * t89 + t92 * t87 + (t33 * t89 - t87 * t97) * qJD(4)) * pkin(3), -t5 * t7 - t6 * t8 + (t1 * t87 + t2 * t89 - t39 * t61 + (-t5 * t87 + t6 * t89) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t129, t130, -t120, t131, 0, t6 * t84 + t127 + t2, t5 * t84 - t1 - t128, 0, 0;];
tauc_reg = t11;

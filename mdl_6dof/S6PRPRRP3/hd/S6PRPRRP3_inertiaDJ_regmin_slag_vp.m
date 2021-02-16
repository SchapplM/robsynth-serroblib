% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:31
% EndTime: 2021-01-16 01:38:38
% DurationCPUTime: 1.35s
% Computational Cost: add. (1683->180), mult. (4286->330), div. (0->0), fcn. (4334->10), ass. (0->106)
t125 = cos(qJ(4));
t60 = sin(pkin(11));
t62 = cos(pkin(11));
t65 = sin(qJ(4));
t74 = t125 * t62 - t65 * t60;
t116 = qJ(3) + pkin(8);
t47 = t116 * t60;
t48 = t116 * t62;
t75 = -t125 * t47 - t65 * t48;
t71 = t74 * qJD(3) + t75 * qJD(4);
t45 = t125 * t60 + t65 * t62;
t52 = -t62 * pkin(3) - pkin(2);
t73 = -pkin(4) * t74 - pkin(9) * t45 + t52;
t132 = -qJD(5) * t73 - t71;
t61 = sin(pkin(6));
t68 = cos(qJ(2));
t118 = t61 * t68;
t64 = sin(qJ(5));
t104 = t64 * t118;
t66 = sin(qJ(2));
t119 = t61 * t66;
t63 = cos(pkin(6));
t38 = -t60 * t119 + t62 * t63;
t39 = t62 * t119 + t63 * t60;
t22 = t125 * t39 + t65 * t38;
t67 = cos(qJ(5));
t17 = t67 * t22 - t104;
t108 = qJD(2) * t66;
t100 = t61 * t108;
t76 = t125 * t38 - t65 * t39;
t109 = qJD(2) * t61;
t99 = t68 * t109;
t70 = t76 * qJD(4) + t74 * t99;
t79 = t67 * t118 + t64 * t22;
t7 = t79 * qJD(5) - t64 * t100 - t67 * t70;
t55 = qJD(5) * t67;
t8 = qJD(5) * t104 + t67 * t100 - t22 * t55 - t64 * t70;
t131 = qJD(5) * (-t17 * t67 - t64 * t79) + t7 * t64 - t8 * t67;
t110 = qJ(6) * t67;
t34 = t125 * t48 - t65 * t47;
t92 = -t34 * t64 + t67 * t73;
t11 = -pkin(5) * t74 - t45 * t110 + t92;
t111 = qJ(6) * t45;
t114 = t67 * t34 + t64 * t73;
t12 = -t64 * t111 + t114;
t41 = t45 * qJD(4);
t127 = t41 * pkin(5);
t106 = t67 * qJD(6);
t40 = t74 * qJD(4);
t89 = pkin(4) * t41 - pkin(9) * t40;
t6 = t132 * t64 - t34 * t55 + t67 * t89;
t96 = qJD(5) * t111;
t69 = -t45 * t106 - t40 * t110 + t64 * t96 + t6;
t3 = t69 + t127;
t103 = t132 * t67 - t64 * t89;
t4 = t67 * t96 + (qJ(6) * t40 + qJD(5) * t34 + qJD(6) * t45) * t64 + t103;
t130 = -t3 * t67 + t4 * t64 + (t11 * t64 - t12 * t67) * qJD(5);
t129 = 0.2e1 * t52;
t128 = 0.2e1 * qJD(5);
t126 = t67 * pkin(5);
t124 = t40 * t64;
t123 = t45 * t40;
t122 = t45 * t64;
t121 = t45 * t67;
t53 = -pkin(4) - t126;
t120 = t53 * t67;
t117 = t64 * t67;
t115 = -qJ(6) - pkin(9);
t113 = t60 ^ 2 + t62 ^ 2;
t58 = t64 ^ 2;
t59 = t67 ^ 2;
t112 = t58 - t59;
t107 = qJD(5) * t64;
t105 = -0.2e1 * pkin(4) * qJD(5);
t102 = pkin(5) * t107;
t101 = t76 * t107;
t98 = t64 * t55;
t95 = -t53 + t126;
t94 = t113 * t68;
t93 = -0.4e1 * t45 * t117;
t91 = t112 * qJD(5);
t90 = 0.2e1 * t113 * qJD(3);
t88 = -pkin(4) * t40 - pkin(9) * t41;
t87 = pkin(4) * t45 - pkin(9) * t74;
t84 = pkin(5) * t58 + t120;
t82 = -t17 * t64 + t67 * t79;
t80 = -t38 * t60 + t39 * t62;
t31 = t45 * t55 + t124;
t78 = t45 * t107 - t40 * t67;
t77 = t41 * t64 - t55 * t74;
t19 = t45 * qJD(3) + t34 * qJD(4);
t50 = t115 * t67;
t49 = t115 * t64;
t42 = t45 ^ 2;
t37 = -t64 * qJD(6) + t115 * t55;
t36 = -t115 * t107 - t106;
t28 = t107 * t74 + t41 * t67;
t20 = pkin(5) * t122 - t75;
t15 = t22 * qJD(4) + t45 * t99;
t13 = t31 * pkin(5) + t19;
t10 = -t15 * t67 - t101;
t9 = t15 * t64 - t55 * t76;
t5 = t34 * t107 + t103;
t2 = t15 * t122 - t31 * t76 - t41 * t79 - t74 * t8;
t1 = t15 * t121 - t17 * t41 - t7 * t74 + t76 * t78;
t14 = [0, 0, 0, 0, 0, 0, 0.2e1 * (t80 - t119) * t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t15 * t76 - 0.2e1 * t17 * t7 - 0.2e1 * t79 * t8; 0, 0, -t100, -t99, -t62 * t100, t94 * t109, t80 * qJD(3) + (-pkin(2) * t66 + qJ(3) * t94) * t109, 0, 0, 0, 0, 0, (-t108 * t74 - t41 * t68) * t61, (t45 * t108 - t40 * t68) * t61, 0, 0, 0, 0, 0, t2, t1, t2, t1, t131 * t45 + t82 * t40, t11 * t8 - t12 * t7 - t13 * t76 + t15 * t20 - t17 * t4 - t3 * t79; 0, 0, 0, 0, 0, t90, qJ(3) * t90, 0.2e1 * t123, 0.2e1 * t40 * t74 - 0.2e1 * t41 * t45, 0, 0, 0, t41 * t129, t40 * t129, 0.2e1 * t59 * t123 - 0.2e1 * t42 * t98, t112 * t42 * t128 + t40 * t93, 0.2e1 * t41 * t121 + 0.2e1 * t74 * t78, -0.2e1 * t41 * t122 + 0.2e1 * t31 * t74, -0.2e1 * t74 * t41, 0.2e1 * t19 * t122 - 0.2e1 * t31 * t75 + 0.2e1 * t41 * t92 - 0.2e1 * t6 * t74, -0.2e1 * t114 * t41 + 0.2e1 * t19 * t121 - 0.2e1 * t5 * t74 + 0.2e1 * t75 * t78, 0.2e1 * t11 * t41 + 0.2e1 * t13 * t122 + 0.2e1 * t20 * t31 - 0.2e1 * t3 * t74, -0.2e1 * t12 * t41 + 0.2e1 * t13 * t121 - 0.2e1 * t20 * t78 - 0.2e1 * t4 * t74, 0.2e1 * (-t11 * t67 - t12 * t64) * t40 + 0.2e1 * t130 * t45, 0.2e1 * t11 * t3 - 0.2e1 * t12 * t4 + 0.2e1 * t13 * t20; 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, 0, 0, 0, 0, t28, -t77, t28, -t77, (-t58 - t59) * t40, -t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t70, 0, 0, 0, 0, 0, t10, t9, t10, t9, qJD(5) * t82 - t8 * t64 - t7 * t67, -pkin(5) * t101 + t15 * t53 - t17 * t36 - t37 * t79 + t49 * t8 + t50 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, 0, -t19, -t71, t40 * t117 - t45 * t91, qJD(5) * t93 - t112 * t40, t77, t28, 0, -t19 * t67 + t88 * t64 + (-t64 * t75 - t67 * t87) * qJD(5), t19 * t64 + t88 * t67 + (t64 * t87 - t67 * t75) * qJD(5), t53 * t124 - t13 * t67 - t37 * t74 + t49 * t41 + (t20 * t64 + t45 * t84) * qJD(5), t40 * t120 + t13 * t64 - t36 * t74 + t50 * t41 + (t95 * t122 + t20 * t67) * qJD(5), (-t37 * t45 - t40 * t49 - t4 + (t45 * t50 - t11) * qJD(5)) * t67 + (t36 * t45 + t40 * t50 - t3 + (t45 * t49 - t12) * qJD(5)) * t64, t102 * t20 + t11 * t37 - t12 * t36 + t13 * t53 + t3 * t49 + t4 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t36 + t67 * t37 + (-t49 * t64 - t50 * t67) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t98, -0.2e1 * t91, 0, 0, 0, t64 * t105, t67 * t105, -0.2e1 * t95 * t107, t84 * t128, -0.2e1 * t36 * t67 - 0.2e1 * t37 * t64 + 0.2e1 * (-t49 * t67 + t50 * t64) * qJD(5), 0.2e1 * t102 * t53 + 0.2e1 * t36 * t50 + 0.2e1 * t49 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, t8, t7, 0, t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t31, t41, t6, t5, t69 + 0.2e1 * t127, t4, t78 * pkin(5), t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t55, -t107, -t55, 0, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t107, 0, -pkin(9) * t55, pkin(9) * t107, t37, t36, -pkin(5) * t55, t37 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t78, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t55, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;

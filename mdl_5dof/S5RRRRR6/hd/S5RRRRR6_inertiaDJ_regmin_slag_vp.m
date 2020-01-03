% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:19
% EndTime: 2020-01-03 12:15:21
% DurationCPUTime: 0.63s
% Computational Cost: add. (1216->111), mult. (2805->167), div. (0->0), fcn. (2518->8), ass. (0->104)
t130 = qJD(3) + qJD(4);
t129 = -pkin(7) - pkin(8);
t102 = qJD(3) * t129;
t84 = sin(qJ(3));
t61 = t84 * t102;
t88 = cos(qJ(3));
t62 = t88 * t102;
t83 = sin(qJ(4));
t87 = cos(qJ(4));
t67 = t129 * t84;
t81 = t88 * pkin(8);
t68 = t88 * pkin(7) + t81;
t93 = t87 * t67 - t83 * t68;
t24 = -qJD(4) * t93 - t87 * t61 - t83 * t62;
t85 = sin(qJ(2));
t74 = t85 * pkin(1) + pkin(7);
t125 = -pkin(8) - t74;
t101 = qJD(3) * t125;
t114 = pkin(1) * qJD(2);
t89 = cos(qJ(2));
t107 = t89 * t114;
t98 = t88 * t107;
t48 = t101 * t84 + t98;
t99 = t84 * t107;
t49 = t101 * t88 - t99;
t56 = t125 * t84;
t57 = t88 * t74 + t81;
t96 = t87 * t56 - t83 * t57;
t15 = -qJD(4) * t96 - t87 * t48 - t83 * t49;
t59 = t83 * t84 - t87 * t88;
t43 = t130 * t59;
t128 = t43 * pkin(9);
t60 = t83 * t88 + t87 * t84;
t127 = t60 * pkin(9);
t126 = t89 * pkin(1);
t82 = sin(qJ(5));
t86 = cos(qJ(5));
t41 = -t82 * t59 + t86 * t60;
t44 = t130 * t60;
t14 = qJD(5) * t41 - t82 * t43 + t86 * t44;
t110 = t84 * qJD(3);
t78 = pkin(3) * t110;
t34 = t44 * pkin(4) + t78;
t79 = t85 * t114;
t29 = t34 + t79;
t40 = t86 * t59 + t82 * t60;
t77 = -t88 * pkin(3) - pkin(2);
t51 = t59 * pkin(4) + t77;
t50 = t51 - t126;
t124 = t50 * t14 + t29 * t40;
t13 = -qJD(5) * t40 - t86 * t43 - t82 * t44;
t123 = t50 * t13 + t29 * t41;
t122 = t51 * t14 + t34 * t40;
t121 = t51 * t13 + t34 * t41;
t63 = t79 + t78;
t66 = t77 - t126;
t120 = t66 * t44 + t63 * t59;
t119 = -t66 * t43 + t63 * t60;
t118 = t77 * t44 + t59 * t78;
t117 = -t77 * t43 + t60 * t78;
t113 = pkin(3) * qJD(4);
t106 = t83 * t113;
t112 = qJD(5) * t82;
t116 = t83 * pkin(3) * t112 + t82 * t106;
t76 = -pkin(2) - t126;
t80 = t88 * qJD(3);
t115 = t76 * t80 + t84 * t79;
t111 = qJD(5) * t86;
t109 = pkin(2) * t110;
t108 = pkin(2) * t80;
t105 = t87 * t113;
t104 = pkin(4) * t112;
t103 = pkin(4) * t111;
t75 = t87 * pkin(3) + pkin(4);
t100 = (-pkin(4) - t75) * qJD(5);
t95 = -t83 * t56 - t87 * t57;
t92 = -t83 * t67 - t87 * t68;
t91 = t76 * t110 - t79 * t88;
t16 = qJD(4) * t95 - t83 * t48 + t87 * t49;
t25 = qJD(4) * t92 - t83 * t61 + t87 * t62;
t90 = (-t83 * t111 + (-t82 * t87 - t83 * t86) * qJD(4)) * pkin(3);
t72 = 0.2e1 * t84 * t80;
t58 = 0.2e1 * (-t84 ^ 2 + t88 ^ 2) * qJD(3);
t55 = t59 * pkin(9);
t42 = t44 * pkin(9);
t38 = -t112 * t75 + t90;
t37 = (-qJD(5) * t75 - t105) * t86 + t116;
t33 = -t55 - t92;
t32 = t93 - t127;
t28 = -0.2e1 * t60 * t43;
t27 = -t55 - t95;
t26 = t96 - t127;
t19 = 0.2e1 * t43 * t59 - 0.2e1 * t60 * t44;
t18 = t25 + t128;
t17 = -t24 - t42;
t10 = t16 + t128;
t9 = -t15 - t42;
t6 = 0.2e1 * t41 * t13;
t5 = -0.2e1 * t13 * t40 - 0.2e1 * t41 * t14;
t4 = -t82 * t17 + t86 * t18 + (-t32 * t82 - t33 * t86) * qJD(5);
t3 = -t86 * t17 - t82 * t18 + (-t32 * t86 + t33 * t82) * qJD(5);
t2 = t86 * t10 - t82 * t9 + (-t26 * t82 - t27 * t86) * qJD(5);
t1 = -t82 * t10 - t86 * t9 + (-t26 * t86 + t27 * t82) * qJD(5);
t7 = [0, 0, 0, 0, -0.2e1 * t79, -0.2e1 * t107, t72, t58, 0, 0, 0, 0.2e1 * t91, 0.2e1 * t115, t28, t19, 0, 0, 0, 0.2e1 * t120, 0.2e1 * t119, t6, t5, 0, 0, 0, 0.2e1 * t124, 0.2e1 * t123; 0, 0, 0, 0, -t79, -t107, t72, t58, 0, 0, 0, t91 - t109, -t108 + t115, t28, t19, 0, 0, 0, t118 + t120, t117 + t119, t6, t5, 0, 0, 0, t122 + t124, t121 + t123; 0, 0, 0, 0, 0, 0, t72, t58, 0, 0, 0, -0.2e1 * t109, -0.2e1 * t108, t28, t19, 0, 0, 0, 0.2e1 * t118, 0.2e1 * t117, t6, t5, 0, 0, 0, 0.2e1 * t122, 0.2e1 * t121; 0, 0, 0, 0, 0, 0, 0, 0, t80, -t110, 0, -t74 * t80 - t99, t110 * t74 - t98, 0, 0, -t43, -t44, 0, t16, t15, 0, 0, t13, -t14, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, t80, -t110, 0, -pkin(7) * t80, pkin(7) * t110, 0, 0, -t43, -t44, 0, t25, t24, 0, 0, t13, -t14, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t106, -0.2e1 * t105, 0, 0, 0, 0, 0, 0.2e1 * t38, 0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, t16, t15, 0, 0, t13, -t14, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, t25, t24, 0, 0, t13, -t14, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, 0, 0, 0, 0, t100 * t82 + t90, (t100 - t105) * t86 + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t104, -0.2e1 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;

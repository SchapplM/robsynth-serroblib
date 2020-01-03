% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:33
% DurationCPUTime: 1.49s
% Computational Cost: add. (1643->185), mult. (4168->339), div. (0->0), fcn. (3597->6), ass. (0->108)
t130 = cos(qJ(4));
t101 = qJD(4) * t130;
t75 = sin(qJ(4));
t123 = qJD(4) * t75;
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t43 = -t74 * t101 + t73 * t123;
t105 = t130 * t74;
t137 = -t75 * t73 + t105;
t76 = sin(qJ(2));
t77 = cos(qJ(2));
t136 = (t76 ^ 2 - t77 ^ 2) * qJD(2);
t135 = -0.2e1 * t43;
t51 = t130 * t73 + t75 * t74;
t44 = t51 * qJD(4);
t134 = 0.2e1 * t44;
t133 = 2 * qJD(5);
t132 = pkin(2) * t76;
t131 = t73 * pkin(6);
t129 = t73 * t76;
t128 = t74 * t76;
t126 = pkin(7) + qJ(3);
t106 = pkin(3) + t131;
t94 = -t77 * pkin(2) - t76 * qJ(3);
t53 = -pkin(1) + t94;
t48 = t74 * t53;
t31 = -pkin(7) * t128 - t106 * t77 + t48;
t63 = t74 * t77 * pkin(6);
t40 = t73 * t53 + t63;
t36 = -pkin(7) * t129 + t40;
t13 = t130 * t36 + t75 * t31;
t120 = t77 * qJD(2);
t111 = t73 * t120;
t66 = pkin(6) * t120;
t45 = pkin(3) * t111 + t66;
t52 = pkin(3) * t129 + t76 * pkin(6);
t124 = qJD(3) * t73;
t122 = t74 * qJD(3);
t67 = t76 * qJD(2);
t121 = t76 * qJD(3);
t119 = t77 * qJD(5);
t25 = t51 * t120 - t43 * t76;
t41 = t51 * t76;
t118 = 0.2e1 * t41 * t25;
t117 = t137 * t134;
t116 = t77 * t131;
t115 = pkin(6) * t128;
t114 = -0.2e1 * pkin(1) * qJD(2);
t113 = pkin(4) * t67;
t110 = t73 * t121;
t109 = t74 * t120;
t108 = t74 * t121;
t107 = t76 * t120;
t65 = -t74 * pkin(3) - pkin(2);
t104 = qJ(5) * t67;
t100 = t130 * qJD(3);
t56 = t126 * t74;
t102 = t126 * t73;
t88 = t130 * t102;
t21 = qJD(4) * t88 - t74 * t100 + (qJD(4) * t56 + t124) * t75;
t22 = t56 * t101 + t75 * t122 + (-t126 * t123 + t100) * t73;
t37 = t75 * t56 + t88;
t38 = -t75 * t102 + t130 * t56;
t103 = -t38 * t21 + t37 * t22;
t99 = 0.2e1 * t107;
t98 = t73 * t109;
t69 = t73 ^ 2;
t70 = t74 ^ 2;
t97 = 0.2e1 * (t69 + t70) * qJD(3);
t96 = -0.2e1 * t136;
t93 = -qJ(3) * t77 + t132;
t24 = -t105 * t120 + t75 * t111 + t76 * t44;
t42 = t137 * t76;
t92 = t24 * t41 - t42 * t25;
t91 = -t137 * t25 + t41 * t44;
t34 = -t108 + (pkin(6) * t129 + t74 * t93) * qJD(2);
t35 = -t110 + (t73 * t93 - t115) * qJD(2);
t90 = -t34 * t73 + t35 * t74;
t89 = -t137 * t43 - t51 * t44;
t87 = -t21 * t77 - t38 * t67;
t86 = t22 * t77 - t37 * t67;
t85 = -t77 * t25 + t41 * t67;
t84 = -t137 * t67 - t77 * t44;
t83 = -t126 * t77 + t132;
t12 = t130 * t31 - t75 * t36;
t78 = -t108 + (t106 * t76 + t83 * t74) * qJD(2);
t79 = -t110 + (t83 * t73 - t115) * qJD(2);
t3 = -t31 * t101 + t36 * t123 - t130 * t79 - t75 * t78;
t82 = t21 * t41 + t22 * t42 - t37 * t24 - t38 * t25;
t81 = -0.2e1 * t137 * t21 + 0.2e1 * t22 * t51 - 0.2e1 * t37 * t43 - 0.2e1 * t38 * t44;
t80 = -t137 * t24 - t51 * t25 + t43 * t41 - t42 * t44;
t4 = -t36 * t101 - t31 * t123 + t130 * t78 - t75 * t79;
t59 = -0.2e1 * t107;
t39 = t48 - t116;
t33 = t51 * t135;
t30 = -pkin(4) * t137 - t51 * qJ(5) + t65;
t27 = t43 * t77 + t51 * t67;
t18 = t44 * pkin(4) + t43 * qJ(5) - t51 * qJD(5);
t17 = t41 * pkin(4) - t42 * qJ(5) + t52;
t15 = -0.2e1 * t42 * t24;
t10 = 0.2e1 * t24 * t77 + 0.2e1 * t42 * t67;
t9 = t77 * pkin(4) - t12;
t7 = -t77 * qJ(5) + t13;
t6 = -t24 * t51 - t42 * t43;
t5 = t25 * pkin(4) + t24 * qJ(5) - t42 * qJD(5) + t45;
t2 = -t4 - t113;
t1 = t104 - t3 - t119;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t96, 0, t59, 0, 0, t76 * t114, t77 * t114, 0, 0, t70 * t99, -0.4e1 * t76 * t98, 0.2e1 * t74 * t136, t69 * t99, t73 * t96, t59, -0.2e1 * t34 * t77 + 0.2e1 * (t39 + 0.2e1 * t116) * t67, 0.2e1 * t35 * t77 + 0.2e1 * (-t40 + 0.2e1 * t63) * t67, 0.2e1 * (-t34 * t74 - t35 * t73) * t76 + 0.2e1 * (-t39 * t74 - t40 * t73) * t120, 0.2e1 * pkin(6) ^ 2 * t107 + 0.2e1 * t39 * t34 + 0.2e1 * t40 * t35, t15, 0.2e1 * t92, t10, t118, -0.2e1 * t85, t59, 0.2e1 * t12 * t67 + 0.2e1 * t52 * t25 - 0.2e1 * t4 * t77 + 0.2e1 * t45 * t41, -0.2e1 * t13 * t67 - 0.2e1 * t52 * t24 - 0.2e1 * t3 * t77 + 0.2e1 * t45 * t42, 0.2e1 * t12 * t24 - 0.2e1 * t13 * t25 + 0.2e1 * t3 * t41 - 0.2e1 * t4 * t42, 0.2e1 * t12 * t4 - 0.2e1 * t13 * t3 + 0.2e1 * t52 * t45, t15, t10, -0.2e1 * t92, t59, 0.2e1 * t85, t118, 0.2e1 * t17 * t25 + 0.2e1 * t2 * t77 + 0.2e1 * t5 * t41 - 0.2e1 * t9 * t67, -0.2e1 * t1 * t41 + 0.2e1 * t2 * t42 - 0.2e1 * t9 * t24 - 0.2e1 * t7 * t25, -0.2e1 * t1 * t77 + 0.2e1 * t17 * t24 - 0.2e1 * t5 * t42 + 0.2e1 * t7 * t67, 0.2e1 * t7 * t1 + 0.2e1 * t17 * t5 + 0.2e1 * t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, -t67, 0, -t66, pkin(6) * t67, 0, 0, t98, (-t69 + t70) * t120, t73 * t67, -t98, t74 * t67, 0, t77 * t124 + (t73 * t94 - t63) * qJD(2), t77 * t122 + (t74 * t94 + t116) * qJD(2), t90, -pkin(2) * t66 + (-t39 * t73 + t40 * t74) * qJD(3) + t90 * qJ(3), t6, t80, t27, t91, -t84, 0, -t137 * t45 + t65 * t25 + t52 * t44 + t86, -t65 * t24 - t52 * t43 + t45 * t51 + t87, t12 * t43 - t13 * t44 - t137 * t3 - t4 * t51 + t82, -t12 * t22 - t13 * t21 - t3 * t38 - t4 * t37 + t45 * t65, t6, t27, -t80, 0, t84, t91, -t137 * t5 + t17 * t44 + t18 * t41 + t30 * t25 + t86, t1 * t137 + t2 * t51 - t9 * t43 - t7 * t44 + t82, t17 * t43 - t18 * t42 + t30 * t24 - t5 * t51 - t87, t1 * t38 + t17 * t18 + t2 * t37 - t7 * t21 + t9 * t22 + t5 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJ(3) * t97, t33, 0.2e1 * t89, 0, -t117, 0, 0, t65 * t134, t65 * t135, t81, 0.2e1 * t103, t33, 0, -0.2e1 * t89, 0, 0, -t117, -0.2e1 * t137 * t18 + 0.2e1 * t30 * t44, t81, -0.2e1 * t18 * t51 + 0.2e1 * t30 * t43, 0.2e1 * t30 * t18 + 0.2e1 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t109, 0, t66, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t45, 0, 0, 0, 0, 0, 0, t25, 0, t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t43, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t25, t67, t4, t3, 0, 0, 0, -t24, 0, t67, t25, 0, t4 + 0.2e1 * t113, pkin(4) * t24 - t25 * qJ(5) - t41 * qJD(5), 0.2e1 * t104 - t3 - 0.2e1 * t119, -t2 * pkin(4) + t1 * qJ(5) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, -t44, 0, -t22, t21, 0, 0, 0, -t43, 0, 0, t44, 0, -t22, pkin(4) * t43 - t44 * qJ(5) + qJD(5) * t137, -t21, -t22 * pkin(4) - t21 * qJ(5) + t38 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, qJ(5) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;

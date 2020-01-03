% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:37
% EndTime: 2019-12-31 18:45:41
% DurationCPUTime: 1.15s
% Computational Cost: add. (767->143), mult. (1806->245), div. (0->0), fcn. (1310->6), ass. (0->97)
t58 = cos(qJ(4));
t56 = sin(qJ(4));
t68 = t58 * pkin(4) + t56 * qJ(5);
t113 = t68 * qJD(4) - t58 * qJD(5);
t57 = sin(qJ(3));
t116 = t113 * t57;
t53 = t58 ^ 2;
t59 = cos(qJ(3));
t94 = t59 * qJD(3);
t45 = t53 * t94;
t51 = t56 ^ 2;
t84 = t51 * t94;
t22 = t45 + t84;
t48 = -cos(pkin(8)) * pkin(1) - pkin(2);
t110 = t57 * pkin(7);
t75 = -t59 * pkin(3) - t110;
t29 = t48 + t75;
t47 = sin(pkin(8)) * pkin(1) + pkin(6);
t107 = t47 * t59;
t31 = t58 * t107;
t115 = t56 * t29 + t31;
t102 = t51 - t53;
t33 = t102 * qJD(4);
t67 = pkin(4) * t56 - qJ(5) * t58;
t64 = t47 + t67;
t16 = t64 * t57;
t21 = t67 * qJD(4) - t56 * qJD(5);
t32 = -pkin(3) - t68;
t114 = (-t32 * t59 + t110) * qJD(3) - qJD(4) * t16 - t21 * t57;
t111 = pkin(7) * t59;
t74 = pkin(3) * t57 - t111;
t63 = t74 * qJD(3);
t49 = t57 * qJD(3);
t81 = t47 * t49;
t4 = -qJD(4) * t115 + t56 * t81 + t58 * t63;
t112 = 0.2e1 * qJD(5);
t109 = t47 * t56;
t108 = t47 * t58;
t106 = t58 * t29;
t50 = qJD(4) * t58;
t105 = t29 * t50 + t56 * t63;
t103 = t22 * pkin(7);
t52 = t57 ^ 2;
t101 = -t59 ^ 2 + t52;
t100 = qJD(3) * t16;
t99 = qJD(3) * t58;
t98 = qJD(4) * t56;
t97 = qJD(4) * t57;
t96 = qJD(4) * t59;
t93 = -0.2e1 * pkin(3) * qJD(4);
t92 = t56 * t107;
t91 = 0.2e1 * qJD(3) * t48;
t90 = pkin(4) * t49;
t89 = pkin(7) * t98;
t88 = pkin(7) * t50;
t87 = t47 * t98;
t86 = t56 * t96;
t85 = t58 * t96;
t83 = t56 * t50;
t82 = t57 * t94;
t80 = t58 * t94;
t79 = t47 * t94;
t78 = t101 * qJD(3);
t77 = t56 * t80;
t76 = t52 * t83;
t6 = -t59 * qJ(5) + t115;
t7 = -t106 + (pkin(4) + t109) * t59;
t73 = t56 * t7 + t58 * t6;
t72 = -t56 * t6 + t58 * t7;
t8 = -t92 + t106;
t71 = -t115 * t56 - t58 * t8;
t70 = -t115 * t58 + t56 * t8;
t24 = t58 * t49 + t86;
t1 = (-qJD(5) - t87) * t59 + (qJ(5) - t108) * t49 + t105;
t2 = -t4 - t90;
t61 = t72 * qJD(4) + t1 * t58 + t2 * t56;
t3 = t24 * t47 - t105;
t60 = t71 * qJD(4) - t3 * t58 - t4 * t56;
t42 = -0.2e1 * t82;
t41 = -0.2e1 * t83;
t40 = 0.2e1 * t83;
t39 = pkin(7) * t85;
t35 = t53 * t82;
t34 = t51 * t82;
t26 = -t56 * t49 + t85;
t25 = t57 * t50 + t56 * t94;
t23 = t56 * t97 - t80;
t18 = 0.2e1 * t35 - 0.2e1 * t76;
t17 = 0.2e1 * t34 + 0.2e1 * t76;
t15 = t102 * t97 - t77;
t14 = -t56 * t78 + t57 * t85;
t13 = 0.4e1 * t57 * t83 - t45 + t84;
t12 = 0.2e1 * t101 * t99 + 0.2e1 * t57 * t86;
t11 = t52 * t33 - 0.2e1 * t57 * t77;
t10 = 0.2e1 * t34 + 0.2e1 * t35 - 0.2e1 * t82;
t5 = t64 * t94 + t116;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t82, -0.2e1 * t78, 0, t42, 0, 0, t57 * t91, t59 * t91, 0, 0, t18, 0.2e1 * t11, t12, t17, 0.2e1 * t14, t42, 0.2e1 * t52 * t47 * t50 - 0.2e1 * t4 * t59 + 0.2e1 * (t8 + 0.2e1 * t92) * t49, -0.2e1 * t52 * t87 - 0.2e1 * t3 * t59 + 0.2e1 * (-t115 + 0.2e1 * t31) * t49, 0.2e1 * t71 * t94 + 0.2e1 * (t70 * qJD(4) + t3 * t56 - t4 * t58) * t57, 0.2e1 * t47 ^ 2 * t82 - 0.2e1 * t115 * t3 + 0.2e1 * t8 * t4, t18, t12, -0.2e1 * t11, t42, -0.2e1 * t14, t17, 0.2e1 * (t56 * t100 + t2) * t59 + 0.2e1 * (-qJD(3) * t7 + t16 * t50 + t5 * t56) * t57, 0.2e1 * t72 * t94 + 0.2e1 * (-t73 * qJD(4) - t1 * t56 + t2 * t58) * t57, 0.2e1 * (-t16 * t99 - t1) * t59 + 0.2e1 * (qJD(3) * t6 + t16 * t98 - t5 * t58) * t57, 0.2e1 * t6 * t1 + 0.2e1 * t16 * t5 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t57 + (t101 * t47 - t70 * t59) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, (t73 * qJD(3) - t5) * t59 + (t61 + t100) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, -t49, 0, -t79, t81, 0, 0, -t15, -t13, -t26, t15, t24, 0, t39 + (-pkin(3) * t58 + t109) * t97 + (t75 * t56 - t31) * qJD(3), (t57 * t108 + t74 * t56) * qJD(4) + (t75 * t58 + t92) * qJD(3), t60, -pkin(3) * t79 + t60 * pkin(7), -t15, -t26, t13, 0, -t24, t15, t39 + (t32 * t97 - t5) * t58 - t114 * t56, t61, (-t5 + (t32 * t57 + t111) * qJD(4)) * t56 + t114 * t58, t61 * pkin(7) + t16 * t21 + t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t94, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, t22, -pkin(3) * t49 + t103, 0, 0, 0, 0, 0, 0, -t24, t22, t26, -t59 * t21 + t32 * t49 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -0.2e1 * t33, 0, t41, 0, 0, t56 * t93, t58 * t93, 0, 0, t40, 0, 0.2e1 * t33, 0, 0, t41, -0.2e1 * t21 * t58 + 0.2e1 * t32 * t98, 0, -0.2e1 * t21 * t56 - 0.2e1 * t32 * t50, 0.2e1 * t32 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t25, t49, t4, t3, 0, 0, 0, -t23, 0, t49, t25, 0, t4 + 0.2e1 * t90, (-pkin(4) * t94 - qJ(5) * t97) * t58 + (-qJ(5) * t94 + (pkin(4) * qJD(4) - qJD(5)) * t57) * t56, (-0.2e1 * qJD(5) - t87) * t59 + (0.2e1 * qJ(5) - t108) * t49 + t105, -t2 * pkin(4) + t1 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t23, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, -t23, -t67 * t94 - t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t98, 0, -t88, t89, 0, 0, 0, t50, 0, 0, t98, 0, -t88, -t113, -t89, -t113 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, qJ(5) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t23, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;

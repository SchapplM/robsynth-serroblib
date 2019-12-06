% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:32
% EndTime: 2019-12-05 16:13:37
% DurationCPUTime: 0.93s
% Computational Cost: add. (412->127), mult. (1319->243), div. (0->0), fcn. (1037->6), ass. (0->95)
t55 = sin(pkin(8));
t51 = t55 ^ 2;
t56 = cos(pkin(8));
t52 = t56 ^ 2;
t109 = (t51 + t52) * qJD(4);
t59 = cos(qJ(3));
t58 = sin(qJ(2));
t57 = sin(qJ(3));
t88 = t57 * qJD(3);
t77 = t58 * t88;
t60 = cos(qJ(2));
t85 = t60 * qJD(2);
t108 = t59 * t85 - t77;
t53 = t57 ^ 2;
t54 = t59 ^ 2;
t72 = qJD(3) * (t53 - t54);
t31 = -pkin(4) * t56 - qJ(5) * t55 - pkin(3);
t94 = t57 * qJ(4);
t107 = (-t31 * t59 + t94) * qJD(3);
t106 = 0.2e1 * t56;
t86 = t59 * qJD(3);
t50 = pkin(6) * t86;
t66 = pkin(4) * t55 - qJ(5) * t56;
t92 = qJD(5) * t57;
t7 = -t56 * t92 + t66 * t86 + t50;
t105 = t7 * t55;
t104 = t7 * t56;
t98 = t58 * t59;
t81 = t56 * t98;
t26 = -t55 * t60 + t81;
t97 = t60 * t56;
t11 = -t56 * t77 + (t55 * t58 + t59 * t97) * qJD(2);
t6 = t11 * t56;
t103 = t26 * t56 * qJD(4) + qJ(4) * t6;
t102 = t55 * t59;
t22 = -t57 * qJD(4) + (pkin(3) * t57 - qJ(4) * t59) * qJD(3);
t101 = t56 * t22;
t67 = -pkin(3) * t59 - t94;
t35 = -pkin(2) + t67;
t100 = t56 * t35;
t99 = t56 * t59;
t47 = pkin(6) * t99;
t15 = t55 * t35 + t47;
t96 = qJ(4) * t109;
t93 = qJD(4) * t59;
t91 = t51 * qJD(5);
t90 = t55 * qJD(4);
t89 = t55 * qJD(5);
t87 = t58 * qJD(2);
t83 = pkin(6) * t102;
t82 = -0.2e1 * pkin(2) * qJD(3);
t80 = t53 * t85;
t43 = t55 * t88;
t79 = t56 * t88;
t45 = t56 * t86;
t78 = t57 * t86;
t76 = t58 * t85;
t74 = pkin(6) * t55 + pkin(4);
t10 = t108 * t55 - t56 * t87;
t73 = t10 * t55 + t6;
t71 = 0.2e1 * t78;
t34 = t55 * t45;
t70 = t58 * t71;
t8 = pkin(6) * t43 + t101;
t19 = t55 * t22;
t9 = -pkin(6) * t79 + t19;
t69 = -t55 * t8 + t56 * t9;
t68 = t57 * t34;
t25 = t55 * t98 + t97;
t64 = qJ(4) * t10 + qJD(4) * t25;
t38 = t53 * t76;
t63 = 0.2e1 * t58 ^ 2 * t78 + 0.2e1 * t10 * t25 + 0.2e1 * t26 * t11 + 0.2e1 * t38;
t27 = t57 * t85 + t58 * t86;
t62 = t10 * t59 - t25 * t88 + (t70 + t80) * t55;
t61 = (t10 * t56 - t11 * t55) * t57 + (t25 * t56 - t26 * t55) * t86;
t44 = t55 * t86;
t42 = t59 * t90;
t41 = -0.2e1 * t78;
t40 = pkin(6) * t80;
t33 = t52 * t71;
t32 = t51 * t71;
t30 = 0.2e1 * t109;
t24 = t55 * t72;
t23 = (t51 - t52) * t86;
t21 = t72 * t106;
t20 = (pkin(6) + t66) * t57;
t17 = t27 * t56;
t16 = t27 * t55;
t14 = -t83 + t100;
t13 = t59 * t74 - t100;
t12 = -qJ(5) * t59 + t15;
t4 = -t74 * t88 - t101;
t3 = -t59 * qJD(5) + t19 + (-pkin(6) * t56 + qJ(5)) * t88;
t1 = t56 * t80 + t11 * t59 + (-t26 + 0.2e1 * t81) * t88;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t38 + 0.2e1 * (t54 - 0.1e1) * t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t85, 0, 0, 0, 0, 0, 0, 0, 0, -t59 * t87 - t60 * t88, t57 * t87 - t60 * t86, (t53 + t54) * t85, t40 + (pkin(6) * t54 * t60 - pkin(2) * t58) * qJD(2), 0, 0, 0, 0, 0, 0, t62, t1, t61, pkin(6) * t70 - t10 * t14 + t11 * t15 - t25 * t8 + t26 * t9 + t40, 0, 0, 0, 0, 0, 0, t62, t61, -t1, t57 * t58 * t7 + t10 * t13 + t11 * t12 + t20 * t27 + t25 * t4 + t26 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -0.2e1 * t72, 0, t41, 0, 0, t57 * t82, t59 * t82, 0, 0, t33, -0.4e1 * t68, t21, t32, -0.2e1 * t24, t41, -0.2e1 * t8 * t59 + 0.2e1 * (t14 + 0.2e1 * t83) * t88, 0.2e1 * t9 * t59 + 0.2e1 * (-t15 + 0.2e1 * t47) * t88, 0.2e1 * (-t55 * t9 - t56 * t8) * t57 + 0.2e1 * (-t14 * t56 - t15 * t55) * t86, 0.2e1 * pkin(6) ^ 2 * t78 + 0.2e1 * t14 * t8 + 0.2e1 * t15 * t9, t33, t21, 0.4e1 * t68, t41, 0.2e1 * t24, t32, 0.2e1 * t57 * t105 + 0.2e1 * t4 * t59 + 0.2e1 * (t102 * t20 - t13 * t57) * qJD(3), 0.2e1 * (-t3 * t55 + t4 * t56) * t57 + 0.2e1 * (-t12 * t55 + t13 * t56) * t86, -0.2e1 * t57 * t104 - 0.2e1 * t3 * t59 + 0.2e1 * (t12 * t57 - t20 * t99) * qJD(3), 0.2e1 * t12 * t3 + 0.2e1 * t13 * t4 + 0.2e1 * t20 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t108, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, t73, -pkin(3) * t27 + t55 * t64 + t103, 0, 0, 0, 0, 0, 0, -t17, t73, -t16, t27 * t31 + (-t58 * t92 + t64) * t55 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, -t88, 0, -t50, pkin(6) * t88, 0, 0, t34, -t23, t43, -t34, t79, 0, t42 + (t55 * t67 - t47) * qJD(3), t56 * t93 + (t56 * t67 + t83) * qJD(3), t69, -pkin(3) * t50 + (-t14 * t55 + t15 * t56) * qJD(4) + t69 * qJ(4), t34, t43, t23, 0, -t79, -t34, -t107 * t55 - t57 * t91 - t104 + t42, t3 * t56 + t4 * t55, -t105 + (t57 * t89 + t107 - t93) * t56, t7 * t31 + (qJ(4) * t3 + qJD(4) * t12) * t56 + (qJ(4) * t4 + qJD(4) * t13 - qJD(5) * t20) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0.2e1 * t96, 0, 0, 0, 0, 0, 0, t89 * t106, t30, 0.2e1 * t91, -0.2e1 * t31 * t89 + 0.2e1 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, t50, 0, 0, 0, 0, 0, 0, t44, 0, -t45, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t45, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;

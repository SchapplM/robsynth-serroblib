% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:12
% EndTime: 2019-12-31 18:22:16
% DurationCPUTime: 1.05s
% Computational Cost: add. (1098->143), mult. (2623->277), div. (0->0), fcn. (2270->8), ass. (0->98)
t55 = sin(pkin(9));
t56 = cos(pkin(9));
t113 = sin(qJ(5));
t84 = qJD(5) * t113;
t114 = cos(qJ(5));
t85 = qJD(5) * t114;
t30 = t55 * t84 - t56 * t85;
t33 = t113 * t55 - t114 * t56;
t120 = -0.2e1 * t30;
t90 = t114 * t55;
t34 = t113 * t56 + t90;
t31 = t34 * qJD(5);
t119 = 0.2e1 * t31;
t118 = 2 * qJD(3);
t58 = sin(qJ(3));
t117 = pkin(3) * t58;
t59 = cos(qJ(3));
t116 = t59 * pkin(3);
t101 = t59 * qJD(3);
t13 = t33 * t101 + t58 * t31;
t27 = t33 * t58;
t115 = t13 * t33 + t27 * t31;
t14 = t34 * t101 - t30 * t58;
t26 = t34 * t58;
t112 = t26 * t14;
t111 = t27 * t13;
t47 = sin(pkin(8)) * pkin(1) + pkin(6);
t110 = t58 * t47;
t109 = t59 * t47;
t108 = pkin(7) + qJ(4);
t35 = t56 * t109;
t104 = t58 * qJ(4);
t48 = -cos(pkin(8)) * pkin(1) - pkin(2);
t69 = t48 - t116;
t65 = t69 - t104;
t21 = t55 * t65 + t35;
t51 = t55 ^ 2;
t52 = t56 ^ 2;
t107 = t51 + t52;
t106 = t58 ^ 2 - t59 ^ 2;
t105 = qJ(4) * t59;
t103 = qJD(4) * t59;
t50 = t58 * qJD(3);
t102 = t58 * qJD(4);
t100 = t55 * t109;
t99 = t56 * t110;
t98 = t48 * t118;
t97 = t51 * t101;
t96 = t55 * t101;
t95 = t55 * t102;
t94 = t56 * t50;
t93 = t56 * t101;
t92 = t56 * t102;
t91 = t58 * t101;
t36 = t47 * t101;
t87 = t55 * t47 + pkin(4);
t86 = qJD(4) * t113;
t83 = t114 * qJD(4);
t82 = t107 * qJD(4);
t81 = 0.2e1 * t91;
t80 = t55 * t93;
t79 = 0.2e1 * t82;
t78 = -0.2e1 * t106 * qJD(3);
t75 = t108 * t113;
t74 = -t104 - t116;
t73 = -t105 + t117;
t72 = -t34 * t14 + t26 * t30;
t18 = -t92 + (t55 * t110 + t56 * t73) * qJD(3);
t19 = -t95 + (t55 * t73 - t99) * qJD(3);
t71 = -t18 * t55 + t19 * t56;
t20 = t56 * t65 - t100;
t70 = -t20 * t55 + t21 * t56;
t68 = t108 * t90;
t67 = -t59 * t31 + t33 * t50;
t66 = -t108 * t59 + t117;
t64 = -t87 * t59 + (-t108 * t58 + t69) * t56;
t63 = t95 - (t66 * t55 - t99) * qJD(3);
t62 = t114 * t64;
t61 = t113 * t64;
t60 = -t92 + (t66 * t56 + t87 * t58) * qJD(3);
t49 = -t56 * pkin(4) - pkin(3);
t43 = t52 * t101;
t42 = t55 * t50;
t41 = -0.2e1 * t91;
t39 = t108 * t56;
t29 = (pkin(4) * t55 + t47) * t58;
t25 = pkin(4) * t96 + t36;
t23 = t114 * t39 - t55 * t75;
t22 = -t113 * t39 - t68;
t16 = t30 * t59 + t34 * t50;
t15 = -t55 * t58 * pkin(7) + t21;
t12 = -t39 * t85 - t56 * t86 + (qJD(5) * t75 - t83) * t55;
t11 = qJD(5) * t68 + t39 * t84 + t55 * t86 - t56 * t83;
t4 = t114 * t15 + t61;
t3 = -t113 * t15 + t62;
t2 = -qJD(5) * t61 + t113 * t63 + t114 * t60 - t15 * t85;
t1 = -qJD(5) * t62 - t113 * t60 + t114 * t63 + t15 * t84;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t78, 0, t41, 0, 0, t58 * t98, t59 * t98, 0, 0, t52 * t81, -0.4e1 * t58 * t80, t106 * t56 * t118, t51 * t81, t55 * t78, t41, -0.2e1 * t18 * t59 + 0.2e1 * (t20 + 0.2e1 * t100) * t50, 0.2e1 * t19 * t59 + 0.2e1 * (-t21 + 0.2e1 * t35) * t50, 0.2e1 * (-t18 * t56 - t19 * t55) * t58 + 0.2e1 * (-t20 * t56 - t21 * t55) * t101, 0.2e1 * t47 ^ 2 * t91 + 0.2e1 * t20 * t18 + 0.2e1 * t21 * t19, 0.2e1 * t111, 0.2e1 * t26 * t13 + 0.2e1 * t27 * t14, 0.2e1 * t13 * t59 - 0.2e1 * t27 * t50, 0.2e1 * t112, 0.2e1 * t59 * t14 - 0.2e1 * t26 * t50, t41, 0.2e1 * t29 * t14 - 0.2e1 * t2 * t59 + 0.2e1 * t25 * t26 + 0.2e1 * t3 * t50, -0.2e1 * t1 * t59 - 0.2e1 * t29 * t13 - 0.2e1 * t25 * t27 - 0.2e1 * t4 * t50, 0.2e1 * t1 * t26 + 0.2e1 * t3 * t13 - 0.2e1 * t4 * t14 + 0.2e1 * t2 * t27, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t29 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t58 + (t106 * t47 + t70 * t59) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t27 - t4 * t13 - t3 * t14 - t2 * t26 - t25 * t59 + t29 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t107) * t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t111 - 0.2e1 * t91 + 0.2e1 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, -t50, 0, -t36, t47 * t50, 0, 0, t80, t43 - t97, t42, -t80, t94, 0, t55 * t103 + (t74 * t55 - t35) * qJD(3), t56 * t103 + (t74 * t56 + t100) * qJD(3), t71, -pkin(3) * t36 + t71 * qJ(4) + t70 * qJD(4), -t13 * t34 + t27 * t30, t72 + t115, t16, t14 * t33 + t26 * t31, -t67, 0, -t12 * t59 + t49 * t14 + t22 * t50 + t25 * t33 + t29 * t31, -t11 * t59 - t49 * t13 - t23 * t50 + t25 * t34 - t29 * t30, t1 * t33 + t11 * t26 + t12 * t27 + t22 * t13 - t23 * t14 - t2 * t34 + t3 * t30 - t4 * t31, -t1 * t23 - t4 * t11 + t3 * t12 + t2 * t22 + t25 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t101, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t42, t43 + t97, t58 * t82 + (t107 * t105 - t117) * qJD(3), 0, 0, 0, 0, 0, 0, t67, t16, -t72 + t115, t27 * t11 - t26 * t12 - t13 * t23 - t14 * t22 + t49 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, qJ(4) * t79, t34 * t120, 0.2e1 * t30 * t33 - 0.2e1 * t34 * t31, 0, t33 * t119, 0, 0, t49 * t119, t49 * t120, 0.2e1 * t11 * t33 - 0.2e1 * t12 * t34 + 0.2e1 * t22 * t30 - 0.2e1 * t23 * t31, -0.2e1 * t23 * t11 + 0.2e1 * t22 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t93, 0, t36, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, -t14, t50, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, -t31, 0, t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;

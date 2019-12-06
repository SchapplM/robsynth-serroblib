% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:23
% EndTime: 2019-12-05 16:23:27
% DurationCPUTime: 0.92s
% Computational Cost: add. (1017->118), mult. (2614->235), div. (0->0), fcn. (2481->8), ass. (0->69)
t60 = sin(qJ(3));
t88 = -qJ(4) - pkin(6);
t47 = t88 * t60;
t57 = sin(pkin(9));
t58 = cos(pkin(9));
t62 = cos(qJ(3));
t91 = t88 * t62;
t92 = -t57 * t47 + t58 * t91;
t23 = t58 * t47 + t57 * t91;
t59 = sin(qJ(5));
t73 = t58 * pkin(3) + pkin(4);
t89 = cos(qJ(5));
t90 = pkin(3) * t57;
t36 = -t59 * t90 + t89 * t73;
t82 = t62 * qJD(3);
t83 = t60 * qJD(3);
t87 = t57 * t82 + t58 * t83;
t55 = t60 ^ 2;
t56 = t62 ^ 2;
t86 = t55 + t56;
t61 = sin(qJ(2));
t85 = qJD(3) * t61;
t84 = qJD(5) * t59;
t54 = t61 * qJD(2);
t63 = cos(qJ(2));
t81 = t63 * qJD(2);
t79 = -0.2e1 * pkin(2) * qJD(3);
t53 = pkin(3) * t83;
t78 = t60 * t82;
t77 = t61 * t81;
t76 = t63 * t83;
t75 = t60 * t81;
t74 = t62 * t81;
t52 = -t62 * pkin(3) - pkin(2);
t71 = t86 * t63;
t70 = qJD(5) * t89;
t41 = t57 * t62 + t58 * t60;
t40 = t57 * t60 - t58 * t62;
t69 = -t41 * pkin(7) + t23;
t33 = t41 * t61;
t34 = t40 * t61;
t15 = -t59 * t33 - t89 * t34;
t22 = -t59 * t40 + t89 * t41;
t67 = t59 * t69;
t66 = t89 * t69;
t37 = t59 * t73 + t89 * t90;
t18 = t92 * qJD(3) - t41 * qJD(4);
t19 = t23 * qJD(3) - t40 * qJD(4);
t39 = -t57 * t83 + t58 * t82;
t65 = -t39 * pkin(7) + t18;
t64 = -t87 * pkin(7) + t19;
t29 = t37 * qJD(5);
t28 = t36 * qJD(5);
t27 = t40 * pkin(4) + t52;
t25 = t87 * pkin(4) + t53;
t21 = t89 * t40 + t59 * t41;
t20 = -t40 * pkin(7) - t92;
t17 = t41 * t85 + t57 * t75 - t58 * t74;
t16 = t40 * t85 - t41 * t81;
t14 = -t89 * t33 + t59 * t34;
t8 = t22 * qJD(5) + t59 * t39 + t89 * t87;
t7 = -t89 * t39 + t40 * t70 + t41 * t84 + t59 * t87;
t6 = t89 * t20 + t67;
t5 = -t59 * t20 + t66;
t4 = -t15 * qJD(5) + t89 * t16 + t59 * t17;
t3 = -t59 * t16 + t89 * t17 + t33 * t70 - t34 * t84;
t2 = -qJD(5) * t67 - t20 * t70 - t59 * t64 + t89 * t65;
t1 = -qJD(5) * t66 + t20 * t84 - t59 * t65 - t89 * t64;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t86) * t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t33 * t16 + 0.2e1 * t34 * t17 - 0.2e1 * t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t14 * t4 - 0.2e1 * t15 * t3 - 0.2e1 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t81, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t54 - t76, t60 * t54 - t63 * t82, qJD(2) * t71, (-pkin(2) * t61 + pkin(6) * t71) * qJD(2), 0, 0, 0, 0, 0, 0, t40 * t54 - t63 * t87, -t63 * t39 + t41 * t54, -t16 * t41 + t17 * t40 + t33 * t39 + t34 * t87, -pkin(3) * t76 + t16 * t23 + t17 * t92 - t33 * t18 - t34 * t19 + t52 * t54, 0, 0, 0, 0, 0, 0, t21 * t54 - t63 * t8, t22 * t54 + t63 * t7, t14 * t7 - t15 * t8 + t3 * t21 - t4 * t22, -t15 * t1 + t14 * t2 - t63 * t25 + t27 * t54 - t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78, 0.2e1 * (-t55 + t56) * qJD(3), 0, -0.2e1 * t78, 0, 0, t60 * t79, t62 * t79, 0, 0, 0.2e1 * t41 * t39, -0.2e1 * t39 * t40 - 0.2e1 * t41 * t87, 0, 0.2e1 * t40 * t87, 0, 0, 0.2e1 * t40 * t53 + 0.2e1 * t52 * t87, 0.2e1 * t52 * t39 + 0.2e1 * t41 * t53, -0.2e1 * t18 * t41 - 0.2e1 * t19 * t40 - 0.2e1 * t23 * t39 + 0.2e1 * t87 * t92, 0.2e1 * t23 * t18 - 0.2e1 * t19 * t92 + 0.2e1 * t52 * t53, -0.2e1 * t22 * t7, 0.2e1 * t7 * t21 - 0.2e1 * t22 * t8, 0, 0.2e1 * t21 * t8, 0, 0, 0.2e1 * t25 * t21 + 0.2e1 * t27 * t8, 0.2e1 * t25 * t22 - 0.2e1 * t27 * t7, 0.2e1 * t1 * t21 - 0.2e1 * t2 * t22 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t27 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t82 - t75, t61 * t83 - t74, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, (t16 * t58 - t17 * t57) * pkin(3), 0, 0, 0, 0, 0, 0, t4, t3, 0, -t14 * t29 + t15 * t28 - t3 * t37 + t4 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, -t83, 0, -pkin(6) * t82, pkin(6) * t83, 0, 0, 0, 0, t39, 0, -t87, 0, t18, -t19, (-t58 * t39 - t57 * t87) * pkin(3), (t18 * t58 + t19 * t57) * pkin(3), 0, 0, -t7, 0, -t8, 0, t2, t1, -t28 * t21 + t29 * t22 + t36 * t7 - t37 * t8, -t1 * t37 + t2 * t36 + t6 * t28 - t5 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t29, -0.2e1 * t28, 0, 0.2e1 * t37 * t28 - 0.2e1 * t36 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t39, 0, t53, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;

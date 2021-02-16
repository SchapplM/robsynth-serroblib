% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP4
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
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:35
% EndTime: 2021-01-15 12:56:39
% DurationCPUTime: 0.56s
% Computational Cost: add. (908->102), mult. (2248->186), div. (0->0), fcn. (1984->6), ass. (0->69)
t52 = sin(qJ(4));
t53 = sin(qJ(3));
t54 = cos(qJ(4));
t55 = cos(qJ(3));
t31 = t52 * t55 + t54 * t53;
t50 = sin(pkin(8));
t25 = t31 * t50;
t70 = qJD(3) + qJD(4);
t51 = cos(pkin(8));
t32 = -t51 * pkin(2) - t50 * pkin(6) - pkin(1);
t84 = pkin(7) * t50;
t60 = -t32 + t84;
t65 = t55 * t51 * qJ(2);
t87 = t60 * t53 - t65;
t86 = 0.2e1 * t51;
t85 = 2 * qJD(3);
t83 = t54 * pkin(3);
t22 = t70 * t31;
t82 = t22 * t51;
t81 = t52 * t53;
t80 = t54 * t55;
t74 = qJD(3) * t55;
t76 = qJD(2) * t55;
t79 = t32 * t74 + t51 * t76;
t63 = t50 * t74;
t28 = pkin(3) * t63 + t50 * qJD(2);
t29 = (pkin(3) * t53 + qJ(2)) * t50;
t78 = qJ(2) * t53;
t77 = qJD(2) * t53;
t75 = qJD(3) * t53;
t73 = qJD(4) * t52;
t72 = qJD(4) * t54;
t71 = qJ(2) * qJD(3);
t69 = t50 * t81;
t68 = t50 * t80;
t67 = pkin(3) * t73;
t66 = pkin(3) * t72;
t64 = t50 * t75;
t62 = t51 * t77;
t61 = t53 * t71;
t59 = t50 * t51 * t85;
t48 = t50 ^ 2;
t58 = 0.2e1 * (t51 ^ 2 + t48) * qJD(2);
t16 = -t60 * t55 + (-pkin(3) - t78) * t51;
t57 = -t52 * t16 + t54 * t87;
t14 = (-t51 * t78 - t55 * t84) * qJD(3) + t79;
t56 = t87 * qJD(3) - t62;
t3 = -t54 * t14 - t16 * t72 - t52 * t56 - t73 * t87;
t13 = -qJD(4) * t69 - t52 * t64 + t70 * t68;
t1 = t13 * qJ(5) + t25 * qJD(5) + t3;
t4 = t57 * qJD(4) - t52 * t14 + t54 * t56;
t12 = t70 * t25;
t26 = t68 - t69;
t2 = t12 * qJ(5) - t26 * qJD(5) + t4;
t45 = pkin(4) + t83;
t44 = -0.2e1 * t66;
t43 = -0.2e1 * t67;
t37 = t51 * t66;
t36 = t51 * t67;
t30 = t80 - t81;
t21 = -t54 * t74 - t55 * t72 + t70 * t81;
t20 = -t62 + (-t53 * t32 - t65) * qJD(3);
t19 = t51 * t61 - t79;
t18 = t25 * pkin(4) + t29;
t17 = t21 * t51;
t7 = t13 * pkin(4) + t28;
t6 = -t25 * qJ(5) - t57;
t5 = -t51 * pkin(4) - t26 * qJ(5) + t54 * t16 + t52 * t87;
t8 = [0, 0, 0, 0, t58, qJ(2) * t58, -0.2e1 * t48 * t53 * t74, (t53 ^ 2 - t55 ^ 2) * t48 * t85, t53 * t59, t55 * t59, 0, -0.2e1 * t20 * t51 + 0.2e1 * (t55 * t71 + t77) * t48, -0.2e1 * t19 * t51 + 0.2e1 * (-t61 + t76) * t48, -0.2e1 * t26 * t12, 0.2e1 * t12 * t25 - 0.2e1 * t26 * t13, t12 * t86, t13 * t86, 0, 0.2e1 * t29 * t13 + 0.2e1 * t28 * t25 - 0.2e1 * t4 * t51, -0.2e1 * t29 * t12 + 0.2e1 * t28 * t26 - 0.2e1 * t3 * t51, 0.2e1 * t18 * t13 - 0.2e1 * t2 * t51 + 0.2e1 * t7 * t25, -0.2e1 * t1 * t51 - 0.2e1 * t18 * t12 + 0.2e1 * t7 * t26, 0.2e1 * t1 * t25 + 0.2e1 * t5 * t12 - 0.2e1 * t6 * t13 - 0.2e1 * t2 * t26, -0.2e1 * t6 * t1 + 0.2e1 * t18 * t7 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t75, t51 * t74, 0, 0, 0, 0, 0, t82, -t17, t82, -t17, t30 * t12 - t31 * t13 + t21 * t25 + t22 * t26, -t1 * t31 + t2 * t30 - t6 * t21 - t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t31 * t21 - 0.2e1 * t30 * t22; 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, t20, t19, 0, 0, -t12, -t13, 0, t36 + t4, t3 + t37, t36 + t2, t1 + t37, t45 * t12 + (-t13 * t52 + (-t25 * t54 + t26 * t52) * qJD(4)) * pkin(3), t2 * t45 + (-t1 * t52 + (-t5 * t52 + t54 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74, 0, 0, 0, 0, 0, -t22, t21, -t22, t21, 0, -t22 * t45 + (-t21 * t52 + (-t30 * t52 + t31 * t54) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t44, t43, t44, 0, 0.2e1 * (-t45 + t83) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, t4, t3, t2, t1, t12 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, -t22, t21, 0, -t22 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t66, -t67, -t66, 0, -pkin(4) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;

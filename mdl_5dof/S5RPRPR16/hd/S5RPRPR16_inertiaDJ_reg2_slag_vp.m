% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR16_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:35
% DurationCPUTime: 0.77s
% Computational Cost: add. (506->94), mult. (985->176), div. (0->0), fcn. (679->4), ass. (0->74)
t38 = cos(qJ(3));
t78 = pkin(3) + pkin(7);
t63 = t78 * t38;
t36 = sin(qJ(3));
t73 = qJ(4) * t36;
t83 = t63 + t73;
t72 = t38 * qJ(4);
t82 = -t36 * pkin(3) + t72;
t32 = t36 ^ 2;
t34 = t38 ^ 2;
t53 = (t32 - t34) * qJD(3);
t35 = sin(qJ(5));
t31 = t35 ^ 2;
t37 = cos(qJ(5));
t33 = t37 ^ 2;
t76 = t31 - t33;
t52 = t76 * qJD(5);
t40 = -pkin(1) - pkin(6);
t77 = pkin(4) - t40;
t21 = t77 * t36;
t66 = t36 * qJD(4);
t81 = (-t36 * t78 + t72) * qJD(3) + qJD(5) * t21 + t66;
t20 = qJ(2) - t82;
t17 = t36 * pkin(7) + t20;
t55 = t77 * t38;
t47 = t37 * t55;
t48 = -t38 * qJD(4) + qJD(2);
t71 = qJD(3) * t21;
t62 = t35 * t71;
t70 = qJD(5) * t35;
t2 = t17 * t70 + t62 - qJD(5) * t47 - t37 * (t83 * qJD(3) + t48);
t5 = t37 * t17 + t35 * t55;
t3 = -t35 * t48 - t5 * qJD(5) + (-t35 * t63 + (-t35 * qJ(4) - t37 * t77) * t36) * qJD(3);
t4 = -t35 * t17 + t47;
t44 = t35 * t4 - t37 * t5;
t1 = -t44 * qJD(5) - t2 * t35 + t3 * t37;
t80 = 0.2e1 * qJD(2);
t79 = 0.2e1 * qJD(4);
t75 = t31 + t33;
t69 = qJD(5) * t37;
t68 = qJD(5) * t38;
t67 = qJD(5) * t78;
t26 = t36 * qJD(3);
t27 = t38 * qJD(3);
t65 = qJ(2) * qJD(3);
t64 = qJ(4) * qJD(5);
t61 = t35 * t68;
t60 = t37 * t68;
t59 = t40 * t26;
t58 = t37 * t27;
t57 = t35 * t69;
t56 = t36 * t27;
t54 = t75 * t36;
t51 = qJD(5) * (t32 + t34);
t22 = 0.2e1 * t56;
t50 = t35 * t58;
t49 = t32 * t57;
t45 = t35 * t5 + t37 * t4;
t25 = t40 * t27;
t18 = -pkin(4) * t27 + t25;
t41 = t83 * qJD(5) + t18;
t9 = t82 * qJD(3) + t66;
t30 = qJ(2) * t80;
t29 = qJ(4) * t79;
t23 = -0.2e1 * t56;
t19 = 0.2e1 * t53;
t15 = -t35 * t26 + t60;
t14 = t35 * t27 + t36 * t69;
t13 = t37 * t26 + t61;
t12 = -t36 * t70 + t58;
t11 = qJD(3) * t54;
t8 = (pkin(3) * t38 + t73) * qJD(3) + t48;
t7 = -t36 * t52 + t50;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t30, t23, t19, 0, t22, 0, 0, 0.2e1 * qJD(2) * t36 + 0.2e1 * t38 * t65, 0.2e1 * qJD(2) * t38 - 0.2e1 * t36 * t65, 0, t30, 0, 0, 0, t23, t19, t22, 0, -0.2e1 * t20 * t27 - 0.2e1 * t8 * t36, 0.2e1 * t20 * t26 - 0.2e1 * t8 * t38, 0.2e1 * t20 * t8, 0.2e1 * t31 * t56 + 0.2e1 * t49, -0.2e1 * t32 * t52 + 0.4e1 * t36 * t50, -0.2e1 * t35 * t53 + 0.2e1 * t36 * t60, 0.2e1 * t33 * t56 - 0.2e1 * t49, -0.2e1 * t36 * t61 - 0.2e1 * t37 * t53, t23, 0.2e1 * (t37 * t71 + t3) * t38 + 0.2e1 * (-qJD(3) * t4 - t18 * t37 - t21 * t70) * t36, 0.2e1 * (t2 - t62) * t38 + 0.2e1 * (qJD(3) * t5 + t18 * t35 - t21 * t69) * t36, -0.2e1 * t44 * t27 + 0.2e1 * (-t45 * qJD(5) - t2 * t37 - t3 * t35) * t36, -0.2e1 * t21 * t18 - 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t51, t37 * t51, 0, (t45 * qJD(3) + t18) * t36 + (-t1 - t71) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (0.1e1 - t75) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, 0, -t59, -t25, 0, 0, 0, t26, t27, 0, 0, 0, -t9, t59, t25, t9 * t40, t7, -t76 * t27 - 0.4e1 * t36 * t57, -t13, -t7, -t15, 0, t41 * t35 - t81 * t37, t81 * t35 + t41 * t37, -t1, t18 * qJ(4) - t21 * qJD(4) - t1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, t9, 0, 0, 0, 0, 0, 0, t14, t12, -t11, t66 + (-t54 * t78 + t72) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t29, -0.2e1 * t57, 0.2e1 * t52, 0, 0.2e1 * t57, 0, 0, 0.2e1 * qJD(4) * t35 + 0.2e1 * t37 * t64, 0.2e1 * qJD(4) * t37 - 0.2e1 * t35 * t64, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, t59, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t12, -t26, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, -t69, 0, t35 * t67, t37 * t67, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;

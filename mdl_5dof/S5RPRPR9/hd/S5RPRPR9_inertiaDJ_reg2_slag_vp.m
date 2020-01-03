% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:37
% EndTime: 2019-12-31 18:24:40
% DurationCPUTime: 0.74s
% Computational Cost: add. (520->88), mult. (1077->168), div. (0->0), fcn. (782->6), ass. (0->73)
t34 = sin(qJ(3));
t78 = pkin(3) + pkin(7);
t61 = t78 * t34;
t36 = cos(qJ(3));
t72 = qJ(4) * t36;
t82 = t61 - t72;
t81 = t78 * t36;
t31 = t36 ^ 2;
t50 = qJD(3) * (t34 ^ 2 - t31);
t33 = sin(qJ(5));
t28 = t33 ^ 2;
t35 = cos(qJ(5));
t30 = t35 ^ 2;
t75 = t28 - t30;
t49 = t75 * qJD(5);
t23 = sin(pkin(8)) * pkin(1) + pkin(6);
t76 = pkin(4) + t23;
t19 = t76 * t36;
t64 = t36 * qJD(4);
t71 = t34 * qJ(4);
t80 = (t71 + t81) * qJD(3) - qJD(5) * t19 - t64;
t52 = t76 * t34;
t46 = t35 * t52;
t70 = qJD(3) * t19;
t60 = t33 * t70;
t65 = t34 * qJD(4);
t69 = qJD(5) * t33;
t24 = -cos(pkin(8)) * pkin(1) - pkin(2);
t40 = t24 - t71;
t9 = t40 - t81;
t2 = t9 * t69 - qJD(5) * t46 - t35 * (t82 * qJD(3) - t65) - t60;
t68 = qJD(5) * t35;
t3 = -t9 * t68 + (-qJD(5) * t76 + qJD(4)) * t33 * t34 + (-t33 * t61 + (t33 * qJ(4) + t35 * t76) * t36) * qJD(3);
t4 = -t33 * t9 + t46;
t5 = t33 * t52 + t35 * t9;
t43 = t33 * t4 - t35 * t5;
t1 = -t43 * qJD(5) - t2 * t33 + t3 * t35;
t79 = 0.2e1 * qJD(4);
t77 = t36 * pkin(3);
t74 = t28 + t30;
t67 = qJD(5) * t36;
t66 = qJD(5) * t78;
t25 = t34 * qJD(3);
t26 = t36 * qJD(3);
t63 = qJ(4) * qJD(5);
t62 = 0.2e1 * qJD(3) * t24;
t59 = t33 * t67;
t58 = t35 * t67;
t57 = t34 * t26;
t56 = t23 * t25;
t55 = t33 * t68;
t54 = t35 * t25;
t53 = t23 * t26;
t51 = t74 * t34;
t21 = 0.2e1 * t57;
t48 = t31 * t55;
t47 = t33 * t54;
t44 = t33 * t5 + t35 * t4;
t11 = t76 * t25;
t39 = t82 * qJD(5) - t11;
t38 = t64 + (-t71 - t77) * qJD(3);
t27 = qJ(4) * t79;
t22 = -0.2e1 * t57;
t20 = -0.2e1 * t50;
t18 = t40 - t77;
t17 = -t33 * t25 + t58;
t16 = t33 * t26 + t34 * t68;
t15 = t54 + t59;
t14 = t35 * t26 - t34 * t69;
t13 = qJD(3) * t51;
t12 = t65 + (-pkin(3) * t34 + t72) * qJD(3);
t7 = -t36 * t49 - t47;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20, 0, t22, 0, 0, t34 * t62, t36 * t62, 0, 0, 0, 0, 0, t21, t20, t22, 0, -0.2e1 * t12 * t36 - 0.2e1 * t18 * t25, 0.2e1 * t12 * t34 - 0.2e1 * t18 * t26, -0.2e1 * t18 * t12, -0.2e1 * t28 * t57 + 0.2e1 * t48, -0.2e1 * t31 * t49 - 0.4e1 * t36 * t47, 0.2e1 * t33 * t50 - 0.2e1 * t34 * t58, -0.2e1 * t30 * t57 - 0.2e1 * t48, 0.2e1 * t34 * t59 + 0.2e1 * t35 * t50, t21, 0.2e1 * (-t35 * t70 + t3) * t34 + 0.2e1 * (qJD(3) * t4 - t11 * t35 - t19 * t69) * t36, 0.2e1 * (t2 + t60) * t34 + 0.2e1 * (-qJD(3) * t5 + t11 * t33 - t19 * t68) * t36, -0.2e1 * t43 * t25 + 0.2e1 * (t44 * qJD(5) + t2 * t35 + t3 * t33) * t36, -0.2e1 * t19 * t11 - 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t44 * qJD(3) - t11) * t34 + (-t1 + t70) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (0.1e1 - t74) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t25, 0, -t53, t56, 0, 0, 0, -t26, t25, 0, 0, 0, t38, t53, -t56, t38 * t23, -t7, -t75 * t25 + 0.4e1 * t36 * t55, t14, t7, -t16, 0, t39 * t33 - t80 * t35, t80 * t33 + t39 * t35, -t1, -t11 * qJ(4) + t19 * qJD(4) - t1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, t12, 0, 0, 0, 0, 0, 0, t16, t14, -t13, t65 + (-t51 * t78 + t72) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t27, -0.2e1 * t55, 0.2e1 * t49, 0, 0.2e1 * t55, 0, 0, 0.2e1 * qJD(4) * t33 + 0.2e1 * t35 * t63, 0.2e1 * qJD(4) * t35 - 0.2e1 * t33 * t63, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, t53, 0, 0, 0, 0, 0, 0, t14, -t16, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, t15, t26, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, -t68, 0, t33 * t66, t35 * t66, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t68, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;

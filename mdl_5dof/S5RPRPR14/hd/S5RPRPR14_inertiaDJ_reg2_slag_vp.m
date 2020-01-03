% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR14_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:19
% DurationCPUTime: 0.75s
% Computational Cost: add. (1091->111), mult. (2153->208), div. (0->0), fcn. (1922->6), ass. (0->79)
t33 = sin(pkin(8));
t37 = cos(qJ(3));
t73 = cos(pkin(8));
t59 = t73 * t37;
t35 = sin(qJ(3));
t70 = t35 * qJD(3);
t17 = -qJD(3) * t59 + t33 * t70;
t21 = t33 * t37 + t73 * t35;
t18 = t21 * qJD(3);
t88 = (t17 * t33 + t73 * t18) * pkin(3);
t34 = sin(qJ(5));
t31 = t34 ^ 2;
t36 = cos(qJ(5));
t32 = t36 ^ 2;
t76 = t31 - t32;
t58 = qJD(5) * t76;
t38 = -pkin(1) - pkin(6);
t74 = qJ(4) - t38;
t23 = t74 * t35;
t60 = t74 * t37;
t13 = -t73 * t23 - t33 * t60;
t16 = -qJD(3) * t60 - t35 * qJD(4);
t44 = t37 * qJD(4) - t74 * t70;
t40 = t73 * t16 - t33 * t44;
t22 = -t33 * t35 + t59;
t29 = t35 * pkin(3) + qJ(2);
t45 = t21 * pkin(4) - t22 * pkin(7) + t29;
t42 = t36 * t45;
t69 = t37 * qJD(3);
t24 = pkin(3) * t69 + qJD(2);
t43 = -t17 * pkin(4) + t18 * pkin(7) + t24;
t72 = qJD(5) * t34;
t1 = -qJD(5) * t42 + t13 * t72 - t34 * t43 - t36 * t40;
t4 = t36 * t13 + t34 * t45;
t2 = -qJD(5) * t4 - t34 * t40 + t36 * t43;
t3 = -t34 * t13 + t42;
t54 = t3 * t34 - t36 * t4;
t87 = t54 * qJD(5) + t1 * t34 - t2 * t36;
t20 = t22 ^ 2;
t86 = 0.2e1 * qJD(2);
t12 = -t33 * t23 + t74 * t59;
t8 = t33 * t16 + t73 * t44;
t85 = t12 * t8;
t84 = t8 * t34;
t83 = t8 * t36;
t28 = -t73 * pkin(3) - pkin(4);
t82 = t18 * t28;
t81 = t18 * t36;
t80 = t21 * t17;
t79 = t22 * t18;
t78 = t34 * t17;
t77 = t36 * t17;
t75 = t31 + t32;
t71 = qJD(5) * t36;
t68 = qJ(2) * qJD(3);
t67 = -0.2e1 * t80;
t66 = t34 * t81;
t65 = 0.2e1 * qJD(5) * t28;
t64 = t34 * t71;
t63 = t35 * t69;
t62 = t21 ^ 2 + t20;
t61 = t75 * t17;
t57 = t20 * t64;
t55 = t3 * t36 + t34 * t4;
t53 = t12 * t18 - t8 * t22;
t52 = t79 + t80;
t27 = t33 * pkin(3) + pkin(7);
t51 = t17 * t27 - t82;
t50 = t21 * t27 - t22 * t28;
t11 = t21 * t71 - t78;
t49 = -t18 * t34 + t22 * t71;
t48 = -t22 * t72 - t81;
t46 = 0.2e1 * t52;
t41 = -t55 * qJD(5) - t1 * t36 - t2 * t34;
t39 = -t13 * t17 + t40 * t21 + t53;
t30 = qJ(2) * t86;
t10 = t21 * t72 + t77;
t5 = t22 * t58 + t66;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t30, -0.2e1 * t63, 0.2e1 * (t35 ^ 2 - t37 ^ 2) * qJD(3), 0, 0.2e1 * t63, 0, 0, 0.2e1 * qJD(2) * t35 + 0.2e1 * t37 * t68, 0.2e1 * qJD(2) * t37 - 0.2e1 * t35 * t68, 0, t30, -0.2e1 * t79, 0.2e1 * t22 * t17 + 0.2e1 * t18 * t21, 0, t67, 0, 0, -0.2e1 * t29 * t17 + 0.2e1 * t24 * t21, -0.2e1 * t29 * t18 + 0.2e1 * t24 * t22, -0.2e1 * t39, 0.2e1 * t13 * t40 + 0.2e1 * t29 * t24 + 0.2e1 * t85, -0.2e1 * t32 * t79 - 0.2e1 * t57, 0.2e1 * t20 * t58 + 0.4e1 * t22 * t66, 0.2e1 * t48 * t21 - 0.2e1 * t22 * t77, -0.2e1 * t31 * t79 + 0.2e1 * t57, -0.2e1 * t49 * t21 + 0.2e1 * t22 * t78, t67, 0.2e1 * t49 * t12 - 0.2e1 * t3 * t17 + 0.2e1 * t2 * t21 + 0.2e1 * t22 * t84, 0.2e1 * t1 * t21 + 0.2e1 * t48 * t12 + 0.2e1 * t4 * t17 + 0.2e1 * t22 * t83, 0.2e1 * t55 * t18 + 0.2e1 * t87 * t22, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t39, 0, 0, 0, 0, 0, 0, t34 * t46 - t62 * t71, t36 * t46 + t62 * t72, 0, t54 * t17 + t41 * t21 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21 * t61 - 0.2e1 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, -t69, 0, -t38 * t70, -t38 * t69, 0, 0, 0, 0, -t18, 0, t17, 0, -t8, -t40, t88, (t40 * t33 - t8 * t73) * pkin(3), -t5, t76 * t18 - 0.4e1 * t22 * t64, t11, t5, -t10, 0, -t83 + t51 * t34 + (t12 * t34 - t50 * t36) * qJD(5), t84 + t51 * t36 + (t12 * t36 + t50 * t34) * qJD(5), t41, t41 * t27 + t8 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, -t88, 0, 0, 0, 0, 0, 0, t48, -t49, -t61, -t27 * t61 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t64, -0.2e1 * t58, 0, -0.2e1 * t64, 0, 0, t34 * t65, t36 * t65, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t24, 0, 0, 0, 0, 0, 0, -t10, -t11, t75 * t18, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t49, -t17, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, -t72, 0, -t27 * t71, t27 * t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;

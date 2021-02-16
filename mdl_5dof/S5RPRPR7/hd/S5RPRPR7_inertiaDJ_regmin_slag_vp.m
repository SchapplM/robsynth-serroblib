% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR7
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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:27
% EndTime: 2021-01-15 12:06:29
% DurationCPUTime: 0.41s
% Computational Cost: add. (485->86), mult. (1127->172), div. (0->0), fcn. (1001->8), ass. (0->63)
t72 = 2 * qJD(5);
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t29 = sin(pkin(8)) * pkin(1) + pkin(6);
t61 = qJ(4) + t29;
t47 = qJD(3) * t61;
t16 = t40 * qJD(4) - t38 * t47;
t35 = sin(pkin(9));
t41 = -t38 * qJD(4) - t40 * t47;
t60 = cos(pkin(9));
t3 = t35 * t16 - t60 * t41;
t37 = sin(qJ(5));
t71 = t3 * t37;
t50 = t60 * t38;
t24 = t35 * t40 + t50;
t19 = t24 * qJD(3);
t49 = t60 * t40;
t67 = t35 * t38;
t23 = -t49 + t67;
t70 = t23 * t19;
t57 = t38 * qJD(3);
t20 = qJD(3) * t49 - t35 * t57;
t69 = t24 * t20;
t39 = cos(qJ(5));
t68 = t24 * t39;
t66 = t37 * t19;
t65 = t39 * t19;
t64 = t39 * t20;
t63 = t23 * t64 + t24 * t65;
t34 = t39 ^ 2;
t62 = t37 ^ 2 - t34;
t59 = qJD(5) * t37;
t58 = qJD(5) * t39;
t56 = t40 * qJD(3);
t30 = -t60 * pkin(3) - pkin(4);
t55 = t30 * t72;
t54 = 0.2e1 * t56;
t32 = pkin(3) * t57;
t53 = t24 * t59;
t52 = t37 * t58;
t31 = -cos(pkin(8)) * pkin(1) - pkin(2);
t51 = -0.4e1 * t37 * t68;
t48 = t62 * qJD(5);
t21 = t61 * t40;
t12 = t60 * t21 - t61 * t67;
t25 = -t40 * pkin(3) + t31;
t6 = t23 * pkin(4) - t24 * pkin(7) + t25;
t46 = t39 * t12 + t37 * t6;
t45 = t37 * t12 - t39 * t6;
t28 = t35 * pkin(3) + pkin(7);
t44 = -t19 * t28 + t20 * t30;
t43 = t23 * t28 - t24 * t30;
t9 = t23 * t58 + t66;
t42 = t37 * t20 + t24 * t58;
t8 = t53 - t64;
t22 = t24 ^ 2;
t11 = t35 * t21 + t61 * t50;
t7 = t23 * t59 - t65;
t5 = t19 * pkin(4) - t20 * pkin(7) + t32;
t4 = t60 * t16 + t35 * t41;
t2 = -t46 * qJD(5) - t37 * t4 + t39 * t5;
t1 = t45 * qJD(5) - t37 * t5 - t39 * t4;
t10 = [0, 0, 0, 0, t38 * t54, 0.2e1 * (-t38 ^ 2 + t40 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t31 * t57, t31 * t54, 0.2e1 * t25 * t19 + 0.2e1 * t23 * t32, 0.2e1 * t25 * t20 + 0.2e1 * t24 * t32, 0.2e1 * t11 * t20 - 0.2e1 * t12 * t19 - 0.2e1 * t4 * t23 + 0.2e1 * t3 * t24, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t4 + 0.2e1 * t25 * t32, -0.2e1 * t22 * t52 + 0.2e1 * t34 * t69, t62 * t22 * t72 + t20 * t51, -0.2e1 * t23 * t53 + 0.2e1 * t63, -0.2e1 * t42 * t23 - 0.2e1 * t24 * t66, 0.2e1 * t70, 0.2e1 * t42 * t11 - 0.2e1 * t45 * t19 + 0.2e1 * t2 * t23 + 0.2e1 * t24 * t71, 0.2e1 * t1 * t23 - 0.2e1 * t8 * t11 - 0.2e1 * t46 * t19 + 0.2e1 * t3 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t19 + t12 * t20 + t3 * t23 + t4 * t24, 0, 0, 0, 0, 0, 0, (-t24 * t19 - t20 * t23) * t39 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t69 + 0.2e1 * t70, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t56, -t57, 0, -t29 * t56, t29 * t57, -t3, -t4, (-t19 * t35 - t60 * t20) * pkin(3), (-t60 * t3 + t35 * t4) * pkin(3), -t24 * t48 + t37 * t64, qJD(5) * t51 - t62 * t20, t9, -t7, 0, -t3 * t39 + t44 * t37 + (t11 * t37 - t43 * t39) * qJD(5), t71 + t44 * t39 + (t11 * t39 + t43 * t37) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, -t19, -t20, 0, (-t60 * t19 + t20 * t35) * pkin(3), 0, 0, 0, 0, 0, t7, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52, -0.2e1 * t48, 0, 0, 0, t37 * t55, t39 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, t32, 0, 0, 0, 0, 0, -t7, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t42, t19, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t59, 0, -t28 * t58, t28 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;

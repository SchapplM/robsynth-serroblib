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
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:52
% DurationCPUTime: 0.38s
% Computational Cost: add. (457->80), mult. (1055->164), div. (0->0), fcn. (943->8), ass. (0->63)
t72 = 2 * qJD(5);
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t28 = sin(pkin(8)) * pkin(1) + pkin(6);
t61 = qJ(4) + t28;
t47 = qJD(3) * t61;
t16 = t39 * qJD(4) - t37 * t47;
t34 = sin(pkin(9));
t40 = -t37 * qJD(4) - t39 * t47;
t60 = cos(pkin(9));
t3 = t34 * t16 - t60 * t40;
t36 = sin(qJ(5));
t71 = t3 * t36;
t50 = t60 * t37;
t24 = t34 * t39 + t50;
t19 = t24 * qJD(3);
t49 = t60 * t39;
t67 = t34 * t37;
t23 = -t49 + t67;
t70 = t23 * t19;
t57 = t37 * qJD(3);
t20 = qJD(3) * t49 - t34 * t57;
t69 = t24 * t20;
t38 = cos(qJ(5));
t68 = t24 * t38;
t66 = t36 * t19;
t65 = t38 * t19;
t64 = t38 * t20;
t63 = t23 * t64 + t24 * t65;
t33 = t38 ^ 2;
t62 = t36 ^ 2 - t33;
t59 = qJD(5) * t36;
t58 = qJD(5) * t38;
t56 = t39 * qJD(3);
t29 = -t60 * pkin(3) - pkin(4);
t55 = t29 * t72;
t54 = 0.2e1 * t56;
t31 = pkin(3) * t57;
t53 = t24 * t59;
t52 = t36 * t58;
t30 = -cos(pkin(8)) * pkin(1) - pkin(2);
t51 = -0.4e1 * t36 * t68;
t48 = t62 * qJD(5);
t21 = t61 * t39;
t12 = t60 * t21 - t61 * t67;
t42 = -t39 * pkin(3) + t30;
t6 = t23 * pkin(4) - t24 * pkin(7) + t42;
t46 = t38 * t12 + t36 * t6;
t45 = t36 * t12 - t38 * t6;
t27 = t34 * pkin(3) + pkin(7);
t44 = -t19 * t27 + t20 * t29;
t43 = t23 * t27 - t24 * t29;
t9 = t23 * t58 + t66;
t41 = t36 * t20 + t24 * t58;
t8 = t53 - t64;
t22 = t24 ^ 2;
t11 = t34 * t21 + t61 * t50;
t7 = t23 * t59 - t65;
t5 = t19 * pkin(4) - t20 * pkin(7) + t31;
t4 = t60 * t16 + t34 * t40;
t2 = -t46 * qJD(5) - t36 * t4 + t38 * t5;
t1 = t45 * qJD(5) - t36 * t5 - t38 * t4;
t10 = [0, 0, 0, 0, t37 * t54, 0.2e1 * (-t37 ^ 2 + t39 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t30 * t57, t30 * t54, 0.2e1 * t11 * t20 - 0.2e1 * t12 * t19 - 0.2e1 * t4 * t23 + 0.2e1 * t3 * t24, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t4 + 0.2e1 * t42 * t31, -0.2e1 * t22 * t52 + 0.2e1 * t33 * t69, t62 * t22 * t72 + t20 * t51, -0.2e1 * t23 * t53 + 0.2e1 * t63, -0.2e1 * t23 * t41 - 0.2e1 * t24 * t66, 0.2e1 * t70, 0.2e1 * t41 * t11 - 0.2e1 * t45 * t19 + 0.2e1 * t2 * t23 + 0.2e1 * t24 * t71, 0.2e1 * t1 * t23 - 0.2e1 * t8 * t11 - 0.2e1 * t46 * t19 + 0.2e1 * t3 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t19 + t12 * t20 + t3 * t23 + t4 * t24, 0, 0, 0, 0, 0, 0, (-t24 * t19 - t20 * t23) * t38 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t69 + 0.2e1 * t70, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t56, -t57, 0, -t28 * t56, t28 * t57, (-t19 * t34 - t60 * t20) * pkin(3), (-t60 * t3 + t34 * t4) * pkin(3), -t24 * t48 + t36 * t64, qJD(5) * t51 - t62 * t20, t9, -t7, 0, -t3 * t38 + t44 * t36 + (t11 * t36 - t43 * t38) * qJD(5), t71 + t44 * t38 + (t11 * t38 + t43 * t36) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, 0, (-t60 * t19 + t20 * t34) * pkin(3), 0, 0, 0, 0, 0, t7, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52, -0.2e1 * t48, 0, 0, 0, t36 * t55, t38 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, -t7, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t41, t19, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t59, 0, -t27 * t58, t27 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;

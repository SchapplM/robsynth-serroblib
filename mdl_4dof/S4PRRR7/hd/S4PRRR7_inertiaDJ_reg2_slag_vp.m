% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:50
% DurationCPUTime: 0.59s
% Computational Cost: add. (378->96), mult. (1197->211), div. (0->0), fcn. (1050->8), ass. (0->74)
t25 = sin(qJ(3));
t78 = -0.4e1 * t25;
t24 = sin(qJ(4));
t18 = t24 ^ 2;
t27 = cos(qJ(4));
t20 = t27 ^ 2;
t71 = t18 - t20;
t48 = qJD(4) * t71;
t23 = cos(pkin(4));
t28 = cos(qJ(3));
t22 = sin(pkin(4));
t26 = sin(qJ(2));
t74 = t22 * t26;
t12 = -t23 * t28 + t25 * t74;
t13 = t23 * t25 + t28 * t74;
t29 = cos(qJ(2));
t68 = qJD(2) * t29;
t52 = t22 * t68;
t6 = qJD(3) * t13 + t25 * t52;
t77 = t12 * t6;
t76 = t24 * pkin(6);
t75 = t6 * t25;
t73 = t22 * t29;
t72 = t27 * t28;
t19 = t25 ^ 2;
t70 = -t28 ^ 2 + t19;
t69 = qJD(2) * t26;
t67 = qJD(3) * t12;
t66 = qJD(3) * t27;
t65 = qJD(4) * t24;
t64 = qJD(4) * t27;
t63 = qJD(4) * t28;
t62 = t25 * qJD(3);
t61 = t28 * qJD(3);
t60 = t28 * t76;
t59 = pkin(6) * t72;
t58 = -0.2e1 * pkin(2) * qJD(3);
t57 = -0.2e1 * pkin(3) * qJD(4);
t56 = pkin(6) * t61;
t55 = t24 * t63;
t54 = t27 * t63;
t53 = t22 * t69;
t51 = t24 * t64;
t50 = t25 * t61;
t49 = t27 * t61;
t47 = t70 * qJD(3);
t46 = 0.2e1 * t50;
t45 = t24 * t49;
t44 = t19 * t51;
t43 = -t28 * pkin(3) - t25 * pkin(7);
t42 = pkin(3) * t25 - pkin(7) * t28;
t38 = -t13 * t27 + t24 * t73;
t7 = -t13 * t24 - t27 * t73;
t41 = t24 * t38 - t27 * t7;
t39 = pkin(2) - t43;
t34 = t27 * t39;
t10 = -t34 - t60;
t11 = -t24 * t39 + t59;
t40 = -t10 * t27 - t11 * t24;
t37 = t12 * t64 + t6 * t24;
t36 = t12 * t65 - t6 * t27;
t35 = t42 * t24;
t33 = t27 * t62 + t55;
t5 = t28 * t52 - t67;
t1 = t38 * qJD(4) - t5 * t24 + t27 * t53;
t2 = qJD(4) * t7 + t24 * t53 + t5 * t27;
t32 = t41 * qJD(4) - t1 * t24 + t2 * t27;
t3 = pkin(6) * t33 - qJD(3) * t35 + qJD(4) * t34;
t4 = -qJD(4) * t11 + (t25 * t76 + t27 * t42) * qJD(3);
t31 = t40 * qJD(4) - t4 * t24 - t3 * t27;
t30 = t75 + t5 * t28 + (t12 * t28 - t13 * t25) * qJD(3);
t16 = -0.2e1 * t50;
t9 = t25 * t48 - t45;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t22 ^ 2 * t26 * t68 + 0.2e1 * t13 * t5 + 0.2e1 * t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t7 * t1 - 0.2e1 * t2 * t38 + 0.2e1 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, 0, 0, 0, 0, 0, 0, 0, (-t28 * t69 - t29 * t62) * t22, (t25 * t69 - t29 * t61) * t22, t30, -pkin(2) * t53 + pkin(6) * t30, 0, 0, 0, 0, 0, 0, (t24 * t67 - t1) * t28 + (qJD(3) * t7 + t37) * t25, (t12 * t66 + t2) * t28 + (qJD(3) * t38 - t36) * t25, t41 * t61 + (-t1 * t27 - t2 * t24 + (t24 * t7 + t27 * t38) * qJD(4)) * t25, t1 * t10 + t2 * t11 + t38 * t3 + t7 * t4 + (t12 * t61 + t75) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -0.2e1 * t47, 0, t16, 0, 0, t25 * t58, t28 * t58, 0, 0, 0.2e1 * t20 * t50 - 0.2e1 * t44, 0.2e1 * t19 * t48 + t45 * t78, 0.2e1 * t25 * t55 + 0.2e1 * t70 * t66, 0.2e1 * t18 * t50 + 0.2e1 * t44, -0.2e1 * t24 * t47 + 0.2e1 * t25 * t54, t16, 0.2e1 * t10 * t62 - 0.2e1 * t4 * t28 + 0.2e1 * (t19 * t64 + t24 * t46) * pkin(6), -0.2e1 * t11 * t62 - 0.2e1 * t3 * t28 + 0.2e1 * (-t19 * t65 + t27 * t46) * pkin(6), 0.2e1 * t40 * t61 + 0.2e1 * (t24 * t3 - t27 * t4 + (t10 * t24 - t11 * t27) * qJD(4)) * t25, 0.2e1 * pkin(6) ^ 2 * t50 + 0.2e1 * t10 * t4 - 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t5, 0, 0, 0, 0, 0, 0, 0, 0, t36, t37, t32, -t6 * pkin(3) + pkin(7) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t62, 0, -t56, pkin(6) * t62, 0, 0, -t9, t51 * t78 - t71 * t61, t24 * t62 - t54, t9, t33, 0, (pkin(7) * t72 + (-t27 * pkin(3) + t76) * t25) * qJD(4) + (t24 * t43 - t59) * qJD(3), (pkin(6) * t25 * t27 + t35) * qJD(4) + (t27 * t43 + t60) * qJD(3), t31, -pkin(3) * t56 + pkin(7) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t51, -0.2e1 * t48, 0, -0.2e1 * t51, 0, 0, t24 * t57, t27 * t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t65 + t49, 0, -t24 * t61 - t25 * t64, t62, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t65, 0, -pkin(7) * t64, pkin(7) * t65, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;

% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:21
% EndTime: 2020-01-03 11:52:24
% DurationCPUTime: 0.71s
% Computational Cost: add. (611->64), mult. (1416->108), div. (0->0), fcn. (1027->8), ass. (0->55)
t48 = cos(pkin(9)) * pkin(1) + pkin(2);
t50 = sin(pkin(9)) * pkin(1);
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t19 = t64 * t48 + t66 * t50;
t36 = sin(qJ(4));
t45 = t66 * t48;
t46 = t64 * t50;
t41 = t46 - t45;
t40 = pkin(3) - t41;
t65 = cos(qJ(4));
t38 = t65 * t40;
t39 = t19 * qJD(3);
t24 = qJD(3) * t46;
t43 = -qJD(3) * t45 + t24;
t58 = -qJD(4) * t38 + t36 * t39 + t65 * t43;
t60 = qJD(4) * t36;
t5 = t19 * t60 + t58;
t35 = sin(qJ(5));
t33 = t35 ^ 2;
t37 = cos(qJ(5));
t34 = t37 ^ 2;
t68 = t33 + t34;
t51 = t68 * t5;
t32 = t37 * qJD(5);
t18 = t65 * t19;
t12 = t36 * t40 + t18;
t42 = t36 * t43 - t65 * t39;
t6 = t12 * qJD(4) - t42;
t63 = t36 * t19;
t11 = t38 - t63;
t9 = -pkin(4) - t11;
t67 = t9 * t32 + t6 * t35;
t54 = t65 * pkin(3);
t49 = qJD(4) * t54;
t20 = t68 * t49;
t31 = -t54 - pkin(4);
t55 = pkin(3) * t60;
t62 = t31 * t32 + t35 * t55;
t61 = qJD(4) * pkin(3);
t59 = t35 * qJD(5);
t57 = pkin(4) * t59;
t56 = pkin(4) * t32;
t53 = t35 * t32;
t7 = t9 * t59;
t52 = -t6 * t37 + t7;
t22 = t31 * t59;
t47 = -t37 * t55 + t22;
t44 = t68 * t65;
t30 = t36 * pkin(3) + pkin(8);
t29 = -0.2e1 * t53;
t28 = 0.2e1 * t53;
t21 = 0.2e1 * (-t33 + t34) * qJD(5);
t10 = pkin(8) + t12;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t39, 0.2e1 * t43, 0, 0.2e1 * (-t24 + (t41 + t45) * qJD(3)) * t19, 0, 0, 0, 0, 0, 0, -0.2e1 * t6, 0.2e1 * t5, 0, -0.2e1 * t11 * t6 - 0.2e1 * t12 * t5, t28, t21, 0, t29, 0, 0, 0.2e1 * t52, 0.2e1 * t67, -0.2e1 * t51, -0.2e1 * t10 * t51 + 0.2e1 * t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t43, 0, 0, 0, 0, 0, 0, 0, 0, (-t18 + (-0.2e1 * pkin(3) + t41) * t36) * qJD(4) + t42, (-t54 + t63) * qJD(4) + t58, 0, (-t65 * t6 - t36 * t5 + (-t11 * t36 + t65 * t12) * qJD(4)) * pkin(3), t28, t21, 0, t29, 0, 0, t22 + t7 + (-t6 - t55) * t37, t62 + t67, t20 - t51, t6 * t31 - t30 * t51 + (t44 * t10 + t36 * t9) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, -0.2e1 * t49, 0, 0, t28, t21, 0, t29, 0, 0, 0.2e1 * t47, 0.2e1 * t62, 0.2e1 * t20, 0.2e1 * (t44 * t30 + t31 * t36) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, t28, t21, 0, t29, 0, 0, t52 - t57, -t56 + t67, -t51, -t6 * pkin(4) - pkin(8) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t49, 0, 0, t28, t21, 0, t29, 0, 0, t47 - t57, -t56 + t62, t20, (-pkin(4) * t36 + t44 * pkin(8)) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t21, 0, t29, 0, 0, -0.2e1 * t57, -0.2e1 * t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t59, 0, -t10 * t32 + t35 * t5, t10 * t59 + t37 * t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t59, 0, -t30 * t32 - t35 * t49, t30 * t59 - t37 * t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t59, 0, -pkin(8) * t32, pkin(8) * t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;

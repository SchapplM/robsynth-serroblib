% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:08
% EndTime: 2019-12-31 17:58:09
% DurationCPUTime: 0.50s
% Computational Cost: add. (402->55), mult. (751->106), div. (0->0), fcn. (830->8), ass. (0->55)
t32 = cos(pkin(9));
t35 = sin(qJ(4));
t30 = sin(pkin(9));
t53 = cos(qJ(4));
t42 = t53 * t30;
t20 = t35 * t32 + t42;
t63 = -0.2e1 * t20;
t31 = sin(pkin(8));
t57 = t31 * pkin(1);
t23 = qJ(3) + t57;
t54 = pkin(6) + t23;
t14 = t54 * t32;
t4 = t35 * t14 + t54 * t42;
t62 = t4 ^ 2;
t49 = t35 * t30;
t18 = -t53 * t32 + t49;
t61 = t18 ^ 2;
t33 = cos(pkin(8));
t56 = t33 * pkin(1);
t25 = -pkin(2) - t56;
t21 = -t32 * pkin(3) + t25;
t60 = 0.2e1 * t21;
t59 = 0.2e1 * t30;
t58 = t18 * pkin(4);
t55 = t4 * t18;
t34 = sin(qJ(5));
t28 = t34 ^ 2;
t52 = t28 * t20;
t10 = t34 * t18;
t51 = t34 * t20;
t36 = cos(qJ(5));
t50 = t34 * t36;
t48 = t36 * t20;
t26 = t30 ^ 2;
t27 = t32 ^ 2;
t47 = t26 + t27;
t29 = t36 ^ 2;
t46 = t28 + t29;
t45 = t18 * t63;
t44 = t34 * t48;
t43 = t46 * pkin(7);
t41 = -pkin(4) * t20 - pkin(7) * t18;
t3 = -t20 * pkin(7) + t21 + t58;
t6 = t53 * t14 - t54 * t49;
t1 = t36 * t3 - t34 * t6;
t2 = t34 * t3 + t36 * t6;
t40 = t1 * t36 + t2 * t34;
t39 = -t1 * t34 + t2 * t36;
t17 = t20 ^ 2;
t13 = t36 * t18;
t12 = t29 * t20;
t11 = t29 * t17;
t9 = t28 * t17;
t7 = -t12 - t52;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t57, 0, (t31 ^ 2 + t33 ^ 2) * pkin(1) ^ 2, t26, t32 * t59, 0, t27, 0, 0, -0.2e1 * t25 * t32, t25 * t59, 0.2e1 * t47 * t23, t47 * t23 ^ 2 + t25 ^ 2, t17, t45, 0, t61, 0, 0, t18 * t60, t20 * t60, -0.2e1 * t6 * t18 + 0.2e1 * t4 * t20, t21 ^ 2 + t6 ^ 2 + t62, t11, -0.2e1 * t17 * t50, 0.2e1 * t18 * t48, t9, t34 * t45, t61, 0.2e1 * t1 * t18 + 0.2e1 * t4 * t51, -0.2e1 * t2 * t18 + 0.2e1 * t4 * t48, t40 * t63, t1 ^ 2 + t2 ^ 2 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t20 + t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t39 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 + t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t9 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, 0, t25, 0, 0, 0, 0, 0, 0, t18, t20, 0, t21, 0, 0, 0, 0, 0, 0, t13, -t10, t7, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, -t4, -t6, 0, 0, t44, t12 - t52, t10, -t44, t13, 0, t34 * t41 - t4 * t36, t4 * t34 + t36 * t41, t39, -t4 * pkin(4) + pkin(7) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t10, -t7, t20 * t43 - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, 0.2e1 * t50, 0, t29, 0, 0, 0.2e1 * pkin(4) * t36, -0.2e1 * pkin(4) * t34, 0.2e1 * t43, t46 * pkin(7) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t51, t18, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t36, 0, -t34 * pkin(7), -t36 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;

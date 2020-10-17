% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaJ_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:23
% EndTime: 2019-12-05 17:03:25
% DurationCPUTime: 0.59s
% Computational Cost: add. (149->47), mult. (508->107), div. (0->0), fcn. (619->8), ass. (0->59)
t29 = sin(qJ(4));
t30 = sin(qJ(3));
t33 = cos(qJ(4));
t34 = cos(qJ(3));
t16 = t29 * t34 + t30 * t33;
t31 = sin(qJ(2));
t7 = t16 * t31;
t62 = t7 ^ 2;
t14 = t29 * t30 - t33 * t34;
t61 = t14 ^ 2;
t60 = -0.2e1 * t16;
t59 = 0.2e1 * t34;
t58 = pkin(2) * t34;
t57 = t29 * pkin(2);
t56 = t33 * pkin(2);
t55 = t33 * t7;
t32 = cos(qJ(5));
t54 = t7 * t32;
t28 = sin(qJ(5));
t53 = t28 * t16;
t52 = t28 * t32;
t51 = t30 * t31;
t50 = t32 * t16;
t49 = t34 * t31;
t35 = cos(qJ(2));
t48 = t35 * t34;
t21 = t28 ^ 2;
t25 = t32 ^ 2;
t47 = t21 + t25;
t23 = t30 ^ 2;
t26 = t34 ^ 2;
t46 = t23 + t26;
t45 = t14 * t60;
t44 = t28 * t56;
t43 = t32 * t56;
t42 = t14 * t58;
t41 = pkin(2) * t47;
t40 = -0.2e1 * t42;
t36 = pkin(2) ^ 2;
t39 = t47 * t36;
t9 = -t29 * t51 + t33 * t49;
t3 = -t28 * t9 - t32 * t35;
t4 = -t28 * t35 + t32 * t9;
t38 = -t28 * t4 - t3 * t32;
t1 = -t28 * t3 + t32 * t4;
t37 = (-t14 * t29 - t16 * t33) * pkin(2);
t27 = t35 ^ 2;
t24 = t31 ^ 2;
t22 = t29 ^ 2;
t20 = t33 ^ 2 * t36;
t18 = 0.2e1 * t52;
t13 = t16 ^ 2;
t12 = t29 * t41;
t11 = t32 * t14;
t10 = t28 * t14;
t6 = t28 * t50;
t5 = t7 * t28;
t2 = (-t21 + t25) * t16;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t24 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 ^ 2 + t27 + t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t31, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t35 * t30, t46 * t31, 0, 0, 0, 0, 0, 0, 0, -t35 * t14, -t35 * t16, -t14 * t9 + t16 * t7, pkin(2) * t48, 0, 0, 0, 0, 0, 0, t14 * t3 + t7 * t53, -t14 * t4 + t7 * t50, t38 * t16, t38 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t23, t30 * t59, 0, t26, 0, 0, 0, 0, 0, 0, t13, t45, 0, t61, 0, 0, t40, t58 * t60, 0, t26 * t36, t25 * t13, -0.2e1 * t13 * t52, 0.2e1 * t14 * t50, t21 * t13, t28 * t45, t61, t32 * t40, 0.2e1 * t28 * t42, t16 * t41 * t59, t26 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t49, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, (t29 * t9 - t55) * pkin(2), 0, 0, 0, 0, 0, 0, -t54, t5, t1, (t1 * t29 - t55) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t34, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, 0, 0, t37, 0, t6, t2, t10, -t6, t11, 0, t28 * t37, t32 * t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t57, 0, t22 * t36 + t20, t21, t18, 0, t25, 0, 0, 0.2e1 * t43, -0.2e1 * t44, 0.2e1 * t12, t22 * t39 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t5, t1, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, 0, 0, 0, 0, t6, t2, t10, -t6, t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t57, 0, 0, t21, t18, 0, t25, 0, 0, t43, -t44, t12, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t21, t18, 0, t25, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t53, t14, -t32 * t58, t28 * t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t32, 0, -t28 * t57, -t32 * t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t32, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;

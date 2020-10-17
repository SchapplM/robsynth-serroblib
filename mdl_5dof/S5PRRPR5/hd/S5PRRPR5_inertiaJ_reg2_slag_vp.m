% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:27:58
% DurationCPUTime: 0.70s
% Computational Cost: add. (462->84), mult. (1039->167), div. (0->0), fcn. (1241->10), ass. (0->63)
t33 = sin(pkin(10));
t35 = cos(pkin(10));
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t19 = t33 * t41 + t35 * t38;
t71 = -0.2e1 * t19;
t36 = cos(pkin(5));
t34 = sin(pkin(5));
t39 = sin(qJ(2));
t62 = t34 * t39;
t14 = t36 * t41 - t38 * t62;
t15 = t36 * t38 + t41 * t62;
t5 = -t35 * t14 + t15 * t33;
t70 = t5 ^ 2;
t56 = -qJ(4) - pkin(7);
t21 = t56 * t41;
t51 = t56 * t38;
t9 = -t33 * t21 - t35 * t51;
t69 = t9 ^ 2;
t17 = t33 * t38 - t35 * t41;
t68 = t17 ^ 2;
t27 = -pkin(3) * t41 - pkin(2);
t67 = 0.2e1 * t27;
t66 = 0.2e1 * t41;
t65 = t5 * t9;
t64 = t33 * pkin(3);
t63 = t35 * pkin(3);
t42 = cos(qJ(2));
t61 = t34 * t42;
t37 = sin(qJ(5));
t60 = t37 * t17;
t59 = t37 * t19;
t40 = cos(qJ(5));
t58 = t37 * t40;
t57 = t40 * t19;
t29 = t37 ^ 2;
t31 = t40 ^ 2;
t55 = t29 + t31;
t30 = t38 ^ 2;
t32 = t41 ^ 2;
t54 = t30 + t32;
t53 = t17 * t71;
t52 = t37 * t57;
t11 = -t35 * t21 + t33 * t51;
t8 = pkin(4) * t17 - pkin(8) * t19 + t27;
t1 = -t11 * t37 + t40 * t8;
t2 = t11 * t40 + t37 * t8;
t50 = t1 * t40 + t2 * t37;
t49 = -t1 * t37 + t2 * t40;
t7 = t14 * t33 + t15 * t35;
t3 = -t37 * t7 - t40 * t61;
t4 = -t37 * t61 + t40 * t7;
t48 = t3 * t40 + t37 * t4;
t47 = -t3 * t37 + t4 * t40;
t46 = -t14 * t38 + t15 * t41;
t25 = pkin(8) + t64;
t26 = -pkin(4) - t63;
t45 = -t17 * t25 + t19 * t26;
t28 = t34 ^ 2;
t23 = t28 * t42 ^ 2;
t16 = t19 ^ 2;
t13 = t40 * t17;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t39 ^ 2 + t36 ^ 2 + t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2 + t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t23 + t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t62, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t61, -t38 * t61, t46, pkin(2) * t61 + pkin(7) * t46, 0, 0, 0, 0, 0, 0, -t17 * t61, -t19 * t61, -t17 * t7 + t19 * t5, t11 * t7 - t27 * t61 + t65, 0, 0, 0, 0, 0, 0, t17 * t3 + t5 * t59, -t17 * t4 + t5 * t57, -t48 * t19, t1 * t3 + t2 * t4 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t30, t38 * t66, 0, t32, 0, 0, pkin(2) * t66, -0.2e1 * pkin(2) * t38, 0.2e1 * t54 * pkin(7), pkin(7) ^ 2 * t54 + pkin(2) ^ 2, t16, t53, 0, t68, 0, 0, t17 * t67, t19 * t67, -0.2e1 * t11 * t17 + 0.2e1 * t19 * t9, t11 ^ 2 + t27 ^ 2 + t69, t31 * t16, -0.2e1 * t16 * t58, 0.2e1 * t17 * t57, t29 * t16, t37 * t53, t68, 0.2e1 * t1 * t17 + 0.2e1 * t59 * t9, -0.2e1 * t17 * t2 + 0.2e1 * t57 * t9, t50 * t71, t1 ^ 2 + t2 ^ 2 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, (t33 * t7 - t35 * t5) * pkin(3), 0, 0, 0, 0, 0, 0, -t5 * t40, t5 * t37, t47, t25 * t47 + t5 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t41, 0, -t38 * pkin(7), -t41 * pkin(7), 0, 0, 0, 0, t19, 0, -t17, 0, -t9, -t11, (-t17 * t33 - t19 * t35) * pkin(3), (t11 * t33 - t35 * t9) * pkin(3), t52, (-t29 + t31) * t19, t60, -t52, t13, 0, t37 * t45 - t9 * t40, t9 * t37 + t40 * t45, t49, t25 * t49 + t9 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t63, -0.2e1 * t64, 0, (t33 ^ 2 + t35 ^ 2) * pkin(3) ^ 2, t29, 0.2e1 * t58, 0, t31, 0, 0, -0.2e1 * t26 * t40, 0.2e1 * t26 * t37, 0.2e1 * t55 * t25, t25 ^ 2 * t55 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19, 0, t27, 0, 0, 0, 0, 0, 0, t13, -t60, -t55 * t19, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t59, t17, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t40, 0, -t37 * t25, -t40 * t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;

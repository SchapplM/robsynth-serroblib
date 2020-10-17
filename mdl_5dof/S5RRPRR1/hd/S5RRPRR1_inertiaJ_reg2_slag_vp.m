% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaJ_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:33
% EndTime: 2019-12-05 18:25:36
% DurationCPUTime: 0.68s
% Computational Cost: add. (340->61), mult. (713->124), div. (0->0), fcn. (741->6), ass. (0->58)
t36 = sin(qJ(4));
t41 = pkin(2) + pkin(1);
t59 = t36 * t41;
t23 = pkin(4) + t59;
t35 = sin(qJ(5));
t30 = t35 ^ 2;
t38 = cos(qJ(5));
t32 = t38 ^ 2;
t52 = t30 + t32;
t54 = t52 * t23;
t37 = sin(qJ(2));
t39 = cos(qJ(4));
t40 = cos(qJ(2));
t17 = t36 * t40 + t37 * t39;
t67 = -0.2e1 * t17;
t55 = pkin(3) + qJ(3);
t20 = t55 * t40;
t47 = t55 * t37;
t6 = t36 * t20 + t39 * t47;
t66 = t6 ^ 2;
t65 = 2 * pkin(1);
t15 = t36 * t37 - t39 * t40;
t64 = t15 ^ 2;
t21 = t41 * t40;
t63 = -0.2e1 * t21;
t62 = t6 * t38;
t11 = t35 * t15;
t61 = t35 * t17;
t60 = t35 * t38;
t58 = t37 * t40;
t12 = t38 * t15;
t57 = t38 * t17;
t56 = t39 * t41;
t53 = t52 * pkin(4);
t51 = t37 * qJ(3);
t50 = t15 * t67;
t49 = t35 * t56;
t48 = t38 * t56;
t8 = t39 * t20 - t36 * t47;
t9 = -t17 * pkin(4) - t21;
t2 = -t35 * t8 + t38 * t9;
t3 = t35 * t9 + t38 * t8;
t46 = t2 * t38 + t3 * t35;
t1 = -t2 * t35 + t3 * t38;
t45 = -t15 * t23 - t17 * t56;
t44 = pkin(1) ^ 2;
t42 = qJ(3) ^ 2;
t34 = t41 ^ 2;
t33 = t40 ^ 2;
t31 = t37 ^ 2;
t27 = t39 ^ 2 * t34;
t25 = 0.2e1 * t58;
t24 = 0.2e1 * t60;
t13 = t17 ^ 2;
t10 = t35 * t57;
t5 = t6 * t35;
t4 = (-t30 + t32) * t17;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, t25, 0, t33, 0, 0, 0, 0, 0, 0, t31, t25, 0, t33, 0, 0, t33 * t65, -0.2e1 * pkin(1) * t58, 0.2e1 * (t31 + t33) * qJ(3), t31 * t42 + (t42 + t44) * t33, t13, t50, 0, t64, 0, 0, t15 * t63, t17 * t63, -0.2e1 * t15 * t8 + 0.2e1 * t17 * t6, t21 ^ 2 + t8 ^ 2 + t66, t32 * t13, -0.2e1 * t13 * t60, 0.2e1 * t15 * t57, t30 * t13, t35 * t50, t64, 0.2e1 * t15 * t2 + 0.2e1 * t6 * t61, -0.2e1 * t15 * t3 + 0.2e1 * t57 * t6, t46 * t67, t2 ^ 2 + t3 ^ 2 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t40, 0, 0, 0, 0, 0, 0, 0, t37, 0, t40, 0, -t51, -t40 * qJ(3), -t37 * pkin(1), -pkin(1) * t51, 0, 0, t17, 0, -t15, 0, -t6, -t8, (-t15 * t36 - t17 * t39) * t41, (t36 * t8 - t39 * t6) * t41, t10, t4, t11, -t10, t12, 0, t35 * t45 - t62, t38 * t45 + t5, t1, t1 * t23 - t56 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, 0, 0, t44, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t59, 0, t34 * t36 ^ 2 + t27, t30, t24, 0, t32, 0, 0, 0.2e1 * t48, -0.2e1 * t49, 0.2e1 * t54, t23 ^ 2 * t52 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t37, 0, -t40 * pkin(1), 0, 0, 0, 0, 0, 0, t15, t17, 0, -t21, 0, 0, 0, 0, 0, 0, t12, -t11, -t52 * t17, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, -t6, -t8, 0, 0, t10, t4, t11, -t10, t12, 0, -pkin(4) * t11 - t62, -pkin(4) * t12 + t5, t1, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t59, 0, 0, t30, t24, 0, t32, 0, 0, t48, -t49, t53 + t54, pkin(4) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t30, t24, 0, t32, 0, 0, 0, 0, 0.2e1 * t53, t52 * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t61, t15, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * t23, -t38 * t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * pkin(4), -t38 * pkin(4), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;

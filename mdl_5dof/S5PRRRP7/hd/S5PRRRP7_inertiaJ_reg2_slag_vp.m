% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:30
% DurationCPUTime: 0.70s
% Computational Cost: add. (307->98), mult. (772->175), div. (0->0), fcn. (840->8), ass. (0->71)
t41 = cos(pkin(5));
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t40 = sin(pkin(5));
t44 = sin(qJ(2));
t63 = t40 * t44;
t13 = -t41 * t46 + t43 * t63;
t74 = t13 ^ 2;
t73 = -0.2e1 * t43;
t72 = 0.2e1 * t43;
t71 = 0.2e1 * t46;
t45 = cos(qJ(4));
t70 = pkin(3) * t45;
t42 = sin(qJ(4));
t69 = pkin(7) * t42;
t37 = t43 ^ 2;
t68 = t37 * pkin(7);
t67 = t42 * pkin(4);
t66 = t43 * pkin(7);
t65 = t13 * t43;
t64 = t13 * t45;
t47 = cos(qJ(2));
t62 = t40 * t47;
t61 = t42 * t43;
t60 = t42 * t45;
t59 = t42 * t46;
t32 = t45 * t43;
t58 = t45 * t46;
t57 = -qJ(5) - pkin(8);
t36 = t42 ^ 2;
t38 = t45 ^ 2;
t56 = t36 + t38;
t55 = qJ(5) * t43;
t54 = t43 * t71;
t53 = pkin(7) * t58;
t19 = -pkin(3) * t46 - pkin(8) * t43 - pkin(2);
t16 = t45 * t19;
t52 = -t45 * t55 + t16;
t15 = t41 * t43 + t46 * t63;
t7 = -t15 * t42 - t45 * t62;
t8 = t15 * t45 - t42 * t62;
t3 = -t42 * t7 + t45 * t8;
t10 = -pkin(7) * t59 + t16;
t11 = t19 * t42 + t53;
t51 = -t10 * t42 + t11 * t45;
t50 = t15 * t46 + t65;
t49 = pkin(7) ^ 2;
t39 = t46 ^ 2;
t35 = t40 ^ 2;
t34 = t37 * t49;
t33 = -pkin(4) * t45 - pkin(3);
t31 = t38 * t37;
t30 = t36 * t37;
t28 = t35 * t47 ^ 2;
t27 = 0.2e1 * t60;
t26 = t42 * t32;
t24 = t58 * t73;
t23 = -0.2e1 * t37 * t60;
t22 = t42 * t54;
t21 = t57 * t45;
t20 = t57 * t42;
t18 = (pkin(7) + t67) * t43;
t17 = (-t36 + t38) * t43;
t12 = t13 * t42;
t9 = t53 + (t19 - t55) * t42;
t6 = (-pkin(4) - t69) * t46 + t52;
t5 = t13 * t32 + t46 * t8;
t4 = t13 * t61 - t46 * t7;
t2 = (-t42 * t8 - t45 * t7) * t43;
t1 = t7 ^ 2 + t8 ^ 2 + t74;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t44 ^ 2 + t41 ^ 2 + t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t28 + t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t63, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t62, -t43 * t62, t50, pkin(2) * t62 + pkin(7) * t50, 0, 0, 0, 0, 0, 0, t4, t5, t2, pkin(7) * t65 + t10 * t7 + t11 * t8, 0, 0, 0, 0, 0, 0, t4, t5, t2, t13 * t18 + t6 * t7 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t54, 0, t39, 0, 0, pkin(2) * t71, pkin(2) * t73, 0.2e1 * (t37 + t39) * pkin(7), pkin(2) ^ 2 + t39 * t49 + t34, t31, t23, t24, t30, t22, t39, -0.2e1 * t10 * t46 + 0.2e1 * t42 * t68, 0.2e1 * t11 * t46 + 0.2e1 * t45 * t68, (-t10 * t45 - t11 * t42) * t72, t10 ^ 2 + t11 ^ 2 + t34, t31, t23, t24, t30, t22, t39, 0.2e1 * t18 * t61 - 0.2e1 * t46 * t6, 0.2e1 * t18 * t32 + 0.2e1 * t46 * t9, (-t42 * t9 - t45 * t6) * t72, t18 ^ 2 + t6 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t12, t3, -t13 * pkin(3) + pkin(8) * t3, 0, 0, 0, 0, 0, 0, -t64, t12, t3, t13 * t33 + t20 * t7 - t21 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t46, 0, -t66, -t46 * pkin(7), 0, 0, t26, t17, -t59, -t26, -t58, 0, -pkin(7) * t32 + (-pkin(3) * t43 + pkin(8) * t46) * t42, pkin(8) * t58 + (t69 - t70) * t43, t51, -pkin(3) * t66 + pkin(8) * t51, t26, t17, -t59, -t26, -t58, 0, -t18 * t45 - t20 * t46 + t33 * t61, t18 * t42 - t21 * t46 + t32 * t33, (-t20 * t43 + t9) * t45 + (t21 * t43 - t6) * t42, t18 * t33 + t20 * t6 - t21 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t27, 0, t38, 0, 0, 0.2e1 * t70, -0.2e1 * pkin(3) * t42, 0.2e1 * t56 * pkin(8), pkin(8) ^ 2 * t56 + pkin(3) ^ 2, t36, t27, 0, t38, 0, 0, -0.2e1 * t33 * t45, 0.2e1 * t33 * t42, -0.2e1 * t20 * t42 - 0.2e1 * t21 * t45, t20 ^ 2 + t21 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, t7 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t61, -t46, t10, -t11, 0, 0, 0, 0, t32, 0, -t61, -t46, (-0.2e1 * pkin(4) - t69) * t46 + t52, -t9, -pkin(4) * t32, t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t45, 0, -t42 * pkin(8), -t45 * pkin(8), 0, 0, 0, 0, t42, 0, t45, 0, t20, t21, -t67, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t32, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t42, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t14;

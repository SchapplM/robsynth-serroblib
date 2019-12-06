% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t41 = cos(pkin(5));
t44 = sin(qJ(3));
t48 = cos(qJ(3));
t40 = sin(pkin(5));
t45 = sin(qJ(2));
t66 = t40 * t45;
t18 = -t41 * t48 + t44 * t66;
t17 = t18 ^ 2;
t47 = cos(qJ(4));
t33 = -t47 * pkin(4) - pkin(3);
t80 = 0.2e1 * t33;
t79 = -0.2e1 * t44;
t78 = -0.2e1 * t48;
t77 = 0.2e1 * t48;
t76 = -pkin(9) - pkin(8);
t75 = pkin(3) * t47;
t43 = sin(qJ(4));
t74 = pkin(7) * t43;
t37 = t44 ^ 2;
t73 = t37 * pkin(7);
t42 = sin(qJ(5));
t72 = t42 * pkin(4);
t71 = t44 * pkin(7);
t46 = cos(qJ(5));
t70 = t46 * pkin(4);
t26 = -t48 * pkin(3) - t44 * pkin(8) - pkin(2);
t60 = t47 * t48;
t57 = pkin(7) * t60;
t9 = t57 + (-pkin(9) * t44 + t26) * t43;
t69 = t46 * t9;
t68 = t48 * pkin(4);
t67 = t18 * t44;
t49 = cos(qJ(2));
t65 = t40 * t49;
t64 = t43 * t44;
t63 = t43 * t47;
t62 = t43 * t48;
t61 = t47 * t44;
t36 = t43 ^ 2;
t38 = t47 ^ 2;
t59 = t36 + t38;
t58 = t44 * t77;
t56 = t43 * t61;
t21 = t47 * t26;
t6 = -pkin(9) * t61 + t21 + (-pkin(4) - t74) * t48;
t1 = -t42 * t9 + t46 * t6;
t20 = t41 * t44 + t48 * t66;
t7 = -t20 * t43 - t47 * t65;
t8 = t20 * t47 - t43 * t65;
t55 = -t7 * t43 + t8 * t47;
t12 = -pkin(7) * t62 + t21;
t13 = t43 * t26 + t57;
t54 = -t12 * t43 + t13 * t47;
t53 = t20 * t48 + t67;
t24 = t42 * t47 + t46 * t43;
t51 = pkin(7) ^ 2;
t39 = t48 ^ 2;
t35 = t40 ^ 2;
t34 = t37 * t51;
t30 = t35 * t49 ^ 2;
t28 = t76 * t47;
t27 = t76 * t43;
t25 = (pkin(4) * t43 + pkin(7)) * t44;
t22 = t42 * t43 - t46 * t47;
t16 = -t42 * t64 + t46 * t61;
t14 = t24 * t44;
t11 = t42 * t27 - t46 * t28;
t10 = t46 * t27 + t42 * t28;
t4 = t42 * t7 + t46 * t8;
t3 = -t42 * t8 + t46 * t7;
t2 = t42 * t6 + t69;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t45 ^ 2 + t41 ^ 2 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 ^ 2 + t17 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t65, -t44 * t65, t53, pkin(2) * t65 + t53 * pkin(7), 0, 0, 0, 0, 0, 0, t18 * t64 - t7 * t48, t18 * t61 + t8 * t48, (-t43 * t8 - t47 * t7) * t44, pkin(7) * t67 + t7 * t12 + t8 * t13, 0, 0, 0, 0, 0, 0, t18 * t14 - t3 * t48, t18 * t16 + t4 * t48, -t4 * t14 - t3 * t16, t3 * t1 + t18 * t25 + t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t58, 0, t39, 0, 0, pkin(2) * t77, pkin(2) * t79, 0.2e1 * (t37 + t39) * pkin(7), pkin(2) ^ 2 + t39 * t51 + t34, t38 * t37, -0.2e1 * t37 * t63, t60 * t79, t36 * t37, t43 * t58, t39, -0.2e1 * t12 * t48 + 0.2e1 * t43 * t73, 0.2e1 * t13 * t48 + 0.2e1 * t47 * t73, 0.2e1 * (-t12 * t47 - t13 * t43) * t44, t12 ^ 2 + t13 ^ 2 + t34, t16 ^ 2, -0.2e1 * t16 * t14, t16 * t78, t14 ^ 2, -t14 * t78, t39, -0.2e1 * t1 * t48 + 0.2e1 * t25 * t14, 0.2e1 * t25 * t16 + 0.2e1 * t2 * t48, -0.2e1 * t1 * t16 - 0.2e1 * t2 * t14, t1 ^ 2 + t2 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t47, t18 * t43, t55, -t18 * pkin(3) + pkin(8) * t55, 0, 0, 0, 0, 0, 0, t18 * t22, t18 * t24, -t4 * t22 - t3 * t24, t3 * t10 + t4 * t11 + t18 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t48, 0, -t71, -t48 * pkin(7), 0, 0, t56, (-t36 + t38) * t44, -t62, -t56, -t60, 0, -pkin(7) * t61 + (-pkin(3) * t44 + pkin(8) * t48) * t43, pkin(8) * t60 + (t74 - t75) * t44, t54, -pkin(3) * t71 + pkin(8) * t54, t16 * t24, -t24 * t14 - t16 * t22, -t24 * t48, t14 * t22, t22 * t48, 0, -t10 * t48 + t33 * t14 + t25 * t22, t11 * t48 + t33 * t16 + t25 * t24, -t1 * t24 - t10 * t16 - t11 * t14 - t2 * t22, t1 * t10 + t2 * t11 + t25 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, 0.2e1 * t63, 0, t38, 0, 0, 0.2e1 * t75, -0.2e1 * pkin(3) * t43, 0.2e1 * t59 * pkin(8), pkin(8) ^ 2 * t59 + pkin(3) ^ 2, t24 ^ 2, -0.2e1 * t24 * t22, 0, t22 ^ 2, 0, 0, t22 * t80, t24 * t80, -0.2e1 * t10 * t24 - 0.2e1 * t11 * t22, t10 ^ 2 + t11 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, (t3 * t46 + t4 * t42) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t64, -t48, t12, -t13, 0, 0, 0, 0, t16, 0, -t14, -t48, -t46 * t68 + t1, -t69 + (-t6 + t68) * t42, (-t14 * t42 - t16 * t46) * pkin(4), (t1 * t46 + t2 * t42) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t47, 0, -t43 * pkin(8), -t47 * pkin(8), 0, 0, 0, 0, t24, 0, -t22, 0, t10, -t11, (-t22 * t42 - t24 * t46) * pkin(4), (t10 * t46 + t11 * t42) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t72, 0, (t42 ^ 2 + t46 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, -t48, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, -t22, 0, t10, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t70, -t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;

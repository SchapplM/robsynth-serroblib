% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t50 = sin(qJ(2));
t44 = t50 ^ 2;
t52 = cos(qJ(2));
t45 = t52 ^ 2;
t82 = t44 + t45;
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t21 = -t49 * t46 + t51 * t47;
t17 = t21 ^ 2;
t19 = t51 * t46 + t49 * t47;
t80 = t19 ^ 2;
t81 = t17 + t80;
t33 = t46 * pkin(4) + qJ(3);
t79 = 0.2e1 * t33;
t78 = -0.2e1 * t50;
t77 = 0.2e1 * t50;
t76 = 0.2e1 * t52;
t75 = 0.2e1 * qJ(3);
t48 = -pkin(2) - qJ(4);
t74 = -pkin(7) + t48;
t10 = t21 * t52;
t73 = t19 * t10;
t72 = t19 * t50;
t11 = t19 * t52;
t71 = t21 * t11;
t70 = t46 * t50;
t69 = t46 * t52;
t68 = t47 * t46;
t67 = t47 * t52;
t65 = t50 * t52;
t59 = -t50 * qJ(3) - pkin(1);
t16 = t48 * t52 + t59;
t37 = t50 * pkin(6);
t28 = t50 * pkin(3) + t37;
t7 = t47 * t16 + t46 * t28;
t63 = t82 * pkin(6) ^ 2;
t39 = t52 * pkin(6);
t29 = t52 * pkin(3) + t39;
t41 = t46 ^ 2;
t42 = t47 ^ 2;
t31 = t41 + t42;
t62 = t52 * qJ(3);
t61 = -0.2e1 * t65;
t60 = t46 * t67;
t23 = t47 * t28;
t4 = t50 * pkin(4) + t23 + (pkin(7) * t52 - t16) * t46;
t5 = -pkin(7) * t67 + t7;
t1 = t51 * t4 - t49 * t5;
t2 = t49 * t4 + t51 * t5;
t58 = t1 * t21 + t2 * t19;
t24 = t74 * t46;
t25 = t74 * t47;
t8 = -t49 * t24 + t51 * t25;
t9 = t51 * t24 + t49 * t25;
t57 = t9 * t19 + t8 * t21;
t6 = -t46 * t16 + t23;
t3 = t7 * t46 + t6 * t47;
t56 = -t50 * pkin(2) + t62;
t55 = t48 * t50 + t62;
t53 = qJ(3) ^ 2;
t35 = t47 * t50;
t32 = 0.2e1 * t65;
t27 = -t52 * pkin(2) + t59;
t26 = 0.2e1 * t82 * pkin(6);
t18 = t31 * t48;
t15 = pkin(4) * t67 + t29;
t14 = t21 * t50;
t12 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t44, t32, 0, t45, 0, 0, pkin(1) * t76, pkin(1) * t78, t26, pkin(1) ^ 2 + t63, 0, 0, 0, t44, t32, t45, t26, t27 * t76, t27 * t78, t27 ^ 2 + t63, t41 * t45, 0.2e1 * t45 * t68, t46 * t61, t42 * t45, t47 * t61, t44, 0.2e1 * t29 * t67 + 0.2e1 * t6 * t50, -0.2e1 * t29 * t69 - 0.2e1 * t7 * t50, (t46 * t6 - t47 * t7) * t76, t29 ^ 2 + t6 ^ 2 + t7 ^ 2, t11 ^ 2, 0.2e1 * t11 * t10, -t11 * t77, t10 ^ 2, -t10 * t77, t44, 0.2e1 * t1 * t50 + 0.2e1 * t15 * t10, -0.2e1 * t15 * t11 - 0.2e1 * t2 * t50, 0.2e1 * t1 * t11 - 0.2e1 * t2 * t10, t1 ^ 2 + t15 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t52, 0, -t37, -t39, 0, 0, 0, -t50, -t52, 0, 0, 0, t56, t37, t39, t56 * pkin(6), -t60, (t41 - t42) * t52, t35, t60, -t70, 0, t29 * t46 + t55 * t47, t29 * t47 - t55 * t46, -t3, t29 * qJ(3) + t3 * t48, -t71, -t21 * t10 + t11 * t19, t14, t73, -t72, 0, t33 * t10 + t15 * t19 + t8 * t50, -t33 * t11 + t15 * t21 - t9 * t50, -t9 * t10 + t8 * t11 - t58, t1 * t8 + t15 * t33 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t75, pkin(2) ^ 2 + t53, t42, -0.2e1 * t68, 0, t41, 0, 0, t46 * t75, t47 * t75, -0.2e1 * t18, t31 * t48 ^ 2 + t53, t17, -0.2e1 * t21 * t19, 0, t80, 0, 0, t19 * t79, t21 * t79, -0.2e1 * t57, t33 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, t37, 0, 0, 0, 0, 0, 0, t35, -t70, 0, t3, 0, 0, 0, 0, 0, 0, t14, -t72, t71 - t73, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t31, t18, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t69, 0, t29, 0, 0, 0, 0, 0, 0, t10, -t11, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t47, 0, qJ(3), 0, 0, 0, 0, 0, 0, t19, t21, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t10, t50, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t12;

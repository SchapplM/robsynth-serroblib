% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t46 = sin(pkin(9));
t77 = t46 * pkin(3);
t37 = pkin(8) + t77;
t48 = sin(qJ(5));
t42 = t48 ^ 2;
t51 = cos(qJ(5));
t44 = t51 ^ 2;
t67 = t42 + t44;
t68 = t67 * t37;
t49 = sin(qJ(3));
t50 = sin(qJ(2));
t52 = cos(qJ(3));
t53 = cos(qJ(2));
t30 = -t49 * t50 + t52 * t53;
t31 = t49 * t53 + t52 * t50;
t47 = cos(pkin(9));
t18 = t46 * t30 + t47 * t31;
t86 = -0.2e1 * t18;
t78 = -pkin(7) - pkin(6);
t62 = t78 * t50;
t63 = t78 * t53;
t20 = t49 * t62 - t52 * t63;
t11 = t30 * qJ(4) + t20;
t19 = t49 * t63 + t52 * t62;
t57 = -t31 * qJ(4) + t19;
t6 = t46 * t11 - t47 * t57;
t85 = t6 ^ 2;
t16 = -t47 * t30 + t46 * t31;
t84 = t16 ^ 2;
t40 = -t53 * pkin(2) - pkin(1);
t23 = -t30 * pkin(3) + t40;
t83 = 0.2e1 * t23;
t82 = 0.2e1 * t31;
t81 = 0.2e1 * t48;
t80 = -0.2e1 * t51;
t79 = 0.2e1 * t53;
t76 = t47 * pkin(3);
t75 = t49 * pkin(2);
t74 = t6 * t51;
t13 = t48 * t16;
t73 = t48 * t18;
t72 = t48 * t51;
t71 = t51 * t18;
t41 = t52 * pkin(2);
t39 = t41 + pkin(3);
t64 = t47 * t75;
t28 = t46 * t39 + t64;
t26 = pkin(8) + t28;
t70 = t67 * t26;
t61 = -t47 * t39 + t46 * t75;
t25 = -pkin(4) + t61;
t38 = -pkin(4) - t76;
t69 = t25 + t38;
t43 = t50 ^ 2;
t45 = t53 ^ 2;
t66 = t43 + t45;
t65 = t16 * t86;
t5 = t16 * pkin(4) - t18 * pkin(8) + t23;
t8 = t47 * t11 + t46 * t57;
t2 = -t48 * t8 + t51 * t5;
t3 = t48 * t5 + t51 * t8;
t60 = t2 * t51 + t3 * t48;
t1 = -t2 * t48 + t3 * t51;
t59 = -t16 * t26 + t18 * t25;
t58 = -t16 * t37 + t18 * t38;
t36 = 0.2e1 * t72;
t15 = t18 ^ 2;
t14 = t51 * t16;
t12 = t48 * t71;
t9 = (-t42 + t44) * t18;
t4 = t6 * t48;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, t50 * t79, 0, t45, 0, 0, pkin(1) * t79, -0.2e1 * pkin(1) * t50, 0.2e1 * t66 * pkin(6), t66 * pkin(6) ^ 2 + pkin(1) ^ 2, t31 ^ 2, t30 * t82, 0, t30 ^ 2, 0, 0, -0.2e1 * t40 * t30, t40 * t82, -0.2e1 * t19 * t31 + 0.2e1 * t20 * t30, t19 ^ 2 + t20 ^ 2 + t40 ^ 2, t15, t65, 0, t84, 0, 0, t16 * t83, t18 * t83, -0.2e1 * t8 * t16 + 0.2e1 * t6 * t18, t23 ^ 2 + t8 ^ 2 + t85, t44 * t15, -0.2e1 * t15 * t72, 0.2e1 * t16 * t71, t42 * t15, t48 * t65, t84, 0.2e1 * t2 * t16 + 0.2e1 * t6 * t73, -0.2e1 * t3 * t16 + 0.2e1 * t6 * t71, t60 * t86, t2 ^ 2 + t3 ^ 2 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t53, 0, -t50 * pkin(6), -t53 * pkin(6), 0, 0, 0, 0, t31, 0, t30, 0, t19, -t20, (t30 * t49 - t31 * t52) * pkin(2), (t19 * t52 + t20 * t49) * pkin(2), 0, 0, t18, 0, -t16, 0, -t6, -t8, -t28 * t16 + t18 * t61, t8 * t28 + t6 * t61, t12, t9, t13, -t12, t14, 0, t59 * t48 - t74, t59 * t51 + t4, t1, t1 * t26 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t75, 0, (t49 ^ 2 + t52 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t61, -0.2e1 * t28, 0, t28 ^ 2 + t61 ^ 2, t42, t36, 0, t44, 0, 0, t25 * t80, t25 * t81, 0.2e1 * t70, t67 * t26 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t30, 0, t19, -t20, 0, 0, 0, 0, t18, 0, -t16, 0, -t6, -t8, (-t16 * t46 - t18 * t47) * pkin(3), (t46 * t8 - t47 * t6) * pkin(3), t12, t9, t13, -t12, t14, 0, t58 * t48 - t74, t58 * t51 + t4, t1, t1 * t37 + t6 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t75, 0, 0, 0, 0, 0, 0, 0, 1, -t61 + t76, -t64 + (-pkin(3) - t39) * t46, 0, (t28 * t46 - t47 * t61) * pkin(3), t42, t36, 0, t44, 0, 0, -t69 * t51, t69 * t48, t68 + t70, t25 * t38 + t26 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t76, -0.2e1 * t77, 0, (t46 ^ 2 + t47 ^ 2) * pkin(3) ^ 2, t42, t36, 0, t44, 0, 0, t38 * t80, t38 * t81, 0.2e1 * t68, t67 * t37 ^ 2 + t38 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, t23, 0, 0, 0, 0, 0, 0, t14, -t13, -t67 * t18, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, -t73, t16, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, t51, 0, -t48 * t26, -t51 * t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, t51, 0, -t48 * t37, -t51 * t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;

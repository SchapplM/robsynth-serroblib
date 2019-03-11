% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPPRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t41 = sin(pkin(11));
t42 = sin(pkin(6));
t43 = cos(pkin(11));
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t11 = (-t41 * t50 + t43 * t47) * t42;
t44 = cos(pkin(6));
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t5 = t11 * t46 + t44 * t49;
t81 = t5 ^ 2;
t9 = (-t41 * t47 - t43 * t50) * t42;
t80 = t9 ^ 2;
t79 = -0.2e1 * t46;
t78 = 0.2e1 * t49;
t48 = cos(qJ(6));
t77 = pkin(5) * t48;
t76 = t5 * t49;
t75 = t9 * t43;
t51 = -pkin(2) - pkin(3);
t22 = t43 * qJ(3) + t41 * t51;
t19 = -pkin(8) + t22;
t74 = t19 * t41;
t38 = t46 ^ 2;
t45 = sin(qJ(6));
t73 = t38 * t45;
t72 = t38 * t48;
t39 = t48 ^ 2;
t71 = t39 * t46;
t70 = t42 * t47;
t69 = t45 * t46;
t68 = t45 * t48;
t67 = t45 * t49;
t66 = t46 * t19;
t65 = t46 * t41;
t64 = t48 * t46;
t30 = t48 * t49;
t63 = t49 * t19;
t62 = t49 * t41;
t37 = t45 ^ 2;
t61 = t37 + t39;
t40 = t49 ^ 2;
t60 = t38 + t40;
t59 = t46 * t78;
t58 = t45 * t64;
t57 = t61 * pkin(9);
t20 = t41 * qJ(3) - t43 * t51;
t18 = pkin(4) + t20;
t7 = t11 * t49 - t44 * t46;
t1 = -t7 * t45 - t9 * t48;
t2 = -t9 * t45 + t7 * t48;
t56 = -t1 * t45 + t2 * t48;
t32 = t49 * pkin(5);
t8 = t46 * pkin(9) + t18 + t32;
t3 = -t45 * t63 + t48 * t8;
t4 = t45 * t8 + t48 * t63;
t55 = -t3 * t45 + t4 * t48;
t54 = t5 * t46 + t7 * t49;
t14 = -t48 * t43 - t45 * t62;
t15 = -t45 * t43 + t48 * t62;
t53 = -t14 * t45 + t15 * t48;
t36 = t44 ^ 2;
t35 = t43 ^ 2;
t33 = t41 ^ 2;
t29 = t39 * t38;
t28 = t37 * t46;
t27 = t37 * t38;
t24 = t42 * t50;
t23 = t38 * t33;
t17 = t19 ^ 2;
t16 = t38 * t17;
t13 = t36 + (t47 ^ 2 + t50 ^ 2) * t42 ^ 2;
t12 = t38 * t74;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t36 + t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t80 + t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t70, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t70 (pkin(2) * t50 + qJ(3) * t47) * t42, 0, 0, 0, 0, 0, 0, -t9, t11, 0, t11 * t22 - t9 * t20, 0, 0, 0, 0, 0, 0, -t9 * t49, t9 * t46, -t54, -t9 * t18 + t19 * t54, 0, 0, 0, 0, 0, 0, t1 * t49 - t5 * t69, -t2 * t49 - t5 * t64 (t1 * t48 + t2 * t45) * t46, t1 * t3 + t2 * t4 + t5 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, 0.2e1 * t22, 0, t20 ^ 2 + t22 ^ 2, t38, t59, 0, t40, 0, 0, t18 * t78, t18 * t79, -0.2e1 * t60 * t19, t40 * t17 + t18 ^ 2 + t16, t29, -0.2e1 * t38 * t68, t30 * t79, t27, t45 * t59, t40, -0.2e1 * t19 * t73 + 0.2e1 * t3 * t49, -0.2e1 * t19 * t72 - 0.2e1 * t4 * t49, 0.2e1 * (t3 * t48 + t4 * t45) * t46, t3 ^ 2 + t4 ^ 2 + t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t41 + t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t54 + t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t14 + t2 * t15 + t5 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t43, t41, 0, -t20 * t43 + t22 * t41, 0, 0, 0, 0, 0, 0, -t43 * t49, t46 * t43, -t60 * t41, -t18 * t43 + t40 * t74 + t12, 0, 0, 0, 0, 0, 0, t14 * t49 - t41 * t73, -t15 * t49 - t41 * t72 (t14 * t48 + t15 * t45) * t46, t3 * t14 + t4 * t15 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 + t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t33 + t23 + t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t46 - t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t56 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t55 - t63) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t53 - t62) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 + t27 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t48, t5 * t45, t56, -t5 * pkin(5) + pkin(9) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, -t49, 0, -t66, -t63, 0, 0, -t58, t28 - t71, t67, t58, t30, 0, -t19 * t64 + (pkin(5) * t46 - pkin(9) * t49) * t45, -pkin(9) * t30 + (t19 * t45 + t77) * t46, t55, -pkin(5) * t66 + pkin(9) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t62, 0, 0, 0, 0, 0, 0, 0, 0, -t41 * t64, t45 * t65, t53, -pkin(5) * t65 + pkin(9) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t67, t28 + t71, t46 * t57 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, 0.2e1 * t68, 0, t39, 0, 0, 0.2e1 * t77, -0.2e1 * pkin(5) * t45, 0.2e1 * t57, pkin(9) ^ 2 * t61 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, t69, t49, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t48, 0, -t45 * pkin(9), -t48 * pkin(9), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t6;

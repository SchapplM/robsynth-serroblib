% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t46 = sin(qJ(6));
t47 = sin(qJ(4));
t49 = cos(qJ(6));
t50 = cos(qJ(4));
t20 = t46 * t47 + t49 * t50;
t51 = cos(qJ(3));
t16 = t51 * t20;
t82 = -0.2e1 * t16;
t63 = t47 * qJ(5);
t77 = pkin(4) + pkin(5);
t17 = t77 * t50 + pkin(3) + t63;
t81 = 0.2e1 * t17;
t80 = -0.2e1 * t47;
t79 = 0.2e1 * t50;
t78 = 2 * qJ(2);
t48 = sin(qJ(3));
t76 = pkin(8) * t48;
t75 = t48 * pkin(4);
t14 = t20 * t48;
t74 = t47 * t50;
t73 = t47 * t51;
t52 = -pkin(1) - pkin(7);
t72 = t48 * t52;
t71 = t49 * t47;
t70 = t49 * t48;
t35 = t50 * t48;
t36 = t50 * t51;
t69 = t50 * t52;
t55 = -t50 * pkin(4) - t63;
t28 = -pkin(3) + t55;
t68 = t51 * t28;
t67 = t51 * t48;
t38 = t51 * t52;
t27 = t48 * pkin(3) - t51 * pkin(8) + qJ(2);
t66 = -t50 * t27 + t47 * t72;
t11 = t47 * t27 + t48 * t69;
t42 = t47 ^ 2;
t44 = t50 ^ 2;
t65 = t42 + t44;
t43 = t48 ^ 2;
t45 = t51 ^ 2;
t64 = t43 + t45;
t62 = t50 * qJ(5);
t61 = -0.2e1 * t67;
t39 = t48 * qJ(5);
t6 = t39 + t11;
t40 = t47 * pkin(8);
t60 = -t47 * pkin(9) + t40;
t59 = t65 * t48;
t58 = -pkin(3) * t51 - t76;
t3 = -pkin(9) * t36 - t77 * t48 + t66;
t4 = pkin(9) * t73 + t6;
t1 = t49 * t3 - t46 * t4;
t2 = t46 * t3 + t49 * t4;
t7 = t66 - t75;
t57 = t7 * t47 + t6 * t50;
t56 = -t68 + t76;
t54 = pkin(4) * t47 - t62;
t41 = t50 * pkin(8);
t34 = t47 * t48;
t29 = -t50 * pkin(9) + t41;
t26 = t49 * qJ(5) - t46 * t77;
t25 = t46 * qJ(5) + t49 * t77;
t23 = t64 * t50;
t22 = -t46 * t50 + t71;
t21 = t64 * t47;
t15 = t46 * t36 - t51 * t71;
t13 = t46 * t35 - t47 * t70;
t12 = t54 * t51 - t38;
t9 = t49 * t29 + t46 * t60;
t8 = t46 * t29 - t49 * t60;
t5 = t38 + (-t77 * t47 + t62) * t51;
t10 = [1, 0, 0, -2 * pkin(1), t78, pkin(1) ^ 2 + qJ(2) ^ 2, t45, t61, 0, 0, 0, t48 * t78, t51 * t78, t44 * t45, -0.2e1 * t45 * t74, t67 * t79, t47 * t61, t43, -0.2e1 * t45 * t52 * t47 - 0.2e1 * t48 * t66, -0.2e1 * t11 * t48 - 0.2e1 * t45 * t69, 0.2e1 * t12 * t73 - 0.2e1 * t7 * t48, 0.2e1 * (-t47 * t6 + t50 * t7) * t51, -0.2e1 * t12 * t36 + 0.2e1 * t6 * t48, t12 ^ 2 + t6 ^ 2 + t7 ^ 2, t16 ^ 2, t15 * t82, t48 * t82, 0.2e1 * t15 * t48, t43, -0.2e1 * t1 * t48 + 0.2e1 * t5 * t15, 0.2e1 * t5 * t16 + 0.2e1 * t2 * t48; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, -t21, 0, t23, -t12 * t51 + t57 * t48, 0, 0, 0, 0, 0, t13 * t48 + t51 * t15, t14 * t48 + t51 * t16; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t43 + t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, t38, -t72, t47 * t36 (-t42 + t44) * t51, t34, t35, 0, t50 * t38 + t58 * t47, -t47 * t38 + t58 * t50, -t12 * t50 - t56 * t47, t57, -t12 * t47 + t56 * t50, t57 * pkin(8) + t12 * t28, t16 * t22, -t22 * t15 - t16 * t20, -t22 * t48, t14, 0, t17 * t15 + t5 * t20 + t8 * t48, t17 * t16 + t5 * t22 + t9 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, 0, 0, 0, 0, t36, -t73, t36, t59, t73, pkin(8) * t59 - t68, 0, 0, 0, 0, 0, t16, t51 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t42, 0.2e1 * t74, 0, 0, 0, pkin(3) * t79, pkin(3) * t80, -0.2e1 * t28 * t50, 0.2e1 * t65 * pkin(8), t28 * t80, t65 * pkin(8) ^ 2 + t28 ^ 2, t22 ^ 2, -0.2e1 * t22 * t20, 0, 0, 0, t20 * t81, t22 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t73, t48, -t66, -t11, -t66 + 0.2e1 * t75, t55 * t51, 0.2e1 * t39 + t11, -t7 * pkin(4) + t6 * qJ(5), 0, 0, -t16, t15, t48, t25 * t48 - t1, t26 * t48 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t35, -t34, 0, t35, -t54 * t48, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t50, 0, -t40, -t41, -t40, -t54, t41, -t54 * pkin(8), 0, 0, -t22, t20, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t25, 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t36, 0, t7, 0, 0, 0, 0, 0, -t70, t46 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, -t49, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, -t48, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;

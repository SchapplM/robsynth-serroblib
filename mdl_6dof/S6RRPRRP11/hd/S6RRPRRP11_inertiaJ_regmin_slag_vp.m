% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t48 = sin(qJ(5));
t49 = sin(qJ(4));
t51 = cos(qJ(5));
t52 = cos(qJ(4));
t25 = -t48 * t49 + t51 * t52;
t23 = t25 ^ 2;
t24 = t48 * t52 + t51 * t49;
t84 = t24 ^ 2 + t23;
t34 = t49 * pkin(4) + qJ(3);
t82 = 0.2e1 * t34;
t50 = sin(qJ(2));
t81 = -0.2e1 * t50;
t80 = 0.2e1 * t50;
t53 = cos(qJ(2));
t79 = 0.2e1 * t53;
t78 = 0.2e1 * qJ(3);
t54 = -pkin(2) - pkin(8);
t77 = t25 * pkin(5);
t76 = t48 * pkin(4);
t75 = t50 * pkin(4);
t41 = t51 * pkin(4);
t61 = -t50 * qJ(3) - pkin(1);
t21 = t54 * t53 + t61;
t62 = pkin(9) * t53 - t21;
t40 = t50 * pkin(7);
t30 = t50 * pkin(3) + t40;
t71 = t49 * t30;
t9 = -t62 * t52 + t71;
t74 = t51 * t9;
t73 = t24 * t50;
t17 = t24 * t53;
t72 = t25 * t17;
t70 = t49 * t50;
t69 = t49 * t53;
t68 = t50 * t53;
t67 = t52 * t49;
t66 = t52 * t53;
t42 = t53 * pkin(7);
t31 = t53 * pkin(3) + t42;
t45 = t50 ^ 2;
t47 = t53 ^ 2;
t65 = t45 + t47;
t64 = t53 * qJ(3);
t63 = -0.2e1 * t68;
t19 = pkin(4) * t66 + t31;
t26 = t52 * t30;
t8 = t62 * t49 + t26 + t75;
t3 = -t48 * t9 + t51 * t8;
t27 = (-pkin(9) + t54) * t49;
t38 = t52 * t54;
t28 = -t52 * pkin(9) + t38;
t13 = -t48 * t27 + t51 * t28;
t1 = t50 * pkin(5) + t17 * qJ(6) + t3;
t16 = t48 * t69 - t51 * t66;
t4 = t48 * t8 + t74;
t2 = t16 * qJ(6) + t4;
t60 = t1 * t25 + t2 * t24;
t6 = -t25 * qJ(6) + t13;
t14 = t51 * t27 + t48 * t28;
t7 = -t24 * qJ(6) + t14;
t59 = t7 * t24 + t6 * t25;
t58 = -t50 * pkin(2) + t64;
t57 = t50 * t54 + t64;
t37 = t41 + pkin(5);
t56 = t24 * t76 + t25 * t37;
t46 = t52 ^ 2;
t44 = t49 ^ 2;
t36 = t52 * t50;
t29 = -t53 * pkin(2) + t61;
t20 = t50 * t25;
t15 = t24 * pkin(5) + t34;
t12 = t52 * t21 + t71;
t11 = -t49 * t21 + t26;
t10 = -t16 * pkin(5) + t19;
t5 = [1, 0, 0, t45, 0.2e1 * t68, 0, 0, 0, pkin(1) * t79, pkin(1) * t81, 0.2e1 * t65 * pkin(7), t29 * t79, t29 * t81, t65 * pkin(7) ^ 2 + t29 ^ 2, t44 * t47, 0.2e1 * t47 * t67, t49 * t63, t52 * t63, t45, 0.2e1 * t11 * t50 + 0.2e1 * t31 * t66, -0.2e1 * t12 * t50 - 0.2e1 * t31 * t69, t17 ^ 2, -0.2e1 * t17 * t16, -t17 * t80, t16 * t80, t45, -0.2e1 * t19 * t16 + 0.2e1 * t3 * t50, -0.2e1 * t19 * t17 - 0.2e1 * t4 * t50, 0.2e1 * t1 * t17 + 0.2e1 * t2 * t16, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, t50, t53, 0, -t40, -t42, t58, t40, t42, t58 * pkin(7), -t49 * t66 (t44 - t46) * t53, t36, -t70, 0, t31 * t49 + t57 * t52, t31 * t52 - t57 * t49, -t72, t25 * t16 + t17 * t24, t20, -t73, 0, t13 * t50 - t34 * t16 + t19 * t24, -t14 * t50 - t34 * t17 + t19 * t25, t7 * t16 + t6 * t17 - t60, t1 * t6 + t10 * t15 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t78, pkin(2) ^ 2 + qJ(3) ^ 2, t46, -0.2e1 * t67, 0, 0, 0, t49 * t78, t52 * t78, t23, -0.2e1 * t25 * t24, 0, 0, 0, t24 * t82, t25 * t82, -0.2e1 * t59, t15 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, t40, 0, 0, 0, 0, 0, t36, -t70, 0, 0, 0, 0, 0, t20, -t73, t24 * t16 + t72, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t66, t50, t11, -t12, 0, 0, -t17, t16, t50, t50 * t41 + t3, -t74 + (-t8 - t75) * t48, t16 * t76 + t37 * t17, t1 * t37 + t2 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t49, 0, t38, -t49 * t54, 0, 0, t25, -t24, 0, t13, -t14, -t56, t6 * t37 + t7 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t49, 0, 0, 0, 0, 0, t25, -t24, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t76, 0, t48 ^ 2 * pkin(4) ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, t50, t3, -t4, pkin(5) * t17, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t13, -t14, -t77, t6 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t76, 0, t37 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;

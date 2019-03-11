% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = sin(pkin(12));
t47 = cos(pkin(12));
t51 = sin(qJ(4));
t70 = cos(qJ(4));
t29 = t51 * t45 - t70 * t47;
t79 = -0.2e1 * t29;
t78 = 0.2e1 * t29;
t38 = -t47 * pkin(3) - pkin(2);
t77 = 0.2e1 * t38;
t54 = cos(qJ(5));
t40 = -t54 * pkin(5) - pkin(4);
t76 = 0.2e1 * t40;
t75 = pkin(9) + pkin(10);
t74 = t29 * pkin(5);
t49 = sin(qJ(6));
t73 = t49 * pkin(5);
t53 = cos(qJ(6));
t72 = t53 * pkin(5);
t30 = t70 * t45 + t51 * t47;
t18 = t29 * pkin(4) - t30 * pkin(9) + t38;
t50 = sin(qJ(5));
t61 = pkin(8) + qJ(3);
t33 = t61 * t45;
t34 = t61 * t47;
t20 = -t51 * t33 + t70 * t34;
t63 = t54 * t20;
t7 = t63 + (-pkin(10) * t30 + t18) * t50;
t71 = t53 * t7;
t32 = t49 * t54 + t53 * t50;
t69 = t32 * t29;
t46 = sin(pkin(6));
t68 = t46 * sin(qJ(2));
t55 = cos(qJ(2));
t67 = t46 * t55;
t66 = t50 * t29;
t65 = t50 * t30;
t64 = t50 * t54;
t62 = t54 * t30;
t60 = t45 ^ 2 + t47 ^ 2;
t59 = t30 * t79;
t8 = t54 * t18 - t50 * t20;
t6 = -pkin(10) * t62 + t74 + t8;
t1 = -t49 * t7 + t53 * t6;
t58 = -pkin(4) * t30 - pkin(9) * t29;
t48 = cos(pkin(6));
t24 = -t45 * t68 + t48 * t47;
t25 = t48 * t45 + t47 * t68;
t57 = -t24 * t45 + t25 * t47;
t31 = t49 * t50 - t53 * t54;
t19 = t70 * t33 + t51 * t34;
t44 = t54 ^ 2;
t43 = t50 ^ 2;
t36 = t75 * t54;
t35 = t75 * t50;
t28 = t30 ^ 2;
t27 = t29 ^ 2;
t26 = t54 * t29;
t23 = -t49 * t35 + t53 * t36;
t22 = -t53 * t35 - t49 * t36;
t21 = t31 * t29;
t16 = t31 * t30;
t15 = t32 * t30;
t14 = t51 * t24 + t70 * t25;
t13 = -t70 * t24 + t51 * t25;
t12 = pkin(5) * t65 + t19;
t11 = t54 * t14 - t50 * t67;
t10 = -t50 * t14 - t54 * t67;
t9 = t50 * t18 + t63;
t4 = t49 * t10 + t53 * t11;
t3 = t53 * t10 - t49 * t11;
t2 = t49 * t6 + t71;
t5 = [1, 0, 0, 0, 0, 0, 0, t46 ^ 2 * t55 ^ 2 + t24 ^ 2 + t25 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t67, -t68, t47 * t67, -t45 * t67, t57, pkin(2) * t67 + t57 * qJ(3), 0, 0, 0, 0, 0, -t29 * t67, -t30 * t67, 0, 0, 0, 0, 0, t10 * t29 + t13 * t65, -t11 * t29 + t13 * t62, 0, 0, 0, 0, 0, t13 * t15 + t3 * t29, -t13 * t16 - t4 * t29; 0, 1, 0, 0, 0.2e1 * pkin(2) * t47, -0.2e1 * pkin(2) * t45, 0.2e1 * t60 * qJ(3), t60 * qJ(3) ^ 2 + pkin(2) ^ 2, t28, t59, 0, 0, 0, t29 * t77, t30 * t77, t44 * t28, -0.2e1 * t28 * t64, t62 * t78, t50 * t59, t27, 0.2e1 * t19 * t65 + 0.2e1 * t8 * t29, 0.2e1 * t19 * t62 - 0.2e1 * t9 * t29, t16 ^ 2, 0.2e1 * t16 * t15, -t16 * t78, t15 * t79, t27, 0.2e1 * t1 * t29 + 0.2e1 * t12 * t15, -0.2e1 * t12 * t16 - 0.2e1 * t2 * t29; 0, 0, 0, 0, 0, 0, 0, -t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t47, t45, 0, -pkin(2), 0, 0, 0, 0, 0, t29, t30, 0, 0, 0, 0, 0, t26, -t66, 0, 0, 0, 0, 0, -t21, -t69; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, 0, 0, 0, 0, -t13 * t54, t13 * t50, 0, 0, 0, 0, 0, t13 * t31, t13 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, -t19, -t20, t50 * t62 (-t43 + t44) * t30, t66, t26, 0, -t19 * t54 + t58 * t50, t19 * t50 + t58 * t54, -t16 * t32, -t32 * t15 + t16 * t31, t69, -t21, 0, t12 * t31 + t40 * t15 + t22 * t29, t12 * t32 - t40 * t16 - t23 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t43, 0.2e1 * t64, 0, 0, 0, 0.2e1 * pkin(4) * t54, -0.2e1 * pkin(4) * t50, t32 ^ 2, -0.2e1 * t32 * t31, 0, 0, 0, t31 * t76, t32 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t65, t29, t8, -t9, 0, 0, -t16, -t15, t29, t29 * t72 + t1, -t71 + (-t6 - t74) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t50, 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t54, 0, -t50 * pkin(9), -t54 * pkin(9), 0, 0, t32, -t31, 0, t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t72, -0.2e1 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t15, t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t72, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;

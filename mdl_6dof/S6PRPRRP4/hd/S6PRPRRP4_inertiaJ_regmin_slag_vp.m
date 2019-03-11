% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t34 = sin(pkin(11));
t36 = cos(pkin(11));
t39 = sin(qJ(4));
t66 = cos(qJ(4));
t24 = t66 * t34 + t39 * t36;
t75 = -0.2e1 * t24;
t29 = -t36 * pkin(3) - pkin(2);
t74 = 0.2e1 * t29;
t38 = sin(qJ(5));
t73 = -0.2e1 * t38;
t23 = t39 * t34 - t66 * t36;
t72 = pkin(9) * t23;
t71 = t23 * pkin(5);
t70 = t38 * pkin(9);
t41 = cos(qJ(5));
t69 = t41 * pkin(9);
t37 = cos(pkin(6));
t35 = sin(pkin(6));
t65 = t35 * sin(qJ(2));
t17 = -t34 * t65 + t37 * t36;
t18 = t37 * t34 + t36 * t65;
t9 = -t66 * t17 + t39 * t18;
t68 = t9 * t38;
t67 = t9 * t41;
t42 = cos(qJ(2));
t64 = t35 * t42;
t19 = t38 * t23;
t63 = t38 * t24;
t62 = t38 * t41;
t20 = t41 * t23;
t21 = t41 * t24;
t61 = pkin(8) + qJ(3);
t13 = t23 * pkin(4) - t24 * pkin(9) + t29;
t26 = t61 * t34;
t27 = t61 * t36;
t16 = -t39 * t26 + t66 * t27;
t4 = t38 * t13 + t41 * t16;
t60 = t34 ^ 2 + t36 ^ 2;
t32 = t38 ^ 2;
t33 = t41 ^ 2;
t59 = t32 + t33;
t58 = t23 * qJ(6);
t57 = t23 * t75;
t10 = t39 * t17 + t66 * t18;
t7 = t38 * t10 + t41 * t64;
t56 = -t7 * t23 + t9 * t63;
t55 = -t41 * t13 + t38 * t16;
t54 = -pkin(4) * t24 - t72;
t1 = t58 + t4;
t2 = t55 - t71;
t53 = t1 * t41 + t2 * t38;
t52 = t1 * t38 - t2 * t41;
t8 = t41 * t10 - t38 * t64;
t51 = t8 * t38 - t7 * t41;
t50 = t7 * t38 + t8 * t41;
t48 = t41 * pkin(5) + t38 * qJ(6);
t25 = -pkin(4) - t48;
t49 = -t24 * t25 + t72;
t47 = pkin(5) * t38 - t41 * qJ(6);
t46 = -t17 * t34 + t18 * t36;
t45 = t9 * t21 - t8 * t23;
t15 = t66 * t26 + t39 * t27;
t22 = t24 ^ 2;
t5 = t47 * t24 + t15;
t3 = [1, 0, 0, 0, 0, 0, 0, t35 ^ 2 * t42 ^ 2 + t17 ^ 2 + t18 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, t64, -t65, t36 * t64, -t34 * t64, t46, pkin(2) * t64 + t46 * qJ(3), 0, 0, 0, 0, 0, -t23 * t64, -t24 * t64, 0, 0, 0, 0, 0, t56, t45, t56, -t51 * t24, -t45, t8 * t1 + t7 * t2 + t9 * t5; 0, 1, 0, 0, 0.2e1 * pkin(2) * t36, -0.2e1 * pkin(2) * t34, 0.2e1 * t60 * qJ(3), t60 * qJ(3) ^ 2 + pkin(2) ^ 2, t22, t57, 0, 0, 0, t23 * t74, t24 * t74, t33 * t22, -0.2e1 * t22 * t62, 0.2e1 * t23 * t21, t38 * t57, t23 ^ 2, 0.2e1 * t15 * t63 - 0.2e1 * t23 * t55, 0.2e1 * t15 * t21 - 0.2e1 * t4 * t23, -0.2e1 * t2 * t23 + 0.2e1 * t5 * t63, t52 * t75, 0.2e1 * t1 * t23 - 0.2e1 * t5 * t21, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, -t36, t34, 0, -pkin(2), 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, t20, -t19, t20, -t59 * t24, t19, t52; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, 0, 0, 0, 0, -t67, t68, -t67, t50, -t68, t50 * pkin(9) + t9 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t15, -t16, t38 * t21 (-t32 + t33) * t24, t19, t20, 0, -t15 * t41 + t54 * t38, t15 * t38 + t54 * t41, -t49 * t38 - t5 * t41, t53, -t5 * t38 + t49 * t41, t53 * pkin(9) + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, 0.2e1 * t62, 0, 0, 0, 0.2e1 * pkin(4) * t41, pkin(4) * t73, -0.2e1 * t25 * t41, 0.2e1 * t59 * pkin(9), t25 * t73, t59 * pkin(9) ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -t7, 0, t8, -t7 * pkin(5) + t8 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t63, t23, -t55, -t4, -t55 + 0.2e1 * t71, -t48 * t24, 0.2e1 * t58 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38, t41, 0, t38, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t70, -t69, -t70, -t47, t69, -t47 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;

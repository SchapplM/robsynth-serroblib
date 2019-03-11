% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t54 = cos(qJ(3));
t42 = -t54 * pkin(3) - pkin(2);
t45 = sin(pkin(12));
t47 = cos(pkin(12));
t51 = sin(qJ(3));
t58 = t45 * t51 - t47 * t54;
t26 = pkin(4) * t58 + t42;
t74 = 0.2e1 * t26;
t73 = 0.2e1 * t54;
t72 = t45 * pkin(3);
t53 = cos(qJ(6));
t62 = -qJ(4) - pkin(8);
t38 = t62 * t51;
t39 = t62 * t54;
t25 = t38 * t45 - t39 * t47;
t15 = -pkin(9) * t58 + t25;
t50 = sin(qJ(5));
t24 = t38 * t47 + t39 * t45;
t34 = t45 * t54 + t47 * t51;
t57 = -pkin(9) * t34 + t24;
t69 = cos(qJ(5));
t7 = t15 * t50 - t57 * t69;
t71 = t7 * t53;
t41 = pkin(3) * t47 + pkin(4);
t29 = t41 * t69 - t50 * t72;
t27 = -pkin(5) - t29;
t70 = pkin(5) - t27;
t48 = cos(pkin(6));
t46 = sin(pkin(6));
t67 = t46 * sin(qJ(2));
t32 = t48 * t54 - t51 * t67;
t33 = t48 * t51 + t54 * t67;
t16 = t32 * t47 - t33 * t45;
t17 = t32 * t45 + t33 * t47;
t10 = -t16 * t69 + t17 * t50;
t68 = t10 * t53;
t55 = cos(qJ(2));
t66 = t46 * t55;
t22 = t34 * t50 + t58 * t69;
t49 = sin(qJ(6));
t19 = t49 * t22;
t23 = t34 * t69 - t50 * t58;
t65 = t49 * t23;
t64 = t49 * t53;
t63 = t53 * t23;
t61 = -0.2e1 * t23 * t22;
t60 = -pkin(5) * t23 - pkin(10) * t22;
t30 = -t41 * t50 - t69 * t72;
t28 = pkin(10) - t30;
t59 = -t22 * t28 + t23 * t27;
t44 = t53 ^ 2;
t43 = t49 ^ 2;
t40 = 0.2e1 * t64;
t21 = t23 ^ 2;
t20 = t53 * t22;
t18 = t49 * t63;
t12 = (-t43 + t44) * t23;
t11 = t16 * t50 + t17 * t69;
t9 = t10 * t49;
t8 = t15 * t69 + t50 * t57;
t6 = t22 * pkin(5) - t23 * pkin(10) + t26;
t5 = t7 * t49;
t4 = t11 * t53 - t49 * t66;
t3 = -t11 * t49 - t53 * t66;
t2 = t49 * t6 + t53 * t8;
t1 = -t49 * t8 + t53 * t6;
t13 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 ^ 2 * t55 ^ 2 + t16 ^ 2 + t17 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t66, -t67, 0, 0, 0, 0, 0, t54 * t66, -t51 * t66, -t16 * t34 - t17 * t58, t16 * t24 + t17 * t25 - t42 * t66, 0, 0, 0, 0, 0, -t22 * t66, -t23 * t66, 0, 0, 0, 0, 0, t10 * t65 + t22 * t3, t10 * t63 - t22 * t4; 0, 1, 0, 0, t51 ^ 2, t51 * t73, 0, 0, 0, pkin(2) * t73, -0.2e1 * pkin(2) * t51, -0.2e1 * t24 * t34 - 0.2e1 * t25 * t58, t24 ^ 2 + t25 ^ 2 + t42 ^ 2, t21, t61, 0, 0, 0, t22 * t74, t23 * t74, t44 * t21, -0.2e1 * t21 * t64, 0.2e1 * t22 * t63, t49 * t61, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t65 * t7, -0.2e1 * t2 * t22 + 0.2e1 * t63 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33, 0 (t16 * t47 + t17 * t45) * pkin(3), 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t68, t9; 0, 0, 0, 0, 0, 0, t51, t54, 0, -t51 * pkin(8), -t54 * pkin(8) (-t47 * t34 - t45 * t58) * pkin(3) (t24 * t47 + t25 * t45) * pkin(3), 0, 0, t23, -t22, 0, -t7, -t8, t18, t12, t19, t20, 0, t49 * t59 - t71, t53 * t59 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t45 ^ 2 + t47 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t29, 0.2e1 * t30, t43, t40, 0, 0, 0, -0.2e1 * t27 * t53, 0.2e1 * t27 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, t20, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t68, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t7, -t8, t18, t12, t19, t20, 0, t49 * t60 - t71, t53 * t60 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, t30, t43, t40, 0, 0, 0, t70 * t53, -t70 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t43, t40, 0, 0, 0, 0.2e1 * pkin(5) * t53, -0.2e1 * pkin(5) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t65, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t53, 0, -t49 * t28, -t53 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t53, 0, -t49 * pkin(10), -t53 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;

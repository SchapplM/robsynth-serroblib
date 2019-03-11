% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t42 = sin(pkin(10));
t43 = cos(pkin(10));
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t26 = t42 * t46 - t43 * t49;
t27 = t42 * t49 + t43 * t46;
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t13 = -t26 * t48 + t27 * t45;
t74 = 0.2e1 * t13;
t44 = sin(qJ(6));
t73 = -0.2e1 * t44;
t47 = cos(qJ(6));
t72 = 0.2e1 * t47;
t71 = 0.2e1 * t49;
t60 = -qJ(3) - pkin(7);
t31 = t60 * t49;
t55 = t60 * t46;
t16 = -t31 * t42 - t43 * t55;
t8 = -pkin(8) * t27 + t16;
t18 = -t31 * t43 + t42 * t55;
t9 = pkin(8) * t26 + t18;
t4 = t45 * t9 - t48 * t8;
t70 = t4 * t44;
t69 = t4 * t47;
t38 = -pkin(2) * t49 - pkin(1);
t36 = pkin(2) * t43 + pkin(3);
t32 = -pkin(4) - t36;
t34 = pkin(2) * t42 + qJ(4);
t21 = -t32 * t48 + t34 * t45;
t19 = pkin(5) + t21;
t68 = pkin(5) + t19;
t67 = t44 * t13;
t14 = t26 * t45 + t27 * t48;
t66 = t44 * t14;
t65 = t44 * t47;
t64 = t47 * t13;
t63 = t47 * t14;
t62 = t48 * t44;
t61 = t48 * t47;
t59 = -0.2e1 * t14 * t13;
t58 = -0.2e1 * t65;
t57 = t16 ^ 2 + t18 ^ 2;
t56 = t44 * t63;
t11 = pkin(3) * t26 - qJ(4) * t27 + t38;
t54 = -pkin(5) * t14 - pkin(9) * t13;
t22 = t32 * t45 + t34 * t48;
t20 = -pkin(9) + t22;
t53 = -t13 * t20 + t14 * t19;
t52 = -t13 * t45 - t14 * t48;
t51 = 0.2e1 * t16 * t27 - 0.2e1 * t18 * t26;
t7 = -pkin(4) * t26 - t11;
t41 = t47 ^ 2;
t40 = t44 ^ 2;
t33 = 0.2e1 * t65;
t12 = t14 ^ 2;
t6 = (t40 - t41) * t14;
t5 = t45 * t8 + t48 * t9;
t3 = pkin(5) * t13 - pkin(9) * t14 + t7;
t2 = t3 * t44 + t47 * t5;
t1 = t3 * t47 - t44 * t5;
t10 = [1, 0, 0, t46 ^ 2, t46 * t71, 0, 0, 0, pkin(1) * t71, -0.2e1 * pkin(1) * t46, t51, t38 ^ 2 + t57, 0.2e1 * t11 * t26, t51, -0.2e1 * t11 * t27, t11 ^ 2 + t57, t12, t59, 0, 0, 0, t7 * t74, 0.2e1 * t7 * t14, t41 * t12, t12 * t58, t63 * t74, t44 * t59, t13 ^ 2, 0.2e1 * t1 * t13 + 0.2e1 * t4 * t66, -0.2e1 * t13 * t2 + 0.2e1 * t4 * t63; 0, 0, 0, 0, 0, t46, t49, 0, -t46 * pkin(7), -t49 * pkin(7) (-t26 * t42 - t27 * t43) * pkin(2) (-t16 * t43 + t18 * t42) * pkin(2), -t16, -t26 * t34 - t27 * t36, t18, -t16 * t36 + t18 * t34, 0, 0, -t14, t13, 0, t4, t5, -t56, t6, -t67, -t64, 0, t44 * t53 + t69, t47 * t53 - t70; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t42 ^ 2 + t43 ^ 2) * pkin(2) ^ 2, 0.2e1 * t36, 0, 0.2e1 * t34, t34 ^ 2 + t36 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t21, 0.2e1 * t22, t40, t33, 0, 0, 0, t19 * t72, t19 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t26, 0, -t27, t11, 0, 0, 0, 0, 0, -t13, -t14, 0, 0, 0, 0, 0, -t64, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t44, t52 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t36, 0, 0, 0, 0, 0, -t48, t45, 0, 0, 0, 0, 0, -t61, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, -t4, -t5, t56, -t6, t67, t64, 0, t44 * t54 - t69, t47 * t54 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t21, -t22, -t40, t58, 0, 0, 0, -t68 * t47, t68 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t45, 0, 0, 0, 0, 0, t61, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t40, t33, 0, 0, 0, pkin(5) * t72, pkin(5) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t66, t13, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t47, 0, -t44 * t20, -t47 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t45, -t47 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t47, 0, -t44 * pkin(9), -t47 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;

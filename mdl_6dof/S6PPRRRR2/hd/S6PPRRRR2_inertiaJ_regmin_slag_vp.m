% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t47 = sin(qJ(6));
t48 = sin(qJ(5));
t51 = cos(qJ(6));
t52 = cos(qJ(5));
t28 = t47 * t48 - t51 * t52;
t49 = sin(qJ(4));
t23 = t28 * t49;
t75 = 0.2e1 * t23;
t35 = -pkin(5) * t52 - pkin(4);
t74 = 0.2e1 * t35;
t53 = cos(qJ(4));
t73 = 0.2e1 * t53;
t72 = pkin(10) + pkin(11);
t71 = pkin(4) * t52;
t70 = pkin(9) * t48;
t69 = t47 * pkin(5);
t68 = t51 * pkin(5);
t67 = t53 * pkin(5);
t42 = sin(pkin(7));
t50 = sin(qJ(3));
t66 = t42 * t50;
t54 = cos(qJ(3));
t65 = t42 * t54;
t44 = cos(pkin(13));
t45 = cos(pkin(7));
t64 = t44 * t45;
t63 = t48 * t49;
t62 = t48 * t52;
t61 = t48 * t53;
t31 = -pkin(4) * t53 - pkin(10) * t49 - pkin(3);
t58 = t52 * t53;
t55 = pkin(9) * t58;
t17 = t55 + (-pkin(11) * t49 + t31) * t48;
t60 = t51 * t17;
t59 = t52 * t49;
t57 = t53 * t49;
t56 = 0.2e1 * t57;
t27 = t52 * t31;
t12 = -pkin(11) * t59 + t27 + (-pkin(5) - t70) * t53;
t5 = t12 * t51 - t17 * t47;
t29 = t47 * t52 + t48 * t51;
t46 = cos(pkin(6));
t43 = sin(pkin(6));
t41 = sin(pkin(13));
t40 = t53 ^ 2;
t39 = t52 ^ 2;
t38 = t49 ^ 2;
t37 = t48 ^ 2;
t33 = t72 * t52;
t32 = t72 * t48;
t30 = (pkin(5) * t48 + pkin(9)) * t49;
t26 = t45 * t49 + t53 * t66;
t25 = -t45 * t53 + t49 * t66;
t24 = -t42 * t43 * t44 + t45 * t46;
t22 = t29 * t49;
t21 = t31 * t48 + t55;
t20 = -pkin(9) * t61 + t27;
t19 = -t32 * t47 + t33 * t51;
t18 = -t32 * t51 - t33 * t47;
t16 = t26 * t52 - t48 * t65;
t15 = -t26 * t48 - t52 * t65;
t14 = t46 * t66 + (t41 * t54 + t50 * t64) * t43;
t13 = -t46 * t65 + (t41 * t50 - t54 * t64) * t43;
t10 = t14 * t53 + t24 * t49;
t9 = t14 * t49 - t24 * t53;
t8 = t15 * t47 + t16 * t51;
t7 = t15 * t51 - t16 * t47;
t6 = t12 * t47 + t60;
t4 = t10 * t52 + t13 * t48;
t3 = -t10 * t48 + t13 * t52;
t2 = t3 * t47 + t4 * t51;
t1 = t3 * t51 - t4 * t47;
t11 = [1, t46 ^ 2 + (t41 ^ 2 + t44 ^ 2) * t43 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t13, -t14, 0, 0, 0, 0, 0, -t13 * t53, t13 * t49, 0, 0, 0, 0, 0, -t3 * t53 + t63 * t9, t4 * t53 + t59 * t9, 0, 0, 0, 0, 0, -t1 * t53 + t22 * t9, t2 * t53 - t23 * t9; 0, 0, 0, t65, -t66, 0, 0, 0, 0, 0, t53 * t65, -t49 * t65, 0, 0, 0, 0, 0, -t15 * t53 + t25 * t63, t16 * t53 + t25 * t59, 0, 0, 0, 0, 0, t22 * t25 - t53 * t7, -t23 * t25 + t53 * t8; 0, 0, 1, 0, 0, t38, t56, 0, 0, 0, pkin(3) * t73, -0.2e1 * pkin(3) * t49, t39 * t38, -0.2e1 * t38 * t62, -0.2e1 * t52 * t57, t48 * t56, t40, -0.2e1 * t20 * t53 + 0.2e1 * t38 * t70, 0.2e1 * pkin(9) * t38 * t52 + 0.2e1 * t21 * t53, t23 ^ 2, t22 * t75, t53 * t75, t22 * t73, t40, 0.2e1 * t22 * t30 - 0.2e1 * t5 * t53, -0.2e1 * t23 * t30 + 0.2e1 * t53 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, 0, 0, 0, 0, -t9 * t52, t9 * t48, 0, 0, 0, 0, 0, t9 * t28, t9 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, 0, 0, 0, 0, -t25 * t52, t25 * t48, 0, 0, 0, 0, 0, t25 * t28, t25 * t29; 0, 0, 0, 0, 0, 0, 0, t49, t53, 0, -t49 * pkin(9), -t53 * pkin(9), t48 * t59 (-t37 + t39) * t49, -t61, -t58, 0, -pkin(9) * t59 + (-pkin(4) * t49 + pkin(10) * t53) * t48, pkin(10) * t58 + (t70 - t71) * t49, -t23 * t29, -t22 * t29 + t23 * t28, -t29 * t53, t28 * t53, 0, -t18 * t53 + t22 * t35 + t28 * t30, t19 * t53 - t23 * t35 + t29 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t37, 0.2e1 * t62, 0, 0, 0, 0.2e1 * t71, -0.2e1 * pkin(4) * t48, t29 ^ 2, -0.2e1 * t29 * t28, 0, 0, 0, t28 * t74, t29 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, 0, 0, 0, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t63, -t53, t20, -t21, 0, 0, -t23, -t22, -t53, -t51 * t67 + t5, -t60 + (-t12 + t67) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t52, 0, -t48 * pkin(10), -t52 * pkin(10), 0, 0, t29, -t28, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22, -t53, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t68, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;

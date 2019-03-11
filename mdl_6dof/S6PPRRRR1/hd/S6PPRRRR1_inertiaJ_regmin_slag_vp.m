% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRRR1
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t54 = cos(qJ(4));
t39 = -t54 * pkin(4) - pkin(3);
t76 = 0.2e1 * t39;
t75 = 0.2e1 * t54;
t74 = pkin(9) + pkin(10);
t43 = sin(pkin(13));
t45 = sin(pkin(6));
t48 = cos(pkin(6));
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t46 = cos(pkin(13));
t47 = cos(pkin(7));
t65 = t46 * t47;
t44 = sin(pkin(7));
t67 = t44 * t52;
t19 = t48 * t67 + (t43 * t55 + t52 * t65) * t45;
t25 = -t45 * t46 * t44 + t48 * t47;
t51 = sin(qJ(4));
t10 = t19 * t54 + t25 * t51;
t50 = sin(qJ(5));
t70 = cos(qJ(5));
t9 = -t19 * t51 + t25 * t54;
t4 = t50 * t10 - t70 * t9;
t53 = cos(qJ(6));
t73 = t4 * t53;
t72 = t50 * pkin(4);
t59 = t70 * pkin(4);
t38 = -t59 - pkin(5);
t71 = pkin(5) - t38;
t26 = t54 * t47 - t51 * t67;
t27 = t51 * t47 + t54 * t67;
t14 = -t70 * t26 + t50 * t27;
t69 = t14 * t53;
t34 = t74 * t54;
t58 = t70 * t51;
t21 = t50 * t34 + t74 * t58;
t68 = t21 * t53;
t66 = t44 * t55;
t32 = t50 * t54 + t58;
t49 = sin(qJ(6));
t64 = t49 * t32;
t63 = t49 * t53;
t62 = t50 * t51;
t61 = t53 * t32;
t31 = -t70 * t54 + t62;
t60 = -0.2e1 * t32 * t31;
t57 = -pkin(5) * t32 - pkin(11) * t31;
t37 = pkin(11) + t72;
t56 = -t31 * t37 + t32 * t38;
t42 = t53 ^ 2;
t41 = t49 ^ 2;
t35 = 0.2e1 * t63;
t30 = t32 ^ 2;
t29 = t53 * t31;
t28 = t49 * t31;
t24 = t49 * t61;
t22 = t70 * t34 - t74 * t62;
t20 = t21 * t49;
t18 = -t48 * t66 + (t43 * t52 - t55 * t65) * t45;
t17 = (-t41 + t42) * t32;
t16 = t31 * pkin(5) - t32 * pkin(11) + t39;
t15 = t50 * t26 + t70 * t27;
t13 = t14 * t49;
t12 = t53 * t15 - t49 * t66;
t11 = -t49 * t15 - t53 * t66;
t7 = t49 * t16 + t53 * t22;
t6 = t53 * t16 - t49 * t22;
t5 = t70 * t10 + t50 * t9;
t3 = t4 * t49;
t2 = t18 * t49 + t53 * t5;
t1 = t18 * t53 - t49 * t5;
t8 = [1, t48 ^ 2 + (t43 ^ 2 + t46 ^ 2) * t45 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t18, -t19, 0, 0, 0, 0, 0, -t18 * t54, t18 * t51, 0, 0, 0, 0, 0, t18 * t31, t18 * t32, 0, 0, 0, 0, 0, t1 * t31 + t4 * t64, -t2 * t31 + t4 * t61; 0, 0, 0, t66, -t67, 0, 0, 0, 0, 0, t54 * t66, -t51 * t66, 0, 0, 0, 0, 0, -t31 * t66, -t32 * t66, 0, 0, 0, 0, 0, t11 * t31 + t14 * t64, -t12 * t31 + t14 * t61; 0, 0, 1, 0, 0, t51 ^ 2, t51 * t75, 0, 0, 0, pkin(3) * t75, -0.2e1 * pkin(3) * t51, t30, t60, 0, 0, 0, t31 * t76, t32 * t76, t42 * t30, -0.2e1 * t30 * t63, 0.2e1 * t31 * t61, t49 * t60, t31 ^ 2, 0.2e1 * t21 * t64 + 0.2e1 * t6 * t31, 0.2e1 * t21 * t61 - 0.2e1 * t7 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, 0, 0, 0, 0, 0, -t4, -t5, 0, 0, 0, 0, 0, -t73, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t69, t13; 0, 0, 0, 0, 0, 0, 0, t51, t54, 0, -t51 * pkin(9), -t54 * pkin(9), 0, 0, t32, -t31, 0, -t21, -t22, t24, t17, t28, t29, 0, t56 * t49 - t68, t56 * t53 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t72, t41, t35, 0, 0, 0, -0.2e1 * t38 * t53, 0.2e1 * t38 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, 0, 0, 0, 0, 0, -t73, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t69, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, -t21, -t22, t24, t17, t28, t29, 0, t57 * t49 - t68, t57 * t53 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t72, t41, t35, 0, 0, 0, t71 * t53, -t71 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t41, t35, 0, 0, 0, 0.2e1 * pkin(5) * t53, -0.2e1 * pkin(5) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t64, t31, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t53, 0, -t49 * t37, -t53 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t53, 0, -t49 * pkin(11), -t53 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;

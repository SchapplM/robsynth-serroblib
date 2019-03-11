% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t60 = cos(qJ(3));
t45 = -t60 * pkin(3) - pkin(2);
t55 = sin(qJ(4));
t56 = sin(qJ(3));
t59 = cos(qJ(4));
t68 = -t55 * t56 + t59 * t60;
t27 = -t68 * pkin(4) + t45;
t86 = 0.2e1 * t27;
t35 = t55 * t60 + t59 * t56;
t85 = 0.2e1 * t35;
t84 = 0.2e1 * t60;
t83 = pkin(8) + pkin(9);
t53 = sin(qJ(6));
t82 = pkin(5) * t53;
t54 = sin(qJ(5));
t81 = t54 * pkin(4);
t80 = t55 * pkin(3);
t58 = cos(qJ(6));
t36 = t83 * t56;
t37 = t83 * t60;
t25 = t55 * t36 - t59 * t37;
t15 = t68 * pkin(10) - t25;
t24 = -t59 * t36 - t55 * t37;
t62 = -t35 * pkin(10) + t24;
t78 = cos(qJ(5));
t8 = t54 * t15 - t78 * t62;
t79 = t8 * t58;
t52 = cos(pkin(6));
t51 = sin(pkin(6));
t74 = t51 * sin(qJ(2));
t32 = t52 * t60 - t56 * t74;
t33 = t52 * t56 + t60 * t74;
t16 = t59 * t32 - t55 * t33;
t17 = t55 * t32 + t59 * t33;
t10 = -t78 * t16 + t54 * t17;
t77 = t10 * t58;
t47 = t59 * pkin(3);
t44 = t47 + pkin(4);
t66 = -t78 * t44 + t54 * t80;
t28 = -pkin(5) + t66;
t76 = t28 * t58;
t46 = t78 * pkin(4);
t43 = -t46 - pkin(5);
t75 = t43 * t58;
t73 = t51 * cos(qJ(2));
t23 = t78 * t35 + t54 * t68;
t72 = t53 * t23;
t71 = t53 * t58;
t70 = t58 * t23;
t22 = t54 * t35 - t78 * t68;
t69 = -0.2e1 * t23 * t22;
t67 = t78 * t80;
t65 = -pkin(5) * t23 - pkin(11) * t22;
t31 = -t54 * t44 - t67;
t29 = pkin(11) - t31;
t64 = -t22 * t29 + t23 * t28;
t42 = pkin(11) + t81;
t63 = -t22 * t42 + t23 * t43;
t50 = t58 ^ 2;
t49 = t53 ^ 2;
t48 = pkin(5) * t58;
t40 = 0.2e1 * t71;
t39 = t43 * t53;
t26 = t28 * t53;
t21 = t23 ^ 2;
t20 = t58 * t22;
t19 = t53 * t22;
t18 = t53 * t70;
t12 = (-t49 + t50) * t23;
t11 = t54 * t16 + t78 * t17;
t9 = t78 * t15 + t54 * t62;
t7 = t22 * pkin(5) - t23 * pkin(11) + t27;
t6 = t10 * t53;
t5 = t8 * t53;
t4 = t58 * t11 - t53 * t73;
t3 = -t53 * t11 - t58 * t73;
t2 = t53 * t7 + t58 * t9;
t1 = -t53 * t9 + t58 * t7;
t13 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t73, -t74, 0, 0, 0, 0, 0, t60 * t73, -t56 * t73, 0, 0, 0, 0, 0, t68 * t73, -t35 * t73, 0, 0, 0, 0, 0, -t22 * t73, -t23 * t73, 0, 0, 0, 0, 0, t10 * t72 + t3 * t22, t10 * t70 - t4 * t22; 0, 1, 0, 0, t56 ^ 2, t56 * t84, 0, 0, 0, pkin(2) * t84, -0.2e1 * pkin(2) * t56, t35 ^ 2, t68 * t85, 0, 0, 0, -0.2e1 * t45 * t68, t45 * t85, t21, t69, 0, 0, 0, t22 * t86, t23 * t86, t50 * t21, -0.2e1 * t21 * t71, 0.2e1 * t22 * t70, t53 * t69, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t8 * t72, -0.2e1 * t2 * t22 + 0.2e1 * t8 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33, 0, 0, 0, 0, 0, t16, -t17, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t77, t6; 0, 0, 0, 0, 0, 0, t56, t60, 0, -t56 * pkin(8), -t60 * pkin(8), 0, 0, t35, t68, 0, t24, t25, 0, 0, t23, -t22, 0, -t8, -t9, t18, t12, t19, t20, 0, t64 * t53 - t79, t64 * t58 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t47, -0.2e1 * t80, 0, 0, 0, 0, 1, -0.2e1 * t66, 0.2e1 * t31, t49, t40, 0, 0, 0, -0.2e1 * t76, 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t77, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t68, 0, t24, t25, 0, 0, t23, -t22, 0, -t8, -t9, t18, t12, t19, t20, 0, t63 * t53 - t79, t63 * t58 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t47, -t80, 0, 0, 0, 0, 1, t46 - t66, -t67 + (-pkin(4) - t44) * t54, t49, t40, 0, 0, 0 (-t28 - t43) * t58, t39 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t81, t49, t40, 0, 0, 0, -0.2e1 * t75, 0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t77, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t8, -t9, t18, t12, t19, t20, 0, t65 * t53 - t79, t65 * t58 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t66, t31, t49, t40, 0, 0, 0, t48 - t76, t26 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, -t81, t49, t40, 0, 0, 0, t48 - t75, t39 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t49, t40, 0, 0, 0, 0.2e1 * t48, -0.2e1 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t72, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t58, 0, -t53 * t29, -t58 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t58, 0, -t53 * t42, -t58 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t58, 0, -t53 * pkin(11), -t58 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;

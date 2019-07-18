% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t67 = cos(qJ(3));
t68 = cos(qJ(2));
t28 = t44 * t45 - t67 * t68;
t77 = -0.2e1 * t28;
t42 = sin(qJ(5));
t43 = sin(qJ(4));
t46 = cos(qJ(5));
t47 = cos(qJ(4));
t78 = -t42 * t43 + t46 * t47;
t51 = t67 * pkin(1);
t38 = -t51 - pkin(2);
t70 = t47 * pkin(3);
t32 = t38 - t70;
t76 = 0.2e1 * t32;
t39 = -pkin(2) - t70;
t75 = 0.2e1 * t39;
t74 = t28 * pkin(3);
t73 = t42 * pkin(3);
t72 = t44 * pkin(1);
t71 = t46 * pkin(3);
t69 = pkin(2) - t38;
t30 = t44 * t68 + t67 * t45;
t52 = t68 * pkin(1);
t10 = t28 * pkin(2) - t30 * pkin(5) - t52;
t65 = t43 * t10;
t20 = t43 * t28;
t64 = t43 * t30;
t37 = pkin(5) + t72;
t63 = t43 * t37;
t62 = t43 * t47;
t61 = t46 * t43;
t21 = t47 * t28;
t59 = t47 * t30;
t58 = t47 * t37;
t57 = t32 + t39;
t56 = -0.2e1 * t20;
t55 = 0.2e1 * t21;
t54 = pkin(3) * t64;
t53 = t10 * t61;
t50 = -0.2e1 * t52;
t9 = t47 * t10;
t5 = t9 + t74;
t2 = -t42 * t65 + t46 * t5;
t49 = -pkin(2) * t30 - pkin(5) * t28;
t48 = -t28 * t37 + t30 * t38;
t29 = t42 * t47 + t61;
t41 = t47 ^ 2;
t40 = t43 ^ 2;
t34 = 0.2e1 * t62;
t26 = t30 ^ 2;
t25 = t29 ^ 2;
t24 = t28 ^ 2;
t23 = t78 * pkin(5);
t22 = t29 * pkin(5);
t19 = t43 * t59;
t18 = -t42 * t63 + t46 * t58;
t17 = t29 * t37;
t16 = t29 * t28;
t15 = t78 * t28;
t14 = 0.2e1 * t29 * t78;
t13 = t29 * t54;
t12 = t78 * t54;
t11 = (-t40 + t41) * t30;
t8 = t78 * t30;
t7 = t29 * t30;
t6 = t8 * t29;
t3 = t42 * t5 + t53;
t1 = -t29 * t7 + t78 * t8;
t4 = [1, 0, 0, t45 ^ 2, 0.2e1 * t45 * t68, 0, 0, 0, 0, 0, t26, t30 * t77, 0, 0, 0, t28 * t50, t30 * t50, t41 * t26, -0.2e1 * t26 * t62, t30 * t55, t30 * t56, t24, t10 * t55, t10 * t56, t8 ^ 2, -0.2e1 * t8 * t7, -t8 * t77, t7 * t77, t24, 0.2e1 * t2 * t28 + 0.2e1 * t7 * t54, -0.2e1 * t3 * t28 + 0.2e1 * t8 * t54; 0, 0, 0, 0, 0, t45, t68, 0, 0, 0, 0, 0, t30, -t28, 0, 0, 0, t19, t11, t20, t21, 0, t48 * t43, t48 * t47, t6, t1, t16, t15, 0, -t17 * t28 + t32 * t7 - t12, -t18 * t28 + t32 * t8 + t13; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t72, t40, t34, 0, 0, 0, -0.2e1 * t38 * t47, 0.2e1 * t38 * t43, t25, t14, 0, 0, 0, -t78 * t76, t29 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, 0, 0, t19, t11, t20, t21, 0, t49 * t43, t49 * t47, t6, t1, t16, t15, 0, -t22 * t28 + t39 * t7 - t12, -t23 * t28 + t39 * t8 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t72, t40, t34, 0, 0, 0, t69 * t47, -t69 * t43, t25, t14, 0, 0, 0, -t57 * t78, t57 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t40, t34, 0, 0, 0, 0.2e1 * pkin(2) * t47, -0.2e1 * pkin(2) * t43, t25, t14, 0, 0, 0, -t78 * t75, t29 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t64, t28, t9, -t65, 0, 0, t8, -t7, t28, t28 * t71 + t2, -t53 + (-t5 - t74) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t47, 0, -t63, -t58, 0, 0, t29, t78, 0, -t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t47, 0, -t43 * pkin(5), -t47 * pkin(5), 0, 0, t29, t78, 0, -t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, t28, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t78, 0, -t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t78, 0, -t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t4;

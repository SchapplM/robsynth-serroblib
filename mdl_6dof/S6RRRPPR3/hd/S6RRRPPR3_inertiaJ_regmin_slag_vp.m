% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = sin(qJ(3));
t68 = t39 * pkin(2);
t26 = qJ(4) + t68;
t76 = t26 ^ 2;
t40 = sin(qJ(2));
t60 = t39 * t40;
t65 = cos(qJ(3));
t66 = cos(qJ(2));
t17 = -t65 * t66 + t60;
t75 = 0.2e1 * t17;
t49 = t65 * t40;
t18 = t39 * t66 + t49;
t74 = -0.2e1 * t18;
t32 = -t66 * pkin(2) - pkin(1);
t73 = 0.2e1 * t32;
t38 = sin(qJ(6));
t72 = -0.2e1 * t38;
t41 = cos(qJ(6));
t71 = 0.2e1 * t41;
t70 = -pkin(4) - pkin(9);
t69 = -pkin(8) - pkin(7);
t51 = t66 * pkin(7);
t20 = t66 * pkin(8) + t51;
t12 = t65 * t20 + t69 * t60;
t7 = t17 * qJ(5) + t12;
t67 = t7 * t38;
t64 = t26 * t17;
t63 = t38 * t17;
t62 = t38 * t18;
t61 = t38 * t41;
t59 = t41 * t17;
t58 = t41 * t18;
t23 = pkin(5) + t26;
t37 = qJ(4) + pkin(5);
t57 = t23 + t37;
t56 = t17 * qJ(4);
t55 = t26 * qJ(4);
t54 = 0.2e1 * t66;
t53 = t17 * t74;
t52 = t38 * t59;
t50 = t65 * pkin(2);
t11 = t39 * t20 - t69 * t49;
t48 = -t17 * pkin(3) + t18 * qJ(4) - t32;
t30 = t50 + pkin(3);
t24 = -pkin(4) - t30;
t22 = -pkin(9) + t24;
t47 = t17 * t23 - t18 * t22;
t34 = -pkin(3) + t70;
t46 = t17 * t37 - t18 * t34;
t45 = 0.2e1 * pkin(3);
t44 = qJ(4) ^ 2;
t43 = 0.2e1 * qJ(4);
t42 = -pkin(3) - pkin(4);
t36 = t41 ^ 2;
t35 = t38 ^ 2;
t28 = t43 + t68;
t25 = 0.2e1 * t61;
t21 = 0.2e1 * t26;
t16 = t18 ^ 2;
t15 = t17 ^ 2;
t10 = (t35 - t36) * t17;
t6 = -t18 * qJ(5) + t11;
t5 = t7 * t41;
t4 = -t17 * pkin(4) + t48;
t3 = t18 * pkin(5) + t70 * t17 + t48;
t2 = t38 * t3 + t41 * t6;
t1 = t41 * t3 - t38 * t6;
t8 = [1, 0, 0, t40 ^ 2, t40 * t54, 0, 0, 0, pkin(1) * t54, -0.2e1 * pkin(1) * t40, t16, t53, 0, 0, 0, t17 * t73, t18 * t73, -t48 * t75, 0.2e1 * t11 * t18 - 0.2e1 * t12 * t17, -t48 * t74, t11 ^ 2 + t12 ^ 2 + t48 ^ 2, 0.2e1 * t4 * t18, t4 * t75, 0.2e1 * t7 * t17 - 0.2e1 * t6 * t18, t4 ^ 2 + t6 ^ 2 + t7 ^ 2, t36 * t15, -0.2e1 * t15 * t61, t58 * t75, t38 * t53, t16, 0.2e1 * t1 * t18 + 0.2e1 * t7 * t63, -0.2e1 * t2 * t18 + 0.2e1 * t7 * t59; 0, 0, 0, 0, 0, t40, t66, 0, -t40 * pkin(7), -t51, 0, 0, t18, -t17, 0, -t11, -t12, -t11, -t30 * t18 - t64, t12, -t11 * t30 + t12 * t26, t7, t6, -t24 * t18 + t64, t6 * t24 + t7 * t26, -t52, t10, -t62, -t58, 0, t47 * t38 + t5, t47 * t41 - t67; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t68, 0.2e1 * t30, 0, t21, t30 ^ 2 + t76, t21, 0.2e1 * t24, 0, t24 ^ 2 + t76, t35, t25, 0, 0, 0, t23 * t71, t23 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, -t11, -t12, -t11, -pkin(3) * t18 - t56, t12, -t11 * pkin(3) + t12 * qJ(4), t7, t6, -t42 * t18 + t56, t7 * qJ(4) + t6 * t42, -t52, t10, -t62, -t58, 0, t46 * t38 + t5, t46 * t41 - t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t68, t45 + t50, 0, t28, t30 * pkin(3) + t55, t28, -0.2e1 * pkin(3) - 0.2e1 * pkin(4) - t50, 0, t24 * t42 + t55, t35, t25, 0, 0, 0, t57 * t41, -t57 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, 0, t43, pkin(3) ^ 2 + t44, t43, 0.2e1 * t42, 0, t42 ^ 2 + t44, t35, t25, 0, 0, 0, t37 * t71, t37 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t11, 0, 0, -t18, t6, 0, 0, 0, 0, 0, -t62, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t30, 0, 1, 0, t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 1, 0, t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, 0, t4, 0, 0, 0, 0, 0, t58, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t63, t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t41, 0, -t38 * t22, -t41 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t41, 0, -t38 * t34, -t41 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;

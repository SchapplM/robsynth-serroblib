% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t40 = sin(pkin(10));
t41 = cos(pkin(10));
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t25 = t44 * t40 - t47 * t41;
t26 = t47 * t40 + t44 * t41;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t20 = t46 * t25 + t43 * t26;
t21 = -t43 * t25 + t46 * t26;
t34 = -t41 * pkin(2) - pkin(1);
t22 = t25 * pkin(3) + t34;
t52 = -t21 * qJ(5) + t22;
t8 = t20 * pkin(4) + t52;
t74 = -0.2e1 * t8;
t73 = 0.2e1 * t22;
t69 = t43 * pkin(3);
t32 = qJ(5) + t69;
t72 = 0.2e1 * t32;
t71 = 0.2e1 * t34;
t49 = 0.2e1 * qJ(5);
t70 = pkin(4) + pkin(9);
t68 = t46 * pkin(3);
t67 = t21 * t20;
t66 = t32 * t20;
t42 = sin(qJ(6));
t65 = t42 * t20;
t64 = t42 * t21;
t45 = cos(qJ(6));
t63 = t45 * t20;
t17 = t45 * t21;
t62 = t45 * t42;
t61 = pkin(7) + qJ(2);
t60 = t40 ^ 2 + t41 ^ 2;
t59 = qJ(5) * t20;
t58 = qJ(5) + t32;
t57 = 0.2e1 * t67;
t35 = -pkin(4) - t68;
t27 = t61 * t40;
t28 = t61 * t41;
t56 = -t47 * t27 - t44 * t28;
t14 = -t26 * pkin(8) + t56;
t54 = t44 * t27 - t47 * t28;
t15 = -t25 * pkin(8) - t54;
t9 = -t46 * t14 + t43 * t15;
t10 = t43 * t14 + t46 * t15;
t29 = -pkin(9) + t35;
t55 = -t21 * t29 + t66;
t53 = t21 * t70 + t59;
t51 = -0.2e1 * pkin(4);
t39 = t45 ^ 2;
t38 = t42 ^ 2;
t30 = -0.2e1 * t62;
t19 = t21 ^ 2;
t18 = t20 ^ 2;
t16 = t20 * t62;
t12 = (-t38 + t39) * t20;
t7 = -t20 * pkin(5) + t10;
t6 = t21 * pkin(5) + t9;
t5 = t70 * t20 + t52;
t4 = t7 * t45;
t3 = t7 * t42;
t2 = t42 * t6 + t45 * t5;
t1 = -t42 * t5 + t45 * t6;
t11 = [1, 0, 0, 0.2e1 * pkin(1) * t41, -0.2e1 * pkin(1) * t40, 0.2e1 * t60 * qJ(2), t60 * qJ(2) ^ 2 + pkin(1) ^ 2, t26 ^ 2, -0.2e1 * t26 * t25, 0, 0, 0, t25 * t71, t26 * t71, t19, -0.2e1 * t67, 0, 0, 0, t20 * t73, t21 * t73, -0.2e1 * t10 * t20 + 0.2e1 * t9 * t21, t20 * t74, t21 * t74, t10 ^ 2 + t8 ^ 2 + t9 ^ 2, t38 * t18, 0.2e1 * t18 * t62, t42 * t57, t45 * t57, t19, 0.2e1 * t1 * t21 - 0.2e1 * t7 * t63, -0.2e1 * t2 * t21 + 0.2e1 * t7 * t65; 0, 0, 0, -t41, t40, 0, -pkin(1), 0, 0, 0, 0, 0, t25, t26, 0, 0, 0, 0, 0, t20, t21, 0, -t20, -t21, t8, 0, 0, 0, 0, 0, -t64, -t17; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t56, t54, 0, 0, t21, -t20, 0, -t9, -t10, t35 * t21 - t66, t9, t10, t10 * t32 + t9 * t35, t16, t12, t17, -t64, 0, -t55 * t45 + t3, t55 * t42 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t69, 0, 0.2e1 * t35, t72, t32 ^ 2 + t35 ^ 2, t39, t30, 0, 0, 0, t42 * t72, t45 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t9, -t10, -pkin(4) * t21 - t59, t9, t10, -t9 * pkin(4) + t10 * qJ(5), t16, t12, t17, -t64, 0, -t53 * t45 + t3, t53 * t42 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t68, -t69, 0, t51 - t68, t49 + t69, -t35 * pkin(4) + t32 * qJ(5), t39, t30, 0, 0, 0, t58 * t42, t58 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t51, t49, pkin(4) ^ 2 + qJ(5) ^ 2, t39, t30, 0, 0, 0, t42 * t49, t45 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, t9, 0, 0, 0, 0, 0, t17, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t63, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, t45 * t29, -t42 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, -t45 * t70, t42 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;

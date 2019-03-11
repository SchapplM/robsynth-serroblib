% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = sin(pkin(10));
t46 = cos(pkin(10));
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t32 = -t45 * t49 + t46 * t52;
t33 = t45 * t52 + t46 * t49;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t20 = -t51 * t32 + t48 * t33;
t21 = t48 * t32 + t51 * t33;
t42 = -t52 * pkin(2) - pkin(1);
t24 = -t32 * pkin(3) + t42;
t57 = -t21 * qJ(5) + t24;
t8 = t20 * pkin(4) + t57;
t76 = -0.2e1 * t8;
t75 = 0.2e1 * t24;
t41 = t46 * pkin(2) + pkin(3);
t71 = pkin(2) * t45;
t58 = t48 * t41 + t51 * t71;
t26 = qJ(5) + t58;
t74 = 0.2e1 * t26;
t73 = 0.2e1 * t52;
t54 = 0.2e1 * qJ(5);
t72 = pkin(4) + pkin(9);
t70 = t21 * t20;
t69 = t26 * t20;
t47 = sin(qJ(6));
t68 = t47 * t20;
t67 = t47 * t21;
t50 = cos(qJ(6));
t66 = t50 * t20;
t17 = t50 * t21;
t65 = t50 * t47;
t64 = -qJ(3) - pkin(7);
t37 = t64 * t49;
t38 = t64 * t52;
t23 = t45 * t37 - t46 * t38;
t63 = qJ(5) * t20;
t62 = qJ(5) + t26;
t61 = 0.2e1 * t70;
t22 = t46 * t37 + t45 * t38;
t14 = -t33 * pkin(8) + t22;
t15 = t32 * pkin(8) + t23;
t9 = -t51 * t14 + t48 * t15;
t29 = t51 * t41 - t48 * t71;
t28 = -pkin(4) - t29;
t10 = t48 * t14 + t51 * t15;
t25 = -pkin(9) + t28;
t60 = -t21 * t25 + t69;
t59 = t21 * t72 + t63;
t55 = -0.2e1 * pkin(4);
t44 = t50 ^ 2;
t43 = t47 ^ 2;
t40 = -0.2e1 * t65;
t19 = t21 ^ 2;
t18 = t20 ^ 2;
t16 = t20 * t65;
t12 = (-t43 + t44) * t20;
t7 = -t20 * pkin(5) + t10;
t6 = t21 * pkin(5) + t9;
t5 = t72 * t20 + t57;
t4 = t7 * t50;
t3 = t7 * t47;
t2 = t47 * t6 + t50 * t5;
t1 = -t47 * t5 + t50 * t6;
t11 = [1, 0, 0, t49 ^ 2, t49 * t73, 0, 0, 0, pkin(1) * t73, -0.2e1 * pkin(1) * t49, -0.2e1 * t22 * t33 + 0.2e1 * t23 * t32, t22 ^ 2 + t23 ^ 2 + t42 ^ 2, t19, -0.2e1 * t70, 0, 0, 0, t20 * t75, t21 * t75, -0.2e1 * t10 * t20 + 0.2e1 * t9 * t21, t20 * t76, t21 * t76, t10 ^ 2 + t8 ^ 2 + t9 ^ 2, t43 * t18, 0.2e1 * t18 * t65, t47 * t61, t50 * t61, t19, 0.2e1 * t1 * t21 - 0.2e1 * t7 * t66, -0.2e1 * t2 * t21 + 0.2e1 * t7 * t68; 0, 0, 0, 0, 0, t49, t52, 0, -t49 * pkin(7), -t52 * pkin(7) (t32 * t45 - t33 * t46) * pkin(2) (t22 * t46 + t23 * t45) * pkin(2), 0, 0, t21, -t20, 0, -t9, -t10, t28 * t21 - t69, t9, t10, t10 * t26 + t9 * t28, t16, t12, t17, -t67, 0, -t60 * t50 + t3, t60 * t47 + t4; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t45 ^ 2 + t46 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t58, 0, 0.2e1 * t28, t74, t26 ^ 2 + t28 ^ 2, t44, t40, 0, 0, 0, t47 * t74, t50 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t20, t21, 0, -t20, -t21, t8, 0, 0, 0, 0, 0, -t67, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, -t9, -t10, -pkin(4) * t21 - t63, t9, t10, -t9 * pkin(4) + t10 * qJ(5), t16, t12, t17, -t67, 0, -t59 * t50 + t3, t59 * t47 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t58, 0, -t29 + t55, t54 + t58, -t28 * pkin(4) + t26 * qJ(5), t44, t40, 0, 0, 0, t62 * t47, t62 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t55, t54, pkin(4) ^ 2 + qJ(5) ^ 2, t44, t40, 0, 0, 0, t47 * t54, t50 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, t9, 0, 0, 0, 0, 0, t17, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t66, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, t50 * t25, -t47 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, -t50 * t72, t47 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;

% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = sin(pkin(11));
t47 = cos(pkin(11));
t50 = sin(qJ(6));
t53 = cos(qJ(6));
t26 = -t50 * t45 + t53 * t47;
t48 = cos(pkin(6));
t51 = sin(qJ(3));
t54 = cos(qJ(3));
t46 = sin(pkin(6));
t71 = t46 * sin(qJ(2));
t21 = t48 * t51 + t54 * t71;
t80 = t21 ^ 2;
t36 = t45 * pkin(5) + qJ(4);
t79 = 0.2e1 * t36;
t78 = -0.2e1 * t51;
t77 = 0.2e1 * t51;
t76 = 0.2e1 * t54;
t75 = 0.2e1 * qJ(4);
t49 = -pkin(3) - qJ(5);
t74 = -pkin(9) + t49;
t25 = t53 * t45 + t50 * t47;
t73 = t25 * t51;
t72 = t45 * t54;
t55 = cos(qJ(2));
t70 = t46 * t55;
t69 = t47 * t54;
t61 = -t51 * qJ(4) - pkin(2);
t23 = t49 * t54 + t61;
t37 = t51 * pkin(8);
t32 = t51 * pkin(4) + t37;
t10 = t47 * t23 + t45 * t32;
t38 = t54 * pkin(8);
t33 = t54 * pkin(4) + t38;
t35 = t45 ^ 2 + t47 ^ 2;
t43 = t51 ^ 2;
t66 = t54 ^ 2 + t43;
t65 = qJ(4) * t54;
t64 = t21 * qJ(4);
t63 = t51 * t70;
t62 = t54 * t70;
t28 = t47 * t32;
t9 = -t45 * t23 + t28;
t3 = t10 * t45 + t9 * t47;
t60 = -pkin(3) * t51 + t65;
t20 = -t48 * t54 + t51 * t71;
t13 = t20 * t47 + t45 * t70;
t14 = t20 * t45 - t47 * t70;
t4 = t13 * t47 + t14 * t45;
t59 = t20 * t51 + t21 * t54;
t58 = t49 * t51 + t65;
t56 = qJ(4) ^ 2;
t31 = -t54 * pkin(3) + t61;
t30 = t74 * t47;
t29 = t74 * t45;
t24 = t35 * t49;
t19 = pkin(5) * t69 + t33;
t18 = t26 * t51;
t16 = t25 * t54;
t15 = t26 * t54;
t12 = t53 * t29 + t50 * t30;
t11 = -t50 * t29 + t53 * t30;
t8 = -pkin(9) * t69 + t10;
t7 = t51 * pkin(5) + t28 + (pkin(9) * t54 - t23) * t45;
t6 = t50 * t13 + t53 * t14;
t5 = t53 * t13 - t50 * t14;
t2 = t50 * t7 + t53 * t8;
t1 = -t50 * t8 + t53 * t7;
t17 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 ^ 2 * t55 ^ 2 + t20 ^ 2 + t80, 0, 0, 0, t13 ^ 2 + t14 ^ 2 + t80, 0, 0, 0, 0, 0, 0, 0; 0, 0, t70, -t71, 0, 0, 0, 0, 0, t62, -t63, t59, -t62, t63, t59 * pkin(8) - t31 * t70, t13 * t51 + t21 * t69, -t14 * t51 - t21 * t72 (t13 * t45 - t14 * t47) * t54, t14 * t10 + t13 * t9 + t21 * t33, 0, 0, 0, 0, 0, t21 * t15 + t5 * t51, -t21 * t16 - t6 * t51; 0, 1, 0, 0, t43, t51 * t76, 0, 0, 0, pkin(2) * t76, pkin(2) * t78, 0.2e1 * t66 * pkin(8), t31 * t76, t31 * t78, t66 * pkin(8) ^ 2 + t31 ^ 2, 0.2e1 * t33 * t69 + 0.2e1 * t9 * t51, -0.2e1 * t10 * t51 - 0.2e1 * t33 * t72 (-t10 * t47 + t45 * t9) * t76, t10 ^ 2 + t33 ^ 2 + t9 ^ 2, t16 ^ 2, 0.2e1 * t16 * t15, -t16 * t77, -t15 * t77, t43, 0.2e1 * t1 * t51 + 0.2e1 * t19 * t15, -0.2e1 * t19 * t16 - 0.2e1 * t2 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, t20, t21, -t20 * pkin(3) + t64, t21 * t45, t21 * t47, -t4, t4 * t49 + t64, 0, 0, 0, 0, 0, t21 * t25, t21 * t26; 0, 0, 0, 0, 0, 0, t51, t54, 0, -t37, -t38, t60, t37, t38, t60 * pkin(8), t33 * t45 + t58 * t47, t33 * t47 - t58 * t45, -t3, t33 * qJ(4) + t3 * t49, -t16 * t26, -t26 * t15 + t16 * t25, t18, -t73, 0, t11 * t51 + t36 * t15 + t19 * t25, -t12 * t51 - t36 * t16 + t19 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t75, pkin(3) ^ 2 + t56, t45 * t75, t47 * t75, -0.2e1 * t24, t35 * t49 ^ 2 + t56, t26 ^ 2, -0.2e1 * t26 * t25, 0, 0, 0, t25 * t79, t26 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, t37, t47 * t51, -t45 * t51, 0, t3, 0, 0, 0, 0, 0, t18, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, -t35, t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t72, 0, t33, 0, 0, 0, 0, 0, t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t47, 0, qJ(4), 0, 0, 0, 0, 0, t25, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t15, t51, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t17;

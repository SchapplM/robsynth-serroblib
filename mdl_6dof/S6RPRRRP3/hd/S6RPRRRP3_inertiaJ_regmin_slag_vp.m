% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(qJ(5));
t45 = sin(qJ(4));
t47 = cos(qJ(5));
t48 = cos(qJ(4));
t24 = t44 * t48 + t47 * t45;
t79 = -0.2e1 * t24;
t36 = -t48 * pkin(4) - pkin(3);
t78 = 0.2e1 * t36;
t49 = cos(qJ(3));
t77 = -0.2e1 * t49;
t76 = -pkin(9) - pkin(8);
t43 = cos(pkin(10));
t32 = -t43 * pkin(1) - pkin(2);
t46 = sin(qJ(3));
t22 = -t49 * pkin(3) - t46 * pkin(8) + t32;
t42 = sin(pkin(10));
t30 = t42 * pkin(1) + pkin(7);
t60 = t49 * t30;
t56 = t48 * t60;
t10 = t56 + (-pkin(9) * t46 + t22) * t45;
t20 = t48 * t22;
t63 = t48 * t46;
t67 = t30 * t45;
t8 = -pkin(9) * t63 + t20 + (-pkin(4) - t67) * t49;
t4 = t47 * t10 + t44 * t8;
t75 = pkin(3) * t48;
t74 = t47 * pkin(4);
t73 = t49 * pkin(4);
t72 = t49 * pkin(5);
t27 = t76 * t48;
t55 = t76 * t45;
t15 = -t44 * t27 - t47 * t55;
t71 = t15 * t49;
t16 = -t47 * t27 + t44 * t55;
t70 = t16 * t49;
t18 = t24 * t46;
t69 = t18 * t24;
t68 = t24 * t49;
t66 = t45 * t46;
t65 = t45 * t48;
t64 = t45 * t49;
t62 = t48 * t49;
t23 = t44 * t45 - t47 * t48;
t61 = t49 * t23;
t59 = t49 * t46;
t26 = t46 * t30;
t21 = pkin(4) * t66 + t26;
t58 = t49 * qJ(6);
t57 = 0.2e1 * t59;
t54 = t44 * t10 - t47 * t8;
t19 = -t44 * t66 + t47 * t63;
t53 = -t18 * pkin(5) + t19 * qJ(6);
t51 = 0.2e1 * pkin(5);
t50 = 0.2e1 * qJ(6);
t41 = t49 ^ 2;
t40 = t48 ^ 2;
t39 = t46 ^ 2;
t38 = t45 ^ 2;
t37 = t44 * pkin(4);
t34 = pkin(5) + t74;
t31 = t37 + qJ(6);
t17 = t19 ^ 2;
t14 = t45 * t22 + t56;
t13 = -t45 * t60 + t20;
t12 = t19 * t23;
t11 = t23 * pkin(5) - t24 * qJ(6) + t36;
t5 = -t53 + t21;
t2 = t54 + t72;
t1 = -t58 + t4;
t3 = [1, 0, 0 (t42 ^ 2 + t43 ^ 2) * pkin(1) ^ 2, t39, t57, 0, 0, 0, t32 * t77, 0.2e1 * t32 * t46, t40 * t39, -0.2e1 * t39 * t65, -0.2e1 * t48 * t59, t45 * t57, t41, -0.2e1 * t13 * t49 + 0.2e1 * t39 * t67, 0.2e1 * t39 * t30 * t48 + 0.2e1 * t14 * t49, t17, -0.2e1 * t18 * t19, t19 * t77, 0.2e1 * t49 * t18, t41, 0.2e1 * t21 * t18 + 0.2e1 * t49 * t54, 0.2e1 * t21 * t19 + 0.2e1 * t4 * t49, 0.2e1 * t5 * t18 + 0.2e1 * t2 * t49, -0.2e1 * t1 * t18 + 0.2e1 * t2 * t19, -0.2e1 * t1 * t49 - 0.2e1 * t5 * t19, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t19 + t2 * t18 - t5 * t49; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 + t17 + t41; 0, 0, 0, 0, 0, 0, t46, t49, 0, -t26, -t60, t45 * t63 (-t38 + t40) * t46, -t64, -t62, 0, -t30 * t63 + (-pkin(3) * t46 + pkin(8) * t49) * t45, pkin(8) * t62 + (t67 - t75) * t46, t19 * t24, -t12 - t69, -t68, t61, 0, t36 * t18 + t21 * t23 + t71, t36 * t19 + t21 * t24 + t70, t11 * t18 + t5 * t23 + t71, -t1 * t23 + t15 * t19 - t16 * t18 + t2 * t24, -t11 * t19 - t5 * t24 - t70, t1 * t16 + t5 * t11 + t2 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46, 0, 0, 0, 0, 0, t62, -t64, 0, 0, 0, 0, 0, -t61, -t68, -t61, -t12 + t69, t68, -t49 * t11 + t18 * t15 + t19 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, 0.2e1 * t65, 0, 0, 0, 0.2e1 * t75, -0.2e1 * pkin(3) * t45, t24 ^ 2, t23 * t79, 0, 0, 0, t23 * t78, t24 * t78, 0.2e1 * t11 * t23, 0.2e1 * t15 * t24 - 0.2e1 * t16 * t23, t11 * t79, t11 ^ 2 + t15 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t66, -t49, t13, -t14, 0, 0, t19, -t18, -t49, -t47 * t73 - t54, t44 * t73 - t4 (-pkin(5) - t34) * t49 - t54, -t31 * t18 - t34 * t19 (-qJ(6) - t31) * t49 + t4, t1 * t31 - t2 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t63, 0, 0, 0, 0, 0, -t18, -t19, -t18, 0, t19, -t18 * t34 + t19 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t48, 0, -t45 * pkin(8), -t48 * pkin(8), 0, 0, t24, -t23, 0, -t15, -t16, -t15, -t31 * t23 - t34 * t24, t16, -t15 * t34 + t16 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t74, -0.2e1 * t37, 0.2e1 * t34, 0, 0.2e1 * t31, t31 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, -t49, -t54, -t4, -t54 - 0.2e1 * t72, -pkin(5) * t19 - t18 * qJ(6), -0.2e1 * t58 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, -t18, 0, t19, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t15, -t16, -t15, -pkin(5) * t24 - t23 * qJ(6), t16, -t15 * pkin(5) + t16 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t74, -t37, t51 + t74, 0, t50 + t37, t34 * pkin(5) + t31 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t51, 0, t50, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;

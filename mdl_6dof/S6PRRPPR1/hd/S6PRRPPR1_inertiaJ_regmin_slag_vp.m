% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t52 = sin(qJ(3));
t53 = sin(qJ(2));
t48 = sin(pkin(6));
t75 = cos(qJ(3));
t64 = t48 * t75;
t67 = cos(pkin(6));
t28 = t67 * t52 + t53 * t64;
t47 = sin(pkin(11));
t50 = cos(pkin(11));
t71 = t48 * t53;
t57 = -t52 * t71 + t67 * t75;
t14 = t47 * t28 - t50 * t57;
t80 = t14 ^ 2;
t65 = t75 * pkin(8);
t37 = t75 * qJ(4) + t65;
t63 = (-qJ(4) - pkin(8)) * t52;
t24 = t47 * t37 - t50 * t63;
t79 = t24 ^ 2;
t33 = t47 * t75 + t50 * t52;
t46 = sin(pkin(12));
t49 = cos(pkin(12));
t51 = sin(qJ(6));
t54 = cos(qJ(6));
t34 = t54 * t46 + t51 * t49;
t12 = t34 * t33;
t78 = -0.2e1 * t12;
t42 = -t50 * pkin(3) - pkin(4);
t36 = -t49 * pkin(5) + t42;
t77 = 0.2e1 * t36;
t39 = t47 * pkin(3) + qJ(5);
t76 = pkin(9) + t39;
t74 = t14 * t24;
t31 = t47 * t52 - t50 * t75;
t73 = t34 * t31;
t72 = t46 * t33;
t55 = cos(qJ(2));
t70 = t48 * t55;
t69 = t49 * t33;
t43 = -t75 * pkin(3) - pkin(2);
t21 = t31 * pkin(4) - t33 * qJ(5) + t43;
t26 = t50 * t37 + t47 * t63;
t8 = t46 * t21 + t49 * t26;
t68 = t46 ^ 2 + t49 ^ 2;
t66 = 0.2e1 * t75;
t7 = t49 * t21 - t46 * t26;
t62 = t8 * t46 + t7 * t49;
t61 = -t7 * t46 + t8 * t49;
t16 = t50 * t28 + t47 * t57;
t10 = t49 * t16 - t46 * t70;
t9 = -t46 * t16 - t49 * t70;
t60 = t10 * t49 - t9 * t46;
t59 = t10 * t46 + t9 * t49;
t58 = -t31 * t39 + t33 * t42;
t32 = t51 * t46 - t54 * t49;
t30 = t76 * t49;
t29 = t76 * t46;
t23 = t32 * t31;
t20 = -t51 * t29 + t54 * t30;
t19 = -t54 * t29 - t51 * t30;
t13 = t32 * t33;
t11 = pkin(5) * t72 + t24;
t6 = -pkin(9) * t72 + t8;
t5 = t31 * pkin(5) - pkin(9) * t69 + t7;
t4 = t54 * t10 + t51 * t9;
t3 = -t51 * t10 + t54 * t9;
t2 = t51 * t5 + t54 * t6;
t1 = t54 * t5 - t51 * t6;
t15 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 ^ 2 * t55 ^ 2 + t16 ^ 2 + t80, 0, 0, 0, t10 ^ 2 + t9 ^ 2 + t80, 0, 0, 0, 0, 0, 0, 0; 0, 0, t70, -t71, 0, 0, 0, 0, 0, t55 * t64, -t52 * t70, t14 * t33 - t16 * t31, t16 * t26 - t43 * t70 + t74, t14 * t72 + t9 * t31, -t10 * t31 + t14 * t69, -t59 * t33, t10 * t8 + t9 * t7 + t74, 0, 0, 0, 0, 0, t14 * t12 + t3 * t31, -t14 * t13 - t4 * t31; 0, 1, 0, 0, t52 ^ 2, t52 * t66, 0, 0, 0, pkin(2) * t66, -0.2e1 * pkin(2) * t52, 0.2e1 * t24 * t33 - 0.2e1 * t26 * t31, t26 ^ 2 + t43 ^ 2 + t79, 0.2e1 * t24 * t72 + 0.2e1 * t7 * t31, 0.2e1 * t24 * t69 - 0.2e1 * t8 * t31, -0.2e1 * t62 * t33, t7 ^ 2 + t8 ^ 2 + t79, t13 ^ 2, -t13 * t78, -0.2e1 * t13 * t31, t31 * t78, t31 ^ 2, 0.2e1 * t1 * t31 + 0.2e1 * t11 * t12, -0.2e1 * t11 * t13 - 0.2e1 * t2 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t28, 0 (-t14 * t50 + t16 * t47) * pkin(3), -t14 * t49, t14 * t46, t60, t14 * t42 + t60 * t39, 0, 0, 0, 0, 0, t14 * t32, t14 * t34; 0, 0, 0, 0, 0, 0, t52, t75, 0, -t52 * pkin(8), -t65 (-t31 * t47 - t33 * t50) * pkin(3) (-t24 * t50 + t26 * t47) * pkin(3), -t24 * t49 + t58 * t46, t24 * t46 + t58 * t49, t61, t24 * t42 + t61 * t39, -t13 * t34, -t34 * t12 + t13 * t32, t73, -t23, 0, t11 * t32 + t36 * t12 + t19 * t31, t11 * t34 - t36 * t13 - t20 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t47 ^ 2 + t50 ^ 2) * pkin(3) ^ 2, -0.2e1 * t42 * t49, 0.2e1 * t42 * t46, 0.2e1 * t68 * t39, t68 * t39 ^ 2 + t42 ^ 2, t34 ^ 2, -0.2e1 * t34 * t32, 0, 0, 0, t32 * t77, t34 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t49 * t31, -t46 * t31, -t68 * t33, t62, 0, 0, 0, 0, 0, -t23, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t69, 0, t24, 0, 0, 0, 0, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t46, 0, t42, 0, 0, 0, 0, 0, t32, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, t31, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t15;

% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t56 = sin(pkin(10));
t57 = cos(pkin(10));
t60 = sin(qJ(2));
t82 = cos(qJ(2));
t36 = t56 * t60 - t57 * t82;
t87 = -0.2e1 * t36;
t86 = 0.2e1 * t36;
t58 = sin(qJ(5));
t61 = cos(qJ(4));
t59 = sin(qJ(4));
t81 = cos(qJ(5));
t68 = t81 * t59;
t41 = t58 * t61 + t68;
t85 = -0.2e1 * t41;
t48 = -t57 * pkin(2) - pkin(3);
t42 = -t61 * pkin(4) + t48;
t84 = 0.2e1 * t42;
t37 = t56 * t82 + t57 * t60;
t52 = -t82 * pkin(2) - pkin(1);
t21 = t36 * pkin(3) - t37 * pkin(8) + t52;
t43 = (-qJ(3) - pkin(7)) * t60;
t70 = t82 * pkin(7);
t44 = t82 * qJ(3) + t70;
t24 = t56 * t43 + t57 * t44;
t73 = t61 * t24;
t10 = t73 + (-pkin(9) * t37 + t21) * t59;
t11 = t61 * t21 - t59 * t24;
t72 = t61 * t37;
t8 = t36 * pkin(4) - pkin(9) * t72 + t11;
t4 = t81 * t10 + t58 * t8;
t32 = t36 * pkin(5);
t53 = t58 * pkin(4);
t46 = t56 * pkin(2) + pkin(8);
t83 = pkin(9) + t46;
t33 = t83 * t61;
t19 = t58 * t33 + t83 * t68;
t80 = t19 * t36;
t77 = t58 * t59;
t20 = t81 * t33 - t83 * t77;
t79 = t20 * t36;
t67 = t81 * t61;
t75 = t59 * t37;
t16 = t37 * t67 - t58 * t75;
t40 = -t67 + t77;
t78 = t40 * t16;
t25 = t40 * t36;
t26 = t41 * t36;
t76 = t59 * t36;
t74 = t59 * t61;
t71 = 0.2e1 * t82;
t31 = t36 * qJ(6);
t1 = t31 + t4;
t69 = t81 * pkin(4);
t3 = -t58 * t10 + t81 * t8;
t22 = -t57 * t43 + t56 * t44;
t2 = -t3 - t32;
t14 = pkin(4) * t75 + t22;
t66 = -t40 * pkin(5) + t41 * qJ(6);
t65 = -t36 * t46 + t37 * t48;
t63 = 0.2e1 * pkin(5);
t62 = 0.2e1 * qJ(6);
t55 = t61 ^ 2;
t54 = t59 ^ 2;
t50 = t69 + pkin(5);
t47 = t53 + qJ(6);
t38 = t41 ^ 2;
t35 = t37 ^ 2;
t34 = t36 ^ 2;
t30 = t61 * t36;
t17 = t42 - t66;
t15 = t41 * t37;
t13 = t41 * t15;
t12 = t59 * t21 + t73;
t5 = t15 * pkin(5) - t16 * qJ(6) + t14;
t6 = [1, 0, 0, t60 ^ 2, t60 * t71, 0, 0, 0, pkin(1) * t71, -0.2e1 * pkin(1) * t60, 0.2e1 * t22 * t37 - 0.2e1 * t24 * t36, t22 ^ 2 + t24 ^ 2 + t52 ^ 2, t55 * t35, -0.2e1 * t35 * t74, t72 * t86, t75 * t87, t34, 0.2e1 * t11 * t36 + 0.2e1 * t22 * t75, -0.2e1 * t12 * t36 + 0.2e1 * t22 * t72, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t86, t15 * t87, t34, 0.2e1 * t14 * t15 + 0.2e1 * t3 * t36, 0.2e1 * t14 * t16 - 0.2e1 * t4 * t36, 0.2e1 * t5 * t15 - 0.2e1 * t2 * t36, -0.2e1 * t1 * t15 + 0.2e1 * t2 * t16, 0.2e1 * t1 * t36 - 0.2e1 * t5 * t16, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t60, t82, 0, -t60 * pkin(7), -t70 (-t36 * t56 - t37 * t57) * pkin(2) (-t22 * t57 + t24 * t56) * pkin(2), t59 * t72 (-t54 + t55) * t37, t76, t30, 0, -t22 * t61 + t59 * t65, t22 * t59 + t61 * t65, t16 * t41, -t13 - t78, t26, -t25, 0, t14 * t40 + t42 * t15 - t80, t14 * t41 + t42 * t16 - t79, t17 * t15 + t5 * t40 - t80, -t1 * t40 - t20 * t15 + t19 * t16 + t2 * t41, -t17 * t16 - t5 * t41 + t79, t1 * t20 + t5 * t17 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t56 ^ 2 + t57 ^ 2) * pkin(2) ^ 2, t54, 0.2e1 * t74, 0, 0, 0, -0.2e1 * t48 * t61, 0.2e1 * t48 * t59, t38, t40 * t85, 0, 0, 0, t40 * t84, t41 * t84, 0.2e1 * t17 * t40, 0.2e1 * t19 * t41 - 0.2e1 * t20 * t40, t17 * t85, t17 ^ 2 + t19 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, t30, -t76, 0, 0, 0, 0, 0, -t25, -t26, -t25, -t13 + t78, t26, t1 * t41 + t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t40 + t20 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t75, t36, t11, -t12, 0, 0, t16, -t15, t36, t36 * t69 + t3, -t36 * t53 - t4, t50 * t36 - t2, -t47 * t15 - t50 * t16, t47 * t36 + t1, t1 * t47 - t2 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t61, 0, -t59 * t46, -t61 * t46, 0, 0, t41, -t40, 0, -t19, -t20, -t19, -t47 * t40 - t50 * t41, t20, -t19 * t50 + t20 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t59, 0, 0, 0, 0, 0, -t40, -t41, -t40, 0, t41, -t40 * t50 + t41 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t53, 0.2e1 * t50, 0, 0.2e1 * t47, t47 ^ 2 + t50 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t36, t3, -t4, -t2 + t32, -pkin(5) * t16 - t15 * qJ(6), 0.2e1 * t31 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, -t19, -t20, -t19, -pkin(5) * t41 - t40 * qJ(6), t20, -t19 * pkin(5) + t20 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t41, -t40, 0, t41, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t53, t63 + t69, 0, t62 + t53, t50 * pkin(5) + t47 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t63, 0, t62, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;

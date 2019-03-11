% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t67 = cos(qJ(2));
t53 = -t67 * pkin(2) - pkin(1);
t62 = sin(qJ(3));
t63 = sin(qJ(2));
t66 = cos(qJ(3));
t69 = t62 * t63 - t66 * t67;
t28 = t69 * pkin(3) + t53;
t94 = 0.2e1 * t28;
t93 = 0.2e1 * t53;
t92 = 0.2e1 * t67;
t91 = pkin(7) + pkin(8);
t60 = sin(qJ(5));
t90 = pkin(4) * t60;
t89 = t60 * pkin(5);
t61 = sin(qJ(4));
t88 = t61 * pkin(3);
t87 = t62 * pkin(2);
t64 = cos(qJ(5));
t86 = t64 * pkin(10);
t43 = t91 * t63;
t44 = t91 * t67;
t25 = t62 * t43 - t66 * t44;
t14 = -t69 * pkin(9) - t25;
t65 = cos(qJ(4));
t24 = -t66 * t43 - t62 * t44;
t38 = t62 * t67 + t66 * t63;
t68 = -t38 * pkin(9) + t24;
t10 = t61 * t14 - t65 * t68;
t85 = t10 * t64;
t56 = t66 * pkin(2);
t51 = t56 + pkin(3);
t76 = -t65 * t51 + t61 * t87;
t30 = -pkin(4) + t76;
t84 = t30 * t64;
t55 = t65 * pkin(3);
t50 = -t55 - pkin(4);
t83 = t50 * t64;
t20 = t65 * t38 - t61 * t69;
t82 = t60 * t20;
t81 = t60 * t64;
t11 = t65 * t14 + t61 * t68;
t80 = t64 * t11;
t79 = t64 * t20;
t74 = t65 * t87;
t33 = -t61 * t51 - t74;
t31 = pkin(10) - t33;
t78 = t64 * t31;
t49 = pkin(10) + t88;
t77 = t64 * t49;
t54 = t64 * qJ(6);
t19 = t61 * t38 + t65 * t69;
t75 = -0.2e1 * t20 * t19;
t52 = -t64 * pkin(5) - pkin(4);
t9 = t19 * pkin(4) - t20 * pkin(10) + t28;
t4 = -t60 * t11 + t64 * t9;
t1 = t19 * pkin(5) - t20 * t54 + t4;
t3 = t80 + (-qJ(6) * t20 + t9) * t60;
t73 = -t1 * t60 + t3 * t64;
t72 = -pkin(4) * t20 - pkin(10) * t19;
t71 = -t19 * t31 + t20 * t30;
t70 = -t19 * t49 + t20 * t50;
t59 = t64 ^ 2;
t58 = t60 ^ 2;
t57 = pkin(4) * t64;
t48 = 0.2e1 * t81;
t46 = t50 * t60;
t42 = t54 + t86;
t41 = (-qJ(6) - pkin(10)) * t60;
t40 = t52 - t55;
t37 = t42 * t64;
t36 = t54 + t77;
t35 = (-qJ(6) - t49) * t60;
t29 = t36 * t64;
t27 = t30 * t60;
t26 = t52 + t76;
t23 = t54 + t78;
t22 = (-qJ(6) - t31) * t60;
t21 = t23 * t64;
t18 = t20 ^ 2;
t17 = t64 * t19;
t16 = t60 * t19;
t15 = t60 * t79;
t12 = (-t58 + t59) * t20;
t8 = t10 * t60;
t6 = pkin(5) * t82 + t10;
t5 = t60 * t9 + t80;
t2 = [1, 0, 0, t63 ^ 2, t63 * t92, 0, 0, 0, pkin(1) * t92, -0.2e1 * pkin(1) * t63, t38 ^ 2, -0.2e1 * t38 * t69, 0, 0, 0, t69 * t93, t38 * t93, t18, t75, 0, 0, 0, t19 * t94, t20 * t94, t59 * t18, -0.2e1 * t18 * t81, 0.2e1 * t19 * t79, t60 * t75, t19 ^ 2, 0.2e1 * t10 * t82 + 0.2e1 * t4 * t19, 0.2e1 * t10 * t79 - 0.2e1 * t5 * t19, 0.2e1 * (-t1 * t64 - t3 * t60) * t20, t1 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t63, t67, 0, -t63 * pkin(7), -t67 * pkin(7), 0, 0, t38, -t69, 0, t24, t25, 0, 0, t20, -t19, 0, -t10, -t11, t15, t12, t16, t17, 0, t71 * t60 - t85, t71 * t64 + t8 (-t22 * t64 - t23 * t60) * t20 + t73, t1 * t22 + t3 * t23 + t6 * t26; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t87, 0, 0, 0, 0, 1, -0.2e1 * t76, 0.2e1 * t33, t58, t48, 0, 0, 0, -0.2e1 * t84, 0.2e1 * t27, -0.2e1 * t22 * t60 + 0.2e1 * t21, t22 ^ 2 + t23 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t69, 0, t24, t25, 0, 0, t20, -t19, 0, -t10, -t11, t15, t12, t16, t17, 0, t70 * t60 - t85, t70 * t64 + t8 (-t35 * t64 - t36 * t60) * t20 + t73, t1 * t35 + t3 * t36 + t6 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t87, 0, 0, 0, 0, 1, t55 - t76, -t74 + (-pkin(3) - t51) * t61, t58, t48, 0, 0, 0 (-t30 - t50) * t64, t46 + t27, t21 + t29 + (-t22 - t35) * t60, t22 * t35 + t23 * t36 + t26 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t88, t58, t48, 0, 0, 0, -0.2e1 * t83, 0.2e1 * t46, -0.2e1 * t35 * t60 + 0.2e1 * t29, t35 ^ 2 + t36 ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, -t10, -t11, t15, t12, t16, t17, 0, t72 * t60 - t85, t72 * t64 + t8 (-t41 * t64 - t42 * t60) * t20 + t73, t1 * t41 + t3 * t42 + t6 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t76, t33, t58, t48, 0, 0, 0, t57 - t84, t27 - t90, t21 + t37 + (-t22 - t41) * t60, t22 * t41 + t23 * t42 + t26 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t88, t58, t48, 0, 0, 0, t57 - t83, t46 - t90, t29 + t37 + (-t35 - t41) * t60, t35 * t41 + t36 * t42 + t40 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t58, t48, 0, 0, 0, 0.2e1 * t57, -0.2e1 * t90, -0.2e1 * t41 * t60 + 0.2e1 * t37, t41 ^ 2 + t42 ^ 2 + t52 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t82, t19, t4, -t5, -pkin(5) * t79, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t64, 0, -t60 * t31, -t78, -t89, t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t64, 0, -t60 * t49, -t77, -t89, t35 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t64, 0, -t60 * pkin(10), -t86, -t89, t41 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;

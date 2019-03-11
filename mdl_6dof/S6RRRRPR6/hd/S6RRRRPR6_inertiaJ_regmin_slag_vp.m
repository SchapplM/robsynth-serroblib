% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t81 = sin(qJ(4));
t82 = sin(qJ(3));
t85 = cos(qJ(4));
t86 = cos(qJ(3));
t58 = t81 * t82 - t85 * t86;
t59 = t81 * t86 + t82 * t85;
t78 = sin(pkin(11));
t79 = cos(pkin(11));
t39 = -t58 * t79 - t59 * t78;
t71 = -pkin(3) * t86 - pkin(2);
t49 = pkin(4) * t58 + t71;
t26 = -pkin(5) * t39 + t49;
t108 = 0.2e1 * t26;
t107 = 0.2e1 * t71;
t83 = sin(qJ(2));
t106 = -0.2e1 * t83;
t87 = cos(qJ(2));
t105 = -0.2e1 * t87;
t104 = 0.2e1 * t87;
t103 = pkin(8) + pkin(9);
t102 = pkin(2) * t86;
t101 = pkin(4) * t78;
t100 = pkin(7) * t82;
t99 = t81 * pkin(3);
t98 = t87 * pkin(3);
t97 = t82 * t83;
t96 = t82 * t86;
t95 = t82 * t87;
t62 = -pkin(2) * t87 - pkin(8) * t83 - pkin(1);
t92 = t86 * t87;
t90 = pkin(7) * t92;
t43 = t90 + (-pkin(9) * t83 + t62) * t82;
t94 = t85 * t43;
t93 = t86 * t83;
t57 = t86 * t62;
t41 = -pkin(9) * t93 + t57 + (-pkin(3) - t100) * t87;
t24 = t41 * t85 - t43 * t81;
t52 = t58 * t83;
t16 = -pkin(4) * t87 + qJ(5) * t52 + t24;
t25 = t41 * t81 + t94;
t51 = t59 * t83;
t20 = -qJ(5) * t51 + t25;
t9 = t16 * t78 + t20 * t79;
t65 = t103 * t82;
t66 = t103 * t86;
t44 = -t65 * t85 - t66 * t81;
t35 = -qJ(5) * t59 + t44;
t45 = -t65 * t81 + t66 * t85;
t36 = -qJ(5) * t58 + t45;
t19 = t35 * t78 + t36 * t79;
t72 = t83 * pkin(7);
t61 = pkin(3) * t97 + t72;
t91 = t83 * t104;
t29 = -t51 * t78 - t52 * t79;
t8 = t16 * t79 - t20 * t78;
t6 = -pkin(5) * t87 - pkin(10) * t29 + t8;
t28 = -t51 * t79 + t52 * t78;
t7 = pkin(10) * t28 + t9;
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t1 = t6 * t84 - t7 * t80;
t73 = t85 * pkin(3);
t70 = t73 + pkin(4);
t55 = t70 * t78 + t79 * t99;
t89 = -t55 - t101;
t18 = t35 * t79 - t36 * t78;
t53 = t70 * t79 - t78 * t99;
t42 = pkin(4) * t51 + t61;
t2 = t6 * t80 + t7 * t84;
t77 = t87 ^ 2;
t76 = t86 ^ 2;
t75 = t83 ^ 2;
t74 = t82 ^ 2;
t68 = pkin(4) * t79 + pkin(5);
t64 = t84 * t68;
t56 = t101 * t84 + t80 * t68;
t54 = -t101 * t80 + t64;
t50 = pkin(5) + t53;
t48 = t62 * t82 + t90;
t47 = -pkin(7) * t95 + t57;
t46 = t84 * t50;
t40 = -t58 * t78 + t59 * t79;
t34 = t80 * t50 + t84 * t55;
t33 = -t55 * t80 + t46;
t23 = t39 * t80 + t40 * t84;
t22 = -t39 * t84 + t40 * t80;
t21 = -pkin(5) * t28 + t42;
t13 = t28 * t80 + t29 * t84;
t12 = -t28 * t84 + t29 * t80;
t11 = pkin(10) * t39 + t19;
t10 = -pkin(10) * t40 + t18;
t4 = t10 * t80 + t11 * t84;
t3 = t10 * t84 - t11 * t80;
t5 = [1, 0, 0, t75, t91, 0, 0, 0, pkin(1) * t104, pkin(1) * t106, t76 * t75, -0.2e1 * t75 * t96, t92 * t106, t82 * t91, t77, 0.2e1 * t100 * t75 - 0.2e1 * t47 * t87, 0.2e1 * pkin(7) * t75 * t86 + 0.2e1 * t48 * t87, t52 ^ 2, 0.2e1 * t52 * t51, -t52 * t105, t51 * t104, t77, -0.2e1 * t24 * t87 + 0.2e1 * t51 * t61, 0.2e1 * t25 * t87 - 0.2e1 * t52 * t61, 0.2e1 * t28 * t9 - 0.2e1 * t29 * t8, t42 ^ 2 + t8 ^ 2 + t9 ^ 2, t13 ^ 2, -0.2e1 * t13 * t12, t13 * t105, t12 * t104, t77, -0.2e1 * t1 * t87 + 0.2e1 * t12 * t21, 0.2e1 * t13 * t21 + 0.2e1 * t2 * t87; 0, 0, 0, 0, 0, t83, t87, 0, -t72, -t87 * pkin(7), t82 * t93 (-t74 + t76) * t83, -t95, -t92, 0, -pkin(7) * t93 + (-pkin(2) * t83 + pkin(8) * t87) * t82, pkin(8) * t92 + (t100 - t102) * t83, -t52 * t59, -t51 * t59 + t52 * t58, -t59 * t87, t58 * t87, 0, -t44 * t87 + t51 * t71 + t58 * t61, t45 * t87 - t52 * t71 + t59 * t61, -t18 * t29 + t19 * t28 + t39 * t9 - t40 * t8, t18 * t8 + t19 * t9 + t42 * t49, t13 * t23, -t12 * t23 - t13 * t22, -t23 * t87, t22 * t87, 0, t12 * t26 + t21 * t22 - t3 * t87, t13 * t26 + t21 * t23 + t4 * t87; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t74, 0.2e1 * t96, 0, 0, 0, 0.2e1 * t102, -0.2e1 * pkin(2) * t82, t59 ^ 2, -0.2e1 * t59 * t58, 0, 0, 0, t58 * t107, t59 * t107, -0.2e1 * t18 * t40 + 0.2e1 * t19 * t39, t18 ^ 2 + t19 ^ 2 + t49 ^ 2, t23 ^ 2, -0.2e1 * t23 * t22, 0, 0, 0, t22 * t108, t23 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t97, -t87, t47, -t48, 0, 0, -t52, -t51, -t87, -t85 * t98 + t24, -t94 + (-t41 + t98) * t81, t28 * t55 - t29 * t53, t53 * t8 + t55 * t9, 0, 0, t13, -t12, -t87, -t33 * t87 + t1, t34 * t87 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t86, 0, -t82 * pkin(8), -t86 * pkin(8), 0, 0, t59, -t58, 0, t44, -t45, t39 * t55 - t40 * t53, t18 * t53 + t19 * t55, 0, 0, t23, -t22, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t73, -0.2e1 * t99, 0, t53 ^ 2 + t55 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t51, -t87, t24, -t25 (t28 * t78 - t29 * t79) * pkin(4) (t78 * t9 + t79 * t8) * pkin(4), 0, 0, t13, -t12, -t87, -t54 * t87 + t1, t56 * t87 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t58, 0, t44, -t45 (t39 * t78 - t40 * t79) * pkin(4) (t18 * t79 + t19 * t78) * pkin(4), 0, 0, t23, -t22, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t73, -t99, 0 (t53 * t79 + t55 * t78) * pkin(4), 0, 0, 0, 0, 1, t80 * t89 + t46 + t64, t89 * t84 + (-t50 - t68) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t78 ^ 2 + t79 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t54, -0.2e1 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, t22, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t87, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t33, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t54, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;

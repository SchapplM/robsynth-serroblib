% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t110 = -pkin(3) * t78 - qJ(4) * t75;
t76 = sin(qJ(2));
t62 = t78 * t76;
t72 = sin(pkin(10));
t73 = cos(pkin(10));
t96 = t75 * t76;
t34 = t62 * t72 - t73 * t96;
t39 = t72 * t75 + t73 * t78;
t35 = t39 * t76;
t74 = sin(qJ(6));
t77 = cos(qJ(6));
t11 = t34 * t77 + t35 * t74;
t109 = -0.2e1 * t11;
t51 = -pkin(2) + t110;
t38 = pkin(4) * t78 - t51;
t23 = pkin(5) * t39 + t38;
t108 = 0.2e1 * t23;
t107 = 0.2e1 * t38;
t106 = -0.2e1 * t51;
t105 = -0.2e1 * t76;
t79 = cos(qJ(2));
t104 = 0.2e1 * t79;
t103 = -pkin(3) - pkin(4);
t102 = pkin(2) * t75;
t101 = pkin(2) * t78;
t100 = pkin(3) * t75;
t99 = pkin(7) * t75;
t98 = pkin(7) * t78;
t97 = t75 * pkin(8);
t95 = t75 * t78;
t94 = t75 * t79;
t93 = t78 * t79;
t67 = t79 * pkin(3);
t52 = -pkin(2) * t79 - pkin(8) * t76 - pkin(1);
t92 = pkin(7) * t94 - t52 * t78;
t29 = t67 + t92;
t89 = t78 * qJ(5);
t16 = pkin(4) * t79 - t76 * t89 + t29;
t32 = pkin(7) * t93 + t52 * t75;
t88 = t79 * qJ(4);
t28 = -t88 + t32;
t20 = qJ(5) * t96 + t28;
t8 = t16 * t72 + t20 * t73;
t64 = t78 * pkin(8);
t53 = t64 - t89;
t86 = (pkin(8) - qJ(5)) * t75;
t26 = t53 * t73 + t72 * t86;
t68 = t75 ^ 2;
t70 = t78 ^ 2;
t91 = t68 + t70;
t90 = t78 * qJ(4);
t87 = t76 * t104;
t85 = -t16 * t73 + t20 * t72;
t24 = t53 * t72 - t73 * t86;
t48 = qJ(4) * t72 - t103 * t73;
t84 = -pkin(5) - t48;
t3 = pkin(5) * t79 - pkin(9) * t35 - t85;
t4 = -pkin(9) * t34 + t8;
t1 = t3 * t77 - t4 * t74;
t2 = t3 * t74 + t4 * t77;
t83 = t90 - t100;
t82 = t28 * t78 + t29 * t75;
t41 = -t72 * t78 + t73 * t75;
t81 = -pkin(9) * t41 - t24;
t55 = t76 * t90;
t27 = t55 + (t103 * t75 - pkin(7)) * t76;
t71 = t79 ^ 2;
t69 = t76 ^ 2;
t57 = pkin(8) * t94;
t50 = qJ(4) * t73 + t103 * t72;
t42 = t72 * t77 + t73 * t74;
t40 = t72 * t74 - t73 * t77;
t33 = -t55 + (pkin(7) + t100) * t76;
t22 = t50 * t77 + t74 * t84;
t21 = t50 * t74 - t77 * t84;
t19 = -t39 * t74 + t41 * t77;
t18 = t39 * t77 + t41 * t74;
t13 = -pkin(9) * t39 + t26;
t12 = -t34 * t74 + t35 * t77;
t9 = t34 * pkin(5) + t27;
t6 = t13 * t77 + t74 * t81;
t5 = t13 * t74 - t77 * t81;
t7 = [1, 0, 0, t69, t87, 0, 0, 0, pkin(1) * t104, pkin(1) * t105, t70 * t69, -0.2e1 * t69 * t95, t93 * t105, t75 * t87, t71, 0.2e1 * t69 * t99 + 0.2e1 * t79 * t92, 0.2e1 * t32 * t79 + 0.2e1 * t69 * t98, 0.2e1 * t29 * t79 + 0.2e1 * t33 * t96, 0.2e1 * (-t28 * t75 + t29 * t78) * t76, -0.2e1 * t28 * t79 - 0.2e1 * t33 * t62, t28 ^ 2 + t29 ^ 2 + t33 ^ 2, 0.2e1 * t27 * t34 - 0.2e1 * t79 * t85, 0.2e1 * t27 * t35 - 0.2e1 * t79 * t8, -0.2e1 * t34 * t8 + 0.2e1 * t35 * t85, t27 ^ 2 + t8 ^ 2 + t85 ^ 2, t12 ^ 2, t12 * t109, t12 * t104, t79 * t109, t71, 0.2e1 * t1 * t79 + 0.2e1 * t11 * t9, 0.2e1 * t12 * t9 - 0.2e1 * t2 * t79; 0, 0, 0, 0, 0, t76, t79, 0, -t76 * pkin(7), -t79 * pkin(7), t75 * t62 (-t68 + t70) * t76, -t94, -t93, 0, t57 + (-t98 - t102) * t76, pkin(8) * t93 + (t99 - t101) * t76, -t33 * t78 + t51 * t96 + t57, t82, -t33 * t75 + (-pkin(8) * t79 - t51 * t76) * t78, pkin(8) * t82 + t33 * t51, -t24 * t79 + t27 * t39 + t34 * t38, -t26 * t79 + t27 * t41 + t35 * t38, t24 * t35 - t26 * t34 - t39 * t8 + t41 * t85, t24 * t85 + t26 * t8 + t27 * t38, t12 * t19, -t11 * t19 - t12 * t18, t19 * t79, -t18 * t79, 0, t11 * t23 + t18 * t9 - t5 * t79, t12 * t23 + t19 * t9 - t6 * t79; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t68, 0.2e1 * t95, 0, 0, 0, 0.2e1 * t101, -0.2e1 * t102, t78 * t106, 0.2e1 * t91 * pkin(8), t75 * t106, pkin(8) ^ 2 * t91 + t51 ^ 2, t39 * t107, t41 * t107, 0.2e1 * t24 * t41 - 0.2e1 * t26 * t39, t24 ^ 2 + t26 ^ 2 + t38 ^ 2, t19 ^ 2, -0.2e1 * t19 * t18, 0, 0, 0, t18 * t108, t19 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t96, -t79, -t92, -t32, -0.2e1 * t67 - t92, t110 * t76, -0.2e1 * t88 + t32, -pkin(3) * t29 + qJ(4) * t28, -t48 * t79 + t85, -t50 * t79 + t8, -t34 * t50 + t35 * t48, t48 * t85 + t50 * t8, 0, 0, -t12, t11, -t79, -t21 * t79 - t1, -t22 * t79 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t78, 0, -t97, -t64, -t97, t83, t64, t83 * pkin(8), t24, t26, -t39 * t50 + t41 * t48, t24 * t48 + t26 * t50, 0, 0, -t19, t18, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0.2e1 * t48, 0.2e1 * t50, 0, t48 ^ 2 + t50 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t21, 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t62, 0, t29, t73 * t79, -t72 * t79, -t34 * t72 - t35 * t73, t72 * t8 - t73 * t85, 0, 0, 0, 0, 0, -t40 * t79, -t42 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t97, 0, 0, -t39 * t72 - t41 * t73, -t24 * t73 + t26 * t72, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), -t73, t72, 0, -t48 * t73 + t50 * t72, 0, 0, 0, 0, 0, t40, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t72 ^ 2 + t73 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, 0, t27, 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, 0, t38, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t79, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t7;

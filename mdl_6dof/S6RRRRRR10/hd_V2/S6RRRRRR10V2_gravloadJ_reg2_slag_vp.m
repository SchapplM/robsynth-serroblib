% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR10V2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t53 = qJ(2) + qJ(3);
t50 = sin(t53);
t58 = sin(qJ(1));
t100 = t50 * t58;
t51 = cos(t53);
t56 = sin(qJ(4));
t63 = cos(qJ(1));
t88 = t63 * t56;
t61 = cos(qJ(4));
t91 = t58 * t61;
t32 = t51 * t91 - t88;
t55 = sin(qJ(5));
t60 = cos(qJ(5));
t15 = t55 * t100 + t32 * t60;
t87 = t63 * t61;
t92 = t58 * t56;
t31 = t51 * t92 + t87;
t54 = sin(qJ(6));
t59 = cos(qJ(6));
t110 = t15 * t54 - t31 * t59;
t109 = t15 * t59 + t31 * t54;
t105 = g(3) * t56;
t33 = t51 * t88 - t91;
t66 = g(1) * t33 + g(2) * t31 + t50 * t105;
t10 = t66 * t55;
t86 = t51 * pkin(3) + t50 * pkin(5);
t99 = t50 * t60;
t74 = -t32 * t55 + t58 * t99;
t108 = g(1) * t74;
t107 = g(1) * t58;
t106 = g(3) * t50;
t94 = t55 * t61;
t73 = t50 * t94 + t51 * t60;
t22 = t73 * t58;
t104 = t22 * pkin(6);
t24 = t73 * t63;
t103 = t24 * pkin(6);
t98 = t50 * t63;
t97 = t51 * t63;
t96 = t54 * t56;
t95 = t54 * t60;
t93 = t56 * t59;
t90 = t59 * t60;
t89 = t60 * t61;
t85 = t50 * t96;
t84 = t50 * t93;
t83 = t50 * t88;
t62 = cos(qJ(2));
t52 = t62 * pkin(2);
t49 = t52 + pkin(1);
t82 = pkin(3) * t97 + pkin(5) * t98 + t63 * t49;
t38 = t58 * t51 * pkin(5);
t81 = -pkin(3) * t100 + t38;
t40 = pkin(5) * t97;
t80 = -pkin(3) * t98 + t40;
t28 = t51 * t94 - t99;
t79 = t28 * pkin(6) + t86;
t57 = sin(qJ(2));
t78 = -pkin(2) * t57 - pkin(3) * t50;
t34 = t51 * t87 + t92;
t18 = t34 * t55 - t60 * t98;
t77 = g(2) * t18 + t108;
t76 = g(1) * t31 - g(2) * t33;
t36 = g(1) * t63 + g(2) * t58;
t75 = -g(2) * t63 + t107;
t27 = t50 * t89 - t51 * t55;
t72 = t36 * t50;
t71 = g(1) * t18 - g(2) * t74 + g(3) * t73;
t19 = t34 * t60 + t55 * t98;
t70 = g(1) * t19 + g(2) * t15 + g(3) * t27;
t4 = g(1) * t24 + g(2) * t22 - g(3) * t28;
t69 = t78 * t58 + t38;
t68 = t78 * t63 + t40;
t67 = (-t49 - t86) * t107;
t65 = g(1) * t34 + g(2) * t32 + t61 * t106;
t20 = -g(3) * t51 + t72;
t64 = -g(3) * t62 + t36 * t57;
t30 = t75 * t50;
t29 = t50 * t55 + t51 * t89;
t25 = t27 * t63;
t23 = t27 * t58;
t21 = t36 * t51 + t106;
t13 = t20 * t61;
t12 = -t51 * t105 + t56 * t72;
t11 = -g(1) * t80 - g(2) * t81 - g(3) * t86;
t9 = t19 * t59 + t33 * t54;
t8 = -t19 * t54 + t33 * t59;
t7 = -g(2) * t82 - t67;
t6 = -g(1) * t68 - g(2) * t69 - g(3) * (t52 + t86);
t5 = g(1) * t25 + g(2) * t23 - g(3) * t29;
t2 = -g(1) * (-t25 * t59 - t54 * t83) - g(2) * (-t23 * t59 - t58 * t85) - g(3) * (t29 * t59 + t51 * t96);
t1 = -g(1) * (t25 * t54 - t59 * t83) - g(2) * (t23 * t54 - t58 * t84) - g(3) * (-t29 * t54 + t51 * t93);
t3 = [0, 0, 0, 0, 0, 0, t75, t36, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t62, -t75 * t57, -t36, t75 * pkin(1), 0, 0, 0, 0, 0, 0, t75 * t51, -t30, -t36, t75 * t49, 0, 0, 0, 0, 0, 0, g(1) * t32 - g(2) * t34, -t76, t30, t7, 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t19, t77, t76, t7, 0, 0, 0, 0, 0, 0, g(1) * t109 - g(2) * t9, -g(1) * t110 - g(2) * t8, -t77, -pkin(6) * t108 - g(2) * (t18 * pkin(6) + t82) - t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, g(3) * t57 + t36 * t62, 0, 0, 0, 0, 0, 0, 0, 0, t20, t21, 0, t64 * pkin(2), 0, 0, 0, 0, 0, 0, t13, -t12, -t21, t6, 0, 0, 0, 0, 0, 0, t5, -t4, t12, t6, 0, 0, 0, 0, 0, 0, t2, t1, t4, -g(1) * (t68 - t103) - g(2) * (t69 - t104) - g(3) * (t52 + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t21, t11, 0, 0, 0, 0, 0, 0, t5, -t4, t12, t11, 0, 0, 0, 0, 0, 0, t2, t1, t4, -g(1) * (t80 - t103) - g(2) * (t81 - t104) - g(3) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t60, -t10, -t65, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t33 * t90 + t34 * t54) - g(2) * (-t31 * t90 + t32 * t54) - (t54 * t61 - t56 * t90) * t106, -g(1) * (t33 * t95 + t34 * t59) - g(2) * (t31 * t95 + t32 * t59) - (t56 * t95 + t59 * t61) * t106, t10, pkin(6) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t70, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t59, -t71 * t54, -t70, -t70 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t110 - g(3) * (-t27 * t54 + t84) g(1) * t9 + g(2) * t109 - g(3) * (-t27 * t59 - t85) 0, 0;];
taug_reg  = t3;

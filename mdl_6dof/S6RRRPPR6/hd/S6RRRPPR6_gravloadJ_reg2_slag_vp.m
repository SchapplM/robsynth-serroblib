% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t63 = sin(pkin(6));
t72 = cos(qJ(1));
t105 = t63 * t72;
t67 = sin(qJ(2));
t68 = sin(qJ(1));
t71 = cos(qJ(2));
t98 = cos(pkin(6));
t90 = t72 * t98;
t39 = t67 * t90 + t68 * t71;
t62 = qJ(3) + pkin(11);
t59 = sin(t62);
t60 = cos(t62);
t15 = t60 * t105 + t39 * t59;
t38 = t68 * t67 - t71 * t90;
t65 = sin(qJ(6));
t69 = cos(qJ(6));
t120 = t15 * t65 + t38 * t69;
t119 = t15 * t69 - t38 * t65;
t106 = t63 * t71;
t99 = qJ(5) * t59;
t118 = (pkin(4) * t60 + t99) * t106;
t117 = g(3) * t63;
t116 = t38 * t60;
t91 = t68 * t98;
t40 = t72 * t67 + t71 * t91;
t113 = t40 * t60;
t41 = -t67 * t91 + t72 * t71;
t66 = sin(qJ(3));
t112 = t41 * t66;
t111 = t59 * t65;
t110 = t59 * t69;
t109 = t63 * t67;
t108 = t63 * t68;
t70 = cos(qJ(3));
t107 = t63 * t70;
t104 = t65 * t71;
t103 = t69 * t71;
t58 = t70 * pkin(3) + pkin(2);
t64 = -qJ(4) - pkin(9);
t102 = -t38 * t58 - t39 * t64;
t101 = -t40 * t58 - t41 * t64;
t100 = t72 * pkin(1) + pkin(8) * t108;
t97 = t66 * t109;
t96 = t66 * t108;
t95 = t68 * t107;
t52 = t66 * t105;
t94 = t70 * t105;
t93 = -t68 * pkin(1) + pkin(8) * t105;
t16 = -t59 * t105 + t39 * t60;
t92 = t39 * t70 - t52;
t89 = t98 * t70;
t88 = -pkin(4) * t116 - t38 * t99 + t102;
t87 = -pkin(4) * t113 - t40 * t99 + t101;
t44 = t58 * t106;
t86 = -t64 * t109 + t44;
t85 = pkin(3) * t96 - t40 * t64 + t41 * t58 + t100;
t19 = -t60 * t108 + t41 * t59;
t84 = -g(1) * t15 + g(2) * t19;
t20 = t59 * t108 + t41 * t60;
t83 = -g(1) * t16 + g(2) * t20;
t10 = g(1) * t38 - g(2) * t40;
t82 = g(1) * t72 + g(2) * t68;
t81 = t39 * t66 + t94;
t80 = pkin(3) * t52 + t38 * t64 - t39 * t58 + t93;
t32 = t59 * t109 - t98 * t60;
t2 = g(1) * t19 + g(2) * t15 + g(3) * t32;
t33 = t60 * t109 + t98 * t59;
t79 = g(1) * t20 + g(2) * t16 + g(3) * t33;
t49 = pkin(3) * t95;
t78 = -pkin(3) * t112 - t19 * pkin(4) + t20 * qJ(5) + t49;
t8 = -g(1) * t40 - g(2) * t38 + g(3) * t106;
t77 = g(1) * t41 + g(2) * t39 + g(3) * t109;
t76 = t20 * pkin(4) + t19 * qJ(5) + t85;
t57 = pkin(3) * t89;
t75 = -pkin(3) * t97 - t32 * pkin(4) + t33 * qJ(5) + t57;
t74 = -pkin(4) * t16 - qJ(5) * t15 + t80;
t73 = -t81 * pkin(3) - t15 * pkin(4) + t16 * qJ(5);
t22 = t41 * t70 + t96;
t21 = t95 - t112;
t7 = t19 * t65 + t40 * t69;
t6 = t19 * t69 - t40 * t65;
t5 = t8 * t60;
t4 = t8 * t59;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t72, t82, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t39 - g(2) * t41, -t10, -t82 * t63, -g(1) * t93 - g(2) * t100, 0, 0, 0, 0, 0, 0, g(1) * t92 - g(2) * t22, -g(1) * t81 - g(2) * t21, t10, -g(1) * (-t39 * pkin(2) - t38 * pkin(9) + t93) - g(2) * (t41 * pkin(2) + t40 * pkin(9) + t100) 0, 0, 0, 0, 0, 0, -t83, t84, t10, -g(1) * t80 - g(2) * t85, 0, 0, 0, 0, 0, 0, t10, t83, -t84, -g(1) * t74 - g(2) * t76, 0, 0, 0, 0, 0, 0, g(1) * t120 - g(2) * t7, g(1) * t119 - g(2) * t6, -t83, -g(1) * (-t38 * pkin(5) - pkin(10) * t16 + t74) - g(2) * (t40 * pkin(5) + t20 * pkin(10) + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t77, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t70, t8 * t66, -t77, -g(1) * (-t40 * pkin(2) + t41 * pkin(9)) - g(2) * (-t38 * pkin(2) + t39 * pkin(9)) - (pkin(2) * t71 + pkin(9) * t67) * t117, 0, 0, 0, 0, 0, 0, -t5, t4, -t77, -g(1) * t101 - g(2) * t102 - g(3) * t86, 0, 0, 0, 0, 0, 0, -t77, t5, -t4, -g(1) * t87 - g(2) * t88 - g(3) * (t86 + t118) 0, 0, 0, 0, 0, 0, -g(1) * (-t40 * t111 + t41 * t69) - g(2) * (-t38 * t111 + t39 * t69) - (t59 * t104 + t67 * t69) * t117, -g(1) * (-t40 * t110 - t41 * t65) - g(2) * (-t38 * t110 - t39 * t65) - (t59 * t103 - t65 * t67) * t117, -t5, -g(1) * (t41 * pkin(5) - pkin(10) * t113 + t87) - g(2) * (t39 * pkin(5) - pkin(10) * t116 + t88) - g(3) * (t44 + t118) - (pkin(10) * t60 * t71 + (pkin(5) - t64) * t67) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t21 + g(2) * t81 - g(3) * (t89 - t97) g(1) * t22 + g(2) * t92 - g(3) * (-t67 * t107 - t98 * t66) 0, 0, 0, 0, 0, 0, 0, 0, t2, t79, 0, -g(1) * t49 - g(3) * t57 + (g(2) * t94 + t77 * t66) * pkin(3), 0, 0, 0, 0, 0, 0, 0, -t2, -t79, -g(1) * t78 - g(2) * t73 - g(3) * t75, 0, 0, 0, 0, 0, 0, -t79 * t65, -t79 * t69, t2, -g(1) * (-t19 * pkin(10) + t78) - g(2) * (-t15 * pkin(10) + t73) - g(3) * (-t32 * pkin(10) + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t119 - g(3) * (t63 * t104 + t32 * t69) g(1) * t7 + g(2) * t120 - g(3) * (t63 * t103 - t32 * t65) 0, 0;];
taug_reg  = t1;

% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t63 = sin(qJ(2));
t64 = sin(qJ(1));
t66 = cos(qJ(2));
t104 = cos(qJ(1));
t88 = cos(pkin(6));
t73 = t88 * t104;
t37 = t63 * t73 + t64 * t66;
t57 = pkin(11) + qJ(4);
t54 = sin(t57);
t55 = cos(t57);
t59 = sin(pkin(6));
t86 = t59 * t104;
t15 = t37 * t54 + t55 * t86;
t36 = t64 * t63 - t66 * t73;
t62 = sin(qJ(6));
t65 = cos(qJ(6));
t108 = t15 * t62 + t36 * t65;
t107 = t15 * t65 - t36 * t62;
t89 = qJ(5) * t54;
t95 = t59 * t66;
t106 = (pkin(4) * t55 + t89) * t95;
t105 = g(3) * t59;
t103 = t36 * t55;
t84 = t64 * t88;
t38 = t104 * t63 + t66 * t84;
t100 = t38 * t55;
t99 = t54 * t62;
t98 = t54 * t65;
t97 = t59 * t63;
t96 = t59 * t64;
t94 = t62 * t66;
t93 = t65 * t66;
t60 = cos(pkin(11));
t53 = t60 * pkin(3) + pkin(2);
t61 = -pkin(9) - qJ(3);
t92 = -t36 * t53 - t37 * t61;
t39 = t104 * t66 - t63 * t84;
t91 = -t38 * t53 - t39 * t61;
t90 = t104 * pkin(1) + pkin(8) * t96;
t58 = sin(pkin(11));
t87 = t58 * t96;
t85 = -t64 * pkin(1) + pkin(8) * t86;
t16 = t37 * t55 - t54 * t86;
t83 = -t15 * pkin(4) + t16 * qJ(5);
t19 = t39 * t54 - t55 * t96;
t20 = t39 * t55 + t54 * t96;
t82 = -t19 * pkin(4) + t20 * qJ(5);
t30 = t54 * t97 - t88 * t55;
t31 = t88 * t54 + t55 * t97;
t81 = -t30 * pkin(4) + t31 * qJ(5);
t80 = -pkin(4) * t103 - t36 * t89 + t92;
t79 = -pkin(4) * t100 - t38 * t89 + t91;
t78 = t58 * t86;
t41 = t53 * t95;
t77 = -t61 * t97 + t41;
t76 = pkin(3) * t87 - t38 * t61 + t39 * t53 + t90;
t75 = -g(1) * t15 + g(2) * t19;
t74 = -g(1) * t16 + g(2) * t20;
t10 = g(1) * t36 - g(2) * t38;
t72 = g(1) * t104 + g(2) * t64;
t71 = pkin(3) * t78 + t36 * t61 - t37 * t53 + t85;
t2 = g(1) * t19 + g(2) * t15 + g(3) * t30;
t70 = g(1) * t20 + g(2) * t16 + g(3) * t31;
t8 = -g(1) * t38 - g(2) * t36 + g(3) * t95;
t69 = g(1) * t39 + g(2) * t37 + g(3) * t97;
t68 = t20 * pkin(4) + t19 * qJ(5) + t76;
t67 = -pkin(4) * t16 - qJ(5) * t15 + t71;
t7 = t19 * t62 + t38 * t65;
t6 = t19 * t65 - t38 * t62;
t5 = t8 * t55;
t4 = t8 * t54;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t64 - g(2) * t104, t72, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t37 - g(2) * t39, -t10, -t72 * t59, -g(1) * t85 - g(2) * t90, 0, 0, 0, 0, 0, 0, -g(1) * (-t37 * t60 + t78) - g(2) * (t39 * t60 + t87) -g(1) * (t37 * t58 + t60 * t86) - g(2) * (-t39 * t58 + t60 * t96) t10, -g(1) * (-t37 * pkin(2) - t36 * qJ(3) + t85) - g(2) * (t39 * pkin(2) + t38 * qJ(3) + t90) 0, 0, 0, 0, 0, 0, -t74, t75, t10, -g(1) * t71 - g(2) * t76, 0, 0, 0, 0, 0, 0, t10, t74, -t75, -g(1) * t67 - g(2) * t68, 0, 0, 0, 0, 0, 0, g(1) * t108 - g(2) * t7, g(1) * t107 - g(2) * t6, -t74, -g(1) * (-t36 * pkin(5) - pkin(10) * t16 + t67) - g(2) * (t38 * pkin(5) + t20 * pkin(10) + t68); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t69, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t60, t8 * t58, -t69, -g(1) * (-t38 * pkin(2) + t39 * qJ(3)) - g(2) * (-t36 * pkin(2) + t37 * qJ(3)) - (pkin(2) * t66 + qJ(3) * t63) * t105, 0, 0, 0, 0, 0, 0, -t5, t4, -t69, -g(1) * t91 - g(2) * t92 - g(3) * t77, 0, 0, 0, 0, 0, 0, -t69, t5, -t4, -g(1) * t79 - g(2) * t80 - g(3) * (t77 + t106) 0, 0, 0, 0, 0, 0, -g(1) * (-t38 * t99 + t39 * t65) - g(2) * (-t36 * t99 + t37 * t65) - (t54 * t94 + t63 * t65) * t105, -g(1) * (-t38 * t98 - t39 * t62) - g(2) * (-t36 * t98 - t37 * t62) - (t54 * t93 - t62 * t63) * t105, -t5, -g(1) * (t39 * pkin(5) - pkin(10) * t100 + t79) - g(2) * (t37 * pkin(5) - pkin(10) * t103 + t80) - g(3) * (t41 + t106) - (pkin(10) * t55 * t66 + (pkin(5) - t61) * t63) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t70, -g(1) * t82 - g(2) * t83 - g(3) * t81, 0, 0, 0, 0, 0, 0, -t70 * t62, -t70 * t65, t2, -g(1) * (-t19 * pkin(10) + t82) - g(2) * (-t15 * pkin(10) + t83) - g(3) * (-t30 * pkin(10) + t81); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t107 - g(3) * (t30 * t65 + t59 * t94) g(1) * t7 + g(2) * t108 - g(3) * (-t30 * t62 + t59 * t93) 0, 0;];
taug_reg  = t1;

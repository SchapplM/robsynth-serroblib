% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t63 = cos(pkin(11));
t57 = pkin(4) * t63 + pkin(3);
t67 = cos(qJ(3));
t64 = -pkin(9) - qJ(4);
t65 = sin(qJ(3));
t95 = t64 * t65;
t106 = -t57 * t67 + t95;
t62 = sin(pkin(6));
t105 = g(3) * t62;
t66 = sin(qJ(2));
t68 = cos(qJ(2));
t88 = cos(pkin(10));
t89 = cos(pkin(6));
t79 = t89 * t88;
t87 = sin(pkin(10));
t41 = t66 * t79 + t68 * t87;
t61 = sin(pkin(11));
t104 = t41 * t61;
t78 = t89 * t87;
t43 = -t66 * t78 + t68 * t88;
t103 = t43 * t61;
t60 = pkin(11) + qJ(5);
t58 = sin(t60);
t101 = t58 * t67;
t59 = cos(t60);
t100 = t59 * t67;
t99 = t61 * t67;
t98 = t62 * t66;
t97 = t62 * t68;
t96 = t63 * t67;
t94 = t67 * t68;
t83 = t62 * t88;
t22 = t41 * t65 + t67 * t83;
t23 = t41 * t67 - t65 * t83;
t93 = -t22 * t57 - t23 * t64;
t82 = t62 * t87;
t24 = t43 * t65 - t67 * t82;
t25 = t43 * t67 + t65 * t82;
t92 = -t24 * t57 - t25 * t64;
t44 = t65 * t98 - t67 * t89;
t45 = t65 * t89 + t67 * t98;
t91 = -t44 * t57 - t45 * t64;
t90 = pkin(2) * t97 + pkin(8) * t98;
t86 = t58 * t97;
t40 = t66 * t87 - t68 * t79;
t85 = -pkin(2) * t40 + t41 * pkin(8);
t42 = t66 * t88 + t68 * t78;
t84 = -pkin(2) * t42 + t43 * pkin(8);
t81 = pkin(3) * t67 + qJ(4) * t65;
t80 = -pkin(5) * t59 - qJ(6) * t58;
t18 = t45 * t58 + t59 * t97;
t7 = t23 * t58 - t40 * t59;
t9 = t25 * t58 - t42 * t59;
t1 = g(1) * t9 + g(2) * t7 + g(3) * t18;
t10 = t25 * t59 + t42 * t58;
t19 = t45 * t59 - t86;
t8 = t23 * t59 + t40 * t58;
t77 = g(1) * t10 + g(2) * t8 + g(3) * t19;
t12 = -t101 * t40 - t41 * t59;
t14 = -t101 * t42 - t43 * t59;
t26 = -t59 * t98 + t67 * t86;
t76 = g(1) * t14 + g(2) * t12 + g(3) * t26;
t75 = g(1) * t24 + g(2) * t22 + g(3) * t44;
t74 = g(1) * t25 + g(2) * t23 + g(3) * t45;
t73 = -g(1) * t42 - g(2) * t40 + g(3) * t97;
t72 = g(1) * t43 + g(2) * t41 + g(3) * t98;
t71 = pkin(4) * t61 * t98 + t57 * t62 * t94 - t95 * t97 + t90;
t70 = pkin(4) * t104 + t106 * t40 + t85;
t69 = pkin(4) * t103 + t106 * t42 + t84;
t27 = (t58 * t66 + t59 * t94) * t62;
t15 = -t100 * t42 + t43 * t58;
t13 = -t100 * t40 + t41 * t58;
t11 = t73 * t65;
t4 = t75 * t59;
t3 = t75 * t58;
t2 = -g(1) * t15 - g(2) * t13 - g(3) * t27;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t72, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * t67, t11, -t72, -g(1) * t84 - g(2) * t85 - g(3) * t90, 0, 0, 0, 0, 0, 0, -g(1) * (-t42 * t96 + t103) - g(2) * (-t40 * t96 + t104) - (t61 * t66 + t63 * t94) * t105, -g(1) * (t42 * t99 + t43 * t63) - g(2) * (t40 * t99 + t41 * t63) - (-t61 * t94 + t63 * t66) * t105, -t11, -g(1) * (-t42 * t81 + t84) - g(2) * (-t40 * t81 + t85) - g(3) * (t81 * t97 + t90) 0, 0, 0, 0, 0, 0, t2, t76, -t11, -g(1) * t69 - g(2) * t70 - g(3) * t71, 0, 0, 0, 0, 0, 0, t2, -t11, -t76, -g(1) * (pkin(5) * t15 + qJ(6) * t14 + t69) - g(2) * (pkin(5) * t13 + qJ(6) * t12 + t70) - g(3) * (pkin(5) * t27 + qJ(6) * t26 + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t63, -t75 * t61, -t74, -g(1) * (-pkin(3) * t24 + qJ(4) * t25) - g(2) * (-pkin(3) * t22 + qJ(4) * t23) - g(3) * (-pkin(3) * t44 + qJ(4) * t45) 0, 0, 0, 0, 0, 0, t4, -t3, -t74, -g(1) * t92 - g(2) * t93 - g(3) * t91, 0, 0, 0, 0, 0, 0, t4, -t74, t3, -g(1) * (t24 * t80 + t92) - g(2) * (t22 * t80 + t93) - g(3) * (t44 * t80 + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t77, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t77, -g(1) * (-pkin(5) * t9 + qJ(6) * t10) - g(2) * (-pkin(5) * t7 + qJ(6) * t8) - g(3) * (-pkin(5) * t18 + qJ(6) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;

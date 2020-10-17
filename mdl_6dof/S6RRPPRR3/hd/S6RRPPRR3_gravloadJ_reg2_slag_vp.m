% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:03:57
% EndTime: 2019-05-06 10:04:00
% DurationCPUTime: 0.80s
% Computational Cost: add. (738->136), mult. (1583->216), div. (0->0), fcn. (1998->14), ass. (0->80)
t61 = sin(pkin(11));
t67 = sin(qJ(2));
t70 = cos(qJ(2));
t96 = cos(pkin(11));
t44 = -t70 * t61 - t67 * t96;
t68 = sin(qJ(1));
t71 = cos(qJ(1));
t64 = cos(pkin(6));
t78 = -t67 * t61 + t70 * t96;
t73 = t78 * t64;
t22 = t68 * t44 + t71 * t73;
t66 = sin(qJ(6));
t69 = cos(qJ(6));
t62 = sin(pkin(6));
t102 = t62 * t71;
t36 = t44 * t64;
t23 = -t71 * t36 + t68 * t78;
t59 = pkin(12) + qJ(5);
t57 = sin(t59);
t58 = cos(t59);
t8 = -t57 * t102 + t23 * t58;
t112 = t22 * t69 + t8 * t66;
t111 = -t22 * t66 + t8 * t69;
t103 = t62 * t70;
t98 = t71 * t67;
t99 = t68 * t70;
t40 = -t64 * t99 - t98;
t110 = -g(1) * t40 - g(3) * t103;
t106 = t58 * t66;
t105 = t58 * t69;
t104 = t62 * t68;
t100 = t68 * t67;
t97 = t71 * t70;
t60 = sin(pkin(12));
t94 = t60 * t104;
t93 = t60 * t102;
t92 = t64 * t97;
t34 = t78 * t62;
t35 = t44 * t62;
t53 = pkin(2) * t103;
t63 = cos(pkin(12));
t55 = t63 * pkin(4) + pkin(3);
t65 = -pkin(9) - qJ(4);
t91 = t34 * t55 + t35 * t65 + t53;
t90 = -t58 * t102 - t23 * t57;
t37 = t64 * t67 * pkin(2) + (-pkin(8) - qJ(3)) * t62;
t56 = t70 * pkin(2) + pkin(1);
t89 = -t68 * t37 + t71 * t56;
t50 = pkin(2) * t92;
t87 = -pkin(2) * t100 + t50;
t24 = -t68 * t36 - t71 * t78;
t11 = -t58 * t104 - t24 * t57;
t86 = g(1) * t90 + g(2) * t11;
t85 = pkin(5) * t58 + pkin(10) * t57;
t25 = t71 * t44 - t68 * t73;
t84 = g(1) * t22 - g(2) * t25;
t83 = g(1) * t71 + g(2) * t68;
t82 = g(1) * t68 - g(2) * t71;
t81 = -t71 * t37 - t68 * t56;
t80 = t22 * t55 - t23 * t65 + t87;
t79 = pkin(4) * t94 - t24 * t55 + t25 * t65 + t89;
t27 = t35 * t57 + t64 * t58;
t77 = g(1) * t11 - g(2) * t90 - g(3) * t27;
t12 = t57 * t104 - t24 * t58;
t28 = -t35 * t58 + t64 * t57;
t76 = g(1) * t12 + g(2) * t8 + g(3) * t28;
t4 = g(1) * t24 - g(2) * t23 + g(3) * t35;
t5 = g(1) * t25 + g(2) * t22 + g(3) * t34;
t75 = t40 * pkin(2);
t74 = pkin(4) * t93 - t22 * t65 - t23 * t55 + t81;
t72 = t24 * t65 + t25 * t55 + t75;
t42 = t83 * t62;
t41 = -t64 * t100 + t97;
t39 = -t64 * t98 - t99;
t38 = -t92 + t100;
t33 = -g(3) * t64 - t82 * t62;
t3 = t12 * t69 - t25 * t66;
t2 = -t12 * t66 - t25 * t69;
t1 = t5 * t57;
t6 = [0, 0, 0, 0, 0, 0, t82, t83, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t39 - g(2) * t41, -g(1) * t38 - g(2) * t40, -t42, -g(1) * (-t68 * pkin(1) + pkin(8) * t102) - g(2) * (t71 * pkin(1) + pkin(8) * t104) 0, 0, 0, 0, 0, 0, g(1) * t23 + g(2) * t24, t84, -t42, -g(1) * t81 - g(2) * t89, 0, 0, 0, 0, 0, 0, -g(1) * (-t23 * t63 + t93) - g(2) * (-t24 * t63 + t94) -g(1) * (t63 * t102 + t23 * t60) - g(2) * (t63 * t104 + t24 * t60) -t84, -g(1) * (-t23 * pkin(3) + t22 * qJ(4) + t81) - g(2) * (-pkin(3) * t24 - t25 * qJ(4) + t89) 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t12, t86, -t84, -g(1) * t74 - g(2) * t79, 0, 0, 0, 0, 0, 0, g(1) * t111 - g(2) * t3, -g(1) * t112 - g(2) * t2, -t86, -g(1) * (-pkin(5) * t8 + pkin(10) * t90 + t74) - g(2) * (t12 * pkin(5) + t11 * pkin(10) + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t38 + t110, g(3) * t62 * t67 + g(1) * t41 - g(2) * t39, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t4, 0, -g(2) * t50 + (g(2) * t100 + t110) * pkin(2), 0, 0, 0, 0, 0, 0, -t5 * t63, t5 * t60, t4, -g(1) * (t25 * pkin(3) - t24 * qJ(4) + t75) - g(2) * (t22 * pkin(3) + qJ(4) * t23 + t87) - g(3) * (t34 * pkin(3) - t35 * qJ(4) + t53) 0, 0, 0, 0, 0, 0, -t5 * t58, t1, t4, -g(1) * t72 - g(2) * t80 - g(3) * t91, 0, 0, 0, 0, 0, 0, -g(1) * (t25 * t105 - t24 * t66) - g(2) * (t22 * t105 + t23 * t66) - g(3) * (t34 * t105 - t35 * t66) -g(1) * (-t25 * t106 - t24 * t69) - g(2) * (-t22 * t106 + t23 * t69) - g(3) * (-t34 * t106 - t35 * t69) -t1, -g(1) * (t85 * t25 + t72) - g(2) * (t85 * t22 + t80) - g(3) * (t85 * t34 + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t69, -t77 * t66, -t76, -g(1) * (-t11 * pkin(5) + t12 * pkin(10)) - g(2) * (pkin(5) * t90 + t8 * pkin(10)) - g(3) * (t27 * pkin(5) + t28 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t112 - g(3) * (-t28 * t66 - t34 * t69) g(1) * t3 + g(2) * t111 - g(3) * (-t28 * t69 + t34 * t66) 0, 0;];
taug_reg  = t6;

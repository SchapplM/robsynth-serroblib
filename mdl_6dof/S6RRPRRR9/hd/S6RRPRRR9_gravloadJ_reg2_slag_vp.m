% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:07:04
% EndTime: 2019-05-06 23:07:06
% DurationCPUTime: 0.83s
% Computational Cost: add. (813->150), mult. (1137->235), div. (0->0), fcn. (1354->14), ass. (0->81)
t67 = sin(qJ(2));
t68 = sin(qJ(1));
t70 = cos(qJ(2));
t71 = cos(qJ(1));
t93 = cos(pkin(6));
t80 = t71 * t93;
t37 = t67 * t80 + t68 * t70;
t61 = pkin(12) + qJ(4);
t57 = qJ(5) + t61;
t52 = sin(t57);
t53 = cos(t57);
t63 = sin(pkin(6));
t99 = t63 * t71;
t14 = t37 * t53 - t52 * t99;
t36 = t68 * t67 - t70 * t80;
t66 = sin(qJ(6));
t69 = cos(qJ(6));
t110 = t14 * t66 - t36 * t69;
t109 = t14 * t69 + t36 * t66;
t108 = g(3) * t63;
t64 = cos(pkin(12));
t54 = t64 * pkin(3) + pkin(2);
t81 = t68 * t93;
t39 = -t67 * t81 + t71 * t70;
t55 = sin(t61);
t105 = t39 * t55;
t104 = t53 * t66;
t103 = t53 * t69;
t102 = t63 * t67;
t101 = t63 * t68;
t100 = t63 * t70;
t98 = t66 * t70;
t97 = t69 * t70;
t65 = -pkin(9) - qJ(3);
t56 = cos(t61);
t40 = pkin(4) * t56 + t54;
t60 = -pkin(10) + t65;
t96 = -t36 * t40 - t37 * t60;
t38 = t71 * t67 + t70 * t81;
t95 = -t38 * t40 - t39 * t60;
t94 = t71 * pkin(1) + pkin(8) * t101;
t92 = t55 * t102;
t91 = t56 * t101;
t62 = sin(pkin(12));
t90 = t62 * t101;
t89 = t56 * t99;
t88 = t62 * t99;
t87 = -t68 * pkin(1) + pkin(8) * t99;
t83 = -t37 * t52 - t53 * t99;
t86 = pkin(5) * t83 + t14 * pkin(11);
t17 = -t53 * t101 + t39 * t52;
t18 = t52 * t101 + t39 * t53;
t85 = -t17 * pkin(5) + t18 * pkin(11);
t27 = -t52 * t102 + t93 * t53;
t28 = t53 * t102 + t93 * t52;
t84 = t27 * pkin(5) + t28 * pkin(11);
t82 = t37 * t56 - t55 * t99;
t79 = t93 * t56;
t45 = t62 * pkin(3) + pkin(4) * t55;
t78 = t45 * t101 - t38 * t60 + t39 * t40 + t94;
t77 = pkin(5) * t53 + pkin(11) * t52;
t76 = g(1) * t83 + g(2) * t17;
t19 = g(1) * t36 - g(2) * t38;
t75 = g(1) * t71 + g(2) * t68;
t74 = t37 * t55 + t89;
t73 = t36 * t60 - t37 * t40 + t45 * t99 + t87;
t3 = g(1) * t17 - g(2) * t83 - g(3) * t27;
t5 = g(1) * t18 + g(2) * t14 + g(3) * t28;
t9 = -g(1) * t38 - g(2) * t36 + g(3) * t100;
t72 = g(1) * t39 + g(2) * t37 + g(3) * t102;
t47 = pkin(4) * t79;
t44 = pkin(4) * t91;
t33 = t40 * t100;
t21 = t55 * t101 + t39 * t56;
t20 = t91 - t105;
t8 = t18 * t69 + t38 * t66;
t7 = -t18 * t66 + t38 * t69;
t6 = t9 * t52;
t2 = t3 * t69;
t1 = t3 * t66;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t71, t75, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t37 - g(2) * t39, -t19, -t75 * t63, -g(1) * t87 - g(2) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (-t37 * t64 + t88) - g(2) * (t39 * t64 + t90) -g(1) * (t37 * t62 + t64 * t99) - g(2) * (t64 * t101 - t39 * t62) t19, -g(1) * (-t37 * pkin(2) - t36 * qJ(3) + t87) - g(2) * (t39 * pkin(2) + t38 * qJ(3) + t94) 0, 0, 0, 0, 0, 0, g(1) * t82 - g(2) * t21, -g(1) * t74 - g(2) * t20, t19, -g(1) * (pkin(3) * t88 + t36 * t65 - t37 * t54 + t87) - g(2) * (pkin(3) * t90 - t38 * t65 + t39 * t54 + t94) 0, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t76, t19, -g(1) * t73 - g(2) * t78, 0, 0, 0, 0, 0, 0, g(1) * t109 - g(2) * t8, -g(1) * t110 - g(2) * t7, -t76, -g(1) * (-pkin(5) * t14 + pkin(11) * t83 + t73) - g(2) * (t18 * pkin(5) + t17 * pkin(11) + t78); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t72, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t64, t9 * t62, -t72, -g(1) * (-t38 * pkin(2) + t39 * qJ(3)) - g(2) * (-t36 * pkin(2) + t37 * qJ(3)) - (pkin(2) * t70 + qJ(3) * t67) * t108, 0, 0, 0, 0, 0, 0, -t9 * t56, t9 * t55, -t72, -g(1) * (-t38 * t54 - t39 * t65) - g(2) * (-t36 * t54 - t37 * t65) - (t54 * t70 - t65 * t67) * t108, 0, 0, 0, 0, 0, 0, -t9 * t53, t6, -t72, -g(1) * t95 - g(2) * t96 - g(3) * (-t60 * t102 + t33) 0, 0, 0, 0, 0, 0, -g(1) * (-t38 * t103 + t39 * t66) - g(2) * (-t36 * t103 + t37 * t66) - (t53 * t97 + t66 * t67) * t108, -g(1) * (t38 * t104 + t39 * t69) - g(2) * (t36 * t104 + t37 * t69) - (-t53 * t98 + t67 * t69) * t108, -t6, -g(1) * (-t77 * t38 + t95) - g(2) * (-t77 * t36 + t96) - g(3) * t33 - (-t60 * t67 + t77 * t70) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + g(2) * t74 - g(3) * (t79 - t92) g(1) * t21 + g(2) * t82 - g(3) * (-t56 * t102 - t93 * t55) 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t44 - g(3) * t47 + (g(2) * t89 + t72 * t55) * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (-pkin(4) * t105 + t44 + t85) - g(2) * (-t74 * pkin(4) + t86) - g(3) * (-pkin(4) * t92 + t47 + t84); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t85 - g(2) * t86 - g(3) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t110 - g(3) * (-t28 * t66 - t63 * t97) g(1) * t8 + g(2) * t109 - g(3) * (-t28 * t69 + t63 * t98) 0, 0;];
taug_reg  = t4;

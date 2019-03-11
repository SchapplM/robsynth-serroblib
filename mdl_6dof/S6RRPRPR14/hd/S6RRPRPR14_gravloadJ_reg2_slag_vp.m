% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t53 = sin(qJ(2));
t54 = sin(qJ(1));
t57 = cos(qJ(2));
t58 = cos(qJ(1));
t88 = cos(pkin(6));
t84 = t58 * t88;
t34 = t53 * t84 + t54 * t57;
t85 = t54 * t88;
t36 = -t53 * t85 + t58 * t57;
t111 = -g(1) * t36 - g(2) * t34;
t51 = sin(qJ(6));
t55 = cos(qJ(6));
t33 = t54 * t53 - t57 * t84;
t52 = sin(qJ(4));
t56 = cos(qJ(4));
t50 = sin(pkin(6));
t94 = t50 * t58;
t67 = t33 * t56 + t52 * t94;
t110 = -t34 * t55 + t51 * t67;
t109 = t34 * t51 + t55 * t67;
t108 = pkin(5) + pkin(9);
t107 = pkin(4) * t52;
t104 = g(3) * t50;
t103 = t33 * pkin(9);
t102 = t34 * pkin(9);
t35 = t58 * t53 + t57 * t85;
t101 = t35 * pkin(9);
t100 = t36 * pkin(9);
t97 = t50 * t53;
t96 = t50 * t54;
t95 = t50 * t57;
t93 = t51 * t56;
t92 = t55 * t56;
t91 = pkin(2) * t95 + qJ(3) * t97;
t90 = t58 * pkin(1) + pkin(8) * t96;
t89 = qJ(5) * t56;
t87 = pkin(9) * t95 + t91;
t86 = -t54 * pkin(1) + pkin(8) * t94;
t27 = t33 * pkin(2);
t83 = t34 * qJ(3) - t27;
t29 = t35 * pkin(2);
t82 = t36 * qJ(3) - t29;
t15 = -t35 * t56 + t52 * t96;
t16 = t35 * t52 + t56 * t96;
t81 = -t15 * pkin(4) + t16 * qJ(5);
t68 = -t33 * t52 + t56 * t94;
t80 = pkin(4) * t67 - qJ(5) * t68;
t31 = t88 * t52 + t56 * t95;
t32 = -t52 * t95 + t88 * t56;
t79 = -t31 * pkin(4) + t32 * qJ(5);
t78 = t97 * t107 + t87;
t77 = qJ(3) - t89;
t76 = t34 * t107 - t103 - t27;
t75 = t36 * t107 - t101 - t29;
t74 = g(1) * t67 + g(2) * t15;
t73 = g(1) * t68 + g(2) * t16;
t72 = g(1) * t33 - g(2) * t35;
t10 = g(1) * t34 - g(2) * t36;
t71 = g(1) * t58 + g(2) * t54;
t70 = pkin(10) * t52 - t89;
t69 = t36 * pkin(2) + t35 * qJ(3) + t90;
t65 = pkin(3) * t96 + t69;
t64 = -t34 * pkin(2) - t33 * qJ(3) + t86;
t2 = g(1) * t15 - g(2) * t67 + g(3) * t31;
t63 = g(1) * t16 - g(2) * t68 + g(3) * t32;
t62 = pkin(3) * t94 + t64;
t8 = -g(1) * t35 - g(2) * t33 + g(3) * t95;
t61 = g(3) * t97 - t111;
t60 = t16 * pkin(4) + t15 * qJ(5) + t65;
t59 = pkin(4) * t68 + qJ(5) * t67 + t62;
t37 = t71 * t50;
t7 = t61 * t56;
t6 = t61 * t52;
t5 = t15 * t51 + t36 * t55;
t4 = t15 * t55 - t36 * t51;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t58, t71, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t72, -t37, -g(1) * t86 - g(2) * t90, 0, 0, 0, 0, 0, 0, -t37, -t10, t72, -g(1) * t64 - g(2) * t69, 0, 0, 0, 0, 0, 0, -t73, t74, t10, -g(1) * (t62 - t102) - g(2) * (t65 + t100) 0, 0, 0, 0, 0, 0, t10, t73, -t74, -g(1) * (t59 - t102) - g(2) * (t60 + t100) 0, 0, 0, 0, 0, 0, -g(1) * t110 - g(2) * t5, -g(1) * t109 - g(2) * t4, -t73, -g(1) * (pkin(10) * t68 - t108 * t34 + t59) - g(2) * (t16 * pkin(10) + t108 * t36 + t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t61, -g(1) * t82 - g(2) * t83 - g(3) * t91, 0, 0, 0, 0, 0, 0, -t6, -t7, -t8, -g(1) * (t82 - t101) - g(2) * (t83 - t103) - g(3) * t87, 0, 0, 0, 0, 0, 0, -t8, t6, t7, -g(1) * (t77 * t36 + t75) - g(2) * (t77 * t34 + t76) - g(3) * (-t89 * t97 + t78) 0, 0, 0, 0, 0, 0, -g(1) * (-t35 * t55 - t36 * t93) - g(2) * (-t33 * t55 - t34 * t93) - (-t53 * t93 + t55 * t57) * t104, -g(1) * (t35 * t51 - t36 * t92) - g(2) * (t33 * t51 - t34 * t92) - (-t51 * t57 - t53 * t92) * t104, -t6, -g(1) * (-t35 * pkin(5) + t75) - g(2) * (-t33 * pkin(5) + t76) - g(3) * t78 - (pkin(5) * t57 + t70 * t53) * t104 + t111 * (qJ(3) + t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t63, -g(1) * t81 - g(2) * t80 - g(3) * t79, 0, 0, 0, 0, 0, 0, -t63 * t51, -t63 * t55, t2, -g(1) * (-t15 * pkin(10) + t81) - g(2) * (pkin(10) * t67 + t80) - g(3) * (-t31 * pkin(10) + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t109 - g(3) * (t31 * t55 - t51 * t97) g(1) * t5 - g(2) * t110 - g(3) * (-t31 * t51 - t55 * t97) 0, 0;];
taug_reg  = t1;

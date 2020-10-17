% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:13:51
% EndTime: 2019-05-06 19:13:53
% DurationCPUTime: 0.61s
% Computational Cost: add. (501->127), mult. (1264->185), div. (0->0), fcn. (1529->10), ass. (0->82)
t53 = sin(qJ(2));
t57 = cos(qJ(2));
t58 = cos(qJ(1));
t54 = sin(qJ(1));
t88 = cos(pkin(6));
t82 = t54 * t88;
t34 = t58 * t53 + t57 * t82;
t52 = sin(qJ(4));
t56 = cos(qJ(4));
t49 = sin(pkin(6));
t98 = t49 * t54;
t18 = t34 * t52 + t56 * t98;
t35 = -t53 * t82 + t58 * t57;
t51 = sin(qJ(5));
t55 = cos(qJ(5));
t11 = -t18 * t51 + t35 * t55;
t81 = t58 * t88;
t33 = t53 * t81 + t54 * t57;
t32 = t54 * t53 - t57 * t81;
t96 = t49 * t58;
t70 = -t32 * t52 + t56 * t96;
t111 = t33 * t55 + t51 * t70;
t97 = t49 * t57;
t31 = -t52 * t97 + t88 * t56;
t91 = t53 * t55;
t1 = -g(2) * t111 - g(3) * (-t31 * t51 + t49 * t91) - g(1) * t11;
t112 = -g(1) * t35 - g(2) * t33;
t110 = -t33 * t51 + t55 * t70;
t106 = g(3) * t49;
t105 = t32 * pkin(9);
t104 = t34 * pkin(9);
t103 = t32 * t51;
t100 = t34 * t51;
t99 = t49 * t53;
t95 = t51 * t52;
t94 = t51 * t53;
t93 = t51 * t57;
t92 = t52 * t55;
t90 = pkin(2) * t97 + qJ(3) * t99;
t89 = t58 * pkin(1) + pkin(8) * t98;
t87 = pkin(9) * t97 + t90;
t86 = pkin(5) * t51 + pkin(9);
t85 = -t54 * pkin(1) + pkin(8) * t96;
t26 = t32 * pkin(2);
t84 = -t26 - t105;
t28 = t34 * pkin(2);
t83 = -t28 - t104;
t80 = t33 * qJ(3) - t26;
t79 = t35 * qJ(3) - t28;
t78 = g(3) * t87;
t77 = pkin(4) * t52 - pkin(10) * t56;
t17 = -t34 * t56 + t52 * t98;
t69 = t32 * t56 + t52 * t96;
t76 = g(1) * t69 + g(2) * t17;
t75 = g(1) * t32 - g(2) * t34;
t16 = g(1) * t33 - g(2) * t35;
t74 = g(1) * t58 + g(2) * t54;
t47 = t55 * pkin(5) + pkin(4);
t50 = -qJ(6) - pkin(10);
t73 = t47 * t52 + t50 * t56;
t72 = t35 * pkin(2) + t34 * qJ(3) + t89;
t67 = pkin(3) * t98 + t72;
t66 = -t33 * pkin(2) - t32 * qJ(3) + t85;
t30 = t88 * t52 + t56 * t97;
t65 = g(1) * t17 - g(2) * t69 + g(3) * t30;
t64 = g(1) * t18 - g(2) * t70 + g(3) * t31;
t63 = pkin(3) * t96 + t66;
t14 = -g(1) * t34 - g(2) * t32 + g(3) * t97;
t62 = g(3) * t99 - t112;
t61 = t35 * pkin(9) + t67;
t60 = -t33 * pkin(9) + t63;
t36 = t74 * t49;
t13 = t62 * t56;
t12 = t18 * t55 + t35 * t51;
t8 = t65 * t55;
t7 = t65 * t51;
t6 = -g(1) * t110 - g(2) * t12;
t5 = g(1) * t111 - g(2) * t11;
t4 = -g(1) * (t35 * t92 - t100) - g(2) * (t33 * t92 - t103) - (t52 * t91 + t93) * t106;
t3 = -g(1) * (-t34 * t55 - t35 * t95) - g(2) * (-t32 * t55 - t33 * t95) - (-t52 * t94 + t55 * t57) * t106;
t2 = g(1) * t12 - g(2) * t110 - g(3) * (-t31 * t55 - t49 * t94);
t9 = [0, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t58, t74, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t75, -t36, -g(1) * t85 - g(2) * t89, 0, 0, 0, 0, 0, 0, -t36, -t16, t75, -g(1) * t66 - g(2) * t72, 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t18, t76, t16, -g(1) * t60 - g(2) * t61, 0, 0, 0, 0, 0, 0, t6, t5, -t76, -g(1) * (pkin(4) * t70 + pkin(10) * t69 + t60) - g(2) * (t18 * pkin(4) + t17 * pkin(10) + t61) 0, 0, 0, 0, 0, 0, t6, t5, -t76, -g(1) * (-t86 * t33 + t47 * t70 - t50 * t69 + t63) - g(2) * (-t17 * t50 + t18 * t47 + t86 * t35 + t67); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t62, -g(1) * t79 - g(2) * t80 - g(3) * t90, 0, 0, 0, 0, 0, 0, -t62 * t52, -t13, -t14, -g(1) * (t79 - t104) - g(2) * (t80 - t105) - t78, 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * t83 - g(2) * t84 - g(3) * (t77 * t99 + t87) + t112 * (qJ(3) + t77) 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * (-pkin(5) * t100 + t83) - g(2) * (-pkin(5) * t103 + t84) - t78 - (pkin(5) * t93 + t73 * t53) * t106 + t112 * (qJ(3) + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t64, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t64, -g(1) * (-t17 * pkin(4) + t18 * pkin(10)) - g(2) * (pkin(4) * t69 - pkin(10) * t70) - g(3) * (-t30 * pkin(4) + t31 * pkin(10)) 0, 0, 0, 0, 0, 0, t8, -t7, -t64, -g(1) * (-t17 * t47 - t18 * t50) - g(2) * (t47 * t69 + t50 * t70) - g(3) * (-t30 * t47 - t31 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65;];
taug_reg  = t9;

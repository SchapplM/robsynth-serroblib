% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:57:05
% EndTime: 2019-05-05 20:57:08
% DurationCPUTime: 0.87s
% Computational Cost: add. (864->125), mult. (2393->196), div. (0->0), fcn. (3073->14), ass. (0->82)
t55 = sin(pkin(12));
t63 = sin(qJ(1));
t102 = t63 * t55;
t59 = cos(pkin(6));
t58 = cos(pkin(12));
t66 = cos(qJ(1));
t99 = t66 * t58;
t43 = -t59 * t99 + t102;
t56 = sin(pkin(7));
t57 = sin(pkin(6));
t96 = cos(pkin(7));
t93 = t57 * t96;
t116 = -t43 * t56 + t66 * t93;
t100 = t66 * t55;
t101 = t63 * t58;
t44 = t59 * t100 + t101;
t62 = sin(qJ(3));
t108 = cos(qJ(3));
t80 = t96 * t108;
t95 = t56 * t108;
t88 = t57 * t95;
t22 = t43 * t80 + t44 * t62 + t66 * t88;
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t10 = -t116 * t65 + t22 * t61;
t105 = t57 * t62;
t92 = t62 * t96;
t23 = -t66 * t56 * t105 + t44 * t108 - t43 * t92;
t60 = sin(qJ(6));
t64 = cos(qJ(6));
t120 = t10 * t60 - t23 * t64;
t119 = t10 * t64 + t23 * t60;
t78 = t59 * t101 + t100;
t35 = t78 * t56 + t63 * t93;
t115 = t116 * t61 + t22 * t65;
t45 = -t59 * t102 + t99;
t72 = t78 * t96;
t27 = t45 * t108 + (t63 * t57 * t56 - t72) * t62;
t30 = t59 * t56 * t62 + (t108 * t55 + t58 * t92) * t57;
t75 = g(1) * t27 + g(2) * t23 + g(3) * t30;
t111 = t22 * pkin(10);
t26 = t108 * t72 + t45 * t62 - t63 * t88;
t110 = t26 * pkin(10);
t106 = t57 * t58;
t29 = t55 * t105 - t80 * t106 - t59 * t95;
t109 = t29 * pkin(10);
t104 = t60 * t61;
t103 = t61 * t64;
t97 = qJ(2) * t57;
t98 = t66 * pkin(1) + t63 * t97;
t94 = -t63 * pkin(1) + t66 * t97;
t16 = t22 * pkin(3);
t91 = t23 * qJ(4) - t16;
t18 = t26 * pkin(3);
t90 = t27 * qJ(4) - t18;
t28 = t29 * pkin(3);
t89 = t30 * qJ(4) - t28;
t11 = -t26 * t65 + t35 * t61;
t85 = g(1) * t115 + g(2) * t11;
t84 = -g(1) * t22 + g(2) * t26;
t83 = -g(1) * t23 + g(2) * t27;
t82 = g(1) * t66 + g(2) * t63;
t81 = g(1) * t63 - g(2) * t66;
t42 = -t56 * t106 + t59 * t96;
t20 = t29 * t65 - t42 * t61;
t77 = g(1) * t11 - g(2) * t115 - g(3) * t20;
t12 = t26 * t61 + t35 * t65;
t21 = t29 * t61 + t42 * t65;
t76 = g(1) * t12 + g(2) * t10 + g(3) * t21;
t5 = g(1) * t26 + g(2) * t22 + g(3) * t29;
t73 = -t44 * pkin(2) + t116 * pkin(9) + t94;
t71 = -pkin(3) * t23 - qJ(4) * t22 + t73;
t70 = t45 * pkin(2) + t35 * pkin(9) + t98;
t69 = pkin(4) * t116 - pkin(10) * t23 + t71;
t68 = t27 * pkin(3) + t26 * qJ(4) + t70;
t67 = t35 * pkin(4) + t27 * pkin(10) + t68;
t39 = -g(3) * t59 - t81 * t57;
t13 = -g(1) * t116 - g(2) * t35;
t3 = t12 * t64 + t27 * t60;
t2 = -t12 * t60 + t27 * t64;
t1 = t75 * t65;
t4 = [0, 0, 0, 0, 0, 0, t81, t82, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t45, -g(1) * t43 + g(2) * t78, -t82 * t57, -g(1) * t94 - g(2) * t98, 0, 0, 0, 0, 0, 0, -t83, t84, t13, -g(1) * t73 - g(2) * t70, 0, 0, 0, 0, 0, 0, t13, t83, -t84, -g(1) * t71 - g(2) * t68, 0, 0, 0, 0, 0, 0, g(1) * t10 - g(2) * t12, t85, -t83, -g(1) * t69 - g(2) * t67, 0, 0, 0, 0, 0, 0, g(1) * t119 - g(2) * t3, -g(1) * t120 - g(2) * t2, -t85, -g(1) * (-pkin(5) * t10 + pkin(11) * t115 + t69) - g(2) * (t12 * pkin(5) + t11 * pkin(11) + t67); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t75, -g(1) * t90 - g(2) * t91 - g(3) * t89, 0, 0, 0, 0, 0, 0, -t75 * t61, -t1, t5, -g(1) * (t90 - t110) - g(2) * (t91 - t111) - g(3) * (t89 - t109) 0, 0, 0, 0, 0, 0, -g(1) * (t27 * t103 - t26 * t60) - g(2) * (t23 * t103 - t22 * t60) - g(3) * (t30 * t103 - t29 * t60) -g(1) * (-t27 * t104 - t26 * t64) - g(2) * (-t23 * t104 - t22 * t64) - g(3) * (-t30 * t104 - t29 * t64) t1, -g(1) * (-t18 - t110) - g(2) * (-t16 - t111) - g(3) * (-t28 - t109) - t75 * (pkin(5) * t61 - pkin(11) * t65 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t64, -t77 * t60, -t76, -g(1) * (-t11 * pkin(5) + t12 * pkin(11)) - g(2) * (pkin(5) * t115 + t10 * pkin(11)) - g(3) * (t20 * pkin(5) + t21 * pkin(11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t120 - g(3) * (-t21 * t60 + t30 * t64) g(1) * t3 + g(2) * t119 - g(3) * (-t21 * t64 - t30 * t60) 0, 0;];
taug_reg  = t4;

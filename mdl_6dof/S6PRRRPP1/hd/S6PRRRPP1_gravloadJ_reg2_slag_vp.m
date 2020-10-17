% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:37:38
% EndTime: 2019-05-05 06:37:40
% DurationCPUTime: 0.68s
% Computational Cost: add. (594->141), mult. (1275->208), div. (0->0), fcn. (1561->12), ass. (0->86)
t64 = -qJ(5) - pkin(9);
t66 = sin(qJ(3));
t102 = t64 * t66;
t68 = cos(qJ(4));
t59 = t68 * pkin(4) + pkin(3);
t69 = cos(qJ(3));
t115 = -t59 * t69 + t102;
t63 = sin(pkin(6));
t114 = g(3) * t63;
t67 = sin(qJ(2));
t70 = cos(qJ(2));
t92 = cos(pkin(10));
t93 = cos(pkin(6));
t82 = t93 * t92;
t91 = sin(pkin(10));
t43 = t67 * t82 + t91 * t70;
t86 = t63 * t92;
t23 = t43 * t69 - t66 * t86;
t65 = sin(qJ(4));
t113 = t23 * t65;
t81 = t93 * t91;
t45 = -t67 * t81 + t92 * t70;
t85 = t63 * t91;
t25 = t45 * t69 + t66 * t85;
t112 = t25 * t65;
t42 = t91 * t67 - t70 * t82;
t111 = t42 * t68;
t110 = t43 * t65;
t44 = t92 * t67 + t70 * t81;
t109 = t44 * t68;
t108 = t45 * t65;
t62 = qJ(4) + pkin(11);
t60 = sin(t62);
t106 = t60 * t69;
t61 = cos(t62);
t105 = t61 * t69;
t104 = t63 * t67;
t103 = t63 * t70;
t101 = t65 * t67;
t100 = t65 * t69;
t99 = t68 * t69;
t98 = t69 * t70;
t22 = t43 * t66 + t69 * t86;
t97 = -t22 * t59 - t23 * t64;
t24 = t45 * t66 - t69 * t85;
t96 = -t24 * t59 - t25 * t64;
t46 = t66 * t104 - t93 * t69;
t47 = t69 * t104 + t93 * t66;
t95 = -t46 * t59 - t47 * t64;
t94 = pkin(2) * t103 + pkin(8) * t104;
t90 = t60 * t103;
t89 = t68 * t103;
t88 = -t42 * pkin(2) + t43 * pkin(8);
t87 = -t44 * pkin(2) + t45 * pkin(8);
t84 = pkin(3) * t69 + pkin(9) * t66;
t83 = -pkin(5) * t61 - qJ(6) * t60;
t80 = -t47 * t65 - t89;
t18 = t61 * t103 + t47 * t60;
t7 = t23 * t60 - t42 * t61;
t9 = t25 * t60 - t44 * t61;
t1 = g(1) * t9 + g(2) * t7 + g(3) * t18;
t10 = t25 * t61 + t44 * t60;
t19 = t47 * t61 - t90;
t8 = t23 * t61 + t42 * t60;
t79 = g(1) * t10 + g(2) * t8 + g(3) * t19;
t12 = -t42 * t106 - t43 * t61;
t14 = -t44 * t106 - t45 * t61;
t26 = -t61 * t104 + t69 * t90;
t78 = g(1) * t14 + g(2) * t12 + g(3) * t26;
t77 = g(1) * t24 + g(2) * t22 + g(3) * t46;
t76 = g(1) * t25 + g(2) * t23 + g(3) * t47;
t75 = -g(1) * t44 - g(2) * t42 + g(3) * t103;
t74 = g(1) * t45 + g(2) * t43 + g(3) * t104;
t73 = -t102 * t103 + t94 + (pkin(4) * t101 + t59 * t98) * t63;
t72 = pkin(4) * t110 + t115 * t42 + t88;
t71 = pkin(4) * t108 + t115 * t44 + t87;
t37 = pkin(4) * t109;
t35 = pkin(4) * t111;
t27 = (t60 * t67 + t61 * t98) * t63;
t15 = -t44 * t105 + t45 * t60;
t13 = -t42 * t105 + t43 * t60;
t11 = t75 * t66;
t4 = t77 * t61;
t3 = t77 * t60;
t2 = -g(1) * t15 - g(2) * t13 - g(3) * t27;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t74, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t69, t11, -t74, -g(1) * t87 - g(2) * t88 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (-t44 * t99 + t108) - g(2) * (-t42 * t99 + t110) - (t68 * t98 + t101) * t114, -g(1) * (t44 * t100 + t45 * t68) - g(2) * (t42 * t100 + t43 * t68) - (-t65 * t98 + t67 * t68) * t114, -t11, -g(1) * (-t84 * t44 + t87) - g(2) * (-t84 * t42 + t88) - g(3) * (t84 * t103 + t94) 0, 0, 0, 0, 0, 0, t2, t78, -t11, -g(1) * t71 - g(2) * t72 - g(3) * t73, 0, 0, 0, 0, 0, 0, t2, -t11, -t78, -g(1) * (t15 * pkin(5) + t14 * qJ(6) + t71) - g(2) * (t13 * pkin(5) + t12 * qJ(6) + t72) - g(3) * (t27 * pkin(5) + t26 * qJ(6) + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t68, -t77 * t65, -t76, -g(1) * (-t24 * pkin(3) + t25 * pkin(9)) - g(2) * (-t22 * pkin(3) + t23 * pkin(9)) - g(3) * (-t46 * pkin(3) + t47 * pkin(9)) 0, 0, 0, 0, 0, 0, t4, -t3, -t76, -g(1) * t96 - g(2) * t97 - g(3) * t95, 0, 0, 0, 0, 0, 0, t4, -t76, t3, -g(1) * (t83 * t24 + t96) - g(2) * (t83 * t22 + t97) - g(3) * (t83 * t46 + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t109 - t112) - g(2) * (t111 - t113) - g(3) * t80, -g(1) * (-t25 * t68 - t44 * t65) - g(2) * (-t23 * t68 - t42 * t65) - g(3) * (t65 * t103 - t47 * t68) 0, 0, 0, 0, 0, 0, 0, 0, t1, t79, 0, -g(1) * t37 - g(2) * t35 + (g(3) * t89 + t76 * t65) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t79, -g(1) * (-pkin(4) * t112 - t9 * pkin(5) + t10 * qJ(6) + t37) - g(2) * (-pkin(4) * t113 - t7 * pkin(5) + t8 * qJ(6) + t35) - g(3) * (t80 * pkin(4) - t18 * pkin(5) + t19 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;

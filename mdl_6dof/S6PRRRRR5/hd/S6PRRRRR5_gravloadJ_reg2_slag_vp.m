% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:01:40
% EndTime: 2019-05-05 12:01:42
% DurationCPUTime: 1.11s
% Computational Cost: add. (1232->198), mult. (3301->322), div. (0->0), fcn. (4249->16), ass. (0->101)
t112 = cos(pkin(13));
t113 = cos(pkin(7));
t109 = sin(pkin(13));
t120 = sin(qJ(2));
t122 = cos(qJ(2));
t114 = cos(pkin(6));
t92 = t114 * t112;
t68 = t109 * t120 - t122 * t92;
t110 = sin(pkin(7));
t111 = sin(pkin(6));
t88 = t111 * t110;
t125 = t112 * t88 + t68 * t113;
t90 = t114 * t109;
t69 = t112 * t120 + t122 * t90;
t87 = t111 * t109;
t124 = -t110 * t87 + t69 * t113;
t89 = t113 * t111;
t123 = t114 * t110 + t122 * t89;
t121 = cos(qJ(3));
t58 = qJ(5) + qJ(6);
t56 = sin(t58);
t63 = cos(qJ(4));
t119 = t56 * t63;
t57 = cos(t58);
t118 = t57 * t63;
t59 = sin(qJ(5));
t117 = t59 * t63;
t62 = cos(qJ(5));
t116 = t62 * t63;
t74 = t120 * t88;
t96 = t122 * t111;
t115 = pkin(2) * t96 + pkin(9) * t74;
t61 = sin(qJ(3));
t76 = t120 * t89;
t41 = t121 * t96 - t61 * t76;
t108 = t41 * pkin(3) + t115;
t107 = pkin(5) * t59 + pkin(10);
t45 = t109 * t122 + t120 * t92;
t16 = t125 * t121 + t45 * t61;
t14 = t16 * pkin(3);
t17 = t45 * t121 - t125 * t61;
t106 = t17 * pkin(10) - t14;
t46 = t112 * t122 - t120 * t90;
t18 = t124 * t121 + t46 * t61;
t15 = t18 * pkin(3);
t19 = t46 * t121 - t124 * t61;
t105 = t19 * pkin(10) - t15;
t95 = t111 * t120;
t31 = -t123 * t121 + t61 * t95;
t30 = t31 * pkin(3);
t32 = t121 * t95 + t123 * t61;
t104 = t32 * pkin(10) - t30;
t103 = t45 * t110;
t102 = t46 * t110;
t60 = sin(qJ(4));
t101 = t60 * t110;
t100 = t61 * t113;
t99 = t63 * t110;
t98 = -pkin(4) * t63 - pkin(11) * t60;
t97 = t113 * t121;
t55 = t62 * pkin(5) + pkin(4);
t64 = -pkin(12) - pkin(11);
t94 = -t55 * t63 + t60 * t64;
t40 = t121 * t76 + t61 * t96;
t93 = t40 * pkin(10) + t108;
t86 = -t68 * pkin(2) + pkin(9) * t103;
t85 = -t69 * pkin(2) + pkin(9) * t102;
t44 = t114 * t113 - t122 * t88;
t20 = -t32 * t60 + t44 * t63;
t33 = t68 * t110 - t112 * t89;
t6 = -t17 * t60 + t33 * t63;
t34 = t69 * t110 + t113 * t87;
t8 = -t19 * t60 + t34 * t63;
t84 = g(1) * t8 + g(2) * t6 + g(3) * t20;
t21 = t32 * t63 + t44 * t60;
t7 = t17 * t63 + t33 * t60;
t9 = t19 * t63 + t34 * t60;
t83 = g(1) * t9 + g(2) * t7 + g(3) * t21;
t25 = -t45 * t100 - t68 * t121;
t82 = t25 * pkin(3) + t86;
t27 = -t46 * t100 - t69 * t121;
t81 = t27 * pkin(3) + t85;
t10 = t25 * t60 - t45 * t99;
t12 = t27 * t60 - t46 * t99;
t28 = t41 * t60 - t63 * t74;
t80 = g(1) * t12 + g(2) * t10 + g(3) * t28;
t79 = g(1) * t18 + g(2) * t16 + g(3) * t31;
t78 = g(1) * t19 + g(2) * t17 + g(3) * t32;
t24 = t45 * t97 - t68 * t61;
t26 = t46 * t97 - t69 * t61;
t77 = g(1) * t26 + g(2) * t24 + g(3) * t40;
t71 = t24 * pkin(10) + t82;
t70 = t26 * pkin(10) + t81;
t65 = -g(1) * (t18 * t62 - t9 * t59) - g(2) * (t16 * t62 - t7 * t59) - g(3) * (-t21 * t59 + t31 * t62);
t29 = t41 * t63 + t60 * t74;
t13 = t46 * t101 + t27 * t63;
t11 = t45 * t101 + t25 * t63;
t5 = t79 * t60;
t2 = -g(1) * (-t18 * t56 - t9 * t57) - g(2) * (-t16 * t56 - t7 * t57) - g(3) * (-t21 * t57 - t31 * t56);
t1 = -g(1) * (t18 * t57 - t9 * t56) - g(2) * (t16 * t57 - t7 * t56) - g(3) * (-t21 * t56 + t31 * t57);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t69 + g(2) * t68 - g(3) * t96, g(1) * t46 + g(2) * t45 + g(3) * t95, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t25 - g(3) * t41, t77, -g(1) * t102 - g(2) * t103 - g(3) * t74, -g(1) * t85 - g(2) * t86 - g(3) * t115, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t29, t80, -t77, -g(1) * t70 - g(2) * t71 - g(3) * t93, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t62 + t26 * t59) - g(2) * (t11 * t62 + t24 * t59) - g(3) * (t29 * t62 + t40 * t59) -g(1) * (-t13 * t59 + t26 * t62) - g(2) * (-t11 * t59 + t24 * t62) - g(3) * (-t29 * t59 + t40 * t62) -t80, -g(1) * (t13 * pkin(4) + t12 * pkin(11) + t70) - g(2) * (t11 * pkin(4) + t10 * pkin(11) + t71) - g(3) * (t29 * pkin(4) + t28 * pkin(11) + t93) 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t57 + t26 * t56) - g(2) * (t11 * t57 + t24 * t56) - g(3) * (t29 * t57 + t40 * t56) -g(1) * (-t13 * t56 + t26 * t57) - g(2) * (-t11 * t56 + t24 * t57) - g(3) * (-t29 * t56 + t40 * t57) -t80, -g(1) * (t107 * t26 - t12 * t64 + t13 * t55 + t81) - g(2) * (-t10 * t64 + t107 * t24 + t11 * t55 + t82) - g(3) * (t107 * t40 - t28 * t64 + t29 * t55 + t108); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t78, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t63, -t5, -t78, -g(1) * t105 - g(2) * t106 - g(3) * t104, 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t116 + t19 * t59) - g(2) * (-t16 * t116 + t17 * t59) - g(3) * (-t31 * t116 + t32 * t59) -g(1) * (t18 * t117 + t19 * t62) - g(2) * (t16 * t117 + t17 * t62) - g(3) * (t31 * t117 + t32 * t62) t5, -g(1) * (t98 * t18 + t105) - g(2) * (t98 * t16 + t106) - g(3) * (t98 * t31 + t104) 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t118 + t19 * t56) - g(2) * (-t16 * t118 + t17 * t56) - g(3) * (-t31 * t118 + t32 * t56) -g(1) * (t18 * t119 + t19 * t57) - g(2) * (t16 * t119 + t17 * t57) - g(3) * (t31 * t119 + t32 * t57) t5, -g(1) * (t107 * t19 + t94 * t18 - t15) - g(2) * (t107 * t17 + t94 * t16 - t14) - g(3) * (t107 * t32 + t94 * t31 - t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t83, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t62, t84 * t59, -t83, -g(1) * (t8 * pkin(4) + t9 * pkin(11)) - g(2) * (t6 * pkin(4) + t7 * pkin(11)) - g(3) * (t20 * pkin(4) + t21 * pkin(11)) 0, 0, 0, 0, 0, 0, -t84 * t57, t84 * t56, -t83, -g(1) * (t8 * t55 - t9 * t64) - g(2) * (t6 * t55 - t7 * t64) - g(3) * (t20 * t55 - t21 * t64); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -g(1) * (-t18 * t59 - t9 * t62) - g(2) * (-t16 * t59 - t7 * t62) - g(3) * (-t21 * t62 - t31 * t59) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t65 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;

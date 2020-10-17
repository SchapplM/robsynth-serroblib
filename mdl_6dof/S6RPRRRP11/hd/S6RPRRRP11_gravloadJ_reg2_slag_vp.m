% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:08:30
% EndTime: 2019-05-06 02:08:34
% DurationCPUTime: 1.23s
% Computational Cost: add. (1167->132), mult. (3228->213), div. (0->0), fcn. (4170->14), ass. (0->85)
t111 = sin(pkin(7));
t121 = cos(qJ(1));
t110 = sin(pkin(12));
t119 = sin(qJ(1));
t113 = cos(pkin(12));
t115 = cos(pkin(6));
t97 = t115 * t113;
t82 = t119 * t110 - t121 * t97;
t112 = sin(pkin(6));
t114 = cos(pkin(7));
t94 = t114 * t112;
t125 = t82 * t111 - t121 * t94;
t120 = cos(qJ(3));
t93 = t112 * t111;
t130 = t82 * t114 + t121 * t93;
t95 = t115 * t110;
t45 = t119 * t113 + t121 * t95;
t61 = sin(qJ(3));
t29 = -t45 * t120 + t130 * t61;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t17 = -t125 * t60 + t29 * t63;
t26 = t130 * t120 + t45 * t61;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t134 = t17 * t59 + t26 * t62;
t133 = t17 * t62 - t26 * t59;
t16 = t125 * t63 + t29 * t60;
t77 = t121 * t110 + t119 * t97;
t129 = t77 * t111 + t119 * t94;
t126 = t77 * t114 - t119 * t93;
t124 = t115 * t111 + t113 * t94;
t46 = t121 * t113 - t119 * t95;
t31 = t46 * t120 - t126 * t61;
t19 = t129 * t60 + t31 * t63;
t30 = t126 * t120 + t46 * t61;
t12 = -t19 * t59 + t30 * t62;
t92 = t112 * t110;
t37 = t120 * t92 + t124 * t61;
t75 = -t113 * t93 + t115 * t114;
t25 = t37 * t63 + t75 * t60;
t36 = -t124 * t120 + t61 * t92;
t1 = -g(1) * t12 - g(2) * t134 - g(3) * (-t25 * t59 + t36 * t62);
t118 = t59 * t63;
t117 = t62 * t63;
t99 = t112 * t119;
t116 = t121 * pkin(1) + qJ(2) * t99;
t109 = t59 * pkin(5) + pkin(10);
t20 = t26 * pkin(3);
t108 = -pkin(10) * t29 - t20;
t22 = t30 * pkin(3);
t107 = t31 * pkin(10) - t22;
t35 = t36 * pkin(3);
t106 = t37 * pkin(10) - t35;
t100 = t121 * t112;
t105 = -t119 * pkin(1) + qJ(2) * t100;
t104 = -pkin(4) * t63 - pkin(11) * t60;
t18 = -t129 * t63 + t31 * t60;
t103 = g(1) * t16 + g(2) * t18;
t102 = -g(1) * t26 + g(2) * t30;
t56 = t62 * pkin(5) + pkin(4);
t58 = -qJ(6) - pkin(11);
t98 = -t56 * t63 + t58 * t60;
t24 = t37 * t60 - t75 * t63;
t91 = g(1) * t18 - g(2) * t16 + g(3) * t24;
t90 = g(1) * t19 - g(2) * t17 + g(3) * t25;
t89 = g(1) * t30 + g(2) * t26 + g(3) * t36;
t88 = g(1) * t31 - g(2) * t29 + g(3) * t37;
t72 = -t45 * pkin(2) - t125 * pkin(9) + t105;
t69 = t29 * pkin(3) + t72;
t68 = t46 * pkin(2) + t129 * pkin(9) + t116;
t67 = t31 * pkin(3) + t68;
t66 = -pkin(10) * t26 + t69;
t65 = t30 * pkin(10) + t67;
t42 = -g(1) * t99 + g(2) * t100 - g(3) * t115;
t13 = t19 * t62 + t30 * t59;
t11 = t89 * t60;
t8 = t91 * t62;
t7 = t91 * t59;
t6 = -g(1) * t133 - g(2) * t13;
t5 = g(1) * t134 - g(2) * t12;
t4 = -g(1) * (-t30 * t117 + t31 * t59) - g(2) * (-t26 * t117 - t29 * t59) - g(3) * (-t36 * t117 + t37 * t59);
t3 = -g(1) * (t30 * t118 + t31 * t62) - g(2) * (t26 * t118 - t29 * t62) - g(3) * (t36 * t118 + t37 * t62);
t2 = g(1) * t13 - g(2) * t133 - g(3) * (-t25 * t62 - t36 * t59);
t9 = [0, 0, 0, 0, 0, 0, g(1) * t119 - g(2) * t121, g(1) * t121 + g(2) * t119, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t45 - g(2) * t46, -g(1) * t82 + g(2) * t77, -g(1) * t100 - g(2) * t99, -g(1) * t105 - g(2) * t116, 0, 0, 0, 0, 0, 0, -g(1) * t29 - g(2) * t31, t102, g(1) * t125 - g(2) * t129, -g(1) * t72 - g(2) * t68, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, t103, -t102, -g(1) * t66 - g(2) * t65, 0, 0, 0, 0, 0, 0, t6, t5, -t103, -g(1) * (t17 * pkin(4) + t16 * pkin(11) + t66) - g(2) * (t19 * pkin(4) + t18 * pkin(11) + t65) 0, 0, 0, 0, 0, 0, t6, t5, -t103, -g(1) * (-t109 * t26 - t16 * t58 + t17 * t56 + t69) - g(2) * (t109 * t30 - t18 * t58 + t19 * t56 + t67); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t63, -t11, -t88, -g(1) * t107 - g(2) * t108 - g(3) * t106, 0, 0, 0, 0, 0, 0, t4, t3, t11, -g(1) * (t104 * t30 + t107) - g(2) * (t104 * t26 + t108) - g(3) * (t104 * t36 + t106) 0, 0, 0, 0, 0, 0, t4, t3, t11, -g(1) * (t109 * t31 + t98 * t30 - t22) - g(2) * (-t109 * t29 + t98 * t26 - t20) - g(3) * (t109 * t37 + t98 * t36 - t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t90, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t90, -g(1) * (-t18 * pkin(4) + t19 * pkin(11)) - g(2) * (pkin(4) * t16 - pkin(11) * t17) - g(3) * (-t24 * pkin(4) + t25 * pkin(11)) 0, 0, 0, 0, 0, 0, t8, -t7, -t90, -g(1) * (-t18 * t56 - t19 * t58) - g(2) * (t16 * t56 + t17 * t58) - g(3) * (-t24 * t56 - t25 * t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91;];
taug_reg  = t9;

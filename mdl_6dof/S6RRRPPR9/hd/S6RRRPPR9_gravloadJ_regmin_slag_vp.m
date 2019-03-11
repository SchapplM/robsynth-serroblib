% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t110 = cos(qJ(3));
t109 = sin(qJ(1));
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t103 = cos(pkin(6));
t111 = cos(qJ(1));
t87 = t103 * t111;
t51 = t109 * t74 + t72 * t87;
t71 = sin(qJ(3));
t68 = sin(pkin(6));
t99 = t68 * t111;
t28 = t51 * t110 - t71 * t99;
t50 = t109 * t72 - t74 * t87;
t67 = sin(pkin(11));
t69 = cos(pkin(11));
t10 = t28 * t67 - t50 * t69;
t11 = t28 * t69 + t50 * t67;
t70 = sin(qJ(6));
t73 = cos(qJ(6));
t117 = t10 * t73 - t11 * t70;
t116 = t10 * t70 + t11 * t73;
t115 = qJ(4) * t71 + pkin(2);
t98 = t68 * t110;
t27 = t111 * t98 + t51 * t71;
t86 = t103 * t109;
t53 = t111 * t74 - t72 * t86;
t97 = t68 * t109;
t31 = -t110 * t97 + t53 * t71;
t106 = t68 * t72;
t48 = -t103 * t110 + t71 * t106;
t79 = g(1) * t31 + g(2) * t27 + g(3) * t48;
t105 = t68 * t74;
t102 = t50 * t110;
t52 = t111 * t72 + t74 * t86;
t101 = t52 * t110;
t100 = t67 * t110;
t96 = t69 * t110;
t95 = t74 * t110;
t94 = -t27 * pkin(3) + t28 * qJ(4);
t32 = t53 * t110 + t71 * t97;
t93 = -t31 * pkin(3) + t32 * qJ(4);
t49 = t103 * t71 + t72 * t98;
t92 = -t48 * pkin(3) + t49 * qJ(4);
t90 = t68 * t95;
t91 = pkin(3) * t90 + pkin(9) * t106 + t115 * t105;
t14 = t32 * t67 - t52 * t69;
t89 = -g(1) * t10 + g(2) * t14;
t88 = -g(1) * t27 + g(2) * t31;
t85 = -pkin(4) * t69 - qJ(5) * t67;
t82 = -pkin(3) * t102 + t51 * pkin(9) - t115 * t50;
t81 = -pkin(3) * t101 + t53 * pkin(9) - t115 * t52;
t17 = -t50 * t100 - t51 * t69;
t19 = -t100 * t52 - t53 * t69;
t33 = -t69 * t106 + t67 * t90;
t80 = g(1) * t19 + g(2) * t17 + g(3) * t33;
t78 = g(1) * t32 + g(2) * t28 + g(3) * t49;
t77 = -g(1) * t52 - g(2) * t50 + g(3) * t105;
t76 = t111 * pkin(1) + t53 * pkin(2) + t32 * pkin(3) + pkin(8) * t97 + t52 * pkin(9) + t31 * qJ(4);
t75 = -t109 * pkin(1) - t51 * pkin(2) - pkin(3) * t28 + pkin(8) * t99 - t50 * pkin(9) - qJ(4) * t27;
t34 = (t67 * t72 + t69 * t95) * t68;
t26 = -t67 * t105 + t49 * t69;
t25 = t69 * t105 + t49 * t67;
t20 = -t52 * t96 + t53 * t67;
t18 = -t50 * t96 + t51 * t67;
t16 = t77 * t71;
t15 = t32 * t69 + t52 * t67;
t6 = t79 * t69;
t5 = t79 * t67;
t4 = g(1) * t11 - g(2) * t15;
t3 = -g(1) * t20 - g(2) * t18 - g(3) * t34;
t2 = t14 * t70 + t15 * t73;
t1 = t14 * t73 - t15 * t70;
t7 = [0, g(1) * t109 - g(2) * t111, g(1) * t111 + g(2) * t109, 0, 0, 0, 0, 0, g(1) * t51 - g(2) * t53, -g(1) * t50 + g(2) * t52, 0, 0, 0, 0, 0, g(1) * t28 - g(2) * t32, t88, t4, t89, -t88, -g(1) * t75 - g(2) * t76, t4, -t88, -t89, -g(1) * (-pkin(4) * t11 - qJ(5) * t10 + t75) - g(2) * (t15 * pkin(4) + t14 * qJ(5) + t76) 0, 0, 0, 0, 0, g(1) * t116 - g(2) * t2, g(1) * t117 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t77, g(1) * t53 + g(2) * t51 + g(3) * t106, 0, 0, 0, 0, 0, g(1) * t101 + g(2) * t102 - g(3) * t90, t16, t3, t80, -t16, -g(1) * t81 - g(2) * t82 - g(3) * t91, t3, -t16, -t80, -g(1) * (t20 * pkin(4) + t19 * qJ(5) + t81) - g(2) * (t18 * pkin(4) + t17 * qJ(5) + t82) - g(3) * (t34 * pkin(4) + t33 * qJ(5) + t91) 0, 0, 0, 0, 0, -g(1) * (t19 * t70 + t20 * t73) - g(2) * (t17 * t70 + t18 * t73) - g(3) * (t33 * t70 + t34 * t73) -g(1) * (t19 * t73 - t20 * t70) - g(2) * (t17 * t73 - t18 * t70) - g(3) * (t33 * t73 - t34 * t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t78, t6, -t5, -t78, -g(1) * t93 - g(2) * t94 - g(3) * t92, t6, -t78, t5, -g(1) * (t85 * t31 + t93) - g(2) * (t85 * t27 + t94) - g(3) * (t85 * t48 + t92) 0, 0, 0, 0, 0, t79 * (t67 * t70 + t69 * t73) t79 * (t67 * t73 - t69 * t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, 0, 0, -t79, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t10 - g(3) * t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t117 - g(3) * (t25 * t73 - t26 * t70) g(1) * t2 + g(2) * t116 - g(3) * (-t25 * t70 - t26 * t73);];
taug_reg  = t7;

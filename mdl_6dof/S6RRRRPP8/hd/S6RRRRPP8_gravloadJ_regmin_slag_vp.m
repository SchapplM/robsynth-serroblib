% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t87 = sin(qJ(4));
t134 = qJ(5) * t87 + pkin(3);
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t133 = -pkin(3) * t91 - pkin(10) * t88 - pkin(2);
t130 = cos(qJ(1));
t86 = sin(pkin(6));
t109 = t86 * t130;
t115 = cos(pkin(6));
t102 = t115 * t130;
t129 = sin(qJ(1));
t89 = sin(qJ(2));
t92 = cos(qJ(2));
t69 = t89 * t102 + t129 * t92;
t42 = -t88 * t109 + t69 * t91;
t68 = -t92 * t102 + t129 * t89;
t90 = cos(qJ(4));
t17 = t42 * t87 - t68 * t90;
t18 = t42 * t90 + t68 * t87;
t107 = -t91 * t109 - t69 * t88;
t128 = t107 * t90;
t108 = t86 * t129;
t101 = t115 * t129;
t71 = -t89 * t101 + t130 * t92;
t45 = -t91 * t108 + t71 * t88;
t127 = t45 * t90;
t123 = t86 * t89;
t66 = t115 * t91 - t88 * t123;
t126 = t66 * t90;
t122 = t86 * t92;
t121 = t87 * t91;
t120 = t90 * t91;
t119 = t91 * t92;
t118 = pkin(10) - qJ(6);
t116 = qJ(6) * t88;
t114 = t88 * t122;
t113 = t87 * t122;
t112 = pkin(4) * t128 + t134 * t107;
t111 = -pkin(4) * t127 - t134 * t45;
t110 = pkin(4) * t126 + t134 * t66;
t106 = -t17 * pkin(4) + t18 * qJ(5);
t46 = t88 * t108 + t71 * t91;
t70 = t92 * t101 + t130 * t89;
t21 = t46 * t87 - t70 * t90;
t22 = t46 * t90 + t70 * t87;
t105 = -t21 * pkin(4) + t22 * qJ(5);
t67 = t115 * t88 + t91 * t123;
t39 = t90 * t122 + t67 * t87;
t40 = t67 * t90 - t113;
t104 = -t39 * pkin(4) + t40 * qJ(5);
t103 = -g(1) * t17 + g(2) * t21;
t12 = g(1) * t107 + g(2) * t45;
t2 = g(1) * t21 + g(2) * t17 + g(3) * t39;
t100 = g(1) * t22 + g(2) * t18 + g(3) * t40;
t26 = -t68 * t121 - t69 * t90;
t28 = -t70 * t121 - t71 * t90;
t48 = t91 * t113 - t90 * t123;
t99 = g(1) * t28 + g(2) * t26 + g(3) * t48;
t10 = g(1) * t45 - g(2) * t107 - g(3) * t66;
t11 = g(1) * t46 + g(2) * t42 + g(3) * t67;
t49 = (t90 * t119 + t87 * t89) * t86;
t98 = t86 * pkin(3) * t119 + pkin(2) * t122 + t49 * pkin(4) + pkin(9) * t123 + pkin(10) * t114 + t48 * qJ(5);
t97 = -g(1) * t70 - g(2) * t68 + g(3) * t122;
t27 = -t68 * t120 + t69 * t87;
t96 = t27 * pkin(4) + t69 * pkin(9) + t26 * qJ(5) + t133 * t68;
t29 = -t70 * t120 + t71 * t87;
t95 = t29 * pkin(4) + t71 * pkin(9) + t28 * qJ(5) + t133 * t70;
t94 = t130 * pkin(1) + t71 * pkin(2) + t46 * pkin(3) + t22 * pkin(4) + pkin(8) * t108 + t70 * pkin(9) + t21 * qJ(5);
t93 = -t129 * pkin(1) - t69 * pkin(2) - pkin(3) * t42 - pkin(4) * t18 + pkin(8) * t109 - t68 * pkin(9) - qJ(5) * t17;
t23 = t97 * t88;
t9 = t10 * t90;
t8 = t10 * t87;
t7 = g(1) * t18 - g(2) * t22;
t5 = -g(1) * t29 - g(2) * t27 - g(3) * t49;
t1 = [0, g(1) * t129 - g(2) * t130, g(1) * t130 + g(2) * t129, 0, 0, 0, 0, 0, g(1) * t69 - g(2) * t71, -g(1) * t68 + g(2) * t70, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t46, t12, 0, 0, 0, 0, 0, t7, t103, t7, -t12, -t103, -g(1) * (pkin(10) * t107 + t93) - g(2) * (t45 * pkin(10) + t94) t7, -t103, t12, -g(1) * (-pkin(5) * t18 + t107 * t118 + t93) - g(2) * (t22 * pkin(5) + t118 * t45 + t94); 0, 0, 0, 0, 0, 0, 0, 0, -t97, g(1) * t71 + g(2) * t69 + g(3) * t123, 0, 0, 0, 0, 0, -t97 * t91, t23, 0, 0, 0, 0, 0, t5, t99, t5, -t23, -t99, -g(1) * t95 - g(2) * t96 - g(3) * t98, t5, -t99, t23, -g(1) * (t29 * pkin(5) + t70 * t116 + t95) - g(2) * (t27 * pkin(5) + t68 * t116 + t96) - g(3) * (t49 * pkin(5) - qJ(6) * t114 + t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, t9, -t8, t9, -t11, t8, -g(1) * (t46 * pkin(10) + t111) - g(2) * (t42 * pkin(10) + t112) - g(3) * (t67 * pkin(10) + t110) t9, t8, t11, -g(1) * (-pkin(5) * t127 + t118 * t46 + t111) - g(2) * (pkin(5) * t128 + t118 * t42 + t112) - g(3) * (pkin(5) * t126 + t118 * t67 + t110); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t100, t2, 0, -t100, -g(1) * t105 - g(2) * t106 - g(3) * t104, t2, -t100, 0, -g(1) * (-t21 * pkin(5) + t105) - g(2) * (-t17 * pkin(5) + t106) - g(3) * (-t39 * pkin(5) + t104); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10;];
taug_reg  = t1;

% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP9
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:20:50
% EndTime: 2019-05-07 19:20:52
% DurationCPUTime: 0.62s
% Computational Cost: add. (674->140), mult. (1766->191), div. (0->0), fcn. (2227->10), ass. (0->74)
t90 = cos(qJ(3));
t137 = -pkin(3) * t90 - pkin(2);
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t136 = pkin(4) * t89 + qJ(5) * t86 + pkin(3);
t132 = cos(qJ(1));
t85 = sin(pkin(6));
t113 = t85 * t132;
t119 = cos(pkin(6));
t104 = t119 * t132;
t131 = sin(qJ(1));
t88 = sin(qJ(2));
t91 = cos(qJ(2));
t68 = t88 * t104 + t131 * t91;
t87 = sin(qJ(3));
t41 = -t87 * t113 + t68 * t90;
t67 = -t91 * t104 + t131 * t88;
t16 = t41 * t86 - t67 * t89;
t17 = t41 * t89 + t67 * t86;
t135 = pkin(5) + pkin(10);
t129 = t67 * t87;
t103 = t119 * t131;
t69 = t91 * t103 + t132 * t88;
t127 = t69 * t87;
t126 = t85 * t88;
t125 = t85 * t91;
t124 = t86 * t90;
t123 = t89 * t90;
t122 = t90 * t91;
t120 = qJ(6) * t89;
t118 = t87 * t125;
t117 = t86 * t125;
t111 = -t90 * t113 - t68 * t87;
t116 = t136 * t111;
t112 = t85 * t131;
t70 = -t88 * t103 + t132 * t91;
t44 = -t112 * t90 + t70 * t87;
t115 = t136 * t44;
t65 = t119 * t90 - t126 * t87;
t114 = t136 * t65;
t110 = -t16 * pkin(4) + t17 * qJ(5);
t45 = t112 * t87 + t70 * t90;
t20 = t45 * t86 - t69 * t89;
t21 = t45 * t89 + t69 * t86;
t109 = -t20 * pkin(4) + t21 * qJ(5);
t66 = t119 * t87 + t126 * t90;
t38 = t125 * t89 + t66 * t86;
t39 = t66 * t89 - t117;
t108 = -t38 * pkin(4) + t39 * qJ(5);
t107 = -g(1) * t16 + g(2) * t20;
t106 = -g(1) * t17 + g(2) * t21;
t105 = g(1) * t111 + g(2) * t44;
t2 = g(1) * t20 + g(2) * t16 + g(3) * t38;
t102 = g(1) * t21 + g(2) * t17 + g(3) * t39;
t25 = -t67 * t124 - t68 * t89;
t27 = -t124 * t69 - t70 * t89;
t47 = t117 * t90 - t126 * t89;
t101 = g(1) * t27 + g(2) * t25 + g(3) * t47;
t26 = -t67 * t123 + t68 * t86;
t28 = -t123 * t69 + t70 * t86;
t48 = (t122 * t89 + t86 * t88) * t85;
t100 = g(1) * t28 + g(2) * t26 + g(3) * t48;
t99 = g(1) * t44 - g(2) * t111 - g(3) * t65;
t98 = g(1) * t45 + g(2) * t41 + g(3) * t66;
t97 = t85 * pkin(3) * t122 + pkin(2) * t125 + t48 * pkin(4) + pkin(9) * t126 + pkin(10) * t118 + t47 * qJ(5);
t96 = -g(1) * t69 - g(2) * t67 + g(3) * t125;
t95 = t26 * pkin(4) + t68 * pkin(9) - pkin(10) * t129 + t25 * qJ(5) + t137 * t67;
t94 = t28 * pkin(4) + t70 * pkin(9) - pkin(10) * t127 + t27 * qJ(5) + t137 * t69;
t93 = t132 * pkin(1) + t70 * pkin(2) + t45 * pkin(3) + t21 * pkin(4) + pkin(8) * t112 + t69 * pkin(9) + t20 * qJ(5);
t92 = -t131 * pkin(1) - t68 * pkin(2) - pkin(3) * t41 - pkin(4) * t17 + pkin(8) * t113 - t67 * pkin(9) - qJ(5) * t16;
t22 = t96 * t87;
t9 = t99 * t89;
t8 = t99 * t86;
t1 = [0, g(1) * t131 - g(2) * t132, g(1) * t132 + g(2) * t131, 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t70, -g(1) * t67 + g(2) * t69, 0, 0, 0, 0, 0, g(1) * t41 - g(2) * t45, t105, 0, 0, 0, 0, 0, -t106, t107, -t105, t106, -t107, -g(1) * (pkin(10) * t111 + t92) - g(2) * (t44 * pkin(10) + t93) -t105, -t107, -t106, -g(1) * (-qJ(6) * t17 + t111 * t135 + t92) - g(2) * (t21 * qJ(6) + t135 * t44 + t93); 0, 0, 0, 0, 0, 0, 0, 0, -t96, g(1) * t70 + g(2) * t68 + g(3) * t126, 0, 0, 0, 0, 0, -t96 * t90, t22, 0, 0, 0, 0, 0, -t100, t101, -t22, t100, -t101, -g(1) * t94 - g(2) * t95 - g(3) * t97, -t22, -t101, -t100, -g(1) * (-pkin(5) * t127 + t28 * qJ(6) + t94) - g(2) * (-pkin(5) * t129 + t26 * qJ(6) + t95) - g(3) * (pkin(5) * t118 + t48 * qJ(6) + t97); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t98, 0, 0, 0, 0, 0, t9, -t8, -t98, -t9, t8, -g(1) * (t45 * pkin(10) - t115) - g(2) * (t41 * pkin(10) + t116) - g(3) * (t66 * pkin(10) + t114) -t98, t8, t9, -g(1) * (-t44 * t120 + t135 * t45 - t115) - g(2) * (t111 * t120 + t135 * t41 + t116) - g(3) * (t65 * t120 + t135 * t66 + t114); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t102, 0, -t2, -t102, -g(1) * t109 - g(2) * t110 - g(3) * t108, 0, -t102, t2, -g(1) * (-t20 * qJ(6) + t109) - g(2) * (-t16 * qJ(6) + t110) - g(3) * (-t38 * qJ(6) + t108); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102;];
taug_reg  = t1;

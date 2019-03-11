% Calculate inertial parameters regressor of potential energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:17
% EndTime: 2019-03-08 20:12:18
% DurationCPUTime: 0.23s
% Computational Cost: add. (346->91), mult. (579->131), div. (0->0), fcn. (697->12), ass. (0->55)
t108 = cos(pkin(10));
t105 = sin(pkin(10));
t106 = sin(pkin(6));
t132 = t105 * t106;
t134 = t108 * pkin(1) + pkin(7) * t132;
t104 = sin(pkin(11));
t109 = cos(pkin(6));
t133 = t104 * t109;
t131 = t106 * t108;
t112 = sin(qJ(2));
t130 = t106 * t112;
t114 = cos(qJ(2));
t129 = t106 * t114;
t128 = t109 * t112;
t127 = t109 * t114;
t126 = t109 * pkin(7) + qJ(1);
t125 = t104 * t132;
t100 = t105 * pkin(1);
t124 = -pkin(7) * t131 + t100;
t110 = -pkin(8) - qJ(3);
t87 = t105 * t127 + t108 * t112;
t88 = -t105 * t128 + t108 * t114;
t107 = cos(pkin(11));
t97 = pkin(3) * t107 + pkin(2);
t123 = pkin(3) * t125 - t87 * t110 + t88 * t97 + t134;
t122 = pkin(3) * t133 + t110 * t129 + t97 * t130 + t126;
t121 = g(1) * t105 - g(2) * t108;
t111 = sin(qJ(5));
t113 = cos(qJ(5));
t86 = t105 * t114 + t108 * t128;
t103 = pkin(11) + qJ(4);
t98 = sin(t103);
t99 = cos(t103);
t71 = -t98 * t131 + t86 * t99;
t85 = t105 * t112 - t108 * t127;
t63 = t111 * t71 - t85 * t113;
t73 = t98 * t132 + t88 * t99;
t65 = t111 * t73 - t87 * t113;
t80 = t109 * t98 + t99 * t130;
t74 = t80 * t111 + t113 * t129;
t120 = g(1) * t65 + g(2) * t63 + g(3) * t74;
t70 = t99 * t131 + t86 * t98;
t72 = -t99 * t132 + t88 * t98;
t79 = -t109 * t99 + t98 * t130;
t119 = g(1) * t72 + g(2) * t70 + g(3) * t79;
t118 = t73 * pkin(4) + pkin(9) * t72 + t123;
t67 = -g(1) * t87 - g(2) * t85 + g(3) * t129;
t117 = t80 * pkin(4) + pkin(9) * t79 + t122;
t116 = t100 + t86 * t97 - t85 * t110 + (-pkin(3) * t104 - pkin(7)) * t131;
t115 = t71 * pkin(4) + t70 * pkin(9) + t116;
t75 = -t111 * t129 + t80 * t113;
t66 = t111 * t87 + t113 * t73;
t64 = t111 * t85 + t113 * t71;
t61 = -g(1) * t66 - g(2) * t64 - g(3) * t75;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t108 - g(2) * t105, t121, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t88 - g(2) * t86 - g(3) * t130, -t67, -g(3) * t109 - t121 * t106, -g(1) * t134 - g(2) * t124 - g(3) * t126, 0, 0, 0, 0, 0, 0, -g(1) * (t107 * t88 + t125) - g(2) * (-t104 * t131 + t107 * t86) - g(3) * (t107 * t130 + t133) -g(1) * (-t104 * t88 + t107 * t132) - g(2) * (-t104 * t86 - t107 * t131) - g(3) * (-t104 * t130 + t107 * t109) t67, -g(1) * (pkin(2) * t88 + qJ(3) * t87 + t134) - g(2) * (pkin(2) * t86 + qJ(3) * t85 + t124) - g(3) * ((pkin(2) * t112 - qJ(3) * t114) * t106 + t126) 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - g(3) * t80, t119, t67, -g(1) * t123 - g(2) * t116 - g(3) * t122, 0, 0, 0, 0, 0, 0, t61, t120, -t119, -g(1) * t118 - g(2) * t115 - g(3) * t117, 0, 0, 0, 0, 0, 0, t61, -t119, -t120, -g(1) * (pkin(5) * t66 + qJ(6) * t65 + t118) - g(2) * (t64 * pkin(5) + t63 * qJ(6) + t115) - g(3) * (pkin(5) * t75 + qJ(6) * t74 + t117);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:53:18
% EndTime: 2019-03-09 14:53:18
% DurationCPUTime: 0.22s
% Computational Cost: add. (248->97), mult. (524->134), div. (0->0), fcn. (622->12), ass. (0->56)
t106 = cos(pkin(6));
t105 = sin(pkin(6));
t109 = sin(qJ(2));
t139 = t105 * t109;
t143 = t106 * pkin(3) + pkin(9) * t139;
t142 = t106 * pkin(8) + pkin(7);
t107 = sin(qJ(5));
t114 = cos(qJ(1));
t132 = t114 * t109;
t110 = sin(qJ(1));
t113 = cos(qJ(2));
t133 = t110 * t113;
t88 = t106 * t132 + t133;
t141 = t88 * t107;
t138 = t105 * t110;
t140 = t114 * pkin(1) + pkin(8) * t138;
t137 = t105 * t113;
t136 = t105 * t114;
t135 = t107 * t109;
t134 = t110 * t109;
t131 = t114 * t113;
t130 = pkin(2) * t139 + t142;
t129 = pkin(8) * t136;
t128 = (-pkin(3) - pkin(8)) * t114;
t102 = t110 * pkin(1);
t87 = -t106 * t131 + t134;
t127 = t88 * pkin(2) + t87 * qJ(3) + t102;
t126 = g(1) * t110 - g(2) * t114;
t89 = t106 * t133 + t132;
t90 = -t106 * t134 + t131;
t125 = t90 * pkin(2) + t89 * qJ(3) + t140;
t124 = pkin(3) * t138 + t125;
t123 = -qJ(3) * t137 + t130;
t108 = sin(qJ(4));
t112 = cos(qJ(4));
t77 = t108 * t138 - t89 * t112;
t79 = t108 * t136 + t87 * t112;
t85 = t106 * t108 + t112 * t137;
t122 = g(1) * t77 - g(2) * t79 + g(3) * t85;
t121 = t88 * pkin(9) + t127;
t120 = g(1) * t90 + g(2) * t88 + g(3) * t139;
t119 = -g(1) * t89 - g(2) * t87 + g(3) * t137;
t118 = t123 + t143;
t117 = t90 * pkin(9) + t124;
t116 = t105 * t128 + t121;
t115 = -pkin(11) - pkin(10);
t111 = cos(qJ(5));
t104 = qJ(5) + qJ(6);
t99 = cos(t104);
t98 = sin(t104);
t97 = t111 * pkin(5) + pkin(4);
t86 = t106 * t112 - t108 * t137;
t81 = -g(3) * t106 - t126 * t105;
t80 = t87 * t108 - t112 * t136;
t78 = t89 * t108 + t112 * t138;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t114 - g(2) * t110, t126, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t120, -t119, t81, -g(1) * t140 - g(2) * (t102 - t129) - g(3) * t142, 0, 0, 0, 0, 0, 0, t81, t120, t119, -g(1) * t125 - g(2) * (t127 - t129) - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t80 - g(3) * t86, t122, -t120, -g(1) * t117 - g(2) * t116 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t107 + t78 * t111) - g(2) * (t80 * t111 + t141) - g(3) * (t105 * t135 + t86 * t111) -g(1) * (-t78 * t107 + t90 * t111) - g(2) * (-t80 * t107 + t88 * t111) - g(3) * (-t86 * t107 + t111 * t139) -t122, -g(1) * (t78 * pkin(4) + t77 * pkin(10) + t117) - g(2) * (t80 * pkin(4) - t79 * pkin(10) + t116) - g(3) * (t86 * pkin(4) + t85 * pkin(10) + t118) 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t99 + t90 * t98) - g(2) * (t80 * t99 + t88 * t98) - g(3) * (t98 * t139 + t86 * t99) -g(1) * (-t78 * t98 + t90 * t99) - g(2) * (-t80 * t98 + t88 * t99) - g(3) * (t99 * t139 - t86 * t98) -t122, -g(1) * (-t77 * t115 + t78 * t97 + (pkin(5) * t107 + pkin(9)) * t90 + t124) - g(2) * (pkin(5) * t141 + t79 * t115 + t80 * t97 + t121) - g(3) * (-t85 * t115 + t86 * t97 + t130 + t143) + (-g(3) * (pkin(5) * t135 - qJ(3) * t113) - g(2) * t128) * t105;];
U_reg  = t1;

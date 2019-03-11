% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:38:28
% EndTime: 2019-03-09 21:38:28
% DurationCPUTime: 0.19s
% Computational Cost: add. (311->79), mult. (727->108), div. (0->0), fcn. (903->10), ass. (0->55)
t110 = sin(qJ(3));
t115 = cos(qJ(1));
t108 = sin(pkin(6));
t140 = cos(qJ(3));
t131 = t108 * t140;
t111 = sin(qJ(2));
t112 = sin(qJ(1));
t114 = cos(qJ(2));
t137 = cos(pkin(6));
t129 = t115 * t137;
t97 = t111 * t129 + t112 * t114;
t85 = t97 * t110 + t115 * t131;
t143 = t85 * pkin(10);
t130 = t112 * t137;
t99 = -t111 * t130 + t115 * t114;
t87 = t99 * t110 - t112 * t131;
t142 = t87 * pkin(10);
t136 = t108 * t111;
t94 = t110 * t136 - t137 * t140;
t141 = t94 * pkin(10);
t139 = pkin(10) - qJ(6);
t138 = t137 * pkin(8) + pkin(7);
t135 = t108 * t112;
t134 = t108 * t114;
t133 = t108 * t115;
t132 = t115 * pkin(1) + pkin(8) * t135;
t128 = t112 * pkin(1) - pkin(8) * t133;
t127 = g(1) * t112 - g(2) * t115;
t98 = t115 * t111 + t114 * t130;
t126 = t99 * pkin(2) + t98 * pkin(9) + t132;
t88 = t110 * t135 + t99 * t140;
t125 = t88 * pkin(3) + t126;
t124 = pkin(2) * t136 - pkin(9) * t134 + t138;
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t86 = -t110 * t133 + t97 * t140;
t96 = t112 * t111 - t114 * t129;
t76 = t86 * t109 - t96 * t113;
t78 = t88 * t109 - t98 * t113;
t95 = t137 * t110 + t111 * t131;
t83 = t95 * t109 + t113 * t134;
t123 = g(1) * t78 + g(2) * t76 + g(3) * t83;
t73 = g(1) * t87 + g(2) * t85 + g(3) * t94;
t122 = t95 * pkin(3) + t124;
t121 = t97 * pkin(2) + t96 * pkin(9) + t128;
t120 = -g(1) * t98 - g(2) * t96 + g(3) * t134;
t119 = t86 * pkin(3) + t121;
t79 = t98 * t109 + t88 * t113;
t118 = t79 * pkin(4) + t78 * qJ(5) + t125;
t84 = -t109 * t134 + t95 * t113;
t117 = t84 * pkin(4) + t83 * qJ(5) + t122;
t77 = t96 * t109 + t86 * t113;
t116 = t77 * pkin(4) + t76 * qJ(5) + t119;
t71 = -g(1) * t79 - g(2) * t77 - g(3) * t84;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t115 - g(2) * t112, t127, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t99 - g(2) * t97 - g(3) * t136, -t120, -g(3) * t137 - t127 * t108, -g(1) * t132 - g(2) * t128 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * t88 - g(2) * t86 - g(3) * t95, t73, t120, -g(1) * t126 - g(2) * t121 - g(3) * t124, 0, 0, 0, 0, 0, 0, t71, t123, -t73, -g(1) * (t125 + t142) - g(2) * (t119 + t143) - g(3) * (t122 + t141) 0, 0, 0, 0, 0, 0, t71, -t73, -t123, -g(1) * (t118 + t142) - g(2) * (t116 + t143) - g(3) * (t117 + t141) 0, 0, 0, 0, 0, 0, t71, -t123, t73, -g(1) * (t79 * pkin(5) + t139 * t87 + t118) - g(2) * (t77 * pkin(5) + t139 * t85 + t116) - g(3) * (t84 * pkin(5) + t139 * t94 + t117);];
U_reg  = t1;

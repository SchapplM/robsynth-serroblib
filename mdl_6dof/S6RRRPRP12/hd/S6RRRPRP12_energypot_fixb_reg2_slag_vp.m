% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:59:16
% EndTime: 2019-03-09 17:59:16
% DurationCPUTime: 0.18s
% Computational Cost: add. (277->78), mult. (638->108), div. (0->0), fcn. (780->10), ass. (0->53)
t144 = cos(qJ(3));
t142 = cos(pkin(6));
t143 = pkin(8) * t142 + pkin(7);
t111 = sin(pkin(6));
t114 = sin(qJ(2));
t141 = t111 * t114;
t115 = sin(qJ(1));
t140 = t111 * t115;
t117 = cos(qJ(2));
t139 = t111 * t117;
t118 = cos(qJ(1));
t138 = t111 * t118;
t137 = pkin(1) * t118 + pkin(8) * t140;
t136 = pkin(9) * t139;
t135 = pkin(2) * t141 + t143;
t134 = t111 * t144;
t133 = t115 * t142;
t132 = t118 * t142;
t131 = pkin(1) * t115 - pkin(8) * t138;
t130 = g(1) * t115 - g(2) * t118;
t100 = t114 * t118 + t117 * t133;
t101 = -t114 * t133 + t117 * t118;
t129 = pkin(2) * t101 + pkin(9) * t100 + t137;
t113 = sin(qJ(3));
t96 = t113 * t141 - t142 * t144;
t97 = t113 * t142 + t114 * t134;
t128 = pkin(3) * t97 + t96 * qJ(4) + t135;
t112 = sin(qJ(5));
t116 = cos(qJ(5));
t99 = t114 * t132 + t115 * t117;
t87 = t113 * t99 + t118 * t134;
t98 = t114 * t115 - t117 * t132;
t76 = t112 * t98 - t116 * t87;
t89 = t101 * t113 - t115 * t134;
t78 = t100 * t112 - t116 * t89;
t85 = t112 * t139 + t116 * t96;
t127 = g(1) * t78 + g(2) * t76 - g(3) * t85;
t126 = g(1) * t89 + g(2) * t87 + g(3) * t96;
t88 = -t113 * t138 + t144 * t99;
t90 = t101 * t144 + t113 * t140;
t125 = g(1) * t90 + g(2) * t88 + g(3) * t97;
t124 = pkin(2) * t99 + pkin(9) * t98 + t131;
t80 = -g(1) * t100 - g(2) * t98 + g(3) * t139;
t123 = pkin(3) * t90 + qJ(4) * t89 + t129;
t122 = pkin(3) * t88 + qJ(4) * t87 + t124;
t121 = pkin(4) * t100 + pkin(10) * t90 + t123;
t120 = t97 * pkin(10) + (-pkin(4) - pkin(9)) * t139 + t128;
t119 = pkin(4) * t98 + pkin(10) * t88 + t122;
t86 = t112 * t96 - t116 * t139;
t79 = t100 * t116 + t112 * t89;
t77 = t112 * t87 + t116 * t98;
t74 = -g(1) * t79 - g(2) * t77 - g(3) * t86;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t118 - g(2) * t115, t130, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t99 - g(3) * t141, -t80, -g(3) * t142 - t111 * t130, -g(1) * t137 - g(2) * t131 - g(3) * t143, 0, 0, 0, 0, 0, 0, -t125, t126, t80, -g(1) * t129 - g(2) * t124 - g(3) * (t135 - t136) 0, 0, 0, 0, 0, 0, t80, t125, -t126, -g(1) * t123 - g(2) * t122 - g(3) * (t128 - t136) 0, 0, 0, 0, 0, 0, t74, t127, -t125, -g(1) * t121 - g(2) * t119 - g(3) * t120, 0, 0, 0, 0, 0, 0, t74, -t125, -t127, -g(1) * (pkin(5) * t79 + qJ(6) * t78 + t121) - g(2) * (pkin(5) * t77 + qJ(6) * t76 + t119) - g(3) * (t86 * pkin(5) - t85 * qJ(6) + t120);];
U_reg  = t1;

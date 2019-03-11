% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:24
% EndTime: 2019-03-09 12:05:24
% DurationCPUTime: 0.25s
% Computational Cost: add. (364->90), mult. (841->136), div. (0->0), fcn. (1059->12), ass. (0->57)
t110 = sin(pkin(11));
t116 = sin(qJ(2));
t120 = cos(qJ(2));
t144 = cos(pkin(11));
t99 = -t116 * t110 + t120 * t144;
t112 = cos(pkin(6));
t146 = t112 * pkin(8) + pkin(7);
t107 = t120 * pkin(2) + pkin(1);
t117 = sin(qJ(1));
t121 = cos(qJ(1));
t111 = sin(pkin(6));
t97 = t112 * t116 * pkin(2) + (-pkin(8) - qJ(3)) * t111;
t145 = t117 * t107 + t121 * t97;
t143 = t111 * t116;
t142 = t111 * t117;
t141 = t111 * t121;
t139 = t117 * t116;
t138 = t117 * t120;
t137 = t121 * t116;
t136 = t121 * t120;
t98 = -t120 * t110 - t116 * t144;
t96 = t98 * t112;
t85 = t117 * t99 - t121 * t96;
t135 = t85 * pkin(3) + t145;
t114 = sin(qJ(5));
t134 = pkin(5) * t114 + pkin(9);
t132 = t121 * t107 - t117 * t97;
t131 = pkin(2) * t143 + t112 * qJ(3) + t146;
t87 = t117 * t96 + t121 * t99;
t130 = t87 * pkin(3) + t132;
t95 = t98 * t111;
t129 = -t95 * pkin(3) + t131;
t128 = g(1) * t117 - g(2) * t121;
t122 = t99 * t112;
t84 = t117 * t98 + t121 * t122;
t127 = -t84 * pkin(9) + t135;
t115 = sin(qJ(4));
t119 = cos(qJ(4));
t78 = t85 * t115 + t119 * t141;
t80 = t87 * t115 - t119 * t142;
t88 = -t112 * t119 - t95 * t115;
t126 = g(1) * t80 + g(2) * t78 + g(3) * t88;
t86 = -t117 * t122 + t121 * t98;
t94 = t99 * t111;
t125 = g(1) * t86 + g(2) * t84 + g(3) * t94;
t124 = -t86 * pkin(9) + t130;
t123 = -t94 * pkin(9) + t129;
t118 = cos(qJ(5));
t113 = -qJ(6) - pkin(10);
t106 = t118 * pkin(5) + pkin(4);
t93 = -g(3) * t112 - t128 * t111;
t89 = t112 * t115 - t95 * t119;
t81 = t115 * t142 + t87 * t119;
t79 = -t115 * t141 + t85 * t119;
t76 = -g(1) * (-t86 * t114 + t81 * t118) - g(2) * (-t84 * t114 + t79 * t118) - g(3) * (-t94 * t114 + t89 * t118);
t75 = -g(1) * (-t81 * t114 - t86 * t118) - g(2) * (-t79 * t114 - t84 * t118) - g(3) * (-t89 * t114 - t94 * t118);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t117, t128, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t112 * t139 + t136) - g(2) * (t112 * t137 + t138) - g(3) * t143, -g(1) * (-t112 * t138 - t137) - g(2) * (t112 * t136 - t139) - g(3) * t111 * t120, t93, -g(1) * (t121 * pkin(1) + pkin(8) * t142) - g(2) * (t117 * pkin(1) - pkin(8) * t141) - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t85 + g(3) * t95, -t125, t93, -g(1) * t132 - g(2) * t145 - g(3) * t131, 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t89, t126, t125, -g(1) * t124 - g(2) * t127 - g(3) * t123, 0, 0, 0, 0, 0, 0, t76, t75, -t126, -g(1) * (t81 * pkin(4) + t80 * pkin(10) + t124) - g(2) * (t79 * pkin(4) + t78 * pkin(10) + t127) - g(3) * (t89 * pkin(4) + t88 * pkin(10) + t123) 0, 0, 0, 0, 0, 0, t76, t75, -t126, -g(1) * (t81 * t106 - t80 * t113 - t134 * t86 + t130) - g(2) * (t79 * t106 - t78 * t113 - t134 * t84 + t135) - g(3) * (t89 * t106 - t88 * t113 - t134 * t94 + t129);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:01
% EndTime: 2019-03-09 09:05:02
% DurationCPUTime: 0.23s
% Computational Cost: add. (312->85), mult. (712->131), div. (0->0), fcn. (881->12), ass. (0->53)
t110 = sin(pkin(11));
t112 = cos(pkin(11));
t116 = sin(qJ(2));
t120 = cos(qJ(2));
t100 = -t116 * t110 + t120 * t112;
t113 = cos(pkin(6));
t147 = t113 * pkin(8) + pkin(7);
t106 = t120 * pkin(2) + pkin(1);
t117 = sin(qJ(1));
t121 = cos(qJ(1));
t111 = sin(pkin(6));
t98 = t113 * t116 * pkin(2) + (-pkin(8) - qJ(3)) * t111;
t146 = t117 * t106 + t121 * t98;
t145 = t111 * t116;
t144 = t111 * t117;
t143 = t111 * t121;
t141 = t117 * t116;
t140 = t117 * t120;
t138 = t121 * t116;
t137 = t121 * t120;
t136 = t121 * t106 - t117 * t98;
t135 = pkin(2) * t145 + t113 * qJ(3) + t147;
t134 = g(1) * t117 - g(2) * t121;
t133 = t120 * t110 + t116 * t112;
t127 = t100 * t113;
t84 = -t117 * t133 + t121 * t127;
t128 = t133 * t113;
t85 = t117 * t100 + t121 * t128;
t132 = t85 * pkin(3) - t84 * qJ(4) + t146;
t115 = sin(qJ(5));
t119 = cos(qJ(5));
t86 = -t117 * t127 - t121 * t133;
t77 = t115 * t144 + t86 * t119;
t79 = t115 * t143 - t84 * t119;
t96 = t100 * t111;
t88 = t113 * t115 + t96 * t119;
t131 = g(1) * t77 - g(2) * t79 + g(3) * t88;
t130 = g(1) * t86 + g(2) * t84 + g(3) * t96;
t87 = t121 * t100 - t117 * t128;
t97 = t133 * t111;
t129 = g(1) * t87 + g(2) * t85 + g(3) * t97;
t126 = t87 * pkin(3) - t86 * qJ(4) + t136;
t125 = t97 * pkin(3) - t96 * qJ(4) + t135;
t124 = pkin(4) * t144 + t87 * pkin(9) + t126;
t123 = t113 * pkin(4) + t97 * pkin(9) + t125;
t122 = -pkin(4) * t143 + t85 * pkin(9) + t132;
t118 = cos(qJ(6));
t114 = sin(qJ(6));
t95 = -g(3) * t113 - t134 * t111;
t89 = t113 * t119 - t96 * t115;
t80 = -t84 * t115 - t119 * t143;
t78 = -t86 * t115 + t119 * t144;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t117, t134, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t113 * t141 + t137) - g(2) * (t113 * t138 + t140) - g(3) * t145, -g(1) * (-t113 * t140 - t138) - g(2) * (t113 * t137 - t141) - g(3) * t111 * t120, t95, -g(1) * (t121 * pkin(1) + pkin(8) * t144) - g(2) * (t117 * pkin(1) - pkin(8) * t143) - g(3) * t147, 0, 0, 0, 0, 0, 0, -t129, -t130, t95, -g(1) * t136 - g(2) * t146 - g(3) * t135, 0, 0, 0, 0, 0, 0, t95, t129, t130, -g(1) * t126 - g(2) * t132 - g(3) * t125, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t80 - g(3) * t89, t131, -t129, -g(1) * t124 - g(2) * t122 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t87 * t114 + t78 * t118) - g(2) * (t85 * t114 + t80 * t118) - g(3) * (t97 * t114 + t89 * t118) -g(1) * (-t78 * t114 + t87 * t118) - g(2) * (-t80 * t114 + t85 * t118) - g(3) * (-t89 * t114 + t97 * t118) -t131, -g(1) * (t78 * pkin(5) + t77 * pkin(10) + t124) - g(2) * (t80 * pkin(5) - t79 * pkin(10) + t122) - g(3) * (t89 * pkin(5) + t88 * pkin(10) + t123);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:27:52
% EndTime: 2019-03-09 16:27:52
% DurationCPUTime: 0.23s
% Computational Cost: add. (277->96), mult. (600->131), div. (0->0), fcn. (727->12), ass. (0->54)
t145 = pkin(4) + pkin(9);
t112 = sin(qJ(2));
t113 = sin(qJ(1));
t114 = cos(qJ(2));
t115 = cos(qJ(1));
t138 = cos(pkin(6));
t127 = t115 * t138;
t90 = t113 * t112 - t114 * t127;
t144 = t90 * pkin(9);
t128 = t113 * t138;
t92 = t115 * t112 + t114 * t128;
t143 = t92 * pkin(9);
t142 = cos(qJ(3));
t109 = cos(pkin(11));
t141 = t109 * pkin(5) + t145;
t140 = t138 * pkin(8) + pkin(7);
t108 = sin(pkin(6));
t136 = t108 * t113;
t139 = t115 * pkin(1) + pkin(8) * t136;
t137 = t108 * t112;
t135 = t108 * t114;
t134 = t108 * t115;
t133 = pkin(2) * t137 + t140;
t132 = pkin(9) * t135;
t93 = -t112 * t128 + t115 * t114;
t131 = t93 * pkin(2) + t139;
t130 = t108 * t142;
t107 = sin(pkin(11));
t129 = pkin(5) * t107 + qJ(4);
t111 = sin(qJ(3));
t89 = t138 * t111 + t112 * t130;
t126 = t89 * pkin(3) + t133;
t84 = t111 * t136 + t93 * t142;
t125 = t84 * pkin(3) + t131;
t124 = t113 * pkin(1) - pkin(8) * t134;
t123 = g(1) * t113 - g(2) * t115;
t91 = t112 * t127 + t113 * t114;
t122 = t91 * pkin(2) + t124;
t88 = t111 * t137 - t138 * t142;
t121 = t88 * qJ(4) + t126;
t82 = -t111 * t134 + t91 * t142;
t120 = t82 * pkin(3) + t122;
t83 = t93 * t111 - t113 * t130;
t119 = t83 * qJ(4) + t125;
t81 = t91 * t111 + t115 * t130;
t118 = g(1) * t83 + g(2) * t81 + g(3) * t88;
t117 = g(1) * t84 + g(2) * t82 + g(3) * t89;
t78 = -g(1) * t92 - g(2) * t90 + g(3) * t135;
t116 = t81 * qJ(4) + t120;
t110 = -pkin(10) - qJ(5);
t106 = pkin(11) + qJ(6);
t102 = cos(t106);
t101 = sin(t106);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t115 - g(2) * t113, t123, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t137, -t78, -g(3) * t138 - t123 * t108, -g(1) * t139 - g(2) * t124 - g(3) * t140, 0, 0, 0, 0, 0, 0, -t117, t118, t78, -g(1) * (t131 + t143) - g(2) * (t122 + t144) - g(3) * (-t132 + t133) 0, 0, 0, 0, 0, 0, t78, t117, -t118, -g(1) * (t119 + t143) - g(2) * (t116 + t144) - g(3) * (t121 - t132) 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t107 + t92 * t109) - g(2) * (t81 * t107 + t90 * t109) - g(3) * (t88 * t107 - t109 * t135) -g(1) * (-t92 * t107 + t83 * t109) - g(2) * (-t90 * t107 + t81 * t109) - g(3) * (t107 * t135 + t88 * t109) -t117, -g(1) * (t84 * qJ(5) + t145 * t92 + t119) - g(2) * (t82 * qJ(5) + t145 * t90 + t116) - g(3) * (t89 * qJ(5) - t145 * t135 + t121) 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t101 + t92 * t102) - g(2) * (t81 * t101 + t90 * t102) - g(3) * (t88 * t101 - t102 * t135) -g(1) * (-t92 * t101 + t83 * t102) - g(2) * (-t90 * t101 + t81 * t102) - g(3) * (t101 * t135 + t88 * t102) -t117, -g(1) * (-t84 * t110 + t129 * t83 + t141 * t92 + t125) - g(2) * (-t82 * t110 + t129 * t81 + t141 * t90 + t120) - g(3) * (-t89 * t110 + t129 * t88 - t141 * t135 + t126);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:29
% EndTime: 2019-03-09 16:09:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (254->84), mult. (578->108), div. (0->0), fcn. (697->10), ass. (0->51)
t103 = sin(pkin(6));
t109 = cos(qJ(2));
t128 = t103 * t109;
t105 = sin(qJ(3));
t106 = sin(qJ(2));
t138 = cos(qJ(3));
t125 = t103 * t138;
t131 = cos(pkin(6));
t88 = t131 * t105 + t106 * t125;
t139 = t88 * pkin(4) + qJ(5) * t128;
t137 = pkin(5) + qJ(4);
t136 = pkin(9) - qJ(5);
t135 = t131 * pkin(8) + pkin(7);
t110 = cos(qJ(1));
t107 = sin(qJ(1));
t123 = t110 * t131;
t90 = t106 * t123 + t107 * t109;
t79 = t90 * t105 + t110 * t125;
t134 = t79 * qJ(4);
t124 = t107 * t131;
t92 = -t106 * t124 + t110 * t109;
t81 = t92 * t105 - t107 * t125;
t133 = t81 * qJ(4);
t129 = t103 * t107;
t132 = t110 * pkin(1) + pkin(8) * t129;
t130 = t103 * t106;
t127 = t103 * t110;
t126 = t92 * pkin(2) + t132;
t122 = t107 * pkin(1) - pkin(8) * t127;
t121 = g(1) * t107 - g(2) * t110;
t91 = t110 * t106 + t109 * t124;
t120 = t91 * pkin(9) + t126;
t119 = t90 * pkin(2) + t122;
t118 = pkin(2) * t130 - pkin(9) * t128 + t135;
t87 = t105 * t130 - t131 * t138;
t117 = g(1) * t81 + g(2) * t79 + g(3) * t87;
t80 = -t105 * t127 + t90 * t138;
t82 = t105 * t129 + t92 * t138;
t116 = g(1) * t82 + g(2) * t80 + g(3) * t88;
t115 = t88 * pkin(3) + t118;
t89 = t107 * t106 - t109 * t123;
t114 = t89 * pkin(9) + t119;
t73 = -g(1) * t91 - g(2) * t89 + g(3) * t128;
t78 = t82 * pkin(3);
t113 = t82 * pkin(4) + t136 * t91 + t126 + t78;
t112 = t87 * qJ(4) + t115;
t76 = t80 * pkin(3);
t111 = t80 * pkin(4) + t136 * t89 + t119 + t76;
t108 = cos(qJ(6));
t104 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t110 - g(2) * t107, t121, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t92 - g(2) * t90 - g(3) * t130, -t73, -g(3) * t131 - t121 * t103, -g(1) * t132 - g(2) * t122 - g(3) * t135, 0, 0, 0, 0, 0, 0, -t116, t117, t73, -g(1) * t120 - g(2) * t114 - g(3) * t118, 0, 0, 0, 0, 0, 0, -t116, t73, -t117, -g(1) * (t120 + t78 + t133) - g(2) * (t114 + t76 + t134) - g(3) * t112, 0, 0, 0, 0, 0, 0, -t117, t116, -t73, -g(1) * (t113 + t133) - g(2) * (t111 + t134) - g(3) * (t112 + t139) 0, 0, 0, 0, 0, 0, -g(1) * (-t91 * t104 + t81 * t108) - g(2) * (-t89 * t104 + t79 * t108) - g(3) * (t104 * t128 + t87 * t108) -g(1) * (-t81 * t104 - t91 * t108) - g(2) * (-t79 * t104 - t89 * t108) - g(3) * (-t87 * t104 + t108 * t128) -t116, -g(1) * (t82 * pkin(10) + t137 * t81 + t113) - g(2) * (t80 * pkin(10) + t137 * t79 + t111) - g(3) * (t88 * pkin(10) + t137 * t87 + t115 + t139);];
U_reg  = t1;

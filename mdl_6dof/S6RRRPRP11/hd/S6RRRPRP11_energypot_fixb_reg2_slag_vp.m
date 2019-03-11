% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP11
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:02
% EndTime: 2019-03-09 17:49:02
% DurationCPUTime: 0.20s
% Computational Cost: add. (265->85), mult. (600->113), div. (0->0), fcn. (727->10), ass. (0->53)
t142 = pkin(4) + pkin(9);
t108 = sin(qJ(2));
t109 = sin(qJ(1));
t111 = cos(qJ(2));
t112 = cos(qJ(1));
t135 = cos(pkin(6));
t124 = t112 * t135;
t90 = t109 * t108 - t111 * t124;
t141 = t90 * pkin(9);
t125 = t109 * t135;
t92 = t112 * t108 + t111 * t125;
t140 = t92 * pkin(9);
t139 = cos(qJ(3));
t110 = cos(qJ(5));
t138 = t110 * pkin(5) + t142;
t137 = t135 * pkin(8) + pkin(7);
t104 = sin(pkin(6));
t133 = t104 * t109;
t136 = t112 * pkin(1) + pkin(8) * t133;
t134 = t104 * t108;
t132 = t104 * t111;
t131 = t104 * t112;
t130 = pkin(2) * t134 + t137;
t129 = pkin(9) * t132;
t93 = -t108 * t125 + t112 * t111;
t128 = t93 * pkin(2) + t136;
t127 = t104 * t139;
t106 = sin(qJ(5));
t126 = pkin(5) * t106 + qJ(4);
t107 = sin(qJ(3));
t89 = t135 * t107 + t108 * t127;
t123 = t89 * pkin(3) + t130;
t84 = t107 * t133 + t93 * t139;
t122 = t84 * pkin(3) + t128;
t121 = t109 * pkin(1) - pkin(8) * t131;
t120 = g(1) * t109 - g(2) * t112;
t91 = t108 * t124 + t109 * t111;
t119 = t91 * pkin(2) + t121;
t88 = t107 * t134 - t135 * t139;
t118 = t88 * qJ(4) + t123;
t82 = -t107 * t131 + t91 * t139;
t117 = t82 * pkin(3) + t119;
t83 = t93 * t107 - t109 * t127;
t116 = t83 * qJ(4) + t122;
t81 = t91 * t107 + t112 * t127;
t115 = g(1) * t83 + g(2) * t81 + g(3) * t88;
t114 = g(1) * t84 + g(2) * t82 + g(3) * t89;
t78 = -g(1) * t92 - g(2) * t90 + g(3) * t132;
t113 = t81 * qJ(4) + t117;
t105 = -qJ(6) - pkin(10);
t76 = -g(1) * (t83 * t106 + t92 * t110) - g(2) * (t81 * t106 + t90 * t110) - g(3) * (t88 * t106 - t110 * t132);
t75 = -g(1) * (-t92 * t106 + t83 * t110) - g(2) * (-t90 * t106 + t81 * t110) - g(3) * (t106 * t132 + t88 * t110);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t112 - g(2) * t109, t120, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t134, -t78, -g(3) * t135 - t120 * t104, -g(1) * t136 - g(2) * t121 - g(3) * t137, 0, 0, 0, 0, 0, 0, -t114, t115, t78, -g(1) * (t128 + t140) - g(2) * (t119 + t141) - g(3) * (-t129 + t130) 0, 0, 0, 0, 0, 0, t78, t114, -t115, -g(1) * (t116 + t140) - g(2) * (t113 + t141) - g(3) * (t118 - t129) 0, 0, 0, 0, 0, 0, t76, t75, -t114, -g(1) * (t84 * pkin(10) + t142 * t92 + t116) - g(2) * (t82 * pkin(10) + t142 * t90 + t113) - g(3) * (t89 * pkin(10) - t142 * t132 + t118) 0, 0, 0, 0, 0, 0, t76, t75, -t114, -g(1) * (-t84 * t105 + t126 * t83 + t138 * t92 + t122) - g(2) * (-t82 * t105 + t126 * t81 + t138 * t90 + t117) - g(3) * (-t89 * t105 + t126 * t88 - t138 * t132 + t123);];
U_reg  = t1;

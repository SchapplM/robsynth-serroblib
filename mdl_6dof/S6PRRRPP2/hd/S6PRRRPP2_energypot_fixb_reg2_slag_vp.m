% Calculate inertial parameters regressor of potential energy for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:39
% EndTime: 2019-03-08 22:52:39
% DurationCPUTime: 0.18s
% Computational Cost: add. (311->79), mult. (727->108), div. (0->0), fcn. (903->10), ass. (0->55)
t98 = sin(pkin(6));
t132 = pkin(7) * t98;
t101 = sin(qJ(3));
t128 = cos(qJ(3));
t120 = t98 * t128;
t104 = cos(qJ(2));
t102 = sin(qJ(2));
t121 = cos(pkin(6));
t119 = t102 * t121;
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t84 = t97 * t104 + t99 * t119;
t72 = t84 * t101 + t99 * t120;
t131 = pkin(9) * t72;
t86 = t99 * t104 - t97 * t119;
t74 = t101 * t86 - t97 * t120;
t130 = t74 * pkin(9);
t124 = t102 * t98;
t87 = t101 * t124 - t121 * t128;
t129 = t87 * pkin(9);
t127 = pkin(9) - qJ(6);
t126 = t99 * pkin(1) + t97 * t132;
t125 = t101 * t98;
t123 = t104 * t98;
t122 = t121 * pkin(7) + qJ(1);
t118 = t104 * t121;
t117 = t97 * pkin(1) - t99 * t132;
t116 = g(1) * t97 - g(2) * t99;
t85 = t99 * t102 + t97 * t118;
t115 = t86 * pkin(2) + t85 * pkin(8) + t126;
t75 = t97 * t125 + t86 * t128;
t114 = t75 * pkin(3) + t115;
t113 = pkin(2) * t124 - pkin(8) * t123 + t122;
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t73 = -t99 * t125 + t84 * t128;
t83 = t102 * t97 - t99 * t118;
t65 = t100 * t73 - t83 * t103;
t67 = t100 * t75 - t85 * t103;
t88 = t121 * t101 + t102 * t120;
t76 = t88 * t100 + t103 * t123;
t112 = g(1) * t67 + g(2) * t65 + g(3) * t76;
t62 = g(1) * t74 + g(2) * t72 + g(3) * t87;
t111 = t88 * pkin(3) + t113;
t110 = t84 * pkin(2) + pkin(8) * t83 + t117;
t109 = -g(1) * t85 - g(2) * t83 + g(3) * t123;
t108 = t73 * pkin(3) + t110;
t68 = t100 * t85 + t103 * t75;
t107 = t68 * pkin(4) + qJ(5) * t67 + t114;
t77 = -t100 * t123 + t88 * t103;
t106 = t77 * pkin(4) + t76 * qJ(5) + t111;
t66 = t100 * t83 + t103 * t73;
t105 = t66 * pkin(4) + qJ(5) * t65 + t108;
t60 = -g(1) * t68 - g(2) * t66 - g(3) * t77;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t99 - g(2) * t97, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t86 - g(2) * t84 - g(3) * t124, -t109, -g(3) * t121 - t116 * t98, -g(1) * t126 - g(2) * t117 - g(3) * t122, 0, 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t73 - g(3) * t88, t62, t109, -g(1) * t115 - g(2) * t110 - g(3) * t113, 0, 0, 0, 0, 0, 0, t60, t112, -t62, -g(1) * (t114 + t130) - g(2) * (t108 + t131) - g(3) * (t111 + t129) 0, 0, 0, 0, 0, 0, t60, -t62, -t112, -g(1) * (t107 + t130) - g(2) * (t105 + t131) - g(3) * (t106 + t129) 0, 0, 0, 0, 0, 0, t60, -t112, t62, -g(1) * (pkin(5) * t68 + t127 * t74 + t107) - g(2) * (pkin(5) * t66 + t127 * t72 + t105) - g(3) * (t77 * pkin(5) + t127 * t87 + t106);];
U_reg  = t1;

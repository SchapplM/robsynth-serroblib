% Calculate inertial parameters regressor of potential energy for
% S6PRRRPP3
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:28
% EndTime: 2019-03-08 22:58:29
% DurationCPUTime: 0.18s
% Computational Cost: add. (311->82), mult. (727->108), div. (0->0), fcn. (903->10), ass. (0->55)
t132 = pkin(5) + pkin(9);
t96 = sin(pkin(6));
t131 = pkin(7) * t96;
t127 = cos(qJ(3));
t120 = t96 * t127;
t102 = cos(qJ(2));
t100 = sin(qJ(2));
t121 = cos(pkin(6));
t119 = t100 * t121;
t95 = sin(pkin(10));
t97 = cos(pkin(10));
t81 = t95 * t102 + t97 * t119;
t99 = sin(qJ(3));
t69 = t97 * t120 + t81 * t99;
t130 = t69 * pkin(9);
t83 = t97 * t102 - t95 * t119;
t71 = -t95 * t120 + t83 * t99;
t129 = t71 * pkin(9);
t124 = t100 * t96;
t84 = -t121 * t127 + t99 * t124;
t128 = t84 * pkin(9);
t126 = t96 * t99;
t125 = t97 * pkin(1) + t95 * t131;
t123 = t102 * t96;
t122 = t121 * pkin(7) + qJ(1);
t118 = t102 * t121;
t117 = t95 * pkin(1) - t97 * t131;
t116 = g(1) * t95 - g(2) * t97;
t82 = t97 * t100 + t95 * t118;
t115 = t83 * pkin(2) + t82 * pkin(8) + t125;
t72 = t95 * t126 + t83 * t127;
t114 = t72 * pkin(3) + t115;
t113 = pkin(2) * t124 - pkin(8) * t123 + t122;
t101 = cos(qJ(4));
t70 = -t97 * t126 + t81 * t127;
t80 = t95 * t100 - t97 * t118;
t98 = sin(qJ(4));
t62 = -t80 * t101 + t70 * t98;
t64 = -t82 * t101 + t72 * t98;
t85 = t100 * t120 + t121 * t99;
t73 = t101 * t123 + t85 * t98;
t112 = g(1) * t64 + g(2) * t62 + g(3) * t73;
t63 = t70 * t101 + t80 * t98;
t65 = t72 * t101 + t82 * t98;
t74 = t85 * t101 - t98 * t123;
t111 = g(1) * t65 + g(2) * t63 + g(3) * t74;
t110 = g(1) * t71 + g(2) * t69 + g(3) * t84;
t109 = t85 * pkin(3) + t113;
t108 = t81 * pkin(2) + t80 * pkin(8) + t117;
t107 = -g(1) * t82 - g(2) * t80 + g(3) * t123;
t106 = t70 * pkin(3) + t108;
t105 = t65 * pkin(4) + t64 * qJ(5) + t114;
t104 = t74 * pkin(4) + t73 * qJ(5) + t109;
t103 = t63 * pkin(4) + t62 * qJ(5) + t106;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t97 - g(2) * t95, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t124, -t107, -g(3) * t121 - t116 * t96, -g(1) * t125 - g(2) * t117 - g(3) * t122, 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t85, t110, t107, -g(1) * t115 - g(2) * t108 - g(3) * t113, 0, 0, 0, 0, 0, 0, -t111, t112, -t110, -g(1) * (t114 + t129) - g(2) * (t106 + t130) - g(3) * (t109 + t128) 0, 0, 0, 0, 0, 0, -t110, t111, -t112, -g(1) * (t105 + t129) - g(2) * (t103 + t130) - g(3) * (t104 + t128) 0, 0, 0, 0, 0, 0, -t110, -t112, -t111, -g(1) * (t65 * qJ(6) + t132 * t71 + t105) - g(2) * (t63 * qJ(6) + t132 * t69 + t103) - g(3) * (t74 * qJ(6) + t132 * t84 + t104);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% U_reg [1x(7*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S7RRRRRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 07:42:51
% EndTime: 2019-03-10 07:42:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (267->71), mult. (674->117), div. (0->0), fcn. (857->14), ass. (0->59)
t127 = g(3) * pkin(1);
t110 = cos(qJ(2));
t126 = t110 * pkin(2) + pkin(1);
t102 = sin(qJ(3));
t103 = sin(qJ(2));
t125 = t102 * t103;
t104 = sin(qJ(1));
t124 = t103 * t104;
t109 = cos(qJ(3));
t123 = t103 * t109;
t111 = cos(qJ(1));
t122 = t103 * t111;
t121 = t104 * t110;
t120 = t111 * t102;
t119 = t111 * t109;
t101 = sin(qJ(4));
t108 = cos(qJ(4));
t88 = t101 * t123 + t110 * t108;
t118 = t88 * pkin(3) + t126;
t117 = g(1) * t111 + g(2) * t104;
t91 = t109 * t121 + t120;
t82 = t91 * t101 - t108 * t124;
t116 = -pkin(2) * t124 + t82 * pkin(3);
t93 = -t104 * t102 + t110 * t119;
t84 = t93 * t101 - t108 * t122;
t115 = -pkin(2) * t122 + t84 * pkin(3);
t114 = t117 * t103;
t100 = sin(qJ(5));
t107 = cos(qJ(5));
t83 = t101 * t124 + t91 * t108;
t90 = -t102 * t121 + t119;
t75 = t83 * t100 - t90 * t107;
t85 = t101 * t122 + t93 * t108;
t92 = -t104 * t109 - t110 * t120;
t77 = t85 * t100 - t92 * t107;
t89 = -t110 * t101 + t108 * t123;
t80 = t89 * t100 + t107 * t125;
t113 = g(1) * t77 + g(2) * t75 + g(3) * t80;
t112 = g(1) * t84 + g(2) * t82 + g(3) * t88;
t106 = cos(qJ(6));
t105 = cos(qJ(7));
t99 = sin(qJ(6));
t98 = sin(qJ(7));
t94 = g(1) * t104 - g(2) * t111;
t87 = -g(3) * t110 + t114;
t86 = pkin(2) * t114 - g(3) * t126;
t81 = -t100 * t125 + t89 * t107;
t79 = -g(1) * t92 - g(2) * t90 + g(3) * t125;
t78 = t92 * t100 + t85 * t107;
t76 = t90 * t100 + t83 * t107;
t74 = t81 * t106 + t88 * t99;
t73 = t88 * t106 - t81 * t99;
t72 = t78 * t106 + t84 * t99;
t71 = t84 * t106 - t78 * t99;
t70 = t76 * t106 + t82 * t99;
t69 = t82 * t106 - t76 * t99;
t68 = -g(1) * t115 - g(2) * t116 - g(3) * t118;
t67 = -g(1) * t71 - g(2) * t69 - g(3) * t73;
t1 = [0, 0, 0, 0, 0, 0, -t117, t94, -g(3), -t127, 0, 0, 0, 0, 0, 0, -g(3) * t103 - t117 * t110, t87, -t94, -t127, 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t123, t79, t87, t86, 0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t83 - g(3) * t89, t112, t79, t86, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t81, t113, -t112, t68, 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t74, t67, -t113, t68, 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t105 - t77 * t98) - g(2) * (t70 * t105 - t75 * t98) - g(3) * (t74 * t105 - t80 * t98) -g(1) * (-t77 * t105 - t72 * t98) - g(2) * (-t75 * t105 - t70 * t98) - g(3) * (-t80 * t105 - t74 * t98) t67, -g(1) * (t71 * pkin(4) + t115) - g(2) * (t69 * pkin(4) + t116) - g(3) * (t73 * pkin(4) + t118);];
U_reg  = t1;

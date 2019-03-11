% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:32
% EndTime: 2019-03-09 16:51:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (223->85), mult. (243->108), div. (0->0), fcn. (247->10), ass. (0->42)
t120 = g(3) * pkin(6);
t99 = sin(qJ(3));
t119 = t99 * pkin(3);
t102 = cos(qJ(3));
t87 = t102 * pkin(3) + pkin(2);
t100 = sin(qJ(2));
t118 = g(3) * t100;
t98 = -qJ(4) - pkin(8);
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t117 = t104 * pkin(1) + t101 * pkin(7);
t116 = t101 * t99;
t115 = t100 * t101;
t103 = cos(qJ(2));
t114 = t101 * t103;
t113 = t103 * t104;
t112 = t104 * t102;
t97 = qJ(3) + pkin(10);
t89 = cos(t97);
t80 = pkin(4) * t89 + t87;
t96 = -pkin(9) + t98;
t111 = t100 * t80 + t103 * t96 + pkin(6);
t92 = t101 * pkin(1);
t110 = -t104 * pkin(7) + t92;
t109 = pkin(2) * t103 + pkin(8) * t100;
t108 = g(1) * t104 + g(2) * t101;
t88 = sin(t97);
t81 = pkin(4) * t88 + t119;
t107 = -t104 * t100 * t96 + t101 * t81 + t80 * t113 + t117;
t90 = qJ(5) + t97;
t85 = sin(t90);
t86 = cos(t90);
t71 = t104 * t86 + t85 * t114;
t73 = -t101 * t86 + t85 * t113;
t106 = g(1) * t73 + g(2) * t71 + t85 * t118;
t105 = -t96 * t115 + t80 * t114 + t92 + (-pkin(7) - t81) * t104;
t83 = g(1) * t101 - g(2) * t104;
t75 = -g(3) * t103 + t108 * t100;
t74 = t101 * t85 + t86 * t113;
t72 = -t104 * t85 + t86 * t114;
t70 = -g(1) * t74 - g(2) * t72 - t86 * t118;
t1 = [0, 0, 0, 0, 0, 0, -t108, t83, -g(3), -t120, 0, 0, 0, 0, 0, 0, -t108 * t103 - t118, t75, -t83, -g(1) * t117 - g(2) * t110 - t120, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t112 + t116) - g(2) * (t102 * t114 - t104 * t99) - t102 * t118, -g(1) * (t101 * t102 - t99 * t113) - g(2) * (-t99 * t114 - t112) + t99 * t118, -t75, -g(1) * (t109 * t104 + t117) - g(2) * (t109 * t101 + t110) - g(3) * (t100 * pkin(2) - t103 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t88 + t89 * t113) - g(2) * (-t104 * t88 + t89 * t114) - t89 * t118, -g(1) * (t101 * t89 - t88 * t113) - g(2) * (-t104 * t89 - t88 * t114) + t88 * t118, -t75, -g(1) * (pkin(3) * t116 + t117) - g(2) * (t87 * t114 - t98 * t115 + t92) - g(3) * (t100 * t87 + t103 * t98 + pkin(6)) + (-g(1) * (-t100 * t98 + t103 * t87) - g(2) * (-pkin(7) - t119)) * t104, 0, 0, 0, 0, 0, 0, t70, t106, -t75, -g(1) * t107 - g(2) * t105 - g(3) * t111, 0, 0, 0, 0, 0, 0, t70, -t75, -t106, -g(1) * (t74 * pkin(5) + t73 * qJ(6) + t107) - g(2) * (t72 * pkin(5) + t71 * qJ(6) + t105) - g(3) * ((pkin(5) * t86 + qJ(6) * t85) * t100 + t111);];
U_reg  = t1;

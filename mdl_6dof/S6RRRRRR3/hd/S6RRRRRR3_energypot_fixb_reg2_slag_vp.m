% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:37
% EndTime: 2019-03-10 03:42:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (212->86), mult. (200->110), div. (0->0), fcn. (196->12), ass. (0->45)
t121 = g(3) * pkin(6);
t101 = -pkin(10) - pkin(9);
t94 = qJ(2) + qJ(3);
t85 = sin(t94);
t120 = g(3) * t85;
t95 = sin(qJ(4));
t119 = t95 * pkin(4);
t96 = sin(qJ(2));
t118 = t96 * pkin(2) + pkin(6);
t98 = cos(qJ(4));
t81 = t98 * pkin(4) + pkin(3);
t92 = -pkin(11) + t101;
t117 = t85 * t92;
t93 = qJ(4) + qJ(5);
t88 = qJ(6) + t93;
t79 = sin(t88);
t97 = sin(qJ(1));
t116 = t97 * t79;
t80 = cos(t88);
t115 = t97 * t80;
t84 = sin(t93);
t114 = t97 * t84;
t86 = cos(t93);
t113 = t97 * t86;
t112 = t97 * t95;
t111 = t97 * t98;
t100 = cos(qJ(1));
t102 = -pkin(8) - pkin(7);
t99 = cos(qJ(2));
t82 = t99 * pkin(2) + pkin(1);
t110 = t100 * t102 + t97 * t82;
t87 = cos(t94);
t109 = t100 * t87;
t108 = t100 * t95;
t107 = t100 * t98;
t106 = t101 * t85;
t77 = t100 * t82;
t105 = -t97 * t102 + t77;
t104 = pkin(3) * t87 + pkin(9) * t85;
t103 = g(1) * t100 + g(2) * t97;
t75 = g(1) * t97 - g(2) * t100;
t74 = pkin(5) * t84 + t119;
t73 = pkin(5) * t86 + t81;
t72 = -g(3) * t87 + t103 * t85;
t1 = [0, 0, 0, 0, 0, 0, -t103, t75, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t96 - t103 * t99, -g(3) * t99 + t103 * t96, -t75, -g(1) * (t100 * pkin(1) + t97 * pkin(7)) - g(2) * (t97 * pkin(1) - t100 * pkin(7)) - t121, 0, 0, 0, 0, 0, 0, -t103 * t87 - t120, t72, -t75, -g(1) * t105 - g(2) * t110 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t87 * t107 + t112) - g(2) * (t87 * t111 - t108) - t98 * t120, -g(1) * (-t87 * t108 + t111) - g(2) * (-t87 * t112 - t107) + t95 * t120, -t72, -g(1) * (t104 * t100 + t105) - g(2) * (t104 * t97 + t110) - g(3) * (t85 * pkin(3) - t87 * pkin(9) + t118) 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t86 + t114) - g(2) * (-t100 * t84 + t87 * t113) - t86 * t120, -g(1) * (-t109 * t84 + t113) - g(2) * (-t100 * t86 - t87 * t114) + t84 * t120, -t72, -g(1) * (-t100 * t106 + t109 * t81 + t77) - g(2) * (-pkin(4) * t108 + t110) - g(3) * (t87 * t101 + t85 * t81 + t118) + (-g(1) * (-t102 + t119) - g(2) * (t81 * t87 - t106)) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t80 + t116) - g(2) * (-t100 * t79 + t87 * t115) - t80 * t120, -g(1) * (-t109 * t79 + t115) - g(2) * (-t100 * t80 - t87 * t116) + t79 * t120, -t72, -g(1) * (-t100 * t117 + t109 * t73 + t77) - g(2) * (-t100 * t74 + t110) - g(3) * (t85 * t73 + t87 * t92 + t118) + (-g(1) * (-t102 + t74) - g(2) * (t73 * t87 - t117)) * t97;];
U_reg  = t1;

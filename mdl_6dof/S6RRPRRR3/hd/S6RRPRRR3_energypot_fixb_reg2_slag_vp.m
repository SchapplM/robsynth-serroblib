% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:17
% EndTime: 2019-03-09 13:25:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (212->86), mult. (200->110), div. (0->0), fcn. (196->12), ass. (0->47)
t123 = g(3) * pkin(6);
t102 = -pkin(9) - pkin(8);
t92 = qJ(2) + pkin(11);
t84 = sin(t92);
t122 = g(3) * t84;
t96 = sin(qJ(4));
t121 = t96 * pkin(4);
t97 = sin(qJ(2));
t120 = t97 * pkin(2) + pkin(6);
t99 = cos(qJ(4));
t82 = t99 * pkin(4) + pkin(3);
t93 = -pkin(10) + t102;
t119 = t84 * t93;
t94 = qJ(4) + qJ(5);
t88 = qJ(6) + t94;
t80 = sin(t88);
t98 = sin(qJ(1));
t118 = t98 * t80;
t81 = cos(t88);
t117 = t98 * t81;
t86 = sin(t94);
t116 = t98 * t86;
t87 = cos(t94);
t115 = t98 * t87;
t114 = t98 * t96;
t113 = t98 * t99;
t101 = cos(qJ(1));
t100 = cos(qJ(2));
t83 = t100 * pkin(2) + pkin(1);
t95 = -pkin(7) - qJ(3);
t112 = t101 * t95 + t98 * t83;
t85 = cos(t92);
t111 = t101 * t85;
t110 = t101 * t86;
t109 = t101 * t87;
t108 = t101 * t96;
t107 = t101 * t99;
t106 = t102 * t84;
t77 = t101 * t83;
t105 = -t98 * t95 + t77;
t104 = pkin(3) * t85 + pkin(8) * t84;
t103 = g(1) * t101 + g(2) * t98;
t75 = g(1) * t98 - g(2) * t101;
t74 = pkin(5) * t86 + t121;
t73 = pkin(5) * t87 + t82;
t72 = -g(3) * t85 + t103 * t84;
t1 = [0, 0, 0, 0, 0, 0, -t103, t75, -g(3), -t123, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t103 * t100, -g(3) * t100 + t103 * t97, -t75, -g(1) * (t101 * pkin(1) + t98 * pkin(7)) - g(2) * (t98 * pkin(1) - t101 * pkin(7)) - t123, 0, 0, 0, 0, 0, 0, -t103 * t85 - t122, t72, -t75, -g(1) * t105 - g(2) * t112 - g(3) * t120, 0, 0, 0, 0, 0, 0, -g(1) * (t85 * t107 + t114) - g(2) * (t85 * t113 - t108) - t99 * t122, -g(1) * (-t85 * t108 + t113) - g(2) * (-t85 * t114 - t107) + t96 * t122, -t72, -g(1) * (t104 * t101 + t105) - g(2) * (t104 * t98 + t112) - g(3) * (t84 * pkin(3) - t85 * pkin(8) + t120) 0, 0, 0, 0, 0, 0, -g(1) * (t85 * t109 + t116) - g(2) * (t85 * t115 - t110) - t87 * t122, -g(1) * (-t85 * t110 + t115) - g(2) * (-t85 * t116 - t109) + t86 * t122, -t72, -g(1) * (-t101 * t106 + t82 * t111 + t77) - g(2) * (-pkin(4) * t108 + t112) - g(3) * (t85 * t102 + t84 * t82 + t120) + (-g(1) * (-t95 + t121) - g(2) * (t82 * t85 - t106)) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t111 + t118) - g(2) * (-t101 * t80 + t85 * t117) - t81 * t122, -g(1) * (-t80 * t111 + t117) - g(2) * (-t101 * t81 - t85 * t118) + t80 * t122, -t72, -g(1) * (-t101 * t119 + t73 * t111 + t77) - g(2) * (-t101 * t74 + t112) - g(3) * (t84 * t73 + t85 * t93 + t120) + (-g(1) * (t74 - t95) - g(2) * (t73 * t85 - t119)) * t98;];
U_reg  = t1;

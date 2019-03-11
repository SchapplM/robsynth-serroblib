% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:18
% EndTime: 2019-03-09 22:09:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (212->86), mult. (200->111), div. (0->0), fcn. (196->12), ass. (0->44)
t120 = g(3) * pkin(6);
t94 = qJ(2) + qJ(3);
t87 = sin(t94);
t119 = g(3) * t87;
t96 = sin(qJ(4));
t118 = t96 * pkin(4);
t97 = sin(qJ(2));
t117 = t97 * pkin(2) + pkin(6);
t99 = cos(qJ(4));
t81 = t99 * pkin(4) + pkin(3);
t93 = qJ(4) + pkin(11);
t86 = qJ(6) + t93;
t78 = sin(t86);
t98 = sin(qJ(1));
t116 = t98 * t78;
t79 = cos(t86);
t115 = t98 * t79;
t84 = sin(t93);
t114 = t98 * t84;
t85 = cos(t93);
t113 = t98 * t85;
t112 = t98 * t96;
t111 = t98 * t99;
t95 = -qJ(5) - pkin(9);
t101 = cos(qJ(1));
t102 = -pkin(8) - pkin(7);
t100 = cos(qJ(2));
t82 = t100 * pkin(2) + pkin(1);
t110 = t101 * t102 + t98 * t82;
t109 = t101 * t87;
t88 = cos(t94);
t108 = t101 * t88;
t107 = t101 * t96;
t106 = t101 * t99;
t77 = t101 * t82;
t105 = -t98 * t102 + t77;
t104 = pkin(3) * t88 + pkin(9) * t87;
t103 = g(1) * t101 + g(2) * t98;
t92 = -pkin(10) + t95;
t75 = g(1) * t98 - g(2) * t101;
t74 = pkin(5) * t84 + t118;
t73 = pkin(5) * t85 + t81;
t72 = -g(3) * t88 + t103 * t87;
t1 = [0, 0, 0, 0, 0, 0, -t103, t75, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t103 * t100, -g(3) * t100 + t103 * t97, -t75, -g(1) * (t101 * pkin(1) + t98 * pkin(7)) - g(2) * (t98 * pkin(1) - t101 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -t103 * t88 - t119, t72, -t75, -g(1) * t105 - g(2) * t110 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t88 * t106 + t112) - g(2) * (t88 * t111 - t107) - t99 * t119, -g(1) * (-t88 * t107 + t111) - g(2) * (-t88 * t112 - t106) + t96 * t119, -t72, -g(1) * (t104 * t101 + t105) - g(2) * (t104 * t98 + t110) - g(3) * (t87 * pkin(3) - t88 * pkin(9) + t117) 0, 0, 0, 0, 0, 0, -g(1) * (t85 * t108 + t114) - g(2) * (-t101 * t84 + t88 * t113) - t85 * t119, -g(1) * (-t84 * t108 + t113) - g(2) * (-t101 * t85 - t88 * t114) + t84 * t119, -t72, -g(1) * (t81 * t108 - t95 * t109 + t77) - g(2) * (-pkin(4) * t107 + t110) - g(3) * (t87 * t81 + t88 * t95 + t117) + (-g(1) * (-t102 + t118) - g(2) * (t81 * t88 - t87 * t95)) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t108 + t116) - g(2) * (-t101 * t78 + t88 * t115) - t79 * t119, -g(1) * (-t78 * t108 + t115) - g(2) * (-t101 * t79 - t88 * t116) + t78 * t119, -t72, -g(1) * (t73 * t108 - t92 * t109 + t77) - g(2) * (-t101 * t74 + t110) - g(3) * (t87 * t73 + t88 * t92 + t117) + (-g(1) * (-t102 + t74) - g(2) * (t73 * t88 - t87 * t92)) * t98;];
U_reg  = t1;

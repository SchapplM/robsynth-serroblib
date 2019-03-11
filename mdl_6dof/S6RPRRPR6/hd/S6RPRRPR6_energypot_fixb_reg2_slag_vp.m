% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:42
% EndTime: 2019-03-09 05:17:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (212->86), mult. (200->111), div. (0->0), fcn. (196->12), ass. (0->48)
t120 = g(3) * pkin(6);
t89 = pkin(10) + qJ(3);
t80 = sin(t89);
t119 = g(3) * t80;
t95 = sin(qJ(4));
t118 = t95 * pkin(4);
t91 = sin(pkin(10));
t117 = t91 * pkin(2) + pkin(6);
t97 = cos(qJ(4));
t79 = t97 * pkin(4) + pkin(3);
t98 = cos(qJ(1));
t116 = t80 * t98;
t82 = cos(t89);
t115 = t82 * t98;
t90 = qJ(4) + pkin(11);
t84 = qJ(6) + t90;
t75 = sin(t84);
t96 = sin(qJ(1));
t114 = t96 * t75;
t76 = cos(t84);
t113 = t96 * t76;
t81 = sin(t90);
t112 = t96 * t81;
t83 = cos(t90);
t111 = t96 * t83;
t110 = t96 * t95;
t109 = t96 * t97;
t108 = t98 * t75;
t107 = t98 * t76;
t106 = t98 * t81;
t105 = t98 * t83;
t104 = t98 * t95;
t103 = t98 * t97;
t93 = -qJ(5) - pkin(8);
t92 = cos(pkin(10));
t77 = t92 * pkin(2) + pkin(1);
t94 = -pkin(7) - qJ(2);
t102 = t96 * t77 + t98 * t94;
t72 = t98 * t77;
t101 = -t96 * t94 + t72;
t100 = pkin(3) * t82 + pkin(8) * t80;
t99 = g(1) * t98 + g(2) * t96;
t88 = -pkin(9) + t93;
t73 = g(1) * t96 - g(2) * t98;
t70 = pkin(5) * t81 + t118;
t69 = pkin(5) * t83 + t79;
t68 = -g(3) * t82 + t99 * t80;
t1 = [0, 0, 0, 0, 0, 0, -t99, t73, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t91 - t99 * t92, -g(3) * t92 + t99 * t91, -t73, -g(1) * (t98 * pkin(1) + t96 * qJ(2)) - g(2) * (t96 * pkin(1) - t98 * qJ(2)) - t120, 0, 0, 0, 0, 0, 0, -t99 * t82 - t119, t68, -t73, -g(1) * t101 - g(2) * t102 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t103 + t110) - g(2) * (t82 * t109 - t104) - t97 * t119, -g(1) * (-t82 * t104 + t109) - g(2) * (-t82 * t110 - t103) + t95 * t119, -t68, -g(1) * (t100 * t98 + t101) - g(2) * (t100 * t96 + t102) - g(3) * (t80 * pkin(3) - t82 * pkin(8) + t117) 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t105 + t112) - g(2) * (t82 * t111 - t106) - t83 * t119, -g(1) * (-t82 * t106 + t111) - g(2) * (-t82 * t112 - t105) + t81 * t119, -t68, -g(1) * (t79 * t115 - t93 * t116 + t72) - g(2) * (-pkin(4) * t104 + t102) - g(3) * (t80 * t79 + t82 * t93 + t117) + (-g(1) * (-t94 + t118) - g(2) * (t79 * t82 - t80 * t93)) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t107 + t114) - g(2) * (t82 * t113 - t108) - t76 * t119, -g(1) * (-t82 * t108 + t113) - g(2) * (-t82 * t114 - t107) + t75 * t119, -t68, -g(1) * (t69 * t115 - t88 * t116 + t72) - g(2) * (-t98 * t70 + t102) - g(3) * (t80 * t69 + t82 * t88 + t117) + (-g(1) * (t70 - t94) - g(2) * (t69 * t82 - t80 * t88)) * t96;];
U_reg  = t1;

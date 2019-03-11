% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:32
% EndTime: 2019-03-09 07:14:32
% DurationCPUTime: 0.18s
% Computational Cost: add. (212->86), mult. (200->110), div. (0->0), fcn. (196->12), ass. (0->49)
t121 = g(3) * pkin(6);
t98 = -pkin(9) - pkin(8);
t88 = pkin(11) + qJ(3);
t80 = sin(t88);
t120 = g(3) * t80;
t94 = sin(qJ(4));
t119 = t94 * pkin(4);
t91 = sin(pkin(11));
t118 = t91 * pkin(2) + pkin(6);
t96 = cos(qJ(4));
t79 = t96 * pkin(4) + pkin(3);
t89 = -pkin(10) + t98;
t117 = t80 * t89;
t116 = t80 * t98;
t81 = cos(t88);
t97 = cos(qJ(1));
t115 = t81 * t97;
t90 = qJ(4) + qJ(5);
t85 = qJ(6) + t90;
t77 = sin(t85);
t95 = sin(qJ(1));
t114 = t95 * t77;
t78 = cos(t85);
t113 = t95 * t78;
t82 = sin(t90);
t112 = t95 * t82;
t83 = cos(t90);
t111 = t95 * t83;
t110 = t95 * t94;
t109 = t95 * t96;
t108 = t97 * t77;
t107 = t97 * t78;
t106 = t97 * t82;
t105 = t97 * t83;
t104 = t97 * t94;
t103 = t97 * t96;
t92 = cos(pkin(11));
t75 = t92 * pkin(2) + pkin(1);
t93 = -pkin(7) - qJ(2);
t102 = t95 * t75 + t97 * t93;
t72 = t97 * t75;
t101 = -t95 * t93 + t72;
t100 = pkin(3) * t81 + pkin(8) * t80;
t99 = g(1) * t97 + g(2) * t95;
t73 = g(1) * t95 - g(2) * t97;
t70 = pkin(5) * t82 + t119;
t69 = pkin(5) * t83 + t79;
t68 = -g(3) * t81 + t99 * t80;
t1 = [0, 0, 0, 0, 0, 0, -t99, t73, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t91 - t99 * t92, -g(3) * t92 + t99 * t91, -t73, -g(1) * (t97 * pkin(1) + t95 * qJ(2)) - g(2) * (t95 * pkin(1) - t97 * qJ(2)) - t121, 0, 0, 0, 0, 0, 0, -t99 * t81 - t120, t68, -t73, -g(1) * t101 - g(2) * t102 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t103 + t110) - g(2) * (t81 * t109 - t104) - t96 * t120, -g(1) * (-t81 * t104 + t109) - g(2) * (-t81 * t110 - t103) + t94 * t120, -t68, -g(1) * (t100 * t97 + t101) - g(2) * (t100 * t95 + t102) - g(3) * (t80 * pkin(3) - t81 * pkin(8) + t118) 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t105 + t112) - g(2) * (t81 * t111 - t106) - t83 * t120, -g(1) * (-t81 * t106 + t111) - g(2) * (-t81 * t112 - t105) + t82 * t120, -t68, -g(1) * (t79 * t115 - t97 * t116 + t72) - g(2) * (-pkin(4) * t104 + t102) - g(3) * (t80 * t79 + t81 * t98 + t118) + (-g(1) * (-t93 + t119) - g(2) * (t79 * t81 - t116)) * t95, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t107 + t114) - g(2) * (t81 * t113 - t108) - t78 * t120, -g(1) * (-t81 * t108 + t113) - g(2) * (-t81 * t114 - t107) + t77 * t120, -t68, -g(1) * (t69 * t115 - t97 * t117 + t72) - g(2) * (-t97 * t70 + t102) - g(3) * (t80 * t69 + t81 * t89 + t118) + (-g(1) * (t70 - t93) - g(2) * (t69 * t81 - t117)) * t95;];
U_reg  = t1;

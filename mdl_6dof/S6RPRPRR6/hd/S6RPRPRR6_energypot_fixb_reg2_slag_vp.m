% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:52
% EndTime: 2019-03-09 03:52:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (212->86), mult. (200->111), div. (0->0), fcn. (196->12), ass. (0->48)
t120 = g(3) * pkin(6);
t90 = pkin(10) + qJ(3);
t81 = sin(t90);
t119 = g(3) * t81;
t91 = sin(pkin(11));
t118 = t91 * pkin(4);
t92 = sin(pkin(10));
t117 = t92 * pkin(2) + pkin(6);
t93 = cos(pkin(11));
t77 = t93 * pkin(4) + pkin(3);
t98 = cos(qJ(1));
t116 = t81 * t98;
t83 = cos(t90);
t115 = t83 * t98;
t89 = pkin(11) + qJ(5);
t84 = qJ(6) + t89;
t75 = sin(t84);
t97 = sin(qJ(1));
t114 = t97 * t75;
t76 = cos(t84);
t113 = t97 * t76;
t80 = sin(t89);
t112 = t97 * t80;
t82 = cos(t89);
t111 = t97 * t82;
t110 = t97 * t91;
t109 = t97 * t93;
t108 = t98 * t75;
t107 = t98 * t76;
t106 = t98 * t80;
t105 = t98 * t82;
t104 = t98 * t91;
t103 = t98 * t93;
t95 = -pkin(8) - qJ(4);
t94 = cos(pkin(10));
t78 = t94 * pkin(2) + pkin(1);
t96 = -pkin(7) - qJ(2);
t102 = t97 * t78 + t98 * t96;
t72 = t98 * t78;
t101 = -t97 * t96 + t72;
t100 = g(1) * t98 + g(2) * t97;
t99 = pkin(3) * t83 + qJ(4) * t81;
t88 = -pkin(9) + t95;
t73 = g(1) * t97 - g(2) * t98;
t70 = pkin(5) * t80 + t118;
t69 = pkin(5) * t82 + t77;
t68 = -g(3) * t83 + t100 * t81;
t1 = [0, 0, 0, 0, 0, 0, -t100, t73, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t92 - t100 * t94, -g(3) * t94 + t100 * t92, -t73, -g(1) * (t98 * pkin(1) + t97 * qJ(2)) - g(2) * (t97 * pkin(1) - t98 * qJ(2)) - t120, 0, 0, 0, 0, 0, 0, -t100 * t83 - t119, t68, -t73, -g(1) * t101 - g(2) * t102 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t103 + t110) - g(2) * (t83 * t109 - t104) - t93 * t119, -g(1) * (-t83 * t104 + t109) - g(2) * (-t83 * t110 - t103) + t91 * t119, -t68, -g(1) * (t99 * t98 + t101) - g(2) * (t99 * t97 + t102) - g(3) * (t81 * pkin(3) - t83 * qJ(4) + t117) 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t105 + t112) - g(2) * (t83 * t111 - t106) - t82 * t119, -g(1) * (-t83 * t106 + t111) - g(2) * (-t83 * t112 - t105) + t80 * t119, -t68, -g(1) * (t77 * t115 - t95 * t116 + t72) - g(2) * (-pkin(4) * t104 + t102) - g(3) * (t81 * t77 + t83 * t95 + t117) + (-g(1) * (-t96 + t118) - g(2) * (t77 * t83 - t81 * t95)) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t107 + t114) - g(2) * (t83 * t113 - t108) - t76 * t119, -g(1) * (-t83 * t108 + t113) - g(2) * (-t83 * t114 - t107) + t75 * t119, -t68, -g(1) * (t69 * t115 - t88 * t116 + t72) - g(2) * (-t98 * t70 + t102) - g(3) * (t81 * t69 + t83 * t88 + t117) + (-g(1) * (t70 - t96) - g(2) * (t69 * t83 - t81 * t88)) * t97;];
U_reg  = t1;

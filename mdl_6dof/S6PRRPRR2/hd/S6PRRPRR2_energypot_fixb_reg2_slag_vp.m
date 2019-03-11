% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:15
% EndTime: 2019-03-08 22:00:15
% DurationCPUTime: 0.28s
% Computational Cost: add. (331->106), mult. (533->155), div. (0->0), fcn. (633->14), ass. (0->54)
t105 = cos(pkin(11));
t103 = sin(pkin(11));
t104 = sin(pkin(6));
t132 = t103 * t104;
t135 = t105 * pkin(1) + pkin(7) * t132;
t108 = sin(qJ(5));
t110 = sin(qJ(2));
t106 = cos(pkin(6));
t113 = cos(qJ(2));
t124 = t106 * t113;
t80 = t103 * t110 - t105 * t124;
t134 = t80 * t108;
t82 = t103 * t124 + t105 * t110;
t133 = t82 * t108;
t131 = t104 * t105;
t109 = sin(qJ(3));
t130 = t104 * t109;
t129 = t104 * t110;
t112 = cos(qJ(3));
t128 = t104 * t112;
t127 = t104 * t113;
t126 = t106 * t109;
t125 = t106 * t110;
t123 = t106 * pkin(7) + qJ(1);
t122 = t103 * t130;
t121 = t108 * t127;
t98 = t103 * pkin(1);
t120 = -pkin(7) * t131 + t98;
t107 = -qJ(4) - pkin(8);
t83 = -t103 * t125 + t105 * t113;
t93 = t112 * pkin(3) + pkin(2);
t119 = pkin(3) * t122 - t82 * t107 + t83 * t93 + t135;
t118 = pkin(3) * t126 + t107 * t127 + t93 * t129 + t123;
t117 = g(1) * t103 - g(2) * t105;
t81 = t103 * t113 + t105 * t125;
t101 = qJ(3) + pkin(12);
t94 = sin(t101);
t95 = cos(t101);
t70 = t95 * t131 + t81 * t94;
t72 = -t95 * t132 + t83 * t94;
t76 = -t106 * t95 + t94 * t129;
t116 = g(1) * t72 + g(2) * t70 + g(3) * t76;
t69 = -g(1) * t82 - g(2) * t80 + g(3) * t127;
t115 = t81 * t93 - t80 * t107 + t98 + (-pkin(3) * t109 - pkin(7)) * t131;
t114 = -pkin(10) - pkin(9);
t111 = cos(qJ(5));
t102 = qJ(5) + qJ(6);
t97 = cos(t102);
t96 = sin(t102);
t92 = t111 * pkin(5) + pkin(4);
t77 = t106 * t94 + t95 * t129;
t73 = t94 * t132 + t83 * t95;
t71 = -t94 * t131 + t81 * t95;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t103, t117, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t129, -t69, -g(3) * t106 - t117 * t104, -g(1) * t135 - g(2) * t120 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t112 + t122) - g(2) * (-t105 * t130 + t81 * t112) - g(3) * (t110 * t128 + t126) -g(1) * (t103 * t128 - t83 * t109) - g(2) * (-t105 * t128 - t81 * t109) - g(3) * (t106 * t112 - t109 * t129) t69, -g(1) * (t83 * pkin(2) + t82 * pkin(8) + t135) - g(2) * (t81 * pkin(2) + t80 * pkin(8) + t120) - g(3) * ((pkin(2) * t110 - pkin(8) * t113) * t104 + t123) 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - g(3) * t77, t116, t69, -g(1) * t119 - g(2) * t115 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t111 + t133) - g(2) * (t71 * t111 + t134) - g(3) * (t77 * t111 - t121) -g(1) * (-t73 * t108 + t82 * t111) - g(2) * (-t71 * t108 + t80 * t111) - g(3) * (-t77 * t108 - t111 * t127) -t116, -g(1) * (t73 * pkin(4) + t72 * pkin(9) + t119) - g(2) * (t71 * pkin(4) + t70 * pkin(9) + t115) - g(3) * (t77 * pkin(4) + t76 * pkin(9) + t118) 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t97 + t82 * t96) - g(2) * (t71 * t97 + t80 * t96) - g(3) * (-t96 * t127 + t77 * t97) -g(1) * (-t73 * t96 + t82 * t97) - g(2) * (-t71 * t96 + t80 * t97) - g(3) * (-t97 * t127 - t77 * t96) -t116, -g(1) * (pkin(5) * t133 - t72 * t114 + t73 * t92 + t119) - g(2) * (pkin(5) * t134 - t70 * t114 + t71 * t92 + t115) - g(3) * (-pkin(5) * t121 - t76 * t114 + t77 * t92 + t118);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:11
% EndTime: 2019-03-08 21:02:11
% DurationCPUTime: 0.28s
% Computational Cost: add. (331->106), mult. (533->155), div. (0->0), fcn. (633->14), ass. (0->54)
t107 = cos(pkin(10));
t104 = sin(pkin(10));
t105 = sin(pkin(6));
t132 = t104 * t105;
t135 = t107 * pkin(1) + pkin(7) * t132;
t103 = sin(pkin(12));
t112 = sin(qJ(2));
t108 = cos(pkin(6));
t114 = cos(qJ(2));
t124 = t108 * t114;
t80 = t104 * t112 - t107 * t124;
t134 = t80 * t103;
t82 = t104 * t124 + t107 * t112;
t133 = t82 * t103;
t131 = t105 * t107;
t111 = sin(qJ(3));
t130 = t105 * t111;
t129 = t105 * t112;
t113 = cos(qJ(3));
t128 = t105 * t113;
t127 = t105 * t114;
t126 = t108 * t111;
t125 = t108 * t112;
t123 = t108 * pkin(7) + qJ(1);
t122 = t104 * t130;
t121 = t103 * t127;
t98 = t104 * pkin(1);
t120 = -pkin(7) * t131 + t98;
t109 = -qJ(4) - pkin(8);
t83 = -t104 * t125 + t107 * t114;
t93 = t113 * pkin(3) + pkin(2);
t119 = pkin(3) * t122 - t82 * t109 + t83 * t93 + t135;
t118 = pkin(3) * t126 + t109 * t127 + t93 * t129 + t123;
t117 = g(1) * t104 - g(2) * t107;
t81 = t104 * t114 + t107 * t125;
t102 = qJ(3) + pkin(11);
t95 = sin(t102);
t97 = cos(t102);
t70 = t97 * t131 + t81 * t95;
t72 = -t97 * t132 + t83 * t95;
t76 = -t108 * t97 + t95 * t129;
t116 = g(1) * t72 + g(2) * t70 + g(3) * t76;
t69 = -g(1) * t82 - g(2) * t80 + g(3) * t127;
t115 = t81 * t93 - t80 * t109 + t98 + (-pkin(3) * t111 - pkin(7)) * t131;
t110 = -pkin(9) - qJ(5);
t106 = cos(pkin(12));
t101 = pkin(12) + qJ(6);
t96 = cos(t101);
t94 = sin(t101);
t92 = t106 * pkin(5) + pkin(4);
t77 = t108 * t95 + t97 * t129;
t73 = t95 * t132 + t83 * t97;
t71 = -t95 * t131 + t81 * t97;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t107 - g(2) * t104, t117, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t129, -t69, -g(3) * t108 - t117 * t105, -g(1) * t135 - g(2) * t120 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t113 + t122) - g(2) * (-t107 * t130 + t81 * t113) - g(3) * (t112 * t128 + t126) -g(1) * (t104 * t128 - t83 * t111) - g(2) * (-t107 * t128 - t81 * t111) - g(3) * (t108 * t113 - t111 * t129) t69, -g(1) * (t83 * pkin(2) + t82 * pkin(8) + t135) - g(2) * (t81 * pkin(2) + t80 * pkin(8) + t120) - g(3) * ((pkin(2) * t112 - pkin(8) * t114) * t105 + t123) 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - g(3) * t77, t116, t69, -g(1) * t119 - g(2) * t115 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t106 + t133) - g(2) * (t71 * t106 + t134) - g(3) * (t77 * t106 - t121) -g(1) * (-t73 * t103 + t82 * t106) - g(2) * (-t71 * t103 + t80 * t106) - g(3) * (-t77 * t103 - t106 * t127) -t116, -g(1) * (t73 * pkin(4) + t72 * qJ(5) + t119) - g(2) * (t71 * pkin(4) + t70 * qJ(5) + t115) - g(3) * (t77 * pkin(4) + t76 * qJ(5) + t118) 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t96 + t82 * t94) - g(2) * (t71 * t96 + t80 * t94) - g(3) * (-t94 * t127 + t77 * t96) -g(1) * (-t73 * t94 + t82 * t96) - g(2) * (-t71 * t94 + t80 * t96) - g(3) * (-t96 * t127 - t77 * t94) -t116, -g(1) * (pkin(5) * t133 - t72 * t110 + t73 * t92 + t119) - g(2) * (pkin(5) * t134 - t70 * t110 + t71 * t92 + t115) - g(3) * (-pkin(5) * t121 - t76 * t110 + t77 * t92 + t118);];
U_reg  = t1;

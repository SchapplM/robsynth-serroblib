% Calculate inertial parameters regressor of potential energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:56
% EndTime: 2019-03-08 19:27:56
% DurationCPUTime: 0.32s
% Computational Cost: add. (359->99), mult. (712->156), div. (0->0), fcn. (880->14), ass. (0->57)
t101 = sin(pkin(11));
t104 = cos(pkin(11));
t110 = sin(qJ(2));
t113 = cos(qJ(2));
t137 = t110 * t101 - t104 * t113;
t102 = sin(pkin(10));
t105 = cos(pkin(10));
t103 = sin(pkin(6));
t106 = cos(pkin(6));
t127 = t106 * t110;
t83 = pkin(2) * t127 + (-pkin(7) - qJ(3)) * t103;
t95 = pkin(2) * t113 + pkin(1);
t136 = t102 * t95 + t105 * t83;
t135 = t106 * pkin(7) + qJ(1);
t134 = t102 * t103;
t133 = t103 * t105;
t109 = sin(qJ(4));
t132 = t103 * t109;
t131 = t103 * t110;
t112 = cos(qJ(4));
t130 = t103 * t112;
t128 = t106 * t109;
t126 = t106 * t113;
t124 = t102 * t132;
t123 = t105 * t132;
t122 = -t102 * t83 + t105 * t95;
t121 = pkin(2) * t131 + t106 * qJ(3) + t135;
t120 = g(1) * t102 - g(2) * t105;
t119 = t101 * t113 + t110 * t104;
t107 = -qJ(5) - pkin(8);
t115 = t137 * t106;
t71 = t102 * t115 - t105 * t119;
t82 = t119 * t106;
t72 = -t102 * t82 - t105 * t137;
t94 = pkin(4) * t112 + pkin(3);
t118 = pkin(4) * t124 + t71 * t107 + t72 * t94 + t122;
t80 = t137 * t103;
t81 = t119 * t103;
t117 = pkin(4) * t128 - t80 * t107 + t81 * t94 + t121;
t70 = -t102 * t137 + t105 * t82;
t100 = qJ(4) + pkin(12);
t96 = sin(t100);
t97 = cos(t100);
t61 = t97 * t133 + t70 * t96;
t63 = -t97 * t134 + t72 * t96;
t73 = -t106 * t97 + t81 * t96;
t116 = g(1) * t63 + g(2) * t61 + g(3) * t73;
t69 = -t102 * t119 - t105 * t115;
t60 = g(1) * t71 + g(2) * t69 - g(3) * t80;
t114 = -pkin(4) * t123 + t69 * t107 + t70 * t94 + t136;
t111 = cos(qJ(6));
t108 = sin(qJ(6));
t79 = -g(3) * t106 - t120 * t103;
t74 = t106 * t96 + t81 * t97;
t64 = t96 * t134 + t72 * t97;
t62 = -t96 * t133 + t70 * t97;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t102, t120, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t102 * t127 + t105 * t113) - g(2) * (t102 * t113 + t105 * t127) - g(3) * t131, -g(1) * (-t102 * t126 - t105 * t110) - g(2) * (-t102 * t110 + t105 * t126) - g(3) * t103 * t113, t79, -g(1) * (pkin(1) * t105 + pkin(7) * t134) - g(2) * (pkin(1) * t102 - pkin(7) * t133) - g(3) * t135, 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t81, -t60, t79, -g(1) * t122 - g(2) * t136 - g(3) * t121, 0, 0, 0, 0, 0, 0, -g(1) * (t112 * t72 + t124) - g(2) * (t112 * t70 - t123) - g(3) * (t112 * t81 + t128) -g(1) * (t102 * t130 - t109 * t72) - g(2) * (-t105 * t130 - t109 * t70) - g(3) * (t106 * t112 - t109 * t81) t60, -g(1) * (pkin(3) * t72 - pkin(8) * t71 + t122) - g(2) * (pkin(3) * t70 - pkin(8) * t69 + t136) - g(3) * (pkin(3) * t81 + t80 * pkin(8) + t121) 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t62 - g(3) * t74, t116, t60, -g(1) * t118 - g(2) * t114 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (-t108 * t71 + t111 * t64) - g(2) * (-t108 * t69 + t111 * t62) - g(3) * (t108 * t80 + t111 * t74) -g(1) * (-t108 * t64 - t111 * t71) - g(2) * (-t108 * t62 - t111 * t69) - g(3) * (-t108 * t74 + t111 * t80) -t116, -g(1) * (pkin(5) * t64 + pkin(9) * t63 + t118) - g(2) * (pkin(5) * t62 + pkin(9) * t61 + t114) - g(3) * (pkin(5) * t74 + pkin(9) * t73 + t117);];
U_reg  = t1;

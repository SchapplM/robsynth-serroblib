% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:07
% EndTime: 2019-03-09 09:26:07
% DurationCPUTime: 0.18s
% Computational Cost: add. (182->85), mult. (350->114), div. (0->0), fcn. (388->10), ass. (0->46)
t95 = cos(pkin(10));
t97 = sin(qJ(2));
t121 = t95 * t97;
t94 = sin(pkin(10));
t122 = t94 * t97;
t127 = pkin(3) * t121 + qJ(4) * t122;
t126 = g(3) * pkin(6);
t125 = g(3) * t97;
t124 = t97 * pkin(2) + pkin(6);
t96 = sin(qJ(5));
t123 = t94 * t96;
t98 = sin(qJ(1));
t120 = t97 * t98;
t101 = cos(qJ(1));
t119 = t101 * pkin(1) + t98 * pkin(7);
t100 = cos(qJ(2));
t118 = t100 * t98;
t117 = t101 * t97;
t102 = -pkin(9) - pkin(8);
t116 = t102 * t97;
t115 = t100 * t101;
t114 = pkin(5) * t96 + qJ(4);
t113 = t98 * pkin(1) - t101 * pkin(7);
t112 = t124 + t127;
t111 = pkin(2) * t115 + qJ(3) * t117 + t119;
t110 = -t100 * qJ(3) + t124;
t76 = t95 * t115 + t98 * t94;
t109 = t76 * pkin(3) + t111;
t108 = g(1) * t101 + g(2) * t98;
t107 = pkin(2) * t118 + qJ(3) * t120 + t113;
t74 = -t101 * t94 + t95 * t118;
t106 = t74 * pkin(3) + t107;
t75 = t94 * t115 - t98 * t95;
t105 = t75 * qJ(4) + t109;
t73 = t101 * t95 + t94 * t118;
t104 = g(1) * t75 + g(2) * t73 + g(3) * t122;
t103 = t73 * qJ(4) + t106;
t99 = cos(qJ(5));
t93 = qJ(5) + qJ(6);
t87 = cos(t93);
t86 = sin(t93);
t85 = t99 * pkin(5) + pkin(4);
t77 = g(1) * t98 - g(2) * t101;
t70 = -g(3) * t100 + t108 * t97;
t69 = -g(1) * t76 - g(2) * t74 - g(3) * t121;
t1 = [0, 0, 0, 0, 0, 0, -t108, t77, -g(3), -t126, 0, 0, 0, 0, 0, 0, -t108 * t100 - t125, t70, -t77, -g(1) * t119 - g(2) * t113 - t126, 0, 0, 0, 0, 0, 0, t69, t104, -t70, -g(1) * t111 - g(2) * t107 - g(3) * t110, 0, 0, 0, 0, 0, 0, t69, -t70, -t104, -g(1) * t105 - g(2) * t103 - g(3) * (t110 + t127) 0, 0, 0, 0, 0, 0, -g(1) * (t75 * t96 + t76 * t99) - g(2) * (t73 * t96 + t74 * t99) - (t95 * t99 + t123) * t125, -g(1) * (t75 * t99 - t76 * t96) - g(2) * (t73 * t99 - t74 * t96) - (t94 * t99 - t95 * t96) * t125, t70, -g(1) * (t76 * pkin(4) - pkin(8) * t117 + t105) - g(2) * (t74 * pkin(4) - pkin(8) * t120 + t103) - g(3) * (pkin(4) * t121 + (pkin(8) - qJ(3)) * t100 + t112) 0, 0, 0, 0, 0, 0, -g(1) * (t75 * t86 + t76 * t87) - g(2) * (t73 * t86 + t74 * t87) - (t86 * t94 + t87 * t95) * t125, -g(1) * (t75 * t87 - t76 * t86) - g(2) * (t73 * t87 - t74 * t86) - (-t86 * t95 + t87 * t94) * t125, t70, -g(1) * (t101 * t116 + t114 * t75 + t76 * t85 + t109) - g(2) * (t114 * t73 + t98 * t116 + t74 * t85 + t106) - g(3) * ((pkin(5) * t123 + t85 * t95) * t97 + (-qJ(3) - t102) * t100 + t112);];
U_reg  = t1;

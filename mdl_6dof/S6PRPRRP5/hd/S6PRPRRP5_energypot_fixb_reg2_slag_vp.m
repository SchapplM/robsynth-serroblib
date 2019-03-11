% Calculate inertial parameters regressor of potential energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:42
% EndTime: 2019-03-08 20:16:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (236->86), mult. (524->119), div. (0->0), fcn. (622->10), ass. (0->55)
t86 = sin(pkin(6));
t92 = sin(qJ(2));
t119 = t86 * t92;
t88 = cos(pkin(6));
t123 = t88 * pkin(3) + pkin(8) * t119;
t116 = t88 * t92;
t85 = sin(pkin(10));
t87 = cos(pkin(10));
t95 = cos(qJ(2));
t69 = t87 * t116 + t85 * t95;
t90 = sin(qJ(5));
t122 = t69 * t90;
t121 = t85 * t86;
t91 = sin(qJ(4));
t120 = t86 * t91;
t94 = cos(qJ(4));
t118 = t86 * t94;
t117 = t86 * t95;
t115 = t88 * t95;
t114 = t90 * t92;
t113 = t87 * pkin(1) + pkin(7) * t121;
t112 = qJ(3) * t95;
t111 = t88 * pkin(7) + qJ(1);
t110 = t87 * t86 * pkin(7);
t109 = pkin(2) * t119 + t111;
t108 = (-pkin(3) - pkin(7)) * t87;
t68 = -t87 * t115 + t85 * t92;
t81 = t85 * pkin(1);
t107 = t69 * pkin(2) + t68 * qJ(3) + t81;
t106 = g(1) * t85 - g(2) * t87;
t70 = t85 * t115 + t87 * t92;
t71 = -t85 * t116 + t87 * t95;
t105 = t71 * pkin(2) + t70 * qJ(3) + t113;
t104 = pkin(3) * t121 + t105;
t103 = -t86 * t112 + t109;
t102 = t69 * pkin(8) + t107;
t60 = t85 * t120 - t70 * t94;
t62 = t87 * t120 + t68 * t94;
t72 = t94 * t117 + t88 * t91;
t101 = g(1) * t60 - g(2) * t62 + g(3) * t72;
t100 = -g(1) * t70 - g(2) * t68 + g(3) * t117;
t99 = g(1) * t71 + g(2) * t69 + g(3) * t119;
t98 = t103 + t123;
t97 = t71 * pkin(8) + t104;
t96 = t86 * t108 + t102;
t93 = cos(qJ(5));
t89 = -qJ(6) - pkin(9);
t80 = t93 * pkin(5) + pkin(4);
t73 = -t91 * t117 + t88 * t94;
t64 = -g(3) * t88 - t106 * t86;
t63 = -t87 * t118 + t68 * t91;
t61 = t85 * t118 + t70 * t91;
t57 = -g(1) * (t61 * t93 + t71 * t90) - g(2) * (t63 * t93 + t122) - g(3) * (t86 * t114 + t73 * t93);
t56 = -g(1) * (-t61 * t90 + t71 * t93) - g(2) * (-t63 * t90 + t69 * t93) - g(3) * (t93 * t119 - t73 * t90);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t85, t106, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t99, -t100, t64, -g(1) * t113 - g(2) * (t81 - t110) - g(3) * t111, 0, 0, 0, 0, 0, 0, t64, t99, t100, -g(1) * t105 - g(2) * (t107 - t110) - g(3) * t103, 0, 0, 0, 0, 0, 0, -g(1) * t61 - g(2) * t63 - g(3) * t73, t101, -t99, -g(1) * t97 - g(2) * t96 - g(3) * t98, 0, 0, 0, 0, 0, 0, t57, t56, -t101, -g(1) * (t61 * pkin(4) + t60 * pkin(9) + t97) - g(2) * (t63 * pkin(4) - t62 * pkin(9) + t96) - g(3) * (t73 * pkin(4) + t72 * pkin(9) + t98) 0, 0, 0, 0, 0, 0, t57, t56, -t101, -g(1) * (-t60 * t89 + t61 * t80 + (pkin(5) * t90 + pkin(8)) * t71 + t104) - g(2) * (pkin(5) * t122 + t62 * t89 + t63 * t80 + t102) - g(3) * (-t72 * t89 + t73 * t80 + t109 + t123) + (-g(3) * (pkin(5) * t114 - t112) - g(2) * t108) * t86;];
U_reg  = t1;

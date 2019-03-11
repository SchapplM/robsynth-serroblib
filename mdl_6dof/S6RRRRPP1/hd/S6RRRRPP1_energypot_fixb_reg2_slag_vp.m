% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:46:46
% EndTime: 2019-03-09 20:46:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (215->71), mult. (215->89), div. (0->0), fcn. (215->10), ass. (0->43)
t120 = g(3) * pkin(6);
t93 = qJ(2) + qJ(3);
t88 = sin(t93);
t119 = g(3) * t88;
t96 = sin(qJ(2));
t118 = t96 * pkin(2) + pkin(6);
t94 = -qJ(5) - pkin(9);
t117 = t88 * t94;
t92 = qJ(4) + pkin(10);
t86 = sin(t92);
t97 = sin(qJ(1));
t116 = t97 * t86;
t87 = cos(t92);
t115 = t97 * t87;
t95 = sin(qJ(4));
t114 = t97 * t95;
t98 = cos(qJ(4));
t113 = t97 * t98;
t100 = cos(qJ(1));
t101 = -pkin(8) - pkin(7);
t99 = cos(qJ(2));
t84 = t99 * pkin(2) + pkin(1);
t112 = t100 * t101 + t97 * t84;
t89 = cos(t93);
t111 = t100 * t89;
t110 = t100 * t95;
t109 = t100 * t98;
t83 = t98 * pkin(4) + pkin(3);
t108 = t88 * t83 + t89 * t94 + t118;
t107 = t100 * t84 - t97 * t101;
t106 = pkin(3) * t89 + pkin(9) * t88;
t105 = g(1) * t100 + g(2) * t97;
t69 = t100 * t87 + t89 * t116;
t71 = t86 * t111 - t115;
t104 = g(1) * t71 + g(2) * t69 + t86 * t119;
t103 = pkin(4) * t114 - t100 * t117 + t83 * t111 + t107;
t102 = -pkin(4) * t110 + t112 + (t83 * t89 - t117) * t97;
t76 = g(1) * t97 - g(2) * t100;
t72 = t87 * t111 + t116;
t70 = -t100 * t86 + t89 * t115;
t68 = -g(3) * t89 + t105 * t88;
t67 = -g(1) * t72 - g(2) * t70 - t87 * t119;
t1 = [0, 0, 0, 0, 0, 0, -t105, t76, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t96 - t105 * t99, -g(3) * t99 + t105 * t96, -t76, -g(1) * (t100 * pkin(1) + t97 * pkin(7)) - g(2) * (t97 * pkin(1) - t100 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -t105 * t89 - t119, t68, -t76, -g(1) * t107 - g(2) * t112 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t109 + t114) - g(2) * (t89 * t113 - t110) - t98 * t119, -g(1) * (-t89 * t110 + t113) - g(2) * (-t89 * t114 - t109) + t95 * t119, -t68, -g(1) * (t106 * t100 + t107) - g(2) * (t106 * t97 + t112) - g(3) * (t88 * pkin(3) - t89 * pkin(9) + t118) 0, 0, 0, 0, 0, 0, t67, t104, -t68, -g(1) * t103 - g(2) * t102 - g(3) * t108, 0, 0, 0, 0, 0, 0, t67, -t68, -t104, -g(1) * (t72 * pkin(5) + t71 * qJ(6) + t103) - g(2) * (t70 * pkin(5) + t69 * qJ(6) + t102) - g(3) * ((pkin(5) * t87 + qJ(6) * t86) * t88 + t108);];
U_reg  = t1;

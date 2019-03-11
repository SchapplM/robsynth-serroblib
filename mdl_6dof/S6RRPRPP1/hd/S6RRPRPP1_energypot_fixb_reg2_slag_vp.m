% Calculate inertial parameters regressor of potential energy for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:50
% EndTime: 2019-03-09 09:47:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (215->71), mult. (215->89), div. (0->0), fcn. (215->10), ass. (0->43)
t120 = g(3) * pkin(6);
t93 = qJ(2) + pkin(9);
t87 = sin(t93);
t119 = g(3) * t87;
t97 = sin(qJ(2));
t118 = t97 * pkin(2) + pkin(6);
t94 = -qJ(5) - pkin(8);
t117 = t87 * t94;
t92 = qJ(4) + pkin(10);
t86 = sin(t92);
t98 = sin(qJ(1));
t116 = t98 * t86;
t88 = cos(t92);
t115 = t98 * t88;
t96 = sin(qJ(4));
t114 = t98 * t96;
t99 = cos(qJ(4));
t113 = t98 * t99;
t101 = cos(qJ(1));
t100 = cos(qJ(2));
t85 = t100 * pkin(2) + pkin(1);
t95 = -pkin(7) - qJ(3);
t112 = t101 * t95 + t98 * t85;
t89 = cos(t93);
t111 = t101 * t89;
t110 = t101 * t96;
t109 = t101 * t99;
t108 = t101 * t85 - t98 * t95;
t84 = t99 * pkin(4) + pkin(3);
t107 = t87 * t84 + t89 * t94 + t118;
t106 = pkin(3) * t89 + pkin(8) * t87;
t105 = g(1) * t101 + g(2) * t98;
t69 = t101 * t88 + t89 * t116;
t71 = t86 * t111 - t115;
t104 = g(1) * t71 + g(2) * t69 + t86 * t119;
t103 = pkin(4) * t114 - t101 * t117 + t84 * t111 + t108;
t102 = -pkin(4) * t110 + t112 + (t84 * t89 - t117) * t98;
t76 = g(1) * t98 - g(2) * t101;
t72 = t88 * t111 + t116;
t70 = -t101 * t86 + t89 * t115;
t68 = -g(3) * t89 + t105 * t87;
t67 = -g(1) * t72 - g(2) * t70 - t88 * t119;
t1 = [0, 0, 0, 0, 0, 0, -t105, t76, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t105 * t100, -g(3) * t100 + t105 * t97, -t76, -g(1) * (t101 * pkin(1) + t98 * pkin(7)) - g(2) * (t98 * pkin(1) - t101 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -t105 * t89 - t119, t68, -t76, -g(1) * t108 - g(2) * t112 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t109 + t114) - g(2) * (t89 * t113 - t110) - t99 * t119, -g(1) * (-t89 * t110 + t113) - g(2) * (-t89 * t114 - t109) + t96 * t119, -t68, -g(1) * (t106 * t101 + t108) - g(2) * (t106 * t98 + t112) - g(3) * (t87 * pkin(3) - t89 * pkin(8) + t118) 0, 0, 0, 0, 0, 0, t67, t104, -t68, -g(1) * t103 - g(2) * t102 - g(3) * t107, 0, 0, 0, 0, 0, 0, t67, -t68, -t104, -g(1) * (t72 * pkin(5) + t71 * qJ(6) + t103) - g(2) * (t70 * pkin(5) + t69 * qJ(6) + t102) - g(3) * ((pkin(5) * t88 + qJ(6) * t86) * t87 + t107);];
U_reg  = t1;

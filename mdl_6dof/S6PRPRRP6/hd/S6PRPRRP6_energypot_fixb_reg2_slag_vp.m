% Calculate inertial parameters regressor of potential energy for
% S6PRPRRP6
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:55
% EndTime: 2019-03-08 20:20:55
% DurationCPUTime: 0.20s
% Computational Cost: add. (251->78), mult. (570->107), div. (0->0), fcn. (686->10), ass. (0->57)
t94 = sin(pkin(10));
t95 = sin(pkin(6));
t130 = t94 * t95;
t96 = cos(pkin(10));
t129 = t95 * t96;
t99 = sin(qJ(4));
t128 = t95 * t99;
t127 = t96 * pkin(1) + pkin(7) * t130;
t100 = sin(qJ(2));
t126 = t100 * t95;
t102 = cos(qJ(4));
t125 = t102 * t95;
t103 = cos(qJ(2));
t124 = t103 * t95;
t123 = t94 * t100;
t122 = t94 * t103;
t121 = t96 * t100;
t120 = t96 * t103;
t97 = cos(pkin(6));
t119 = t97 * pkin(7) + qJ(1);
t118 = pkin(7) * t129;
t77 = -t97 * t120 + t123;
t78 = t97 * t121 + t122;
t90 = t94 * pkin(1);
t117 = t78 * pkin(2) + t77 * qJ(3) + t90;
t116 = g(1) * t94 - g(2) * t96;
t79 = t97 * t122 + t121;
t80 = -t97 * t123 + t120;
t115 = t80 * pkin(2) + t79 * qJ(3) + t127;
t114 = pkin(2) * t126 - qJ(3) * t124 + t119;
t101 = cos(qJ(5));
t65 = t94 * t125 + t79 * t99;
t98 = sin(qJ(5));
t57 = -t80 * t101 + t65 * t98;
t67 = -t96 * t125 + t77 * t99;
t59 = -t78 * t101 + t67 * t98;
t82 = t97 * t102 - t99 * t124;
t68 = -t101 * t126 + t82 * t98;
t113 = g(1) * t57 + g(2) * t59 + g(3) * t68;
t64 = -t79 * t102 + t94 * t128;
t66 = t77 * t102 + t96 * t128;
t81 = t102 * t124 + t97 * t99;
t112 = g(1) * t64 - g(2) * t66 + g(3) * t81;
t111 = g(1) * t80 + g(2) * t78 + g(3) * t126;
t110 = -g(1) * t79 - g(2) * t77 + g(3) * t124;
t109 = t97 * pkin(3) + pkin(8) * t126 + t114;
t108 = pkin(3) * t130 + t80 * pkin(8) + t115;
t107 = t78 * pkin(8) + (-pkin(3) - pkin(7)) * t129 + t117;
t106 = t82 * pkin(4) + t81 * pkin(9) + t109;
t105 = t65 * pkin(4) + t64 * pkin(9) + t108;
t104 = t67 * pkin(4) - t66 * pkin(9) + t107;
t70 = -g(3) * t97 - t116 * t95;
t69 = t82 * t101 + t98 * t126;
t60 = t67 * t101 + t78 * t98;
t58 = t65 * t101 + t80 * t98;
t55 = -g(1) * t58 - g(2) * t60 - g(3) * t69;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t111, -t110, t70, -g(1) * t127 - g(2) * (t90 - t118) - g(3) * t119, 0, 0, 0, 0, 0, 0, t70, t111, t110, -g(1) * t115 - g(2) * (t117 - t118) - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t67 - g(3) * t82, t112, -t111, -g(1) * t108 - g(2) * t107 - g(3) * t109, 0, 0, 0, 0, 0, 0, t55, t113, -t112, -g(1) * t105 - g(2) * t104 - g(3) * t106, 0, 0, 0, 0, 0, 0, t55, -t112, -t113, -g(1) * (t58 * pkin(5) + t57 * qJ(6) + t105) - g(2) * (t60 * pkin(5) + t59 * qJ(6) + t104) - g(3) * (t69 * pkin(5) + t68 * qJ(6) + t106);];
U_reg  = t1;

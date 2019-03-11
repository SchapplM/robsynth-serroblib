% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:10
% EndTime: 2019-03-08 22:34:10
% DurationCPUTime: 0.25s
% Computational Cost: add. (277->96), mult. (600->131), div. (0->0), fcn. (727->12), ass. (0->54)
t133 = pkin(4) + pkin(8);
t96 = sin(pkin(6));
t132 = pkin(7) * t96;
t100 = sin(qJ(2));
t102 = cos(qJ(2));
t122 = cos(pkin(6));
t115 = t102 * t122;
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t76 = t100 * t95 - t97 * t115;
t131 = t76 * pkin(8);
t78 = t97 * t100 + t95 * t115;
t130 = t78 * pkin(8);
t101 = cos(qJ(5));
t129 = pkin(5) * t101 + t133;
t128 = cos(qJ(3));
t99 = sin(qJ(3));
t127 = t96 * t99;
t126 = t97 * pkin(1) + t95 * t132;
t125 = t100 * t96;
t124 = t102 * t96;
t123 = t122 * pkin(7) + qJ(1);
t121 = pkin(8) * t124;
t116 = t100 * t122;
t79 = t97 * t102 - t95 * t116;
t120 = t79 * pkin(2) + t126;
t119 = pkin(2) * t125 + t123;
t118 = t96 * t128;
t98 = sin(qJ(5));
t117 = pkin(5) * t98 + qJ(4);
t72 = t95 * t127 + t79 * t128;
t114 = t72 * pkin(3) + t120;
t81 = t100 * t118 + t122 * t99;
t113 = t81 * pkin(3) + t119;
t112 = t95 * pkin(1) - t97 * t132;
t111 = g(1) * t95 - g(2) * t97;
t77 = t95 * t102 + t97 * t116;
t110 = t77 * pkin(2) + t112;
t70 = -t97 * t127 + t77 * t128;
t109 = t70 * pkin(3) + t110;
t71 = -t95 * t118 + t79 * t99;
t108 = t71 * qJ(4) + t114;
t80 = -t122 * t128 + t99 * t125;
t107 = t80 * qJ(4) + t113;
t69 = t97 * t118 + t77 * t99;
t106 = g(1) * t71 + g(2) * t69 + g(3) * t80;
t105 = g(1) * t72 + g(2) * t70 + g(3) * t81;
t66 = -g(1) * t78 - g(2) * t76 + g(3) * t124;
t104 = t69 * qJ(4) + t109;
t103 = -pkin(10) - pkin(9);
t94 = qJ(5) + qJ(6);
t90 = cos(t94);
t89 = sin(t94);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t97 - g(2) * t95, t111, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t79 - g(2) * t77 - g(3) * t125, -t66, -g(3) * t122 - t111 * t96, -g(1) * t126 - g(2) * t112 - g(3) * t123, 0, 0, 0, 0, 0, 0, -t105, t106, t66, -g(1) * (t120 + t130) - g(2) * (t110 + t131) - g(3) * (t119 - t121) 0, 0, 0, 0, 0, 0, t66, t105, -t106, -g(1) * (t108 + t130) - g(2) * (t104 + t131) - g(3) * (t107 - t121) 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t78 + t71 * t98) - g(2) * (t101 * t76 + t69 * t98) - g(3) * (-t101 * t124 + t80 * t98) -g(1) * (t101 * t71 - t78 * t98) - g(2) * (t101 * t69 - t76 * t98) - g(3) * (t101 * t80 + t98 * t124) -t105, -g(1) * (pkin(9) * t72 + t133 * t78 + t108) - g(2) * (t70 * pkin(9) + t133 * t76 + t104) - g(3) * (t81 * pkin(9) - t133 * t124 + t107) 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t89 + t78 * t90) - g(2) * (t69 * t89 + t76 * t90) - g(3) * (-t90 * t124 + t80 * t89) -g(1) * (t71 * t90 - t78 * t89) - g(2) * (t69 * t90 - t76 * t89) - g(3) * (t89 * t124 + t80 * t90) -t105, -g(1) * (-t103 * t72 + t117 * t71 + t129 * t78 + t114) - g(2) * (-t103 * t70 + t117 * t69 + t129 * t76 + t109) - g(3) * (-t81 * t103 + t117 * t80 - t129 * t124 + t113);];
U_reg  = t1;

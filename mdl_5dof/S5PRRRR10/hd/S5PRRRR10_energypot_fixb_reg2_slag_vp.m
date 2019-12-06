% Calculate inertial parameters regressor of potential energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:05
% EndTime: 2019-12-05 17:26:05
% DurationCPUTime: 0.32s
% Computational Cost: add. (336->85), mult. (878->140), div. (0->0), fcn. (1108->14), ass. (0->56)
t103 = sin(pkin(6));
t106 = cos(pkin(6));
t107 = cos(pkin(5));
t104 = sin(pkin(5));
t114 = cos(qJ(2));
t134 = t104 * t114;
t87 = -t103 * t134 + t107 * t106;
t102 = sin(pkin(11));
t136 = t104 * t106;
t105 = cos(pkin(11));
t111 = sin(qJ(2));
t131 = t107 * t114;
t90 = -t102 * t131 - t105 * t111;
t81 = t102 * t136 - t90 * t103;
t141 = cos(qJ(3));
t138 = t102 * t104;
t139 = t105 * pkin(1) + pkin(7) * t138;
t137 = t104 * t105;
t135 = t104 * t111;
t132 = t107 * t111;
t130 = t107 * pkin(7) + qJ(1);
t127 = t103 * t141;
t126 = t106 * t141;
t125 = t104 * t127;
t124 = t102 * pkin(1) - pkin(7) * t137;
t123 = g(1) * t102 - g(2) * t105;
t88 = -t102 * t111 + t105 * t131;
t80 = -t103 * t88 - t105 * t136;
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t110 = sin(qJ(3));
t89 = t102 * t114 + t105 * t132;
t70 = t89 * t141 + (-t103 * t137 + t106 * t88) * t110;
t63 = t109 * t70 - t80 * t113;
t91 = -t102 * t132 + t105 * t114;
t72 = t91 * t141 + (t103 * t138 + t106 * t90) * t110;
t65 = t109 * t72 - t81 * t113;
t79 = t107 * t103 * t110 + (t106 * t110 * t114 + t111 * t141) * t104;
t73 = t109 * t79 - t87 * t113;
t122 = g(1) * t65 + g(2) * t63 + g(3) * t73;
t69 = t105 * t125 + t89 * t110 - t126 * t88;
t71 = -t102 * t125 + t91 * t110 - t126 * t90;
t78 = -t107 * t127 + t110 * t135 - t126 * t134;
t121 = g(1) * t71 + g(2) * t69 + g(3) * t78;
t120 = t91 * pkin(2) + pkin(8) * t81 + t139;
t119 = pkin(2) * t135 + pkin(8) * t87 + t130;
t118 = t72 * pkin(3) + pkin(9) * t71 + t120;
t117 = t79 * pkin(3) + pkin(9) * t78 + t119;
t116 = t89 * pkin(2) + pkin(8) * t80 + t124;
t115 = t70 * pkin(3) + t69 * pkin(9) + t116;
t112 = cos(qJ(5));
t108 = sin(qJ(5));
t74 = t109 * t87 + t113 * t79;
t66 = t109 * t81 + t113 * t72;
t64 = t109 * t80 + t113 * t70;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t102, t123, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89 - g(3) * t135, -g(1) * t90 - g(2) * t88 - g(3) * t134, -g(3) * t107 - t104 * t123, -g(1) * t139 - g(2) * t124 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t79, t121, -g(1) * t81 - g(2) * t80 - g(3) * t87, -g(1) * t120 - g(2) * t116 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * t66 - g(2) * t64 - g(3) * t74, t122, -t121, -g(1) * t118 - g(2) * t115 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t71 + t112 * t66) - g(2) * (t108 * t69 + t112 * t64) - g(3) * (t108 * t78 + t112 * t74), -g(1) * (-t108 * t66 + t112 * t71) - g(2) * (-t108 * t64 + t112 * t69) - g(3) * (-t108 * t74 + t112 * t78), -t122, -g(1) * (pkin(4) * t66 + pkin(10) * t65 + t118) - g(2) * (t64 * pkin(4) + t63 * pkin(10) + t115) - g(3) * (pkin(4) * t74 + pkin(10) * t73 + t117);];
U_reg = t1;

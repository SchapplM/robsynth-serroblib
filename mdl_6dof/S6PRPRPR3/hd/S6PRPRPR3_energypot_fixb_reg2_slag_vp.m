% Calculate inertial parameters regressor of potential energy for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:38
% EndTime: 2019-03-08 19:37:38
% DurationCPUTime: 0.25s
% Computational Cost: add. (348->90), mult. (806->133), div. (0->0), fcn. (1011->12), ass. (0->59)
t131 = pkin(5) + pkin(8);
t100 = sin(qJ(2));
t103 = cos(qJ(2));
t117 = cos(pkin(11));
t114 = t103 * t117;
t93 = sin(pkin(11));
t82 = -t100 * t93 + t114;
t97 = cos(pkin(6));
t104 = t82 * t97;
t81 = -t100 * t117 - t103 * t93;
t94 = sin(pkin(10));
t96 = cos(pkin(10));
t66 = t96 * t104 + t94 * t81;
t130 = t66 * pkin(8);
t68 = -t94 * t104 + t96 * t81;
t129 = t68 * pkin(8);
t95 = sin(pkin(6));
t124 = t100 * t95;
t77 = -t95 * t114 + t93 * t124;
t128 = t77 * pkin(8);
t127 = t95 * pkin(7);
t99 = sin(qJ(4));
t126 = t95 * t99;
t80 = t97 * t100 * pkin(2) + (-pkin(7) - qJ(3)) * t95;
t90 = t103 * pkin(2) + pkin(1);
t125 = t96 * t80 + t94 * t90;
t102 = cos(qJ(4));
t123 = t102 * t95;
t122 = t94 * t100;
t121 = t94 * t103;
t120 = t96 * t100;
t119 = t96 * t103;
t118 = t97 * pkin(7) + qJ(1);
t79 = t81 * t97;
t67 = -t96 * t79 + t94 * t82;
t116 = t67 * pkin(3) + t125;
t115 = -t94 * t80 + t96 * t90;
t113 = pkin(2) * t124 + t97 * qJ(3) + t118;
t69 = t94 * t79 + t96 * t82;
t112 = t69 * pkin(3) + t115;
t111 = g(1) * t94 - g(2) * t96;
t78 = t81 * t95;
t110 = -t78 * pkin(3) + t113;
t60 = t96 * t123 + t67 * t99;
t61 = t67 * t102 - t96 * t126;
t109 = t61 * pkin(4) + t60 * qJ(5) + t116;
t62 = -t94 * t123 + t69 * t99;
t71 = -t97 * t102 - t78 * t99;
t108 = g(1) * t62 + g(2) * t60 + g(3) * t71;
t63 = t69 * t102 + t94 * t126;
t72 = -t78 * t102 + t97 * t99;
t107 = g(1) * t63 + g(2) * t61 + g(3) * t72;
t57 = g(1) * t68 + g(2) * t66 - g(3) * t77;
t106 = t63 * pkin(4) + t62 * qJ(5) + t112;
t105 = t72 * pkin(4) + t71 * qJ(5) + t110;
t101 = cos(qJ(6));
t98 = sin(qJ(6));
t76 = -g(3) * t97 - t111 * t95;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94, t111, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t97 * t122 + t119) - g(2) * (t97 * t120 + t121) - g(3) * t124, -g(1) * (-t97 * t121 - t120) - g(2) * (t97 * t119 - t122) - g(3) * t95 * t103, t76, -g(1) * (t96 * pkin(1) + t94 * t127) - g(2) * (t94 * pkin(1) - t96 * t127) - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * t69 - g(2) * t67 + g(3) * t78, -t57, t76, -g(1) * t115 - g(2) * t125 - g(3) * t113, 0, 0, 0, 0, 0, 0, -t107, t108, t57, -g(1) * (t112 - t129) - g(2) * (t116 - t130) - g(3) * (t110 + t128) 0, 0, 0, 0, 0, 0, t57, t107, -t108, -g(1) * (t106 - t129) - g(2) * (t109 - t130) - g(3) * (t105 + t128) 0, 0, 0, 0, 0, 0, -g(1) * (-t68 * t101 + t62 * t98) - g(2) * (-t66 * t101 + t60 * t98) - g(3) * (t77 * t101 + t71 * t98) -g(1) * (t62 * t101 + t68 * t98) - g(2) * (t60 * t101 + t66 * t98) - g(3) * (t71 * t101 - t77 * t98) -t107, -g(1) * (t63 * pkin(9) - t131 * t68 + t106) - g(2) * (t61 * pkin(9) - t131 * t66 + t109) - g(3) * (t72 * pkin(9) + t131 * t77 + t105);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:05
% EndTime: 2019-03-08 20:03:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (394->85), mult. (926->133), div. (0->0), fcn. (1177->12), ass. (0->59)
t102 = sin(pkin(11));
t105 = cos(pkin(11));
t110 = sin(qJ(2));
t113 = cos(qJ(2));
t138 = t110 * t102 - t105 * t113;
t104 = sin(pkin(6));
t137 = pkin(7) * t104;
t103 = sin(pkin(10));
t106 = cos(pkin(10));
t107 = cos(pkin(6));
t131 = t107 * t110;
t90 = pkin(2) * t131 + (-pkin(7) - qJ(3)) * t104;
t99 = pkin(2) * t113 + pkin(1);
t136 = t103 * t99 + t106 * t90;
t109 = sin(qJ(4));
t135 = t104 * t109;
t134 = t104 * t110;
t112 = cos(qJ(4));
t133 = t104 * t112;
t130 = t107 * t113;
t128 = t107 * pkin(7) + qJ(1);
t127 = -t103 * t90 + t106 * t99;
t126 = pkin(2) * t134 + t107 * qJ(3) + t128;
t125 = g(1) * t103 - g(2) * t106;
t118 = t138 * t107;
t123 = t102 * t113 + t110 * t105;
t75 = -t103 * t123 - t106 * t118;
t89 = t123 * t107;
t76 = -t103 * t138 + t106 * t89;
t124 = t76 * pkin(3) - t75 * pkin(8) + t136;
t77 = t103 * t118 - t106 * t123;
t78 = -t103 * t89 - t106 * t138;
t122 = t78 * pkin(3) - pkin(8) * t77 + t127;
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t68 = -t106 * t135 + t112 * t76;
t59 = t108 * t68 + t75 * t111;
t70 = t103 * t135 + t112 * t78;
t61 = t108 * t70 + t77 * t111;
t88 = t123 * t104;
t81 = t107 * t109 + t112 * t88;
t87 = t138 * t104;
t63 = t108 * t81 - t87 * t111;
t121 = g(1) * t61 + g(2) * t59 + g(3) * t63;
t67 = t106 * t133 + t109 * t76;
t69 = -t103 * t133 + t109 * t78;
t80 = -t107 * t112 + t109 * t88;
t120 = g(1) * t69 + g(2) * t67 + g(3) * t80;
t119 = g(1) * t77 + g(2) * t75 - g(3) * t87;
t117 = t88 * pkin(3) + pkin(8) * t87 + t126;
t116 = t68 * pkin(4) + pkin(9) * t67 + t124;
t115 = t70 * pkin(4) + pkin(9) * t69 + t122;
t114 = t81 * pkin(4) + pkin(9) * t80 + t117;
t86 = -g(3) * t107 - t125 * t104;
t64 = t108 * t87 + t111 * t81;
t62 = -t108 * t77 + t111 * t70;
t60 = -t108 * t75 + t111 * t68;
t57 = -g(1) * t62 - g(2) * t60 - g(3) * t64;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t103, t125, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t103 * t131 + t106 * t113) - g(2) * (t103 * t113 + t106 * t131) - g(3) * t134, -g(1) * (-t103 * t130 - t106 * t110) - g(2) * (-t103 * t110 + t106 * t130) - g(3) * t104 * t113, t86, -g(1) * (pkin(1) * t106 + t103 * t137) - g(2) * (pkin(1) * t103 - t106 * t137) - g(3) * t128, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t88, -t119, t86, -g(1) * t127 - g(2) * t136 - g(3) * t126, 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t81, t120, t119, -g(1) * t122 - g(2) * t124 - g(3) * t117, 0, 0, 0, 0, 0, 0, t57, t121, -t120, -g(1) * t115 - g(2) * t116 - g(3) * t114, 0, 0, 0, 0, 0, 0, t57, -t120, -t121, -g(1) * (pkin(5) * t62 + qJ(6) * t61 + t115) - g(2) * (pkin(5) * t60 + qJ(6) * t59 + t116) - g(3) * (pkin(5) * t64 + qJ(6) * t63 + t114);];
U_reg  = t1;

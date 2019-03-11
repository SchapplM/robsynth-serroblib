% Calculate inertial parameters regressor of potential energy for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:11
% EndTime: 2019-03-08 20:34:11
% DurationCPUTime: 0.25s
% Computational Cost: add. (333->107), mult. (459->156), div. (0->0), fcn. (534->14), ass. (0->52)
t104 = cos(pkin(12));
t91 = pkin(3) * t104 + pkin(2);
t102 = sin(pkin(11));
t129 = g(1) * t102;
t105 = cos(pkin(11));
t128 = g(2) * t105;
t101 = sin(pkin(12));
t127 = t101 * pkin(3);
t107 = -pkin(8) - qJ(3);
t103 = sin(pkin(6));
t124 = t102 * t103;
t126 = pkin(1) * t105 + pkin(7) * t124;
t106 = cos(pkin(6));
t125 = pkin(7) * t106 + qJ(1);
t123 = t103 * t105;
t109 = sin(qJ(2));
t122 = t103 * t109;
t111 = cos(qJ(2));
t121 = t103 * t111;
t120 = t106 * t101;
t119 = t106 * t109;
t118 = t106 * t111;
t100 = pkin(12) + qJ(4);
t95 = t102 * pkin(1);
t117 = -pkin(7) * t123 + t95;
t79 = t102 * t118 + t105 * t109;
t80 = -t102 * t119 + t105 * t111;
t93 = cos(t100);
t83 = pkin(4) * t93 + t91;
t92 = sin(t100);
t84 = pkin(4) * t92 + t127;
t99 = -pkin(9) + t107;
t116 = t124 * t84 - t79 * t99 + t80 * t83 + t126;
t115 = t106 * t84 + t99 * t121 + t83 * t122 + t125;
t114 = -t128 + t129;
t78 = t102 * t111 + t105 * t119;
t94 = qJ(5) + t100;
t89 = sin(t94);
t90 = cos(t94);
t65 = t123 * t90 + t78 * t89;
t67 = -t124 * t90 + t80 * t89;
t71 = -t106 * t90 + t122 * t89;
t113 = g(1) * t67 + g(2) * t65 + g(3) * t71;
t77 = t102 * t109 - t105 * t118;
t112 = t78 * t83 - t77 * t99 + t95 + (-pkin(7) - t84) * t123;
t64 = -g(1) * t79 - g(2) * t77 + g(3) * t121;
t110 = cos(qJ(6));
t108 = sin(qJ(6));
t72 = t106 * t89 + t122 * t90;
t68 = t124 * t89 + t80 * t90;
t66 = -t123 * t89 + t78 * t90;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t102, t114, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t78 - g(3) * t122, -t64, -g(3) * t106 - t103 * t114, -g(1) * t126 - g(2) * t117 - g(3) * t125, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t124 + t104 * t80) - g(2) * (-t101 * t123 + t104 * t78) - g(3) * (t104 * t122 + t120) -g(1) * (-t101 * t80 + t104 * t124) - g(2) * (-t101 * t78 - t104 * t123) - g(3) * (-t101 * t122 + t104 * t106) t64, -g(1) * (pkin(2) * t80 + qJ(3) * t79 + t126) - g(2) * (pkin(2) * t78 + qJ(3) * t77 + t117) - g(3) * ((pkin(2) * t109 - qJ(3) * t111) * t103 + t125) 0, 0, 0, 0, 0, 0, -g(1) * (t124 * t92 + t80 * t93) - g(2) * (-t123 * t92 + t78 * t93) - g(3) * (t106 * t92 + t122 * t93) -g(1) * (t124 * t93 - t80 * t92) - g(2) * (-t123 * t93 - t78 * t92) - g(3) * (t106 * t93 - t122 * t92) t64, -g(1) * (-t107 * t79 + t80 * t91 + t126) - g(2) * (-t77 * t107 + t78 * t91 + t95) - g(3) * (pkin(3) * t120 + t125) + (-t127 * t129 - g(3) * (t107 * t111 + t109 * t91) - (-pkin(7) - t127) * t128) * t103, 0, 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t66 - g(3) * t72, t113, t64, -g(1) * t116 - g(2) * t112 - g(3) * t115, 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t79 + t110 * t68) - g(2) * (t108 * t77 + t110 * t66) - g(3) * (-t108 * t121 + t110 * t72) -g(1) * (-t108 * t68 + t110 * t79) - g(2) * (-t108 * t66 + t110 * t77) - g(3) * (-t108 * t72 - t110 * t121) -t113, -g(1) * (pkin(5) * t68 + pkin(10) * t67 + t116) - g(2) * (t66 * pkin(5) + t65 * pkin(10) + t112) - g(3) * (pkin(5) * t72 + pkin(10) * t71 + t115);];
U_reg  = t1;

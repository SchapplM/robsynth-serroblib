% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:45
% EndTime: 2019-03-09 13:59:45
% DurationCPUTime: 0.20s
% Computational Cost: add. (164->79), mult. (310->101), div. (0->0), fcn. (338->10), ass. (0->43)
t103 = sin(qJ(2));
t104 = sin(qJ(1));
t107 = cos(qJ(2));
t122 = t104 * t107;
t131 = t104 * t103 * qJ(3) + pkin(2) * t122;
t108 = cos(qJ(1));
t130 = pkin(3) * t122 + t108 * pkin(8);
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t129 = t107 * t102 - t103 * t106;
t128 = g(3) * pkin(6);
t127 = g(3) * t129;
t101 = sin(qJ(5));
t126 = pkin(5) * t101;
t125 = t108 * pkin(1) + t104 * pkin(7);
t123 = t103 * t108;
t120 = t107 * t108;
t96 = t104 * pkin(1);
t119 = -t108 * pkin(7) + t96;
t118 = pkin(2) * t120 + qJ(3) * t123 + t125;
t117 = t103 * pkin(2) - t107 * qJ(3) + pkin(6);
t116 = pkin(3) * t120 + t118;
t115 = g(1) * t108 + g(2) * t104;
t114 = t119 + t131;
t113 = t103 * pkin(3) + t117;
t79 = t103 * t102 + t107 * t106;
t75 = t129 * t104;
t77 = t102 * t120 - t106 * t123;
t112 = g(1) * t77 + g(2) * t75 + g(3) * t79;
t111 = t114 + t130;
t110 = -t104 * pkin(8) + t116;
t109 = -pkin(10) - pkin(9);
t105 = cos(qJ(5));
t100 = qJ(5) + qJ(6);
t92 = cos(t100);
t91 = sin(t100);
t90 = t105 * pkin(5) + pkin(4);
t81 = g(1) * t104 - g(2) * t108;
t78 = t79 * t108;
t76 = t79 * t104;
t74 = -g(3) * t103 - t115 * t107;
t73 = -g(3) * t107 + t115 * t103;
t1 = [0, 0, 0, 0, 0, 0, -t115, t81, -g(3), -t128, 0, 0, 0, 0, 0, 0, t74, t73, -t81, -g(1) * t125 - g(2) * t119 - t128, 0, 0, 0, 0, 0, 0, t74, -t81, -t73, -g(1) * t118 - g(2) * t114 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 + t127, t112, t81, -g(1) * t110 - g(2) * t111 - g(3) * t113, 0, 0, 0, 0, 0, 0, -g(1) * (-t104 * t101 + t78 * t105) - g(2) * (t108 * t101 + t76 * t105) + t105 * t127, -g(1) * (-t78 * t101 - t104 * t105) - g(2) * (-t76 * t101 + t108 * t105) - t101 * t127, -t112, -g(1) * (t78 * pkin(4) + t77 * pkin(9) + t110) - g(2) * (t76 * pkin(4) + t75 * pkin(9) + t111) - g(3) * (-pkin(4) * t129 + t79 * pkin(9) + t113) 0, 0, 0, 0, 0, 0, -g(1) * (-t104 * t91 + t78 * t92) - g(2) * (t108 * t91 + t76 * t92) + t92 * t127, -g(1) * (-t104 * t92 - t78 * t91) - g(2) * (t108 * t92 - t76 * t91) - t91 * t127, -t112, -g(1) * (-t77 * t109 + t78 * t90 + (-pkin(8) - t126) * t104 + t116) - g(2) * (-t75 * t109 + t76 * t90 + t96 + (-pkin(7) + t126) * t108 + t130 + t131) - g(3) * (-t79 * t109 - t129 * t90 + t113);];
U_reg  = t1;

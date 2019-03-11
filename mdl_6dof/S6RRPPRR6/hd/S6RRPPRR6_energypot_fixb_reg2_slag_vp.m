% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:14
% EndTime: 2019-03-09 09:15:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (184->75), mult. (278->97), div. (0->0), fcn. (293->10), ass. (0->46)
t105 = sin(qJ(1));
t107 = cos(qJ(2));
t122 = t105 * t107;
t104 = sin(qJ(2));
t124 = t104 * t105;
t132 = pkin(2) * t122 + qJ(3) * t124;
t108 = cos(qJ(1));
t115 = g(1) * t108 + g(2) * t105;
t131 = g(3) * pkin(6);
t99 = pkin(10) + qJ(5);
t92 = sin(t99);
t93 = cos(t99);
t77 = t104 * t93 - t107 * t92;
t130 = g(3) * t77;
t129 = t104 * pkin(2) + pkin(6);
t126 = t108 * pkin(1) + t105 * pkin(7);
t100 = sin(pkin(10));
t125 = t104 * t100;
t123 = t104 * t108;
t121 = t107 * t108;
t96 = t105 * pkin(1);
t120 = t96 + t132;
t119 = pkin(4) * t125;
t118 = -t108 * pkin(7) + t96;
t117 = pkin(2) * t121 + qJ(3) * t123 + t126;
t116 = -t107 * qJ(3) + t129;
t76 = t104 * t92 + t107 * t93;
t101 = cos(pkin(10));
t114 = t107 * t100 - t104 * t101;
t113 = t107 * t101 + t125;
t102 = -pkin(8) - qJ(4);
t90 = t101 * pkin(4) + pkin(3);
t112 = t105 * t102 + t108 * t119 + t90 * t121 + t117;
t70 = t92 * t122 - t93 * t124;
t72 = t92 * t121 - t93 * t123;
t111 = g(1) * t72 + g(2) * t70 + g(3) * t76;
t110 = t104 * t90 + (-pkin(4) * t100 - qJ(3)) * t107 + t129;
t109 = t90 * t122 + t105 * t119 + (-pkin(7) - t102) * t108 + t120;
t106 = cos(qJ(6));
t103 = sin(qJ(6));
t85 = g(1) * t105 - g(2) * t108;
t75 = -g(3) * t104 - t115 * t107;
t74 = -g(3) * t107 + t115 * t104;
t73 = t76 * t108;
t71 = t76 * t105;
t1 = [0, 0, 0, 0, 0, 0, -t115, t85, -g(3), -t131, 0, 0, 0, 0, 0, 0, t75, t74, -t85, -g(1) * t126 - g(2) * t118 - t131, 0, 0, 0, 0, 0, 0, t75, -t85, -t74, -g(1) * t117 - g(2) * (t118 + t132) - g(3) * t116, 0, 0, 0, 0, 0, 0, g(3) * t114 - t115 * t113, g(3) * t113 + t115 * t114, t85, -g(1) * (pkin(3) * t121 - t105 * qJ(4) + t117) - g(2) * (pkin(3) * t122 + (-pkin(7) + qJ(4)) * t108 + t120) - g(3) * (t104 * pkin(3) + t116) 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - t130, t111, t85, -g(1) * t112 - g(2) * t109 - g(3) * t110, 0, 0, 0, 0, 0, 0, -g(1) * (-t105 * t103 + t73 * t106) - g(2) * (t108 * t103 + t71 * t106) - t106 * t130, -g(1) * (-t73 * t103 - t105 * t106) - g(2) * (-t71 * t103 + t108 * t106) + t103 * t130, -t111, -g(1) * (t73 * pkin(5) + t72 * pkin(9) + t112) - g(2) * (t71 * pkin(5) + t70 * pkin(9) + t109) - g(3) * (t77 * pkin(5) + t76 * pkin(9) + t110);];
U_reg  = t1;

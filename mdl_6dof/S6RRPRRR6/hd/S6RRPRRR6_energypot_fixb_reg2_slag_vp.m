% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR6
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:24
% EndTime: 2019-03-09 13:53:24
% DurationCPUTime: 0.19s
% Computational Cost: add. (184->75), mult. (278->97), div. (0->0), fcn. (293->10), ass. (0->46)
t104 = sin(qJ(1));
t107 = cos(qJ(2));
t123 = t104 * t107;
t103 = sin(qJ(2));
t125 = t103 * t104;
t133 = pkin(2) * t123 + qJ(3) * t125;
t108 = cos(qJ(1));
t116 = g(1) * t108 + g(2) * t104;
t132 = g(3) * pkin(6);
t100 = qJ(4) + qJ(5);
t93 = sin(t100);
t94 = cos(t100);
t78 = t103 * t94 - t107 * t93;
t131 = g(3) * t78;
t130 = t103 * pkin(2) + pkin(6);
t127 = t108 * pkin(1) + t104 * pkin(7);
t102 = sin(qJ(4));
t126 = t103 * t102;
t124 = t103 * t108;
t122 = t107 * t108;
t97 = t104 * pkin(1);
t121 = t97 + t133;
t120 = pkin(4) * t126;
t119 = -t108 * pkin(7) + t97;
t118 = pkin(2) * t122 + qJ(3) * t124 + t127;
t117 = -t107 * qJ(3) + t130;
t77 = t103 * t93 + t107 * t94;
t106 = cos(qJ(4));
t115 = t107 * t102 - t103 * t106;
t114 = t107 * t106 + t126;
t109 = -pkin(9) - pkin(8);
t91 = t106 * pkin(4) + pkin(3);
t113 = t104 * t109 + t108 * t120 + t91 * t122 + t118;
t71 = t93 * t123 - t94 * t125;
t73 = t93 * t122 - t94 * t124;
t112 = g(1) * t73 + g(2) * t71 + g(3) * t77;
t111 = t103 * t91 + (-pkin(4) * t102 - qJ(3)) * t107 + t130;
t110 = t91 * t123 + t104 * t120 + (-pkin(7) - t109) * t108 + t121;
t105 = cos(qJ(6));
t101 = sin(qJ(6));
t85 = g(1) * t104 - g(2) * t108;
t76 = -g(3) * t103 - t116 * t107;
t75 = -g(3) * t107 + t116 * t103;
t74 = t77 * t108;
t72 = t77 * t104;
t1 = [0, 0, 0, 0, 0, 0, -t116, t85, -g(3), -t132, 0, 0, 0, 0, 0, 0, t76, t75, -t85, -g(1) * t127 - g(2) * t119 - t132, 0, 0, 0, 0, 0, 0, t76, -t85, -t75, -g(1) * t118 - g(2) * (t119 + t133) - g(3) * t117, 0, 0, 0, 0, 0, 0, g(3) * t115 - t116 * t114, g(3) * t114 + t116 * t115, t85, -g(1) * (pkin(3) * t122 - t104 * pkin(8) + t118) - g(2) * (pkin(3) * t123 + (-pkin(7) + pkin(8)) * t108 + t121) - g(3) * (t103 * pkin(3) + t117) 0, 0, 0, 0, 0, 0, -g(1) * t74 - g(2) * t72 - t131, t112, t85, -g(1) * t113 - g(2) * t110 - g(3) * t111, 0, 0, 0, 0, 0, 0, -g(1) * (-t104 * t101 + t74 * t105) - g(2) * (t108 * t101 + t72 * t105) - t105 * t131, -g(1) * (-t74 * t101 - t104 * t105) - g(2) * (-t72 * t101 + t108 * t105) + t101 * t131, -t112, -g(1) * (t74 * pkin(5) + t73 * pkin(10) + t113) - g(2) * (t72 * pkin(5) + t71 * pkin(10) + t110) - g(3) * (t78 * pkin(5) + t77 * pkin(10) + t111);];
U_reg  = t1;

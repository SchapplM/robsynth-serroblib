% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:02
% EndTime: 2019-03-09 08:48:03
% DurationCPUTime: 0.18s
% Computational Cost: add. (207->65), mult. (235->85), div. (0->0), fcn. (244->10), ass. (0->40)
t108 = cos(qJ(1));
t99 = qJ(2) + pkin(10);
t95 = cos(t99);
t119 = t108 * t95;
t94 = sin(t99);
t121 = qJ(4) * t94;
t127 = pkin(3) * t119 + t108 * t121;
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t126 = t95 * t102 - t94 * t106;
t125 = g(3) * pkin(6);
t124 = g(3) * t126;
t103 = sin(qJ(2));
t123 = t103 * pkin(2) + pkin(6);
t100 = -pkin(7) - qJ(3);
t104 = sin(qJ(1));
t107 = cos(qJ(2));
t93 = t107 * pkin(2) + pkin(1);
t122 = t108 * t100 + t104 * t93;
t120 = t104 * t95;
t116 = pkin(3) * t120 + t104 * t121 + t122;
t89 = t108 * t93;
t115 = -t104 * t100 + t89;
t114 = g(1) * t108 + g(2) * t104;
t77 = t94 * t102 + t95 * t106;
t113 = t94 * pkin(3) - t95 * qJ(4) + t123;
t112 = pkin(4) * t120 + t108 * pkin(8) + t116;
t111 = t94 * pkin(4) + t113;
t73 = t126 * t104;
t75 = t126 * t108;
t110 = g(1) * t75 + g(2) * t73 + g(3) * t77;
t109 = pkin(4) * t119 + t89 + (-pkin(8) - t100) * t104 + t127;
t105 = cos(qJ(6));
t101 = sin(qJ(6));
t83 = g(1) * t104 - g(2) * t108;
t76 = t77 * t108;
t74 = t77 * t104;
t72 = -g(3) * t94 - t114 * t95;
t71 = -g(3) * t95 + t114 * t94;
t1 = [0, 0, 0, 0, 0, 0, -t114, t83, -g(3), -t125, 0, 0, 0, 0, 0, 0, -g(3) * t103 - t114 * t107, -g(3) * t107 + t114 * t103, -t83, -g(1) * (t108 * pkin(1) + t104 * pkin(7)) - g(2) * (t104 * pkin(1) - t108 * pkin(7)) - t125, 0, 0, 0, 0, 0, 0, t72, t71, -t83, -g(1) * t115 - g(2) * t122 - g(3) * t123, 0, 0, 0, 0, 0, 0, t72, -t83, -t71, -g(1) * (t115 + t127) - g(2) * t116 - g(3) * t113, 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t74 + t124, t110, t83, -g(1) * t109 - g(2) * t112 - g(3) * t111, 0, 0, 0, 0, 0, 0, -g(1) * (-t104 * t101 + t76 * t105) - g(2) * (t108 * t101 + t74 * t105) + t105 * t124, -g(1) * (-t76 * t101 - t104 * t105) - g(2) * (-t74 * t101 + t108 * t105) - t101 * t124, -t110, -g(1) * (t76 * pkin(5) + t75 * pkin(9) + t109) - g(2) * (t74 * pkin(5) + t73 * pkin(9) + t112) - g(3) * (-pkin(5) * t126 + t77 * pkin(9) + t111);];
U_reg  = t1;

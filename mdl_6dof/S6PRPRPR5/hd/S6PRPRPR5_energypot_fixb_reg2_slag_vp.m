% Calculate inertial parameters regressor of potential energy for
% S6PRPRPR5
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:02
% EndTime: 2019-03-08 19:45:02
% DurationCPUTime: 0.23s
% Computational Cost: add. (308->91), mult. (511->131), div. (0->0), fcn. (603->12), ass. (0->48)
t103 = cos(pkin(10));
t100 = sin(pkin(10));
t101 = sin(pkin(6));
t126 = t100 * t101;
t129 = t103 * pkin(1) + pkin(7) * t126;
t104 = cos(pkin(6));
t99 = sin(pkin(11));
t128 = t104 * t99;
t127 = t104 * pkin(7) + qJ(1);
t125 = t101 * t103;
t107 = sin(qJ(2));
t124 = t101 * t107;
t109 = cos(qJ(2));
t123 = t101 * t109;
t122 = t104 * t107;
t121 = t104 * t109;
t120 = t99 * t126;
t95 = t100 * pkin(1);
t119 = -pkin(7) * t125 + t95;
t105 = -pkin(8) - qJ(3);
t81 = t100 * t121 + t103 * t107;
t82 = -t100 * t122 + t103 * t109;
t102 = cos(pkin(11));
t92 = t102 * pkin(3) + pkin(2);
t118 = pkin(3) * t120 - t81 * t105 + t82 * t92 + t129;
t117 = pkin(3) * t128 + t105 * t123 + t92 * t124 + t127;
t116 = g(1) * t100 - g(2) * t103;
t80 = t100 * t109 + t103 * t122;
t98 = pkin(11) + qJ(4);
t93 = sin(t98);
t94 = cos(t98);
t68 = t94 * t125 + t80 * t93;
t70 = -t94 * t126 + t82 * t93;
t75 = -t104 * t94 + t93 * t124;
t115 = g(1) * t70 + g(2) * t68 + g(3) * t75;
t69 = -t93 * t125 + t80 * t94;
t71 = t93 * t126 + t82 * t94;
t76 = t104 * t93 + t94 * t124;
t114 = g(1) * t71 + g(2) * t69 + g(3) * t76;
t79 = t100 * t107 - t103 * t121;
t65 = -g(1) * t81 - g(2) * t79 + g(3) * t123;
t113 = t71 * pkin(4) + t70 * qJ(5) + t118;
t112 = t76 * pkin(4) + t75 * qJ(5) + t117;
t111 = t80 * t92 - t79 * t105 + t95 + (-pkin(3) * t99 - pkin(7)) * t125;
t110 = t69 * pkin(4) + t68 * qJ(5) + t111;
t108 = cos(qJ(6));
t106 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t100, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t124, -t65, -g(3) * t104 - t116 * t101, -g(1) * t129 - g(2) * t119 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t102 + t120) - g(2) * (t80 * t102 - t99 * t125) - g(3) * (t102 * t124 + t128) -g(1) * (t102 * t126 - t82 * t99) - g(2) * (-t102 * t125 - t80 * t99) - g(3) * (t104 * t102 - t99 * t124) t65, -g(1) * (t82 * pkin(2) + t81 * qJ(3) + t129) - g(2) * (t80 * pkin(2) + t79 * qJ(3) + t119) - g(3) * ((pkin(2) * t107 - qJ(3) * t109) * t101 + t127) 0, 0, 0, 0, 0, 0, -t114, t115, t65, -g(1) * t118 - g(2) * t111 - g(3) * t117, 0, 0, 0, 0, 0, 0, t65, t114, -t115, -g(1) * t113 - g(2) * t110 - g(3) * t112, 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t106 + t81 * t108) - g(2) * (t68 * t106 + t79 * t108) - g(3) * (t75 * t106 - t108 * t123) -g(1) * (-t81 * t106 + t70 * t108) - g(2) * (-t79 * t106 + t68 * t108) - g(3) * (t106 * t123 + t75 * t108) -t114, -g(1) * (t81 * pkin(5) + t71 * pkin(9) + t113) - g(2) * (t79 * pkin(5) + t69 * pkin(9) + t110) - g(3) * (-pkin(5) * t123 + t76 * pkin(9) + t112);];
U_reg  = t1;

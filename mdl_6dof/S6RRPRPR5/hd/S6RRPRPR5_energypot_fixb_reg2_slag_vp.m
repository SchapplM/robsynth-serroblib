% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:35:10
% EndTime: 2019-03-09 10:35:10
% DurationCPUTime: 0.29s
% Computational Cost: add. (376->101), mult. (841->154), div. (0->0), fcn. (1059->14), ass. (0->58)
t115 = sin(pkin(11));
t121 = sin(qJ(2));
t124 = cos(qJ(2));
t148 = cos(pkin(11));
t100 = -t121 * t115 + t124 * t148;
t118 = cos(pkin(6));
t150 = t118 * pkin(8) + pkin(7);
t108 = t124 * pkin(2) + pkin(1);
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t116 = sin(pkin(6));
t98 = t118 * t121 * pkin(2) + (-pkin(8) - qJ(3)) * t116;
t149 = t122 * t108 + t125 * t98;
t147 = t116 * t121;
t146 = t116 * t122;
t145 = t116 * t125;
t143 = t122 * t121;
t142 = t122 * t124;
t141 = t125 * t121;
t140 = t125 * t124;
t99 = -t124 * t115 - t121 * t148;
t97 = t99 * t118;
t86 = t122 * t100 - t125 * t97;
t139 = t86 * pkin(3) + t149;
t114 = sin(pkin(12));
t138 = pkin(5) * t114 + pkin(9);
t136 = t125 * t108 - t122 * t98;
t135 = pkin(2) * t147 + t118 * qJ(3) + t150;
t88 = t125 * t100 + t122 * t97;
t134 = t88 * pkin(3) + t136;
t96 = t99 * t116;
t133 = -t96 * pkin(3) + t135;
t132 = g(1) * t122 - g(2) * t125;
t126 = t100 * t118;
t85 = t122 * t99 + t125 * t126;
t131 = -t85 * pkin(9) + t139;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t79 = t86 * t120 + t123 * t145;
t81 = t88 * t120 - t123 * t146;
t89 = -t118 * t123 - t96 * t120;
t130 = g(1) * t81 + g(2) * t79 + g(3) * t89;
t87 = -t122 * t126 + t125 * t99;
t95 = t100 * t116;
t129 = g(1) * t87 + g(2) * t85 + g(3) * t95;
t128 = -t87 * pkin(9) + t134;
t127 = -t95 * pkin(9) + t133;
t119 = -pkin(10) - qJ(5);
t117 = cos(pkin(12));
t113 = pkin(12) + qJ(6);
t110 = cos(t113);
t109 = sin(t113);
t107 = t117 * pkin(5) + pkin(4);
t94 = -g(3) * t118 - t116 * t132;
t90 = t118 * t120 - t96 * t123;
t82 = t120 * t146 + t88 * t123;
t80 = -t120 * t145 + t86 * t123;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t122, t132, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t118 * t143 + t140) - g(2) * (t118 * t141 + t142) - g(3) * t147, -g(1) * (-t118 * t142 - t141) - g(2) * (t118 * t140 - t143) - g(3) * t116 * t124, t94, -g(1) * (t125 * pkin(1) + pkin(8) * t146) - g(2) * (t122 * pkin(1) - pkin(8) * t145) - g(3) * t150, 0, 0, 0, 0, 0, 0, -g(1) * t88 - g(2) * t86 + g(3) * t96, -t129, t94, -g(1) * t136 - g(2) * t149 - g(3) * t135, 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t90, t130, t129, -g(1) * t128 - g(2) * t131 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (-t87 * t114 + t82 * t117) - g(2) * (-t85 * t114 + t80 * t117) - g(3) * (-t95 * t114 + t90 * t117) -g(1) * (-t82 * t114 - t87 * t117) - g(2) * (-t80 * t114 - t85 * t117) - g(3) * (-t90 * t114 - t95 * t117) -t130, -g(1) * (t82 * pkin(4) + t81 * qJ(5) + t128) - g(2) * (t80 * pkin(4) + t79 * qJ(5) + t131) - g(3) * (t90 * pkin(4) + t89 * qJ(5) + t127) 0, 0, 0, 0, 0, 0, -g(1) * (-t87 * t109 + t82 * t110) - g(2) * (-t85 * t109 + t80 * t110) - g(3) * (-t95 * t109 + t90 * t110) -g(1) * (-t82 * t109 - t87 * t110) - g(2) * (-t80 * t109 - t85 * t110) - g(3) * (-t90 * t109 - t95 * t110) -t130, -g(1) * (t82 * t107 - t81 * t119 - t138 * t87 + t134) - g(2) * (t80 * t107 - t79 * t119 - t138 * t85 + t139) - g(3) * (t90 * t107 - t89 * t119 - t138 * t95 + t133);];
U_reg  = t1;

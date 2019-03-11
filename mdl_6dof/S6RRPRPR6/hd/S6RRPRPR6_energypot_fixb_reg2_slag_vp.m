% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:42:54
% EndTime: 2019-03-09 10:42:54
% DurationCPUTime: 0.22s
% Computational Cost: add. (348->89), mult. (806->131), div. (0->0), fcn. (1011->12), ass. (0->57)
t113 = sin(pkin(11));
t118 = sin(qJ(2));
t122 = cos(qJ(2));
t145 = cos(pkin(11));
t102 = -t118 * t113 + t122 * t145;
t151 = pkin(5) + pkin(9);
t101 = -t122 * t113 - t118 * t145;
t119 = sin(qJ(1));
t123 = cos(qJ(1));
t115 = cos(pkin(6));
t124 = t102 * t115;
t86 = t119 * t101 + t123 * t124;
t150 = t86 * pkin(9);
t88 = t123 * t101 - t119 * t124;
t149 = t88 * pkin(9);
t114 = sin(pkin(6));
t97 = t102 * t114;
t148 = t97 * pkin(9);
t147 = t115 * pkin(8) + pkin(7);
t100 = t115 * t118 * pkin(2) + (-pkin(8) - qJ(3)) * t114;
t110 = t122 * pkin(2) + pkin(1);
t146 = t123 * t100 + t119 * t110;
t144 = t114 * t118;
t143 = t114 * t119;
t142 = t114 * t123;
t140 = t119 * t118;
t139 = t119 * t122;
t138 = t123 * t118;
t137 = t123 * t122;
t99 = t101 * t115;
t87 = t119 * t102 - t123 * t99;
t136 = t87 * pkin(3) + t146;
t134 = -t119 * t100 + t123 * t110;
t133 = pkin(2) * t144 + t115 * qJ(3) + t147;
t89 = t123 * t102 + t119 * t99;
t132 = t89 * pkin(3) + t134;
t98 = t101 * t114;
t131 = -t98 * pkin(3) + t133;
t130 = g(1) * t119 - g(2) * t123;
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t80 = t87 * t117 + t121 * t142;
t81 = -t117 * t142 + t87 * t121;
t129 = t81 * pkin(4) + t80 * qJ(5) + t136;
t82 = t89 * t117 - t121 * t143;
t91 = -t115 * t121 - t98 * t117;
t128 = g(1) * t82 + g(2) * t80 + g(3) * t91;
t83 = t117 * t143 + t89 * t121;
t92 = t115 * t117 - t98 * t121;
t127 = g(1) * t83 + g(2) * t81 + g(3) * t92;
t77 = g(1) * t88 + g(2) * t86 + g(3) * t97;
t126 = t83 * pkin(4) + t82 * qJ(5) + t132;
t125 = t92 * pkin(4) + t91 * qJ(5) + t131;
t120 = cos(qJ(6));
t116 = sin(qJ(6));
t96 = -g(3) * t115 - t130 * t114;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t123 - g(2) * t119, t130, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t115 * t140 + t137) - g(2) * (t115 * t138 + t139) - g(3) * t144, -g(1) * (-t115 * t139 - t138) - g(2) * (t115 * t137 - t140) - g(3) * t114 * t122, t96, -g(1) * (t123 * pkin(1) + pkin(8) * t143) - g(2) * (t119 * pkin(1) - pkin(8) * t142) - g(3) * t147, 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 + g(3) * t98, -t77, t96, -g(1) * t134 - g(2) * t146 - g(3) * t133, 0, 0, 0, 0, 0, 0, -t127, t128, t77, -g(1) * (t132 - t149) - g(2) * (t136 - t150) - g(3) * (t131 - t148) 0, 0, 0, 0, 0, 0, t77, t127, -t128, -g(1) * (t126 - t149) - g(2) * (t129 - t150) - g(3) * (t125 - t148) 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t116 - t88 * t120) - g(2) * (t80 * t116 - t86 * t120) - g(3) * (t91 * t116 - t97 * t120) -g(1) * (t88 * t116 + t82 * t120) - g(2) * (t86 * t116 + t80 * t120) - g(3) * (t97 * t116 + t91 * t120) -t127, -g(1) * (t83 * pkin(10) - t151 * t88 + t126) - g(2) * (t81 * pkin(10) - t151 * t86 + t129) - g(3) * (t92 * pkin(10) - t151 * t97 + t125);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:52
% EndTime: 2019-03-09 00:25:52
% DurationCPUTime: 0.34s
% Computational Cost: add. (527->100), mult. (1376->154), div. (0->0), fcn. (1753->14), ass. (0->64)
t117 = sin(pkin(7));
t120 = cos(pkin(7));
t121 = cos(pkin(6));
t118 = sin(pkin(6));
t128 = cos(qJ(2));
t156 = t118 * t128;
t140 = t117 * t156 - t121 * t120;
t116 = sin(pkin(12));
t119 = cos(pkin(12));
t126 = sin(qJ(2));
t153 = t121 * t128;
t103 = -t116 * t153 - t119 * t126;
t158 = t118 * t120;
t141 = -t103 * t117 + t116 * t158;
t163 = cos(qJ(3));
t162 = cos(qJ(4));
t160 = t116 * t118;
t159 = t118 * t119;
t157 = t118 * t126;
t154 = t121 * t126;
t152 = t119 * pkin(1) + pkin(8) * t160;
t151 = t121 * pkin(8) + qJ(1);
t123 = sin(qJ(5));
t148 = pkin(5) * t123 + pkin(10);
t147 = t117 * t163;
t146 = t120 * t163;
t145 = t118 * t147;
t144 = t116 * pkin(1) - pkin(8) * t159;
t143 = g(1) * t116 - g(2) * t119;
t101 = -t116 * t126 + t119 * t153;
t142 = t101 * t117 + t119 * t158;
t124 = sin(qJ(4));
t102 = t116 * t128 + t119 * t154;
t125 = sin(qJ(3));
t86 = t102 * t163 + (t101 * t120 - t117 * t159) * t125;
t79 = t86 * t124 + t142 * t162;
t104 = -t116 * t154 + t119 * t128;
t88 = t104 * t163 + (t103 * t120 + t117 * t160) * t125;
t81 = t88 * t124 - t141 * t162;
t95 = t121 * t117 * t125 + (t120 * t125 * t128 + t126 * t163) * t118;
t89 = t95 * t124 + t140 * t162;
t139 = g(1) * t81 + g(2) * t79 + g(3) * t89;
t85 = -t101 * t146 + t102 * t125 + t119 * t145;
t87 = -t103 * t146 + t104 * t125 - t116 * t145;
t94 = -t121 * t147 + t125 * t157 - t146 * t156;
t138 = g(1) * t87 + g(2) * t85 + g(3) * t94;
t137 = t104 * pkin(2) + t141 * pkin(9) + t152;
t136 = t88 * pkin(3) + t137;
t135 = pkin(2) * t157 - t140 * pkin(9) + t151;
t134 = t95 * pkin(3) + t135;
t133 = t87 * pkin(10) + t136;
t132 = t94 * pkin(10) + t134;
t131 = t102 * pkin(2) - pkin(9) * t142 + t144;
t130 = t86 * pkin(3) + t131;
t129 = t85 * pkin(10) + t130;
t127 = cos(qJ(5));
t122 = -qJ(6) - pkin(11);
t112 = t127 * pkin(5) + pkin(4);
t90 = -t124 * t140 + t162 * t95;
t82 = t124 * t141 + t162 * t88;
t80 = -t124 * t142 + t162 * t86;
t77 = -g(1) * (t87 * t123 + t82 * t127) - g(2) * (t85 * t123 + t80 * t127) - g(3) * (t94 * t123 + t90 * t127);
t76 = -g(1) * (-t82 * t123 + t87 * t127) - g(2) * (-t80 * t123 + t85 * t127) - g(3) * (-t90 * t123 + t94 * t127);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t119 - g(2) * t116, t143, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t157, -g(1) * t103 - g(2) * t101 - g(3) * t156, -g(3) * t121 - t118 * t143, -g(1) * t152 - g(2) * t144 - g(3) * t151, 0, 0, 0, 0, 0, 0, -g(1) * t88 - g(2) * t86 - g(3) * t95, t138, -g(1) * t141 + g(2) * t142 + g(3) * t140, -g(1) * t137 - g(2) * t131 - g(3) * t135, 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t90, t139, -t138, -g(1) * t133 - g(2) * t129 - g(3) * t132, 0, 0, 0, 0, 0, 0, t77, t76, -t139, -g(1) * (t82 * pkin(4) + t81 * pkin(11) + t133) - g(2) * (t80 * pkin(4) + t79 * pkin(11) + t129) - g(3) * (t90 * pkin(4) + t89 * pkin(11) + t132) 0, 0, 0, 0, 0, 0, t77, t76, -t139, -g(1) * (t82 * t112 - t81 * t122 + t148 * t87 + t136) - g(2) * (t80 * t112 - t79 * t122 + t148 * t85 + t130) - g(3) * (t90 * t112 - t89 * t122 + t148 * t94 + t134);];
U_reg  = t1;

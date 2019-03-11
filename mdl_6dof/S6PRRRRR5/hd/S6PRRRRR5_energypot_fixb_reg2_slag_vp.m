% Calculate inertial parameters regressor of potential energy for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:00
% EndTime: 2019-03-09 01:10:01
% DurationCPUTime: 0.38s
% Computational Cost: add. (539->111), mult. (1376->172), div. (0->0), fcn. (1753->16), ass. (0->65)
t121 = sin(pkin(7));
t124 = cos(pkin(7));
t125 = cos(pkin(6));
t122 = sin(pkin(6));
t131 = cos(qJ(2));
t160 = t122 * t131;
t144 = t121 * t160 - t125 * t124;
t120 = sin(pkin(13));
t123 = cos(pkin(13));
t129 = sin(qJ(2));
t157 = t125 * t131;
t104 = -t120 * t157 - t123 * t129;
t162 = t122 * t124;
t145 = -t104 * t121 + t120 * t162;
t167 = cos(qJ(3));
t166 = cos(qJ(4));
t164 = t120 * t122;
t163 = t122 * t123;
t161 = t122 * t129;
t158 = t125 * t129;
t156 = t123 * pkin(1) + pkin(8) * t164;
t155 = t125 * pkin(8) + qJ(1);
t126 = sin(qJ(5));
t152 = pkin(5) * t126 + pkin(10);
t151 = t121 * t167;
t150 = t124 * t167;
t149 = t122 * t151;
t148 = t120 * pkin(1) - pkin(8) * t163;
t147 = g(1) * t120 - g(2) * t123;
t102 = -t120 * t129 + t123 * t157;
t146 = t102 * t121 + t123 * t162;
t127 = sin(qJ(4));
t103 = t120 * t131 + t123 * t158;
t128 = sin(qJ(3));
t87 = t103 * t167 + (t102 * t124 - t121 * t163) * t128;
t80 = t87 * t127 + t146 * t166;
t105 = -t120 * t158 + t123 * t131;
t89 = t105 * t167 + (t104 * t124 + t121 * t164) * t128;
t82 = t89 * t127 - t145 * t166;
t96 = t125 * t121 * t128 + (t124 * t128 * t131 + t167 * t129) * t122;
t90 = t96 * t127 + t144 * t166;
t143 = g(1) * t82 + g(2) * t80 + g(3) * t90;
t86 = -t102 * t150 + t103 * t128 + t123 * t149;
t88 = -t104 * t150 + t105 * t128 - t120 * t149;
t95 = -t125 * t151 + t128 * t161 - t150 * t160;
t142 = g(1) * t88 + g(2) * t86 + g(3) * t95;
t141 = t105 * pkin(2) + t145 * pkin(9) + t156;
t140 = t89 * pkin(3) + t141;
t139 = pkin(2) * t161 - t144 * pkin(9) + t155;
t138 = t96 * pkin(3) + t139;
t137 = t88 * pkin(10) + t140;
t136 = t95 * pkin(10) + t138;
t135 = t103 * pkin(2) - t146 * pkin(9) + t148;
t134 = t87 * pkin(3) + t135;
t133 = t86 * pkin(10) + t134;
t132 = -pkin(12) - pkin(11);
t130 = cos(qJ(5));
t119 = qJ(5) + qJ(6);
t115 = cos(t119);
t114 = sin(t119);
t113 = t130 * pkin(5) + pkin(4);
t91 = -t144 * t127 + t96 * t166;
t83 = t145 * t127 + t89 * t166;
t81 = -t146 * t127 + t87 * t166;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t123 - g(2) * t120, t147, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t103 - g(3) * t161, -g(1) * t104 - g(2) * t102 - g(3) * t160, -g(3) * t125 - t147 * t122, -g(1) * t156 - g(2) * t148 - g(3) * t155, 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 - g(3) * t96, t142, -g(1) * t145 + g(2) * t146 + g(3) * t144, -g(1) * t141 - g(2) * t135 - g(3) * t139, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t91, t143, -t142, -g(1) * t137 - g(2) * t133 - g(3) * t136, 0, 0, 0, 0, 0, 0, -g(1) * (t88 * t126 + t83 * t130) - g(2) * (t86 * t126 + t81 * t130) - g(3) * (t95 * t126 + t91 * t130) -g(1) * (-t83 * t126 + t88 * t130) - g(2) * (-t81 * t126 + t86 * t130) - g(3) * (-t91 * t126 + t95 * t130) -t143, -g(1) * (t83 * pkin(4) + t82 * pkin(11) + t137) - g(2) * (t81 * pkin(4) + t80 * pkin(11) + t133) - g(3) * (t91 * pkin(4) + t90 * pkin(11) + t136) 0, 0, 0, 0, 0, 0, -g(1) * (t88 * t114 + t83 * t115) - g(2) * (t86 * t114 + t81 * t115) - g(3) * (t95 * t114 + t91 * t115) -g(1) * (-t83 * t114 + t88 * t115) - g(2) * (-t81 * t114 + t86 * t115) - g(3) * (-t91 * t114 + t95 * t115) -t143, -g(1) * (t83 * t113 - t82 * t132 + t152 * t88 + t140) - g(2) * (t81 * t113 - t80 * t132 + t152 * t86 + t134) - g(3) * (t91 * t113 - t90 * t132 + t152 * t95 + t138);];
U_reg  = t1;

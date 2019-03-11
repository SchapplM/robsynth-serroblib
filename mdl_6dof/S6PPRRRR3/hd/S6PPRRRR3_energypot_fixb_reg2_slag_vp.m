% Calculate inertial parameters regressor of potential energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:39
% EndTime: 2019-03-08 19:11:40
% DurationCPUTime: 0.45s
% Computational Cost: add. (927->112), mult. (2553->183), div. (0->0), fcn. (3324->18), ass. (0->75)
t136 = sin(pkin(7));
t141 = cos(pkin(7));
t142 = cos(pkin(6));
t137 = sin(pkin(6));
t138 = cos(pkin(14));
t176 = t137 * t138;
t121 = -t136 * t176 + t142 * t141;
t133 = sin(pkin(14));
t139 = cos(pkin(13));
t134 = sin(pkin(13));
t178 = t134 * t142;
t124 = -t139 * t133 - t138 * t178;
t174 = t137 * t141;
t116 = -t124 * t136 + t134 * t174;
t146 = sin(qJ(3));
t149 = cos(qJ(3));
t173 = t138 * t141;
t177 = t136 * t142;
t113 = t149 * t177 + (-t133 * t146 + t149 * t173) * t137;
t135 = sin(pkin(8));
t140 = cos(pkin(8));
t106 = -t113 * t135 + t121 * t140;
t125 = -t133 * t178 + t139 * t138;
t179 = t134 * t137;
t161 = t124 * t141 + t136 * t179;
t104 = -t125 * t146 + t161 * t149;
t96 = -t104 * t135 + t116 * t140;
t172 = t139 * t142;
t123 = t133 * t172 + t134 * t138;
t122 = -t134 * t133 + t138 * t172;
t175 = t137 * t139;
t162 = t122 * t141 - t136 * t175;
t102 = -t123 * t146 + t162 * t149;
t115 = -t122 * t136 - t139 * t174;
t95 = -t102 * t135 + t115 * t140;
t188 = cos(qJ(4));
t180 = t133 * t137;
t170 = t142 * qJ(2) + qJ(1);
t169 = t139 * pkin(1) + qJ(2) * t179;
t166 = t135 * t188;
t165 = t140 * t188;
t164 = g(1) * t134 - g(2) * t139;
t163 = t134 * pkin(1) - qJ(2) * t175;
t144 = sin(qJ(5));
t148 = cos(qJ(5));
t103 = t123 * t149 + t162 * t146;
t145 = sin(qJ(4));
t85 = t103 * t188 + (t102 * t140 + t115 * t135) * t145;
t78 = t85 * t144 - t95 * t148;
t105 = t125 * t149 + t161 * t146;
t87 = t105 * t188 + (t104 * t140 + t116 * t135) * t145;
t80 = t87 * t144 - t96 * t148;
t114 = t146 * t177 + (t133 * t149 + t146 * t173) * t137;
t94 = t114 * t188 + (t113 * t140 + t121 * t135) * t145;
t88 = -t106 * t148 + t94 * t144;
t160 = g(1) * t80 + g(2) * t78 + g(3) * t88;
t84 = -t102 * t165 + t103 * t145 - t115 * t166;
t86 = -t104 * t165 + t105 * t145 - t116 * t166;
t93 = -t113 * t165 + t114 * t145 - t121 * t166;
t159 = g(1) * t86 + g(2) * t84 + g(3) * t93;
t158 = t125 * pkin(2) + t116 * pkin(9) + t169;
t157 = pkin(2) * t180 + t121 * pkin(9) + t170;
t156 = t105 * pkin(3) + pkin(10) * t96 + t158;
t155 = t123 * pkin(2) + t115 * pkin(9) + t163;
t154 = t114 * pkin(3) + pkin(10) * t106 + t157;
t153 = t87 * pkin(4) + t86 * pkin(11) + t156;
t152 = t94 * pkin(4) + t93 * pkin(11) + t154;
t151 = t103 * pkin(3) + pkin(10) * t95 + t155;
t150 = t85 * pkin(4) + t84 * pkin(11) + t151;
t147 = cos(qJ(6));
t143 = sin(qJ(6));
t89 = t106 * t144 + t94 * t148;
t81 = t96 * t144 + t87 * t148;
t79 = t95 * t144 + t85 * t148;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t139 - g(2) * t134, t164, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t123 - g(3) * t180, -g(1) * t124 - g(2) * t122 - g(3) * t176, -g(3) * t142 - t164 * t137, -g(1) * t169 - g(2) * t163 - g(3) * t170, 0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t103 - g(3) * t114, -g(1) * t104 - g(2) * t102 - g(3) * t113, -g(1) * t116 - g(2) * t115 - g(3) * t121, -g(1) * t158 - g(2) * t155 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t85 - g(3) * t94, t159, -g(1) * t96 - g(2) * t95 - g(3) * t106, -g(1) * t156 - g(2) * t151 - g(3) * t154, 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t89, t160, -t159, -g(1) * t153 - g(2) * t150 - g(3) * t152, 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t143 + t81 * t147) - g(2) * (t84 * t143 + t79 * t147) - g(3) * (t93 * t143 + t89 * t147) -g(1) * (-t81 * t143 + t86 * t147) - g(2) * (-t79 * t143 + t84 * t147) - g(3) * (-t89 * t143 + t93 * t147) -t160, -g(1) * (t81 * pkin(5) + t80 * pkin(12) + t153) - g(2) * (t79 * pkin(5) + t78 * pkin(12) + t150) - g(3) * (t89 * pkin(5) + t88 * pkin(12) + t152);];
U_reg  = t1;

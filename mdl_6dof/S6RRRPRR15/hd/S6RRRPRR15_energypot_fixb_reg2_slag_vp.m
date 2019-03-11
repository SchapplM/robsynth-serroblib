% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR15_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:38:38
% EndTime: 2019-03-09 20:38:39
% DurationCPUTime: 0.32s
% Computational Cost: add. (453->95), mult. (1177->145), div. (0->0), fcn. (1483->14), ass. (0->64)
t139 = cos(pkin(6));
t143 = sin(qJ(2));
t148 = cos(qJ(1));
t170 = t148 * t143;
t144 = sin(qJ(1));
t147 = cos(qJ(2));
t171 = t144 * t147;
t123 = -t139 * t171 - t170;
t136 = sin(pkin(7));
t138 = cos(pkin(7));
t137 = sin(pkin(6));
t177 = t137 * t144;
t114 = -t123 * t136 + t138 * t177;
t176 = t137 * t147;
t120 = -t136 * t176 + t139 * t138;
t182 = cos(qJ(3));
t181 = t139 * pkin(9) + pkin(8);
t142 = sin(qJ(3));
t179 = t136 * t142;
t178 = t137 * t143;
t175 = t137 * t148;
t174 = t138 * t142;
t172 = t144 * t143;
t169 = t148 * t147;
t168 = t148 * pkin(1) + pkin(9) * t177;
t165 = t136 * t182;
t164 = t138 * t182;
t163 = t137 * t165;
t162 = t144 * pkin(1) - pkin(9) * t175;
t161 = g(1) * t144 - g(2) * t148;
t121 = t139 * t169 - t172;
t113 = -t121 * t136 - t138 * t175;
t109 = -t139 * t165 + t142 * t178 - t164 * t176;
t141 = sin(qJ(5));
t146 = cos(qJ(5));
t101 = -t109 * t146 + t120 * t141;
t122 = t139 * t170 + t171;
t103 = -t121 * t164 + t122 * t142 + t148 * t163;
t92 = -t103 * t146 + t113 * t141;
t124 = -t139 * t172 + t169;
t105 = -t123 * t164 + t124 * t142 - t144 * t163;
t94 = -t105 * t146 + t114 * t141;
t160 = g(1) * t94 + g(2) * t92 + g(3) * t101;
t159 = g(1) * t105 + g(2) * t103 + g(3) * t109;
t104 = t121 * t174 + t122 * t182 - t175 * t179;
t106 = t124 * t182 + (t123 * t138 + t136 * t177) * t142;
t110 = t139 * t179 + (t182 * t143 + t147 * t174) * t137;
t158 = g(1) * t106 + g(2) * t104 + g(3) * t110;
t157 = t124 * pkin(2) + t114 * pkin(10) + t168;
t156 = pkin(2) * t178 + t120 * pkin(10) + t181;
t155 = t106 * pkin(3) + t105 * qJ(4) + t157;
t154 = t110 * pkin(3) + t109 * qJ(4) + t156;
t153 = t122 * pkin(2) + t113 * pkin(10) + t162;
t152 = t114 * pkin(4) + t106 * pkin(11) + t155;
t151 = t120 * pkin(4) + t110 * pkin(11) + t154;
t150 = t104 * pkin(3) + t103 * qJ(4) + t153;
t149 = t113 * pkin(4) + t104 * pkin(11) + t150;
t145 = cos(qJ(6));
t140 = sin(qJ(6));
t102 = t109 * t141 + t120 * t146;
t96 = -g(1) * t114 - g(2) * t113 - g(3) * t120;
t95 = t105 * t141 + t114 * t146;
t93 = t103 * t141 + t113 * t146;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t148 - g(2) * t144, t161, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t124 - g(2) * t122 - g(3) * t178, -g(1) * t123 - g(2) * t121 - g(3) * t176, -g(3) * t139 - t161 * t137, -g(1) * t168 - g(2) * t162 - g(3) * t181, 0, 0, 0, 0, 0, 0, -t158, t159, t96, -g(1) * t157 - g(2) * t153 - g(3) * t156, 0, 0, 0, 0, 0, 0, t96, t158, -t159, -g(1) * t155 - g(2) * t150 - g(3) * t154, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t102, t160, -t158, -g(1) * t152 - g(2) * t149 - g(3) * t151, 0, 0, 0, 0, 0, 0, -g(1) * (t106 * t140 + t95 * t145) - g(2) * (t104 * t140 + t93 * t145) - g(3) * (t102 * t145 + t110 * t140) -g(1) * (t106 * t145 - t95 * t140) - g(2) * (t104 * t145 - t93 * t140) - g(3) * (-t102 * t140 + t110 * t145) -t160, -g(1) * (t95 * pkin(5) + t94 * pkin(12) + t152) - g(2) * (t93 * pkin(5) + t92 * pkin(12) + t149) - g(3) * (t102 * pkin(5) + t101 * pkin(12) + t151);];
U_reg  = t1;

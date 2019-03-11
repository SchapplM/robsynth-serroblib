% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR15_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:50:11
% EndTime: 2019-03-10 00:50:11
% DurationCPUTime: 0.33s
% Computational Cost: add. (508->100), mult. (1330->146), div. (0->0), fcn. (1691->14), ass. (0->65)
t131 = sin(pkin(7));
t133 = cos(pkin(7));
t134 = cos(pkin(6));
t132 = sin(pkin(6));
t141 = cos(qJ(2));
t172 = t132 * t141;
t155 = t131 * t172 - t133 * t134;
t139 = sin(qJ(1));
t167 = t139 * t141;
t138 = sin(qJ(2));
t142 = cos(qJ(1));
t169 = t138 * t142;
t119 = -t134 * t167 - t169;
t173 = t132 * t139;
t156 = -t119 * t131 + t133 * t173;
t182 = pkin(5) + pkin(11);
t181 = cos(qJ(3));
t180 = cos(qJ(4));
t166 = t141 * t142;
t168 = t139 * t138;
t120 = -t134 * t168 + t166;
t137 = sin(qJ(3));
t162 = t131 * t181;
t160 = t132 * t162;
t161 = t133 * t181;
t105 = -t119 * t161 + t120 * t137 - t139 * t160;
t179 = pkin(11) * t105;
t174 = t132 * t138;
t110 = -t134 * t162 + t137 * t174 - t161 * t172;
t178 = pkin(11) * t110;
t117 = t134 * t166 - t168;
t118 = t134 * t169 + t167;
t103 = -t117 * t161 + t118 * t137 + t142 * t160;
t177 = t103 * pkin(11);
t176 = t134 * pkin(9) + pkin(8);
t171 = t132 * t142;
t165 = t142 * pkin(1) + pkin(9) * t173;
t159 = t139 * pkin(1) - pkin(9) * t171;
t158 = g(1) * t139 - g(2) * t142;
t157 = t117 * t131 + t133 * t171;
t111 = t134 * t131 * t137 + (t133 * t137 * t141 + t181 * t138) * t132;
t136 = sin(qJ(4));
t101 = t111 * t136 + t155 * t180;
t104 = t118 * t181 + (t117 * t133 - t131 * t171) * t137;
t94 = t104 * t136 + t157 * t180;
t106 = t120 * t181 + (t119 * t133 + t131 * t173) * t137;
t96 = t106 * t136 - t156 * t180;
t154 = g(1) * t96 + g(2) * t94 + g(3) * t101;
t102 = t111 * t180 - t155 * t136;
t95 = t104 * t180 - t157 * t136;
t97 = t106 * t180 + t156 * t136;
t153 = g(1) * t97 + g(2) * t95 + g(3) * t102;
t152 = g(1) * t105 + g(2) * t103 + g(3) * t110;
t151 = t120 * pkin(2) + t156 * pkin(10) + t165;
t150 = pkin(2) * t174 - t155 * pkin(10) + t176;
t149 = t106 * pkin(3) + t151;
t148 = t111 * pkin(3) + t150;
t147 = t97 * pkin(4) + qJ(5) * t96 + t149;
t146 = t102 * pkin(4) + qJ(5) * t101 + t148;
t145 = t118 * pkin(2) - t157 * pkin(10) + t159;
t144 = t104 * pkin(3) + t145;
t143 = t95 * pkin(4) + t94 * qJ(5) + t144;
t140 = cos(qJ(6));
t135 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t142 - g(2) * t139, t158, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t120 - g(2) * t118 - g(3) * t174, -g(1) * t119 - g(2) * t117 - g(3) * t172, -g(3) * t134 - t158 * t132, -g(1) * t165 - g(2) * t159 - g(3) * t176, 0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t104 - g(3) * t111, t152, -g(1) * t156 + g(2) * t157 + g(3) * t155, -g(1) * t151 - g(2) * t145 - g(3) * t150, 0, 0, 0, 0, 0, 0, -t153, t154, -t152, -g(1) * (t149 + t179) - g(2) * (t144 + t177) - g(3) * (t148 + t178) 0, 0, 0, 0, 0, 0, -t152, t153, -t154, -g(1) * (t147 + t179) - g(2) * (t143 + t177) - g(3) * (t146 + t178) 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t140 + t135 * t96) - g(2) * (t103 * t140 + t135 * t94) - g(3) * (t101 * t135 + t110 * t140) -g(1) * (-t105 * t135 + t140 * t96) - g(2) * (-t103 * t135 + t140 * t94) - g(3) * (t101 * t140 - t110 * t135) -t153, -g(1) * (pkin(12) * t97 + t182 * t105 + t147) - g(2) * (t95 * pkin(12) + t182 * t103 + t143) - g(3) * (pkin(12) * t102 + t182 * t110 + t146);];
U_reg  = t1;

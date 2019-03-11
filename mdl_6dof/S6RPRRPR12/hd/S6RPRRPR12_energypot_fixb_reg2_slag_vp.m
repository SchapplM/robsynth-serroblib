% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:52:56
% EndTime: 2019-03-09 05:52:56
% DurationCPUTime: 0.33s
% Computational Cost: add. (508->100), mult. (1330->146), div. (0->0), fcn. (1691->14), ass. (0->65)
t136 = cos(pkin(6));
t131 = sin(pkin(12));
t142 = cos(qJ(1));
t167 = t142 * t131;
t134 = cos(pkin(12));
t140 = sin(qJ(1));
t168 = t140 * t134;
t119 = -t136 * t168 - t167;
t132 = sin(pkin(7));
t135 = cos(pkin(7));
t133 = sin(pkin(6));
t172 = t133 * t140;
t156 = -t119 * t132 + t135 * t172;
t173 = t133 * t134;
t155 = t132 * t173 - t136 * t135;
t182 = pkin(5) + pkin(10);
t181 = cos(qJ(3));
t180 = cos(qJ(4));
t166 = t142 * t134;
t169 = t140 * t131;
t117 = t136 * t166 - t169;
t118 = t136 * t167 + t168;
t139 = sin(qJ(3));
t162 = t132 * t181;
t160 = t133 * t162;
t161 = t135 * t181;
t103 = -t117 * t161 + t118 * t139 + t142 * t160;
t179 = t103 * pkin(10);
t120 = -t136 * t169 + t166;
t105 = -t119 * t161 + t120 * t139 - t140 * t160;
t178 = t105 * pkin(10);
t174 = t131 * t133;
t110 = -t136 * t162 + t139 * t174 - t161 * t173;
t177 = t110 * pkin(10);
t176 = t136 * qJ(2) + pkin(8);
t171 = t133 * t142;
t165 = t142 * pkin(1) + qJ(2) * t172;
t159 = g(1) * t140 - g(2) * t142;
t158 = t140 * pkin(1) - qJ(2) * t171;
t157 = t117 * t132 + t135 * t171;
t111 = t136 * t132 * t139 + (t134 * t135 * t139 + t181 * t131) * t133;
t138 = sin(qJ(4));
t101 = t111 * t138 + t155 * t180;
t104 = t118 * t181 + (t117 * t135 - t132 * t171) * t139;
t94 = t104 * t138 + t157 * t180;
t106 = t120 * t181 + (t119 * t135 + t132 * t172) * t139;
t96 = t106 * t138 - t156 * t180;
t154 = g(1) * t96 + g(2) * t94 + g(3) * t101;
t102 = t111 * t180 - t155 * t138;
t95 = t104 * t180 - t157 * t138;
t97 = t106 * t180 + t156 * t138;
t153 = g(1) * t97 + g(2) * t95 + g(3) * t102;
t152 = g(1) * t105 + g(2) * t103 + g(3) * t110;
t151 = t120 * pkin(2) + t156 * pkin(9) + t165;
t150 = pkin(2) * t174 - t155 * pkin(9) + t176;
t149 = t106 * pkin(3) + t151;
t148 = t111 * pkin(3) + t150;
t147 = t97 * pkin(4) + t96 * qJ(5) + t149;
t146 = t102 * pkin(4) + t101 * qJ(5) + t148;
t145 = t118 * pkin(2) - t157 * pkin(9) + t158;
t144 = t104 * pkin(3) + t145;
t143 = t95 * pkin(4) + t94 * qJ(5) + t144;
t141 = cos(qJ(6));
t137 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t142 - g(2) * t140, t159, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t120 - g(2) * t118 - g(3) * t174, -g(1) * t119 - g(2) * t117 - g(3) * t173, -g(3) * t136 - t159 * t133, -g(1) * t165 - g(2) * t158 - g(3) * t176, 0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t104 - g(3) * t111, t152, -g(1) * t156 + g(2) * t157 + g(3) * t155, -g(1) * t151 - g(2) * t145 - g(3) * t150, 0, 0, 0, 0, 0, 0, -t153, t154, -t152, -g(1) * (t149 + t178) - g(2) * (t144 + t179) - g(3) * (t148 + t177) 0, 0, 0, 0, 0, 0, -t152, t153, -t154, -g(1) * (t147 + t178) - g(2) * (t143 + t179) - g(3) * (t146 + t177) 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t141 + t96 * t137) - g(2) * (t103 * t141 + t94 * t137) - g(3) * (t101 * t137 + t110 * t141) -g(1) * (-t105 * t137 + t96 * t141) - g(2) * (-t103 * t137 + t94 * t141) - g(3) * (t101 * t141 - t110 * t137) -t153, -g(1) * (t97 * pkin(11) + t182 * t105 + t147) - g(2) * (t95 * pkin(11) + t182 * t103 + t143) - g(3) * (t102 * pkin(11) + t182 * t110 + t146);];
U_reg  = t1;

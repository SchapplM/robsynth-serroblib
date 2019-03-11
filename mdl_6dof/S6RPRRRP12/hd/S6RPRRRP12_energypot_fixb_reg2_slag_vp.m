% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:50:37
% EndTime: 2019-03-09 06:50:37
% DurationCPUTime: 0.33s
% Computational Cost: add. (576->95), mult. (1520->147), div. (0->0), fcn. (1949->14), ass. (0->67)
t138 = sin(pkin(12));
t143 = cos(pkin(6));
t149 = cos(qJ(1));
t141 = cos(pkin(12));
t147 = sin(qJ(1));
t173 = t147 * t141;
t126 = -t138 * t149 - t143 * t173;
t139 = sin(pkin(7));
t142 = cos(pkin(7));
t140 = sin(pkin(6));
t178 = t140 * t147;
t163 = -t126 * t139 + t142 * t178;
t179 = t140 * t141;
t162 = t139 * t179 - t142 * t143;
t184 = cos(qJ(3));
t183 = cos(qJ(4));
t182 = t143 * qJ(2) + pkin(8);
t180 = t138 * t140;
t177 = t140 * t149;
t175 = t143 * t149;
t174 = t147 * t138;
t172 = t149 * pkin(1) + qJ(2) * t178;
t169 = t139 * t184;
t168 = t142 * t184;
t167 = t140 * t169;
t166 = g(1) * t147 - g(2) * t149;
t165 = t147 * pkin(1) - qJ(2) * t177;
t124 = t141 * t175 - t174;
t164 = t124 * t139 + t142 * t177;
t125 = t138 * t175 + t173;
t146 = sin(qJ(3));
t109 = -t124 * t168 + t125 * t146 + t149 * t167;
t144 = sin(qJ(5));
t148 = cos(qJ(5));
t110 = t125 * t184 + (t124 * t142 - t139 * t177) * t146;
t145 = sin(qJ(4));
t99 = t110 * t183 - t164 * t145;
t90 = -t109 * t148 + t144 * t99;
t127 = t141 * t149 - t143 * t174;
t112 = t127 * t184 + (t126 * t142 + t139 * t178) * t146;
t101 = t112 * t183 + t163 * t145;
t111 = -t126 * t168 + t127 * t146 - t147 * t167;
t92 = t101 * t144 - t111 * t148;
t118 = t143 * t139 * t146 + (t141 * t142 * t146 + t184 * t138) * t140;
t108 = t118 * t183 - t162 * t145;
t117 = -t143 * t169 + t146 * t180 - t168 * t179;
t94 = t108 * t144 - t117 * t148;
t161 = g(1) * t92 + g(2) * t90 + g(3) * t94;
t100 = t112 * t145 - t163 * t183;
t107 = t118 * t145 + t162 * t183;
t98 = t110 * t145 + t164 * t183;
t160 = g(1) * t100 + g(2) * t98 + g(3) * t107;
t159 = g(1) * t111 + g(2) * t109 + g(3) * t117;
t158 = t127 * pkin(2) + pkin(9) * t163 + t172;
t157 = pkin(2) * t180 - pkin(9) * t162 + t182;
t156 = t112 * pkin(3) + pkin(10) * t111 + t158;
t155 = t118 * pkin(3) + pkin(10) * t117 + t157;
t154 = t125 * pkin(2) - t164 * pkin(9) + t165;
t153 = t101 * pkin(4) + pkin(11) * t100 + t156;
t152 = t108 * pkin(4) + pkin(11) * t107 + t155;
t151 = t110 * pkin(3) + t109 * pkin(10) + t154;
t150 = t99 * pkin(4) + t98 * pkin(11) + t151;
t95 = t108 * t148 + t117 * t144;
t93 = t101 * t148 + t111 * t144;
t91 = t109 * t144 + t148 * t99;
t88 = -g(1) * t93 - g(2) * t91 - g(3) * t95;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t149 - g(2) * t147, t166, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t127 - g(2) * t125 - g(3) * t180, -g(1) * t126 - g(2) * t124 - g(3) * t179, -g(3) * t143 - t166 * t140, -g(1) * t172 - g(2) * t165 - g(3) * t182, 0, 0, 0, 0, 0, 0, -g(1) * t112 - g(2) * t110 - g(3) * t118, t159, -g(1) * t163 + g(2) * t164 + g(3) * t162, -g(1) * t158 - g(2) * t154 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t99 - g(3) * t108, t160, -t159, -g(1) * t156 - g(2) * t151 - g(3) * t155, 0, 0, 0, 0, 0, 0, t88, t161, -t160, -g(1) * t153 - g(2) * t150 - g(3) * t152, 0, 0, 0, 0, 0, 0, t88, -t160, -t161, -g(1) * (pkin(5) * t93 + qJ(6) * t92 + t153) - g(2) * (t91 * pkin(5) + t90 * qJ(6) + t150) - g(3) * (pkin(5) * t95 + qJ(6) * t94 + t152);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:36
% EndTime: 2019-03-10 03:20:36
% DurationCPUTime: 0.33s
% Computational Cost: add. (576->95), mult. (1520->146), div. (0->0), fcn. (1949->14), ass. (0->68)
t141 = cos(pkin(6));
t145 = sin(qJ(2));
t149 = cos(qJ(1));
t174 = t149 * t145;
t146 = sin(qJ(1));
t148 = cos(qJ(2));
t175 = t146 * t148;
t126 = -t141 * t175 - t174;
t138 = sin(pkin(7));
t140 = cos(pkin(7));
t139 = sin(pkin(6));
t180 = t139 * t146;
t163 = -t126 * t138 + t140 * t180;
t179 = t139 * t148;
t162 = t138 * t179 - t141 * t140;
t185 = cos(qJ(3));
t184 = cos(qJ(4));
t183 = t141 * pkin(9) + pkin(8);
t181 = t139 * t145;
t178 = t139 * t149;
t176 = t146 * t145;
t173 = t149 * t148;
t172 = t149 * pkin(1) + pkin(9) * t180;
t169 = t138 * t185;
t168 = t140 * t185;
t167 = t139 * t169;
t166 = t146 * pkin(1) - pkin(9) * t178;
t165 = g(1) * t146 - g(2) * t149;
t124 = t141 * t173 - t176;
t164 = t124 * t138 + t140 * t178;
t125 = t141 * t174 + t175;
t144 = sin(qJ(3));
t109 = -t124 * t168 + t125 * t144 + t149 * t167;
t142 = sin(qJ(5));
t147 = cos(qJ(5));
t110 = t125 * t185 + (t124 * t140 - t138 * t178) * t144;
t143 = sin(qJ(4));
t99 = t110 * t184 - t164 * t143;
t90 = -t109 * t147 + t99 * t142;
t127 = -t141 * t176 + t173;
t112 = t127 * t185 + (t126 * t140 + t138 * t180) * t144;
t101 = t112 * t184 + t163 * t143;
t111 = -t126 * t168 + t127 * t144 - t146 * t167;
t92 = t101 * t142 - t111 * t147;
t118 = t141 * t138 * t144 + (t140 * t144 * t148 + t185 * t145) * t139;
t108 = t118 * t184 - t162 * t143;
t117 = -t141 * t169 + t144 * t181 - t168 * t179;
t94 = t108 * t142 - t117 * t147;
t161 = g(1) * t92 + g(2) * t90 + g(3) * t94;
t100 = t112 * t143 - t163 * t184;
t107 = t118 * t143 + t162 * t184;
t98 = t110 * t143 + t164 * t184;
t160 = g(1) * t100 + g(2) * t98 + g(3) * t107;
t159 = g(1) * t111 + g(2) * t109 + g(3) * t117;
t158 = t127 * pkin(2) + t163 * pkin(10) + t172;
t157 = pkin(2) * t181 - t162 * pkin(10) + t183;
t156 = t112 * pkin(3) + t111 * pkin(11) + t158;
t155 = t118 * pkin(3) + t117 * pkin(11) + t157;
t154 = t125 * pkin(2) - t164 * pkin(10) + t166;
t153 = t101 * pkin(4) + t100 * pkin(12) + t156;
t152 = t108 * pkin(4) + t107 * pkin(12) + t155;
t151 = t110 * pkin(3) + t109 * pkin(11) + t154;
t150 = t99 * pkin(4) + t98 * pkin(12) + t151;
t95 = t108 * t147 + t117 * t142;
t93 = t101 * t147 + t111 * t142;
t91 = t109 * t142 + t99 * t147;
t88 = -g(1) * t93 - g(2) * t91 - g(3) * t95;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t149 - g(2) * t146, t165, -g(3), -g(3) * pkin(8), 0, 0, 0, 0, 0, 0, -g(1) * t127 - g(2) * t125 - g(3) * t181, -g(1) * t126 - g(2) * t124 - g(3) * t179, -g(3) * t141 - t165 * t139, -g(1) * t172 - g(2) * t166 - g(3) * t183, 0, 0, 0, 0, 0, 0, -g(1) * t112 - g(2) * t110 - g(3) * t118, t159, -g(1) * t163 + g(2) * t164 + g(3) * t162, -g(1) * t158 - g(2) * t154 - g(3) * t157, 0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t99 - g(3) * t108, t160, -t159, -g(1) * t156 - g(2) * t151 - g(3) * t155, 0, 0, 0, 0, 0, 0, t88, t161, -t160, -g(1) * t153 - g(2) * t150 - g(3) * t152, 0, 0, 0, 0, 0, 0, t88, -t160, -t161, -g(1) * (t93 * pkin(5) + t92 * qJ(6) + t153) - g(2) * (t91 * pkin(5) + t90 * qJ(6) + t150) - g(3) * (t95 * pkin(5) + t94 * qJ(6) + t152);];
U_reg  = t1;

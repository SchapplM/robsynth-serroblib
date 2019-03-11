% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:25
% EndTime: 2019-03-10 02:27:25
% DurationCPUTime: 0.20s
% Computational Cost: add. (322->92), mult. (648->129), div. (0->0), fcn. (793->12), ass. (0->58)
t118 = cos(pkin(6));
t152 = t118 * pkin(8) + pkin(7);
t125 = cos(qJ(2));
t126 = cos(qJ(1));
t141 = t126 * t125;
t121 = sin(qJ(2));
t122 = sin(qJ(1));
t144 = t122 * t121;
t101 = -t118 * t141 + t144;
t119 = sin(qJ(4));
t151 = t101 * t119;
t142 = t126 * t121;
t143 = t122 * t125;
t103 = t118 * t143 + t142;
t150 = t103 * t119;
t117 = sin(pkin(6));
t149 = t117 * t121;
t148 = t117 * t122;
t124 = cos(qJ(3));
t147 = t117 * t124;
t146 = t117 * t125;
t145 = t117 * t126;
t140 = t126 * pkin(1) + pkin(8) * t148;
t139 = pkin(2) * t149 + t152;
t138 = t122 * pkin(1) - pkin(8) * t145;
t137 = g(1) * t122 - g(2) * t126;
t104 = -t118 * t144 + t141;
t136 = t104 * pkin(2) + t103 * pkin(9) + t140;
t135 = -pkin(9) * t146 + t139;
t116 = qJ(4) + qJ(5);
t111 = sin(t116);
t112 = cos(t116);
t102 = t118 * t142 + t143;
t120 = sin(qJ(3));
t88 = t102 * t124 - t120 * t145;
t77 = -t101 * t112 + t88 * t111;
t90 = t104 * t124 + t120 * t148;
t79 = -t103 * t112 + t90 * t111;
t100 = t118 * t120 + t121 * t147;
t83 = t100 * t111 + t112 * t146;
t134 = g(1) * t79 + g(2) * t77 + g(3) * t83;
t87 = t102 * t120 + t124 * t145;
t89 = t104 * t120 - t122 * t147;
t99 = -t118 * t124 + t120 * t149;
t133 = g(1) * t89 + g(2) * t87 + g(3) * t99;
t132 = t102 * pkin(2) + t101 * pkin(9) + t138;
t123 = cos(qJ(4));
t110 = t123 * pkin(4) + pkin(3);
t127 = -pkin(11) - pkin(10);
t131 = pkin(4) * t150 + t90 * t110 - t89 * t127 + t136;
t130 = -g(1) * t103 - g(2) * t101 + g(3) * t146;
t129 = pkin(4) * t151 + t88 * t110 - t87 * t127 + t132;
t128 = t100 * t110 - t99 * t127 + (-pkin(4) * t119 - pkin(9)) * t146 + t139;
t84 = t100 * t112 - t111 * t146;
t80 = t103 * t111 + t90 * t112;
t78 = t101 * t111 + t88 * t112;
t75 = -g(1) * t80 - g(2) * t78 - g(3) * t84;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t122, t137, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t149, -t130, -g(3) * t118 - t137 * t117, -g(1) * t140 - g(2) * t138 - g(3) * t152, 0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t88 - g(3) * t100, t133, t130, -g(1) * t136 - g(2) * t132 - g(3) * t135, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t123 + t150) - g(2) * (t88 * t123 + t151) - g(3) * (t100 * t123 - t119 * t146) -g(1) * (t103 * t123 - t90 * t119) - g(2) * (t101 * t123 - t88 * t119) - g(3) * (-t100 * t119 - t123 * t146) -t133, -g(1) * (t90 * pkin(3) + t89 * pkin(10) + t136) - g(2) * (t88 * pkin(3) + t87 * pkin(10) + t132) - g(3) * (t100 * pkin(3) + t99 * pkin(10) + t135) 0, 0, 0, 0, 0, 0, t75, t134, -t133, -g(1) * t131 - g(2) * t129 - g(3) * t128, 0, 0, 0, 0, 0, 0, t75, -t133, -t134, -g(1) * (t80 * pkin(5) + t79 * qJ(6) + t131) - g(2) * (t78 * pkin(5) + t77 * qJ(6) + t129) - g(3) * (t84 * pkin(5) + t83 * qJ(6) + t128);];
U_reg  = t1;

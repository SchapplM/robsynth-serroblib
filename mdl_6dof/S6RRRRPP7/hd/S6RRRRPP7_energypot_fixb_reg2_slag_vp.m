% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:26:30
% EndTime: 2019-03-09 21:26:30
% DurationCPUTime: 0.21s
% Computational Cost: add. (322->92), mult. (648->129), div. (0->0), fcn. (793->12), ass. (0->58)
t117 = cos(pkin(6));
t151 = t117 * pkin(8) + pkin(7);
t125 = cos(qJ(2));
t126 = cos(qJ(1));
t140 = t126 * t125;
t121 = sin(qJ(2));
t122 = sin(qJ(1));
t143 = t122 * t121;
t100 = -t117 * t140 + t143;
t119 = sin(qJ(4));
t150 = t100 * t119;
t141 = t126 * t121;
t142 = t122 * t125;
t102 = t117 * t142 + t141;
t149 = t102 * t119;
t116 = sin(pkin(6));
t148 = t116 * t121;
t147 = t116 * t122;
t124 = cos(qJ(3));
t146 = t116 * t124;
t145 = t116 * t125;
t144 = t116 * t126;
t139 = t126 * pkin(1) + pkin(8) * t147;
t138 = pkin(2) * t148 + t151;
t137 = t122 * pkin(1) - pkin(8) * t144;
t136 = g(1) * t122 - g(2) * t126;
t103 = -t117 * t143 + t140;
t135 = t103 * pkin(2) + t102 * pkin(9) + t139;
t134 = -pkin(9) * t145 + t138;
t115 = qJ(4) + pkin(11);
t110 = sin(t115);
t111 = cos(t115);
t101 = t117 * t141 + t142;
t120 = sin(qJ(3));
t87 = t101 * t124 - t120 * t144;
t76 = -t100 * t111 + t87 * t110;
t89 = t103 * t124 + t120 * t147;
t78 = -t102 * t111 + t89 * t110;
t99 = t117 * t120 + t121 * t146;
t82 = t99 * t110 + t111 * t145;
t133 = g(1) * t78 + g(2) * t76 + g(3) * t82;
t86 = t101 * t120 + t124 * t144;
t88 = t103 * t120 - t122 * t146;
t98 = -t117 * t124 + t120 * t148;
t132 = g(1) * t88 + g(2) * t86 + g(3) * t98;
t131 = t101 * pkin(2) + t100 * pkin(9) + t137;
t123 = cos(qJ(4));
t109 = t123 * pkin(4) + pkin(3);
t118 = -qJ(5) - pkin(10);
t130 = pkin(4) * t149 + t89 * t109 - t88 * t118 + t135;
t129 = -g(1) * t102 - g(2) * t100 + g(3) * t145;
t128 = pkin(4) * t150 + t87 * t109 - t86 * t118 + t131;
t127 = t99 * t109 - t98 * t118 + (-pkin(4) * t119 - pkin(9)) * t145 + t138;
t83 = -t110 * t145 + t99 * t111;
t79 = t102 * t110 + t89 * t111;
t77 = t100 * t110 + t87 * t111;
t74 = -g(1) * t79 - g(2) * t77 - g(3) * t83;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t122, t136, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t101 - g(3) * t148, -t129, -g(3) * t117 - t136 * t116, -g(1) * t139 - g(2) * t137 - g(3) * t151, 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 - g(3) * t99, t132, t129, -g(1) * t135 - g(2) * t131 - g(3) * t134, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t123 + t149) - g(2) * (t87 * t123 + t150) - g(3) * (-t119 * t145 + t99 * t123) -g(1) * (t102 * t123 - t89 * t119) - g(2) * (t100 * t123 - t87 * t119) - g(3) * (-t99 * t119 - t123 * t145) -t132, -g(1) * (t89 * pkin(3) + t88 * pkin(10) + t135) - g(2) * (t87 * pkin(3) + t86 * pkin(10) + t131) - g(3) * (t99 * pkin(3) + t98 * pkin(10) + t134) 0, 0, 0, 0, 0, 0, t74, t133, -t132, -g(1) * t130 - g(2) * t128 - g(3) * t127, 0, 0, 0, 0, 0, 0, t74, -t132, -t133, -g(1) * (t79 * pkin(5) + t78 * qJ(6) + t130) - g(2) * (t77 * pkin(5) + t76 * qJ(6) + t128) - g(3) * (t83 * pkin(5) + t82 * qJ(6) + t127);];
U_reg  = t1;

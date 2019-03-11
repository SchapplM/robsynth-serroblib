% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:13
% EndTime: 2019-03-09 19:34:14
% DurationCPUTime: 0.22s
% Computational Cost: add. (310->87), mult. (725->126), div. (0->0), fcn. (900->12), ass. (0->56)
t153 = pkin(9) - pkin(10);
t152 = cos(qJ(3));
t121 = sin(qJ(2));
t122 = sin(qJ(1));
t125 = cos(qJ(2));
t126 = cos(qJ(1));
t148 = cos(pkin(6));
t139 = t126 * t148;
t104 = t122 * t121 - t125 * t139;
t151 = t104 * pkin(9);
t140 = t122 * t148;
t106 = t126 * t121 + t125 * t140;
t150 = t106 * pkin(9);
t149 = t148 * pkin(8) + pkin(7);
t117 = sin(pkin(6));
t147 = t117 * t121;
t146 = t117 * t122;
t145 = t117 * t125;
t144 = t117 * t126;
t143 = t126 * pkin(1) + pkin(8) * t146;
t107 = -t121 * t140 + t126 * t125;
t142 = t107 * pkin(2) + t143;
t141 = t117 * t152;
t138 = t122 * pkin(1) - pkin(8) * t144;
t137 = g(1) * t122 - g(2) * t126;
t105 = t121 * t139 + t122 * t125;
t136 = t105 * pkin(2) + t138;
t135 = pkin(2) * t147 - pkin(9) * t145 + t149;
t120 = sin(qJ(3));
t95 = t107 * t120 - t122 * t141;
t96 = t107 * t152 + t120 * t146;
t134 = t96 * pkin(3) + t95 * qJ(4) + t142;
t119 = sin(qJ(5));
t124 = cos(qJ(5));
t93 = t105 * t120 + t126 * t141;
t94 = t105 * t152 - t120 * t144;
t79 = t94 * t119 - t93 * t124;
t81 = t96 * t119 - t95 * t124;
t102 = t120 * t147 - t148 * t152;
t103 = t148 * t120 + t121 * t141;
t85 = -t102 * t124 + t103 * t119;
t133 = g(1) * t81 + g(2) * t79 + g(3) * t85;
t132 = g(1) * t95 + g(2) * t93 + g(3) * t102;
t83 = -g(1) * t106 - g(2) * t104 + g(3) * t145;
t131 = t94 * pkin(3) + t93 * qJ(4) + t136;
t130 = t103 * pkin(3) + t102 * qJ(4) + t135;
t129 = t96 * pkin(4) + t153 * t106 + t134;
t128 = t103 * pkin(4) + pkin(10) * t145 + t130;
t127 = t94 * pkin(4) + t153 * t104 + t131;
t123 = cos(qJ(6));
t118 = sin(qJ(6));
t86 = t102 * t119 + t103 * t124;
t82 = t95 * t119 + t96 * t124;
t80 = t93 * t119 + t94 * t124;
t78 = -g(1) * t96 - g(2) * t94 - g(3) * t103;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t126 - g(2) * t122, t137, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t107 - g(2) * t105 - g(3) * t147, -t83, -g(3) * t148 - t137 * t117, -g(1) * t143 - g(2) * t138 - g(3) * t149, 0, 0, 0, 0, 0, 0, t78, t132, t83, -g(1) * (t142 + t150) - g(2) * (t136 + t151) - g(3) * t135, 0, 0, 0, 0, 0, 0, t78, t83, -t132, -g(1) * (t134 + t150) - g(2) * (t131 + t151) - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t86, t133, -t83, -g(1) * t129 - g(2) * t127 - g(3) * t128, 0, 0, 0, 0, 0, 0, -g(1) * (-t106 * t118 + t82 * t123) - g(2) * (-t104 * t118 + t80 * t123) - g(3) * (t118 * t145 + t86 * t123) -g(1) * (-t106 * t123 - t82 * t118) - g(2) * (-t104 * t123 - t80 * t118) - g(3) * (-t86 * t118 + t123 * t145) -t133, -g(1) * (t82 * pkin(5) + t81 * pkin(11) + t129) - g(2) * (t80 * pkin(5) + t79 * pkin(11) + t127) - g(3) * (t86 * pkin(5) + t85 * pkin(11) + t128);];
U_reg  = t1;

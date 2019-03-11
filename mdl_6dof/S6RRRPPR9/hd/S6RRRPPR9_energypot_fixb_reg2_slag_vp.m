% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:17:53
% EndTime: 2019-03-09 16:17:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (337->88), mult. (795->126), div. (0->0), fcn. (997->12), ass. (0->57)
t152 = cos(qJ(3));
t151 = -pkin(10) + qJ(4);
t147 = cos(pkin(6));
t150 = t147 * pkin(8) + pkin(7);
t120 = sin(qJ(2));
t121 = sin(qJ(1));
t123 = cos(qJ(2));
t124 = cos(qJ(1));
t138 = t124 * t147;
t104 = t120 * t138 + t121 * t123;
t119 = sin(qJ(3));
t116 = sin(pkin(6));
t140 = t116 * t152;
t92 = t104 * t119 + t124 * t140;
t149 = t92 * qJ(4);
t139 = t121 * t147;
t106 = -t120 * t139 + t124 * t123;
t94 = t106 * t119 - t121 * t140;
t148 = t94 * qJ(4);
t145 = t116 * t120;
t101 = t119 * t145 - t147 * t152;
t146 = t101 * qJ(4);
t144 = t116 * t121;
t143 = t116 * t123;
t142 = t116 * t124;
t141 = t124 * pkin(1) + pkin(8) * t144;
t137 = t121 * pkin(1) - pkin(8) * t142;
t136 = g(1) * t121 - g(2) * t124;
t105 = t124 * t120 + t123 * t139;
t135 = t106 * pkin(2) + t105 * pkin(9) + t141;
t134 = pkin(2) * t145 - pkin(9) * t143 + t150;
t95 = t106 * t152 + t119 * t144;
t133 = t95 * pkin(3) + t135;
t103 = t121 * t120 - t123 * t138;
t115 = sin(pkin(11));
t117 = cos(pkin(11));
t93 = t104 * t152 - t119 * t142;
t83 = -t103 * t117 + t93 * t115;
t85 = -t105 * t117 + t95 * t115;
t102 = t147 * t119 + t120 * t140;
t90 = t102 * t115 + t117 * t143;
t132 = g(1) * t85 + g(2) * t83 + g(3) * t90;
t80 = g(1) * t94 + g(2) * t92 + g(3) * t101;
t131 = t102 * pkin(3) + t134;
t130 = t104 * pkin(2) + t103 * pkin(9) + t137;
t129 = -g(1) * t105 - g(2) * t103 + g(3) * t143;
t128 = t93 * pkin(3) + t130;
t86 = t105 * t115 + t95 * t117;
t127 = t86 * pkin(4) + t85 * qJ(5) + t133;
t91 = t102 * t117 - t115 * t143;
t126 = t91 * pkin(4) + t90 * qJ(5) + t131;
t84 = t103 * t115 + t93 * t117;
t125 = t84 * pkin(4) + t83 * qJ(5) + t128;
t122 = cos(qJ(6));
t118 = sin(qJ(6));
t78 = -g(1) * t86 - g(2) * t84 - g(3) * t91;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t124 - g(2) * t121, t136, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t104 - g(3) * t145, -t129, -g(3) * t147 - t136 * t116, -g(1) * t141 - g(2) * t137 - g(3) * t150, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t102, t80, t129, -g(1) * t135 - g(2) * t130 - g(3) * t134, 0, 0, 0, 0, 0, 0, t78, t132, -t80, -g(1) * (t133 + t148) - g(2) * (t128 + t149) - g(3) * (t131 + t146) 0, 0, 0, 0, 0, 0, t78, -t80, -t132, -g(1) * (t127 + t148) - g(2) * (t125 + t149) - g(3) * (t126 + t146) 0, 0, 0, 0, 0, 0, -g(1) * (t85 * t118 + t86 * t122) - g(2) * (t83 * t118 + t84 * t122) - g(3) * (t90 * t118 + t91 * t122) -g(1) * (-t86 * t118 + t85 * t122) - g(2) * (-t84 * t118 + t83 * t122) - g(3) * (-t91 * t118 + t90 * t122) t80, -g(1) * (t86 * pkin(5) + t151 * t94 + t127) - g(2) * (t84 * pkin(5) + t151 * t92 + t125) - g(3) * (t91 * pkin(5) + t151 * t101 + t126);];
U_reg  = t1;

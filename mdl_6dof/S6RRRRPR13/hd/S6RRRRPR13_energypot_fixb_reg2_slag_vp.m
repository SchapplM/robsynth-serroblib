% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:03
% EndTime: 2019-03-10 00:04:03
% DurationCPUTime: 0.20s
% Computational Cost: add. (337->88), mult. (795->126), div. (0->0), fcn. (997->12), ass. (0->57)
t153 = pkin(10) - pkin(11);
t120 = sin(qJ(2));
t121 = sin(qJ(1));
t124 = cos(qJ(2));
t125 = cos(qJ(1));
t147 = cos(pkin(6));
t139 = t125 * t147;
t105 = t120 * t139 + t121 * t124;
t119 = sin(qJ(3));
t116 = sin(pkin(6));
t150 = cos(qJ(3));
t141 = t116 * t150;
t93 = t105 * t119 + t125 * t141;
t152 = t93 * pkin(10);
t140 = t121 * t147;
t107 = -t120 * t140 + t125 * t124;
t95 = t107 * t119 - t121 * t141;
t151 = t95 * pkin(10);
t146 = t116 * t120;
t102 = t119 * t146 - t147 * t150;
t149 = t102 * pkin(10);
t148 = t147 * pkin(8) + pkin(7);
t145 = t116 * t121;
t144 = t116 * t124;
t143 = t116 * t125;
t142 = t125 * pkin(1) + pkin(8) * t145;
t138 = t121 * pkin(1) - pkin(8) * t143;
t137 = g(1) * t121 - g(2) * t125;
t106 = t125 * t120 + t124 * t140;
t136 = t107 * pkin(2) + t106 * pkin(9) + t142;
t135 = pkin(2) * t146 - pkin(9) * t144 + t148;
t96 = t107 * t150 + t119 * t145;
t134 = t96 * pkin(3) + t136;
t104 = t121 * t120 - t124 * t139;
t118 = sin(qJ(4));
t123 = cos(qJ(4));
t94 = t105 * t150 - t119 * t143;
t84 = -t104 * t123 + t94 * t118;
t86 = -t106 * t123 + t96 * t118;
t103 = t147 * t119 + t120 * t141;
t91 = t103 * t118 + t123 * t144;
t133 = g(1) * t86 + g(2) * t84 + g(3) * t91;
t81 = g(1) * t95 + g(2) * t93 + g(3) * t102;
t132 = t103 * pkin(3) + t135;
t131 = t105 * pkin(2) + t104 * pkin(9) + t138;
t130 = -g(1) * t106 - g(2) * t104 + g(3) * t144;
t129 = t94 * pkin(3) + t131;
t87 = t106 * t118 + t96 * t123;
t128 = t87 * pkin(4) + t86 * qJ(5) + t134;
t92 = t103 * t123 - t118 * t144;
t127 = t92 * pkin(4) + t91 * qJ(5) + t132;
t85 = t104 * t118 + t94 * t123;
t126 = t85 * pkin(4) + t84 * qJ(5) + t129;
t122 = cos(qJ(6));
t117 = sin(qJ(6));
t79 = -g(1) * t87 - g(2) * t85 - g(3) * t92;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t121, t137, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t107 - g(2) * t105 - g(3) * t146, -t130, -g(3) * t147 - t137 * t116, -g(1) * t142 - g(2) * t138 - g(3) * t148, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t103, t81, t130, -g(1) * t136 - g(2) * t131 - g(3) * t135, 0, 0, 0, 0, 0, 0, t79, t133, -t81, -g(1) * (t134 + t151) - g(2) * (t129 + t152) - g(3) * (t132 + t149) 0, 0, 0, 0, 0, 0, t79, -t81, -t133, -g(1) * (t128 + t151) - g(2) * (t126 + t152) - g(3) * (t127 + t149) 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t117 + t87 * t122) - g(2) * (t84 * t117 + t85 * t122) - g(3) * (t91 * t117 + t92 * t122) -g(1) * (-t87 * t117 + t86 * t122) - g(2) * (-t85 * t117 + t84 * t122) - g(3) * (-t92 * t117 + t91 * t122) t81, -g(1) * (t87 * pkin(5) + t153 * t95 + t128) - g(2) * (t85 * pkin(5) + t153 * t93 + t126) - g(3) * (t92 * pkin(5) + t153 * t102 + t127);];
U_reg  = t1;

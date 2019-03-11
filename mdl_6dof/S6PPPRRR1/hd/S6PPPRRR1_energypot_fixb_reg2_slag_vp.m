% Calculate inertial parameters regressor of potential energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPPRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:33
% EndTime: 2019-03-08 18:40:34
% DurationCPUTime: 0.47s
% Computational Cost: add. (927->112), mult. (2553->182), div. (0->0), fcn. (3324->18), ass. (0->75)
t118 = sin(pkin(7));
t124 = cos(pkin(7));
t125 = cos(pkin(6));
t119 = sin(pkin(6));
t121 = cos(pkin(13));
t157 = t119 * t121;
t102 = -t118 * t157 + t125 * t124;
t115 = sin(pkin(13));
t122 = cos(pkin(12));
t116 = sin(pkin(12));
t159 = t116 * t125;
t105 = -t122 * t115 - t121 * t159;
t155 = t119 * t124;
t97 = -t105 * t118 + t116 * t155;
t117 = sin(pkin(8));
t123 = cos(pkin(8));
t114 = sin(pkin(14));
t120 = cos(pkin(14));
t154 = t121 * t124;
t158 = t118 * t125;
t94 = t120 * t158 + (-t114 * t115 + t120 * t154) * t119;
t87 = t102 * t123 - t94 * t117;
t106 = -t115 * t159 + t122 * t121;
t160 = t116 * t119;
t142 = t105 * t124 + t118 * t160;
t85 = -t106 * t114 + t142 * t120;
t77 = -t85 * t117 + t97 * t123;
t153 = t122 * t125;
t104 = t115 * t153 + t116 * t121;
t103 = -t116 * t115 + t121 * t153;
t156 = t119 * t122;
t143 = t103 * t124 - t118 * t156;
t83 = -t104 * t114 + t143 * t120;
t96 = -t103 * t118 - t122 * t155;
t76 = -t83 * t117 + t96 * t123;
t169 = cos(qJ(4));
t161 = t115 * t119;
t151 = t125 * qJ(2) + qJ(1);
t150 = t122 * pkin(1) + qJ(2) * t160;
t147 = t117 * t169;
t146 = t123 * t169;
t145 = g(1) * t116 - g(2) * t122;
t144 = t116 * pkin(1) - qJ(2) * t156;
t127 = sin(qJ(5));
t130 = cos(qJ(5));
t128 = sin(qJ(4));
t84 = t104 * t120 + t143 * t114;
t68 = t84 * t169 + (t117 * t96 + t123 * t83) * t128;
t59 = t68 * t127 - t76 * t130;
t86 = t106 * t120 + t142 * t114;
t70 = t86 * t169 + (t117 * t97 + t123 * t85) * t128;
t61 = t70 * t127 - t77 * t130;
t95 = t120 * t161 + (t119 * t154 + t158) * t114;
t75 = t95 * t169 + (t102 * t117 + t123 * t94) * t128;
t65 = t75 * t127 - t87 * t130;
t141 = g(1) * t61 + g(2) * t59 + g(3) * t65;
t67 = t84 * t128 - t83 * t146 - t96 * t147;
t69 = t86 * t128 - t85 * t146 - t97 * t147;
t74 = -t102 * t147 + t95 * t128 - t94 * t146;
t140 = g(1) * t69 + g(2) * t67 + g(3) * t74;
t139 = t106 * pkin(2) + t97 * qJ(3) + t150;
t138 = pkin(2) * t161 + t102 * qJ(3) + t151;
t137 = t86 * pkin(3) + t77 * pkin(9) + t139;
t136 = t95 * pkin(3) + t87 * pkin(9) + t138;
t135 = t104 * pkin(2) + t96 * qJ(3) + t144;
t134 = t70 * pkin(4) + t69 * pkin(10) + t137;
t133 = t75 * pkin(4) + t74 * pkin(10) + t136;
t132 = t84 * pkin(3) + t76 * pkin(9) + t135;
t131 = t68 * pkin(4) + t67 * pkin(10) + t132;
t129 = cos(qJ(6));
t126 = sin(qJ(6));
t66 = t87 * t127 + t75 * t130;
t62 = t77 * t127 + t70 * t130;
t60 = t76 * t127 + t68 * t130;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t122 - g(2) * t116, t145, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t104 - g(3) * t161, -g(1) * t105 - g(2) * t103 - g(3) * t157, -g(3) * t125 - t145 * t119, -g(1) * t150 - g(2) * t144 - g(3) * t151, 0, 0, 0, 0, 0, 0, -g(1) * t86 - g(2) * t84 - g(3) * t95, -g(1) * t85 - g(2) * t83 - g(3) * t94, -g(1) * t97 - g(2) * t96 - g(3) * t102, -g(1) * t139 - g(2) * t135 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t75, t140, -g(1) * t77 - g(2) * t76 - g(3) * t87, -g(1) * t137 - g(2) * t132 - g(3) * t136, 0, 0, 0, 0, 0, 0, -g(1) * t62 - g(2) * t60 - g(3) * t66, t141, -t140, -g(1) * t134 - g(2) * t131 - g(3) * t133, 0, 0, 0, 0, 0, 0, -g(1) * (t69 * t126 + t62 * t129) - g(2) * (t67 * t126 + t60 * t129) - g(3) * (t74 * t126 + t66 * t129) -g(1) * (-t62 * t126 + t69 * t129) - g(2) * (-t60 * t126 + t67 * t129) - g(3) * (-t66 * t126 + t74 * t129) -t141, -g(1) * (t62 * pkin(5) + t61 * pkin(11) + t134) - g(2) * (t60 * pkin(5) + t59 * pkin(11) + t131) - g(3) * (t66 * pkin(5) + t65 * pkin(11) + t133);];
U_reg  = t1;

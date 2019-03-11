% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:55:42
% EndTime: 2019-03-09 22:55:43
% DurationCPUTime: 0.27s
% Computational Cost: add. (331->106), mult. (533->152), div. (0->0), fcn. (633->14), ass. (0->55)
t118 = cos(pkin(6));
t148 = t118 * pkin(8) + pkin(7);
t115 = sin(pkin(12));
t124 = cos(qJ(2));
t125 = cos(qJ(1));
t136 = t125 * t124;
t121 = sin(qJ(2));
t122 = sin(qJ(1));
t139 = t122 * t121;
t92 = -t118 * t136 + t139;
t147 = t92 * t115;
t137 = t125 * t121;
t138 = t122 * t124;
t94 = t118 * t138 + t137;
t146 = t94 * t115;
t116 = sin(pkin(6));
t145 = t116 * t121;
t144 = t116 * t122;
t123 = cos(qJ(3));
t143 = t116 * t123;
t142 = t116 * t124;
t141 = t116 * t125;
t120 = sin(qJ(3));
t140 = t118 * t120;
t135 = t125 * pkin(1) + pkin(8) * t144;
t134 = t115 * t142;
t133 = t120 * t144;
t111 = t122 * pkin(1);
t132 = -pkin(8) * t141 + t111;
t105 = t123 * pkin(3) + pkin(2);
t126 = -pkin(10) - pkin(9);
t95 = -t118 * t139 + t136;
t131 = pkin(3) * t133 + t95 * t105 - t94 * t126 + t135;
t130 = pkin(3) * t140 + t105 * t145 + t126 * t142 + t148;
t129 = g(1) * t122 - g(2) * t125;
t114 = qJ(3) + qJ(4);
t108 = sin(t114);
t109 = cos(t114);
t93 = t118 * t137 + t138;
t82 = t93 * t108 + t109 * t141;
t84 = t95 * t108 - t109 * t144;
t88 = t108 * t145 - t118 * t109;
t128 = g(1) * t84 + g(2) * t82 + g(3) * t88;
t81 = -g(1) * t94 - g(2) * t92 + g(3) * t142;
t127 = t111 + t93 * t105 - t92 * t126 + (-pkin(3) * t120 - pkin(8)) * t141;
t119 = -pkin(11) - qJ(5);
t117 = cos(pkin(12));
t113 = pkin(12) + qJ(6);
t107 = cos(t113);
t106 = sin(t113);
t104 = t117 * pkin(5) + pkin(4);
t89 = t118 * t108 + t109 * t145;
t85 = t108 * t144 + t95 * t109;
t83 = -t108 * t141 + t93 * t109;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t122, t129, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t145, -t81, -g(3) * t118 - t116 * t129, -g(1) * t135 - g(2) * t132 - g(3) * t148, 0, 0, 0, 0, 0, 0, -g(1) * (t95 * t123 + t133) - g(2) * (-t120 * t141 + t93 * t123) - g(3) * (t121 * t143 + t140) -g(1) * (-t95 * t120 + t122 * t143) - g(2) * (-t93 * t120 - t123 * t141) - g(3) * (t118 * t123 - t120 * t145) t81, -g(1) * (t95 * pkin(2) + t94 * pkin(9) + t135) - g(2) * (t93 * pkin(2) + t92 * pkin(9) + t132) - g(3) * ((pkin(2) * t121 - pkin(9) * t124) * t116 + t148) 0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t83 - g(3) * t89, t128, t81, -g(1) * t131 - g(2) * t127 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * (t85 * t117 + t146) - g(2) * (t83 * t117 + t147) - g(3) * (t89 * t117 - t134) -g(1) * (-t85 * t115 + t94 * t117) - g(2) * (-t83 * t115 + t92 * t117) - g(3) * (-t89 * t115 - t117 * t142) -t128, -g(1) * (t85 * pkin(4) + t84 * qJ(5) + t131) - g(2) * (t83 * pkin(4) + t82 * qJ(5) + t127) - g(3) * (t89 * pkin(4) + t88 * qJ(5) + t130) 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t106 + t85 * t107) - g(2) * (t92 * t106 + t83 * t107) - g(3) * (-t106 * t142 + t89 * t107) -g(1) * (-t85 * t106 + t94 * t107) - g(2) * (-t83 * t106 + t92 * t107) - g(3) * (-t89 * t106 - t107 * t142) -t128, -g(1) * (pkin(5) * t146 + t85 * t104 - t84 * t119 + t131) - g(2) * (pkin(5) * t147 + t83 * t104 - t82 * t119 + t127) - g(3) * (-pkin(5) * t134 + t89 * t104 - t88 * t119 + t130);];
U_reg  = t1;

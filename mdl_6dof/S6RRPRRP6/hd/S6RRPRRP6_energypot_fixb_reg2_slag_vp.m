% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:14:08
% EndTime: 2019-03-09 12:14:09
% DurationCPUTime: 0.24s
% Computational Cost: add. (394->85), mult. (926->131), div. (0->0), fcn. (1177->12), ass. (0->60)
t119 = sin(pkin(11));
t121 = cos(pkin(11));
t125 = sin(qJ(2));
t129 = cos(qJ(2));
t109 = -t125 * t119 + t129 * t121;
t122 = cos(pkin(6));
t155 = t122 * pkin(8) + pkin(7);
t120 = sin(pkin(6));
t154 = t120 * t125;
t126 = sin(qJ(1));
t153 = t120 * t126;
t130 = cos(qJ(1));
t152 = t120 * t130;
t150 = t126 * t125;
t149 = t126 * t129;
t147 = t130 * t125;
t146 = t130 * t129;
t107 = t122 * t125 * pkin(2) + (-pkin(8) - qJ(3)) * t120;
t116 = t129 * pkin(2) + pkin(1);
t145 = t130 * t107 + t126 * t116;
t144 = -t126 * t107 + t130 * t116;
t143 = pkin(2) * t154 + t122 * qJ(3) + t155;
t142 = g(1) * t126 - g(2) * t130;
t135 = t109 * t122;
t140 = t129 * t119 + t125 * t121;
t92 = -t126 * t140 + t130 * t135;
t106 = t140 * t122;
t93 = t130 * t106 + t126 * t109;
t141 = t93 * pkin(3) - t92 * pkin(9) + t145;
t123 = sin(qJ(5));
t127 = cos(qJ(5));
t124 = sin(qJ(4));
t128 = cos(qJ(4));
t85 = -t124 * t152 + t93 * t128;
t76 = t85 * t123 + t92 * t127;
t95 = -t126 * t106 + t130 * t109;
t87 = t124 * t153 + t95 * t128;
t94 = -t126 * t135 - t130 * t140;
t78 = t87 * t123 + t94 * t127;
t104 = t109 * t120;
t105 = t140 * t120;
t98 = t105 * t128 + t122 * t124;
t80 = t104 * t127 + t98 * t123;
t139 = g(1) * t78 + g(2) * t76 + g(3) * t80;
t84 = t93 * t124 + t128 * t152;
t86 = t95 * t124 - t128 * t153;
t97 = t105 * t124 - t122 * t128;
t138 = g(1) * t86 + g(2) * t84 + g(3) * t97;
t137 = g(1) * t94 + g(2) * t92 + g(3) * t104;
t136 = t95 * pkin(3) - t94 * pkin(9) + t144;
t134 = t105 * pkin(3) - t104 * pkin(9) + t143;
t133 = t85 * pkin(4) + t84 * pkin(10) + t141;
t132 = t87 * pkin(4) + t86 * pkin(10) + t136;
t131 = t98 * pkin(4) + t97 * pkin(10) + t134;
t103 = -g(3) * t122 - t142 * t120;
t81 = -t104 * t123 + t98 * t127;
t79 = -t94 * t123 + t87 * t127;
t77 = -t92 * t123 + t85 * t127;
t74 = -g(1) * t79 - g(2) * t77 - g(3) * t81;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t130 - g(2) * t126, t142, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t122 * t150 + t146) - g(2) * (t122 * t147 + t149) - g(3) * t154, -g(1) * (-t122 * t149 - t147) - g(2) * (t122 * t146 - t150) - g(3) * t120 * t129, t103, -g(1) * (t130 * pkin(1) + pkin(8) * t153) - g(2) * (t126 * pkin(1) - pkin(8) * t152) - g(3) * t155, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t105, -t137, t103, -g(1) * t144 - g(2) * t145 - g(3) * t143, 0, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t85 - g(3) * t98, t138, t137, -g(1) * t136 - g(2) * t141 - g(3) * t134, 0, 0, 0, 0, 0, 0, t74, t139, -t138, -g(1) * t132 - g(2) * t133 - g(3) * t131, 0, 0, 0, 0, 0, 0, t74, -t138, -t139, -g(1) * (t79 * pkin(5) + t78 * qJ(6) + t132) - g(2) * (t77 * pkin(5) + t76 * qJ(6) + t133) - g(3) * (t81 * pkin(5) + t80 * qJ(6) + t131);];
U_reg  = t1;

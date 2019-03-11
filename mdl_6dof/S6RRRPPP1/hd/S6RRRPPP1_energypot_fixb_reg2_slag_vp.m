% Calculate inertial parameters regressor of potential energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:18:28
% EndTime: 2019-03-09 15:18:28
% DurationCPUTime: 0.20s
% Computational Cost: add. (268->82), mult. (656->116), div. (0->0), fcn. (786->10), ass. (0->57)
t153 = g(3) * pkin(7);
t118 = sin(qJ(2));
t152 = t118 * pkin(2) + pkin(7);
t151 = cos(pkin(10));
t115 = sin(pkin(6));
t150 = qJ(4) * t115;
t121 = cos(qJ(2));
t149 = t115 * t121;
t117 = sin(qJ(3));
t148 = t117 * t118;
t119 = sin(qJ(1));
t147 = t118 * t119;
t120 = cos(qJ(3));
t146 = t118 * t120;
t122 = cos(qJ(1));
t145 = t118 * t122;
t144 = t119 * t121;
t143 = t122 * t117;
t142 = t122 * t120;
t141 = t122 * pkin(1) + t119 * pkin(8);
t140 = t115 * t148;
t116 = cos(pkin(6));
t139 = t116 * t147;
t138 = t116 * t145;
t137 = t119 * pkin(1) - t122 * pkin(8);
t136 = t116 * t151;
t135 = t118 * t151;
t134 = t122 * t121 * pkin(2) + pkin(9) * t145 + t141;
t133 = t115 * t135;
t132 = g(1) * t122 + g(2) * t119;
t131 = pkin(2) * t144 + pkin(9) * t147 + t137;
t114 = sin(pkin(10));
t93 = -t117 * t144 - t142;
t94 = t120 * t144 - t143;
t78 = t94 * t114 - t119 * t133 - t93 * t136;
t95 = t119 * t120 - t121 * t143;
t96 = t119 * t117 + t121 * t142;
t80 = t96 * t114 - t122 * t133 - t95 * t136;
t83 = t151 * t149 + (t114 * t120 + t117 * t136) * t118;
t130 = g(1) * t80 + g(2) * t78 + g(3) * t83;
t79 = t94 * t151 + (t115 * t147 + t116 * t93) * t114;
t81 = t96 * t151 + (t115 * t145 + t116 * t95) * t114;
t84 = t120 * t135 + (-t116 * t148 - t149) * t114;
t129 = g(1) * t81 + g(2) * t79 + g(3) * t84;
t128 = qJ(4) * t140 + pkin(3) * t146 + (-qJ(4) * t116 - pkin(9)) * t121 + t152;
t127 = t96 * pkin(3) + qJ(4) * t138 - t95 * t150 + t134;
t126 = t94 * pkin(3) + qJ(4) * t139 - t93 * t150 + t131;
t125 = t84 * pkin(4) + t83 * qJ(5) + t128;
t124 = t81 * pkin(4) + t80 * qJ(5) + t127;
t123 = t79 * pkin(4) + t78 * qJ(5) + t126;
t104 = g(1) * t119 - g(2) * t122;
t92 = -t121 * t116 + t140;
t89 = -g(3) * t121 + t132 * t118;
t86 = -t95 * t115 + t138;
t85 = -t93 * t115 + t139;
t75 = -g(1) * t86 - g(2) * t85 - g(3) * t92;
t1 = [0, 0, 0, 0, 0, 0, -t132, t104, -g(3), -t153, 0, 0, 0, 0, 0, 0, -g(3) * t118 - t132 * t121, t89, -t104, -g(1) * t141 - g(2) * t137 - t153, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t146, -g(1) * t95 - g(2) * t93 + g(3) * t148, -t89, -g(1) * t134 - g(2) * t131 - g(3) * (-t121 * pkin(9) + t152) 0, 0, 0, 0, 0, 0, -t129, t130, t75, -g(1) * t127 - g(2) * t126 - g(3) * t128, 0, 0, 0, 0, 0, 0, t75, t129, -t130, -g(1) * t124 - g(2) * t123 - g(3) * t125, 0, 0, 0, 0, 0, 0, t75, -t130, -t129, -g(1) * (t86 * pkin(5) + t81 * qJ(6) + t124) - g(2) * (t85 * pkin(5) + t79 * qJ(6) + t123) - g(3) * (t92 * pkin(5) + t84 * qJ(6) + t125);];
U_reg  = t1;

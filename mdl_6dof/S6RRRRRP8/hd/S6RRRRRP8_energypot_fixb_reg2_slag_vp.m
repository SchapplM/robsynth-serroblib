% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP8
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:48
% EndTime: 2019-03-10 01:55:48
% DurationCPUTime: 0.22s
% Computational Cost: add. (346->91), mult. (579->130), div. (0->0), fcn. (697->12), ass. (0->58)
t120 = cos(pkin(6));
t152 = t120 * pkin(8) + pkin(7);
t119 = sin(pkin(6));
t123 = sin(qJ(2));
t151 = t119 * t123;
t124 = sin(qJ(1));
t150 = t119 * t124;
t126 = cos(qJ(3));
t149 = t119 * t126;
t127 = cos(qJ(2));
t148 = t119 * t127;
t128 = cos(qJ(1));
t147 = t119 * t128;
t122 = sin(qJ(3));
t146 = t120 * t122;
t145 = t124 * t123;
t144 = t124 * t127;
t143 = t128 * t123;
t142 = t128 * t127;
t141 = t128 * pkin(1) + pkin(8) * t150;
t140 = t122 * t150;
t116 = t124 * pkin(1);
t139 = -pkin(8) * t147 + t116;
t138 = g(1) * t124 - g(2) * t128;
t112 = t126 * pkin(3) + pkin(2);
t129 = -pkin(10) - pkin(9);
t137 = pkin(3) * t146 + t112 * t151 + t129 * t148 + t152;
t102 = t120 * t144 + t143;
t103 = -t120 * t145 + t142;
t136 = pkin(3) * t140 - t102 * t129 + t103 * t112 + t141;
t100 = -t120 * t142 + t145;
t121 = sin(qJ(5));
t125 = cos(qJ(5));
t101 = t120 * t143 + t144;
t118 = qJ(3) + qJ(4);
t113 = sin(t118);
t114 = cos(t118);
t88 = t101 * t114 - t113 * t147;
t78 = -t100 * t125 + t88 * t121;
t90 = t103 * t114 + t113 * t150;
t80 = -t102 * t125 + t90 * t121;
t95 = t120 * t113 + t114 * t151;
t85 = t95 * t121 + t125 * t148;
t135 = g(1) * t80 + g(2) * t78 + g(3) * t85;
t87 = t101 * t113 + t114 * t147;
t89 = t103 * t113 - t114 * t150;
t94 = t113 * t151 - t120 * t114;
t134 = g(1) * t89 + g(2) * t87 + g(3) * t94;
t133 = t95 * pkin(4) + t94 * pkin(11) + t137;
t132 = t90 * pkin(4) + t89 * pkin(11) + t136;
t82 = -g(1) * t102 - g(2) * t100 + g(3) * t148;
t131 = t116 + t101 * t112 - t100 * t129 + (-pkin(3) * t122 - pkin(8)) * t147;
t130 = t88 * pkin(4) + t87 * pkin(11) + t131;
t86 = -t121 * t148 + t95 * t125;
t81 = t102 * t121 + t90 * t125;
t79 = t100 * t121 + t88 * t125;
t76 = -g(1) * t81 - g(2) * t79 - g(3) * t86;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t128 - g(2) * t124, t138, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t101 - g(3) * t151, -t82, -g(3) * t120 - t138 * t119, -g(1) * t141 - g(2) * t139 - g(3) * t152, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t126 + t140) - g(2) * (t101 * t126 - t122 * t147) - g(3) * (t123 * t149 + t146) -g(1) * (-t103 * t122 + t124 * t149) - g(2) * (-t101 * t122 - t126 * t147) - g(3) * (t120 * t126 - t122 * t151) t82, -g(1) * (t103 * pkin(2) + t102 * pkin(9) + t141) - g(2) * (t101 * pkin(2) + t100 * pkin(9) + t139) - g(3) * ((pkin(2) * t123 - pkin(9) * t127) * t119 + t152) 0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t88 - g(3) * t95, t134, t82, -g(1) * t136 - g(2) * t131 - g(3) * t137, 0, 0, 0, 0, 0, 0, t76, t135, -t134, -g(1) * t132 - g(2) * t130 - g(3) * t133, 0, 0, 0, 0, 0, 0, t76, -t134, -t135, -g(1) * (t81 * pkin(5) + t80 * qJ(6) + t132) - g(2) * (t79 * pkin(5) + t78 * qJ(6) + t130) - g(3) * (t86 * pkin(5) + t85 * qJ(6) + t133);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP10
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:42:31
% EndTime: 2019-03-09 12:42:31
% DurationCPUTime: 0.23s
% Computational Cost: add. (346->91), mult. (579->129), div. (0->0), fcn. (697->12), ass. (0->57)
t121 = cos(pkin(6));
t150 = t121 * pkin(8) + pkin(7);
t119 = sin(pkin(6));
t124 = sin(qJ(2));
t149 = t119 * t124;
t125 = sin(qJ(1));
t148 = t119 * t125;
t127 = cos(qJ(2));
t147 = t119 * t127;
t128 = cos(qJ(1));
t146 = t119 * t128;
t118 = sin(pkin(11));
t145 = t121 * t118;
t144 = t125 * t124;
t143 = t125 * t127;
t142 = t128 * t124;
t141 = t128 * t127;
t140 = t128 * pkin(1) + pkin(8) * t148;
t139 = t118 * t148;
t115 = t125 * pkin(1);
t138 = -pkin(8) * t146 + t115;
t137 = g(1) * t125 - g(2) * t128;
t120 = cos(pkin(11));
t111 = t120 * pkin(3) + pkin(2);
t122 = -pkin(9) - qJ(3);
t136 = pkin(3) * t145 + t111 * t149 + t122 * t147 + t150;
t101 = t121 * t143 + t142;
t102 = -t121 * t144 + t141;
t135 = pkin(3) * t139 - t101 * t122 + t102 * t111 + t140;
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t100 = t121 * t142 + t143;
t117 = pkin(11) + qJ(4);
t112 = sin(t117);
t113 = cos(t117);
t87 = t100 * t113 - t112 * t146;
t99 = -t121 * t141 + t144;
t77 = t87 * t123 - t99 * t126;
t89 = t102 * t113 + t112 * t148;
t79 = -t101 * t126 + t89 * t123;
t94 = t121 * t112 + t113 * t149;
t84 = t94 * t123 + t126 * t147;
t134 = g(1) * t79 + g(2) * t77 + g(3) * t84;
t86 = t100 * t112 + t113 * t146;
t88 = t102 * t112 - t113 * t148;
t93 = t112 * t149 - t121 * t113;
t133 = g(1) * t88 + g(2) * t86 + g(3) * t93;
t81 = -g(1) * t101 - g(2) * t99 + g(3) * t147;
t132 = t94 * pkin(4) + t93 * pkin(10) + t136;
t131 = t89 * pkin(4) + t88 * pkin(10) + t135;
t130 = t115 + t100 * t111 - t99 * t122 + (-pkin(3) * t118 - pkin(8)) * t146;
t129 = t87 * pkin(4) + t86 * pkin(10) + t130;
t85 = -t123 * t147 + t94 * t126;
t80 = t101 * t123 + t89 * t126;
t78 = t99 * t123 + t87 * t126;
t75 = -g(1) * t80 - g(2) * t78 - g(3) * t85;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t128 - g(2) * t125, t137, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t102 - g(2) * t100 - g(3) * t149, -t81, -g(3) * t121 - t137 * t119, -g(1) * t140 - g(2) * t138 - g(3) * t150, 0, 0, 0, 0, 0, 0, -g(1) * (t102 * t120 + t139) - g(2) * (t100 * t120 - t118 * t146) - g(3) * (t120 * t149 + t145) -g(1) * (-t102 * t118 + t120 * t148) - g(2) * (-t100 * t118 - t120 * t146) - g(3) * (-t118 * t149 + t121 * t120) t81, -g(1) * (t102 * pkin(2) + t101 * qJ(3) + t140) - g(2) * (t100 * pkin(2) + t99 * qJ(3) + t138) - g(3) * ((pkin(2) * t124 - qJ(3) * t127) * t119 + t150) 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 - g(3) * t94, t133, t81, -g(1) * t135 - g(2) * t130 - g(3) * t136, 0, 0, 0, 0, 0, 0, t75, t134, -t133, -g(1) * t131 - g(2) * t129 - g(3) * t132, 0, 0, 0, 0, 0, 0, t75, -t133, -t134, -g(1) * (t80 * pkin(5) + t79 * qJ(6) + t131) - g(2) * (t78 * pkin(5) + t77 * qJ(6) + t129) - g(3) * (t85 * pkin(5) + t84 * qJ(6) + t132);];
U_reg  = t1;

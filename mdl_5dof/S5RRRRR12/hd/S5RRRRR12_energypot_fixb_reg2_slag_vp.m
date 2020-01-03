% Calculate inertial parameters regressor of potential energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:10
% EndTime: 2019-12-31 22:55:10
% DurationCPUTime: 0.27s
% Computational Cost: add. (336->85), mult. (878->137), div. (0->0), fcn. (1108->14), ass. (0->57)
t118 = cos(pkin(5));
t122 = sin(qJ(2));
t127 = cos(qJ(1));
t145 = t127 * t122;
t123 = sin(qJ(1));
t126 = cos(qJ(2));
t146 = t123 * t126;
t103 = -t118 * t146 - t145;
t115 = sin(pkin(6));
t117 = cos(pkin(6));
t116 = sin(pkin(5));
t151 = t116 * t123;
t94 = -t103 * t115 + t117 * t151;
t150 = t116 * t126;
t100 = -t115 * t150 + t118 * t117;
t155 = cos(qJ(3));
t154 = t118 * pkin(8) + pkin(7);
t152 = t116 * t122;
t149 = t116 * t127;
t147 = t123 * t122;
t144 = t127 * t126;
t143 = t127 * pkin(1) + pkin(8) * t151;
t140 = t115 * t155;
t139 = t117 * t155;
t138 = t116 * t140;
t137 = t123 * pkin(1) - pkin(8) * t149;
t136 = g(1) * t123 - g(2) * t127;
t101 = t118 * t144 - t147;
t93 = -t101 * t115 - t117 * t149;
t120 = sin(qJ(4));
t125 = cos(qJ(4));
t102 = t118 * t145 + t146;
t121 = sin(qJ(3));
t85 = t102 * t155 + (t101 * t117 - t115 * t149) * t121;
t76 = t85 * t120 - t93 * t125;
t104 = -t118 * t147 + t144;
t87 = t104 * t155 + (t103 * t117 + t115 * t151) * t121;
t78 = t87 * t120 - t94 * t125;
t92 = t118 * t115 * t121 + (t117 * t121 * t126 + t155 * t122) * t116;
t82 = -t100 * t125 + t92 * t120;
t135 = g(1) * t78 + g(2) * t76 + g(3) * t82;
t84 = -t101 * t139 + t102 * t121 + t127 * t138;
t86 = -t103 * t139 + t104 * t121 - t123 * t138;
t91 = -t118 * t140 + t121 * t152 - t139 * t150;
t134 = g(1) * t86 + g(2) * t84 + g(3) * t91;
t133 = t104 * pkin(2) + t94 * pkin(9) + t143;
t132 = pkin(2) * t152 + t100 * pkin(9) + t154;
t131 = t87 * pkin(3) + t86 * pkin(10) + t133;
t130 = t92 * pkin(3) + t91 * pkin(10) + t132;
t129 = t102 * pkin(2) + t93 * pkin(9) + t137;
t128 = t85 * pkin(3) + t84 * pkin(10) + t129;
t124 = cos(qJ(5));
t119 = sin(qJ(5));
t83 = t100 * t120 + t92 * t125;
t79 = t94 * t120 + t87 * t125;
t77 = t93 * t120 + t85 * t125;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t127 - g(2) * t123, t136, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102 - g(3) * t152, -g(1) * t103 - g(2) * t101 - g(3) * t150, -g(3) * t118 - t136 * t116, -g(1) * t143 - g(2) * t137 - g(3) * t154, 0, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t85 - g(3) * t92, t134, -g(1) * t94 - g(2) * t93 - g(3) * t100, -g(1) * t133 - g(2) * t129 - g(3) * t132, 0, 0, 0, 0, 0, 0, -g(1) * t79 - g(2) * t77 - g(3) * t83, t135, -t134, -g(1) * t131 - g(2) * t128 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t119 + t79 * t124) - g(2) * (t84 * t119 + t77 * t124) - g(3) * (t91 * t119 + t83 * t124), -g(1) * (-t79 * t119 + t86 * t124) - g(2) * (-t77 * t119 + t84 * t124) - g(3) * (-t83 * t119 + t91 * t124), -t135, -g(1) * (t79 * pkin(4) + t78 * pkin(11) + t131) - g(2) * (t77 * pkin(4) + t76 * pkin(11) + t128) - g(3) * (t83 * pkin(4) + t82 * pkin(11) + t130);];
U_reg = t1;

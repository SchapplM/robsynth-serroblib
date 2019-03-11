% Calculate inertial parameters regressor of potential energy for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:53
% EndTime: 2019-03-08 18:50:54
% DurationCPUTime: 0.32s
% Computational Cost: add. (508->100), mult. (1330->149), div. (0->0), fcn. (1691->14), ass. (0->64)
t106 = sin(pkin(7));
t110 = cos(pkin(7));
t111 = cos(pkin(6));
t107 = sin(pkin(6));
t108 = cos(pkin(12));
t143 = t107 * t108;
t128 = t106 * t143 - t111 * t110;
t105 = sin(pkin(11));
t141 = t107 * t110;
t104 = sin(pkin(12));
t109 = cos(pkin(11));
t144 = t105 * t111;
t92 = -t109 * t104 - t108 * t144;
t130 = t105 * t141 - t92 * t106;
t154 = pkin(5) + pkin(9);
t114 = sin(qJ(3));
t150 = cos(qJ(3));
t135 = t106 * t150;
t133 = t107 * t135;
t134 = t110 * t150;
t140 = t109 * t111;
t90 = -t105 * t104 + t108 * t140;
t91 = t104 * t140 + t105 * t108;
t74 = t109 * t133 + t91 * t114 - t90 * t134;
t153 = t74 * pkin(9);
t93 = -t104 * t144 + t109 * t108;
t76 = -t105 * t133 + t93 * t114 - t92 * t134;
t152 = t76 * pkin(9);
t146 = t104 * t107;
t83 = -t111 * t135 + t114 * t146 - t134 * t143;
t151 = t83 * pkin(9);
t149 = cos(qJ(4));
t145 = t105 * t107;
t147 = t109 * pkin(1) + qJ(2) * t145;
t142 = t107 * t109;
t138 = t111 * qJ(2) + qJ(1);
t132 = g(1) * t105 - g(2) * t109;
t131 = t105 * pkin(1) - qJ(2) * t142;
t129 = t90 * t106 + t109 * t141;
t113 = sin(qJ(4));
t75 = t91 * t150 + (-t106 * t142 + t110 * t90) * t114;
t67 = t75 * t113 + t129 * t149;
t77 = t93 * t150 + (t106 * t145 + t110 * t92) * t114;
t69 = t77 * t113 - t130 * t149;
t84 = t111 * t106 * t114 + (t108 * t110 * t114 + t150 * t104) * t107;
t78 = t84 * t113 + t128 * t149;
t127 = g(1) * t69 + g(2) * t67 + g(3) * t78;
t68 = -t129 * t113 + t75 * t149;
t70 = t130 * t113 + t77 * t149;
t79 = -t128 * t113 + t84 * t149;
t126 = g(1) * t70 + g(2) * t68 + g(3) * t79;
t125 = g(1) * t76 + g(2) * t74 + g(3) * t83;
t124 = t93 * pkin(2) + t130 * pkin(8) + t147;
t123 = t77 * pkin(3) + t124;
t122 = pkin(2) * t146 - t128 * pkin(8) + t138;
t121 = t84 * pkin(3) + t122;
t120 = t70 * pkin(4) + t69 * qJ(5) + t123;
t119 = t79 * pkin(4) + t78 * qJ(5) + t121;
t118 = t91 * pkin(2) - t129 * pkin(8) + t131;
t117 = t75 * pkin(3) + t118;
t116 = t68 * pkin(4) + t67 * qJ(5) + t117;
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t105, t132, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t146, -g(1) * t92 - g(2) * t90 - g(3) * t143, -g(3) * t111 - t132 * t107, -g(1) * t147 - g(2) * t131 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * t77 - g(2) * t75 - g(3) * t84, t125, -g(1) * t130 + g(2) * t129 + g(3) * t128, -g(1) * t124 - g(2) * t118 - g(3) * t122, 0, 0, 0, 0, 0, 0, -t126, t127, -t125, -g(1) * (t123 + t152) - g(2) * (t117 + t153) - g(3) * (t121 + t151) 0, 0, 0, 0, 0, 0, -t125, t126, -t127, -g(1) * (t120 + t152) - g(2) * (t116 + t153) - g(3) * (t119 + t151) 0, 0, 0, 0, 0, 0, -g(1) * (t69 * t112 + t76 * t115) - g(2) * (t67 * t112 + t74 * t115) - g(3) * (t78 * t112 + t83 * t115) -g(1) * (-t76 * t112 + t69 * t115) - g(2) * (-t74 * t112 + t67 * t115) - g(3) * (-t83 * t112 + t78 * t115) -t126, -g(1) * (t70 * pkin(10) + t154 * t76 + t120) - g(2) * (t68 * pkin(10) + t154 * t74 + t116) - g(3) * (t79 * pkin(10) + t154 * t83 + t119);];
U_reg  = t1;

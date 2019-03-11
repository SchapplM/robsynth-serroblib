% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:10:58
% EndTime: 2019-03-09 09:10:58
% DurationCPUTime: 0.18s
% Computational Cost: add. (208->79), mult. (457->107), div. (0->0), fcn. (530->10), ass. (0->49)
t102 = sin(pkin(6));
t111 = cos(qJ(1));
t127 = t102 * t111;
t103 = cos(pkin(6));
t106 = sin(qJ(2));
t124 = t111 * t106;
t107 = sin(qJ(1));
t110 = cos(qJ(2));
t125 = t107 * t110;
t87 = t103 * t124 + t125;
t136 = t87 * pkin(3) + qJ(4) * t127;
t135 = t103 * pkin(8) + pkin(7);
t134 = -pkin(9) + qJ(3);
t88 = t103 * t125 + t124;
t133 = t88 * qJ(3);
t130 = t102 * t107;
t132 = t111 * pkin(1) + pkin(8) * t130;
t131 = t102 * t106;
t109 = cos(qJ(5));
t129 = t102 * t109;
t128 = t102 * t110;
t126 = t107 * t106;
t123 = t111 * t110;
t89 = -t103 * t126 + t123;
t122 = t89 * pkin(2) + t132;
t121 = t107 * pkin(1) - pkin(8) * t127;
t120 = t87 * pkin(2) + t121;
t119 = pkin(2) * t131 - qJ(3) * t128 + t135;
t105 = sin(qJ(5));
t73 = t87 * t105 - t109 * t127;
t75 = t89 * t105 + t107 * t129;
t84 = t103 * t109 + t105 * t131;
t118 = g(1) * t75 + g(2) * t73 + g(3) * t84;
t117 = t89 * pkin(3) - qJ(4) * t130 + t122;
t86 = -t103 * t123 + t126;
t70 = -g(1) * t88 - g(2) * t86 + g(3) * t128;
t116 = t86 * qJ(3) + t120;
t115 = pkin(3) * t131 - t103 * qJ(4) + t119;
t114 = pkin(4) * t131 + pkin(9) * t128 + t115;
t113 = t87 * pkin(4) + t134 * t86 + t120 + t136;
t112 = t89 * pkin(4) + t134 * t88 + t117;
t108 = cos(qJ(6));
t104 = sin(qJ(6));
t85 = -t103 * t105 + t106 * t129;
t77 = g(1) * t130 - g(2) * t127 + g(3) * t103;
t76 = -t105 * t130 + t89 * t109;
t74 = t105 * t127 + t87 * t109;
t72 = -g(1) * t89 - g(2) * t87 - g(3) * t131;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t111 - g(2) * t107, g(1) * t107 - g(2) * t111, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, t72, -t70, -t77, -g(1) * t132 - g(2) * t121 - g(3) * t135, 0, 0, 0, 0, 0, 0, t72, -t77, t70, -g(1) * (t122 + t133) - g(2) * t116 - g(3) * t119, 0, 0, 0, 0, 0, 0, t72, t70, t77, -g(1) * (t117 + t133) - g(2) * (t116 + t136) - g(3) * t115, 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t74 - g(3) * t85, t118, -t70, -g(1) * t112 - g(2) * t113 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * (-t88 * t104 + t76 * t108) - g(2) * (-t86 * t104 + t74 * t108) - g(3) * (t104 * t128 + t85 * t108) -g(1) * (-t76 * t104 - t88 * t108) - g(2) * (-t74 * t104 - t86 * t108) - g(3) * (-t85 * t104 + t108 * t128) -t118, -g(1) * (t76 * pkin(5) + t75 * pkin(10) + t112) - g(2) * (t74 * pkin(5) + t73 * pkin(10) + t113) - g(3) * (t85 * pkin(5) + t84 * pkin(10) + t114);];
U_reg  = t1;

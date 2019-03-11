% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:51:24
% EndTime: 2019-03-10 03:51:24
% DurationCPUTime: 0.21s
% Computational Cost: add. (221->99), mult. (228->128), div. (0->0), fcn. (228->12), ass. (0->41)
t121 = g(3) * pkin(6);
t108 = -pkin(9) - pkin(8);
t105 = cos(qJ(3));
t88 = t105 * pkin(3) + pkin(2);
t103 = sin(qJ(2));
t120 = g(3) * t103;
t102 = sin(qJ(3));
t94 = t102 * pkin(3);
t101 = qJ(3) + qJ(4);
t90 = sin(t101);
t80 = pkin(4) * t90 + t94;
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t119 = t107 * pkin(1) + t104 * pkin(7);
t118 = t103 * t104;
t117 = t103 * t108;
t116 = t104 * t102;
t106 = cos(qJ(2));
t115 = t104 * t106;
t114 = t106 * t107;
t113 = t107 * t102;
t112 = t107 * t105;
t100 = -pkin(10) + t108;
t91 = cos(t101);
t79 = pkin(4) * t91 + t88;
t92 = qJ(5) + t101;
t96 = t104 * pkin(1);
t111 = -t107 * pkin(7) + t96;
t110 = pkin(2) * t106 + pkin(8) * t103;
t109 = g(1) * t107 + g(2) * t104;
t93 = -pkin(11) + t100;
t89 = qJ(6) + t92;
t87 = cos(t92);
t86 = sin(t92);
t83 = cos(t89);
t82 = sin(t89);
t81 = g(1) * t104 - g(2) * t107;
t78 = pkin(5) * t86 + t80;
t77 = -g(3) * t106 + t109 * t103;
t76 = pkin(5) * t87 + t79;
t1 = [0, 0, 0, 0, 0, 0, -t109, t81, -g(3), -t121, 0, 0, 0, 0, 0, 0, -t109 * t106 - t120, t77, -t81, -g(1) * t119 - g(2) * t111 - t121, 0, 0, 0, 0, 0, 0, -g(1) * (t106 * t112 + t116) - g(2) * (t105 * t115 - t113) - t105 * t120, -g(1) * (t104 * t105 - t106 * t113) - g(2) * (-t102 * t115 - t112) + t102 * t120, -t77, -g(1) * (t110 * t107 + t119) - g(2) * (t110 * t104 + t111) - g(3) * (t103 * pkin(2) - t106 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t90 + t91 * t114) - g(2) * (-t107 * t90 + t91 * t115) - t91 * t120, -g(1) * (t104 * t91 - t90 * t114) - g(2) * (-t107 * t91 - t90 * t115) + t90 * t120, -t77, -g(1) * (pkin(3) * t116 + t119) - g(2) * (-t104 * t117 + t88 * t115 + t96) - g(3) * (t103 * t88 + t106 * t108 + pkin(6)) + (-g(1) * (t106 * t88 - t117) - g(2) * (-pkin(7) - t94)) * t107, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t86 + t87 * t114) - g(2) * (-t107 * t86 + t87 * t115) - t87 * t120, -g(1) * (t104 * t87 - t86 * t114) - g(2) * (-t107 * t87 - t86 * t115) + t86 * t120, -t77, -g(1) * (t104 * t80 + t119) - g(2) * (-t100 * t118 + t79 * t115 + t96) - g(3) * (t106 * t100 + t103 * t79 + pkin(6)) + (-g(1) * (-t100 * t103 + t106 * t79) - g(2) * (-pkin(7) - t80)) * t107, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t82 + t83 * t114) - g(2) * (-t107 * t82 + t83 * t115) - t83 * t120, -g(1) * (t104 * t83 - t82 * t114) - g(2) * (-t107 * t83 - t82 * t115) + t82 * t120, -t77, -g(1) * (t104 * t78 + t119) - g(2) * (t76 * t115 - t93 * t118 + t96) - g(3) * (t103 * t76 + t106 * t93 + pkin(6)) + (-g(1) * (-t103 * t93 + t106 * t76) - g(2) * (-pkin(7) - t78)) * t107;];
U_reg  = t1;

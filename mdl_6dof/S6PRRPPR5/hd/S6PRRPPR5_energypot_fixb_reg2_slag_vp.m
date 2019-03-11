% Calculate inertial parameters regressor of potential energy for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:04
% EndTime: 2019-03-08 21:21:04
% DurationCPUTime: 0.25s
% Computational Cost: add. (277->96), mult. (600->131), div. (0->0), fcn. (727->12), ass. (0->54)
t132 = pkin(4) + pkin(8);
t96 = sin(pkin(6));
t131 = pkin(7) * t96;
t101 = sin(qJ(2));
t102 = cos(qJ(2));
t121 = cos(pkin(6));
t114 = t102 * t121;
t95 = sin(pkin(10));
t98 = cos(pkin(10));
t75 = t95 * t101 - t98 * t114;
t130 = t75 * pkin(8);
t77 = t98 * t101 + t95 * t114;
t129 = t77 * pkin(8);
t97 = cos(pkin(11));
t128 = t97 * pkin(5) + t132;
t127 = cos(qJ(3));
t126 = t98 * pkin(1) + t95 * t131;
t100 = sin(qJ(3));
t125 = t100 * t96;
t124 = t101 * t96;
t123 = t102 * t96;
t122 = t121 * pkin(7) + qJ(1);
t120 = pkin(8) * t123;
t115 = t101 * t121;
t78 = t98 * t102 - t95 * t115;
t119 = t78 * pkin(2) + t126;
t118 = pkin(2) * t124 + t122;
t117 = t96 * t127;
t94 = sin(pkin(11));
t116 = pkin(5) * t94 + qJ(4);
t71 = t95 * t125 + t78 * t127;
t113 = t71 * pkin(3) + t119;
t80 = t121 * t100 + t101 * t117;
t112 = t80 * pkin(3) + t118;
t111 = t95 * pkin(1) - t98 * t131;
t110 = g(1) * t95 - g(2) * t98;
t76 = t95 * t102 + t98 * t115;
t109 = t76 * pkin(2) + t111;
t69 = -t98 * t125 + t76 * t127;
t108 = t69 * pkin(3) + t109;
t70 = t78 * t100 - t95 * t117;
t107 = t70 * qJ(4) + t113;
t79 = t100 * t124 - t121 * t127;
t106 = t79 * qJ(4) + t112;
t68 = t76 * t100 + t98 * t117;
t105 = g(1) * t70 + g(2) * t68 + g(3) * t79;
t104 = g(1) * t71 + g(2) * t69 + g(3) * t80;
t65 = -g(1) * t77 - g(2) * t75 + g(3) * t123;
t103 = t68 * qJ(4) + t108;
t99 = -pkin(9) - qJ(5);
t93 = pkin(11) + qJ(6);
t89 = cos(t93);
t88 = sin(t93);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t98 - g(2) * t95, t110, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t124, -t65, -g(3) * t121 - t110 * t96, -g(1) * t126 - g(2) * t111 - g(3) * t122, 0, 0, 0, 0, 0, 0, -t104, t105, t65, -g(1) * (t119 + t129) - g(2) * (t109 + t130) - g(3) * (t118 - t120) 0, 0, 0, 0, 0, 0, t65, t104, -t105, -g(1) * (t107 + t129) - g(2) * (t103 + t130) - g(3) * (t106 - t120) 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t94 + t77 * t97) - g(2) * (t68 * t94 + t75 * t97) - g(3) * (-t97 * t123 + t79 * t94) -g(1) * (t70 * t97 - t77 * t94) - g(2) * (t68 * t97 - t75 * t94) - g(3) * (t94 * t123 + t79 * t97) -t104, -g(1) * (t71 * qJ(5) + t132 * t77 + t107) - g(2) * (t69 * qJ(5) + t132 * t75 + t103) - g(3) * (t80 * qJ(5) - t132 * t123 + t106) 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t88 + t77 * t89) - g(2) * (t68 * t88 + t75 * t89) - g(3) * (-t89 * t123 + t79 * t88) -g(1) * (t70 * t89 - t77 * t88) - g(2) * (t68 * t89 - t75 * t88) - g(3) * (t88 * t123 + t79 * t89) -t104, -g(1) * (t116 * t70 + t128 * t77 - t71 * t99 + t113) - g(2) * (t116 * t68 + t128 * t75 - t69 * t99 + t108) - g(3) * (t116 * t79 - t128 * t123 - t80 * t99 + t112);];
U_reg  = t1;

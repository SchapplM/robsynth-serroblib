% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:06
% EndTime: 2019-03-09 10:19:06
% DurationCPUTime: 0.16s
% Computational Cost: add. (212->86), mult. (200->111), div. (0->0), fcn. (196->12), ass. (0->44)
t120 = g(3) * pkin(6);
t94 = qJ(2) + pkin(10);
t85 = sin(t94);
t119 = g(3) * t85;
t97 = sin(qJ(4));
t118 = t97 * pkin(4);
t98 = sin(qJ(2));
t117 = t98 * pkin(2) + pkin(6);
t100 = cos(qJ(4));
t82 = t100 * pkin(4) + pkin(3);
t93 = qJ(4) + pkin(11);
t88 = qJ(6) + t93;
t79 = sin(t88);
t99 = sin(qJ(1));
t116 = t99 * t79;
t80 = cos(t88);
t115 = t99 * t80;
t84 = sin(t93);
t114 = t99 * t84;
t86 = cos(t93);
t113 = t99 * t86;
t112 = t99 * t97;
t95 = -qJ(5) - pkin(8);
t102 = cos(qJ(1));
t101 = cos(qJ(2));
t83 = t101 * pkin(2) + pkin(1);
t96 = -pkin(7) - qJ(3);
t111 = t102 * t96 + t99 * t83;
t110 = t102 * t85;
t87 = cos(t94);
t109 = t102 * t87;
t108 = t102 * t97;
t107 = t99 * t100;
t106 = t102 * t100;
t77 = t102 * t83;
t105 = -t99 * t96 + t77;
t104 = pkin(3) * t87 + pkin(8) * t85;
t103 = g(1) * t102 + g(2) * t99;
t92 = -pkin(9) + t95;
t75 = g(1) * t99 - g(2) * t102;
t74 = pkin(5) * t84 + t118;
t73 = pkin(5) * t86 + t82;
t72 = -g(3) * t87 + t103 * t85;
t1 = [0, 0, 0, 0, 0, 0, -t103, t75, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t98 - t103 * t101, -g(3) * t101 + t103 * t98, -t75, -g(1) * (t102 * pkin(1) + t99 * pkin(7)) - g(2) * (t99 * pkin(1) - t102 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -t103 * t87 - t119, t72, -t75, -g(1) * t105 - g(2) * t111 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t87 * t106 + t112) - g(2) * (t87 * t107 - t108) - t100 * t119, -g(1) * (-t87 * t108 + t107) - g(2) * (-t87 * t112 - t106) + t97 * t119, -t72, -g(1) * (t104 * t102 + t105) - g(2) * (t104 * t99 + t111) - g(3) * (t85 * pkin(3) - t87 * pkin(8) + t117) 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t109 + t114) - g(2) * (-t102 * t84 + t87 * t113) - t86 * t119, -g(1) * (-t84 * t109 + t113) - g(2) * (-t102 * t86 - t87 * t114) + t84 * t119, -t72, -g(1) * (t82 * t109 - t95 * t110 + t77) - g(2) * (-pkin(4) * t108 + t111) - g(3) * (t85 * t82 + t87 * t95 + t117) + (-g(1) * (-t96 + t118) - g(2) * (t82 * t87 - t85 * t95)) * t99, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t109 + t116) - g(2) * (-t102 * t79 + t87 * t115) - t80 * t119, -g(1) * (-t79 * t109 + t115) - g(2) * (-t102 * t80 - t87 * t116) + t79 * t119, -t72, -g(1) * (t73 * t109 - t92 * t110 + t77) - g(2) * (-t102 * t74 + t111) - g(3) * (t85 * t73 + t87 * t92 + t117) + (-g(1) * (t74 - t96) - g(2) * (t73 * t87 - t85 * t92)) * t99;];
U_reg  = t1;

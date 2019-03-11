% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:35
% EndTime: 2019-03-09 08:52:35
% DurationCPUTime: 0.18s
% Computational Cost: add. (212->86), mult. (200->112), div. (0->0), fcn. (196->12), ass. (0->41)
t117 = g(3) * pkin(6);
t94 = qJ(2) + pkin(10);
t85 = sin(t94);
t116 = g(3) * t85;
t95 = sin(pkin(11));
t115 = t95 * pkin(4);
t96 = cos(pkin(11));
t81 = t96 * pkin(4) + pkin(3);
t99 = sin(qJ(2));
t114 = t99 * pkin(2) + pkin(6);
t97 = -pkin(8) - qJ(4);
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t101 = cos(qJ(2));
t83 = t101 * pkin(2) + pkin(1);
t98 = -pkin(7) - qJ(3);
t113 = t100 * t83 + t102 * t98;
t87 = cos(t94);
t112 = t100 * t87;
t111 = t100 * t95;
t110 = t100 * t96;
t109 = t102 * t85;
t108 = t102 * t87;
t107 = t102 * t95;
t106 = t102 * t96;
t93 = pkin(11) + qJ(5);
t77 = t102 * t83;
t105 = -t100 * t98 + t77;
t104 = pkin(3) * t87 + qJ(4) * t85;
t103 = g(1) * t102 + g(2) * t100;
t92 = -pkin(9) + t97;
t88 = qJ(6) + t93;
t86 = cos(t93);
t84 = sin(t93);
t80 = cos(t88);
t79 = sin(t88);
t75 = g(1) * t100 - g(2) * t102;
t74 = pkin(5) * t84 + t115;
t73 = pkin(5) * t86 + t81;
t72 = -g(3) * t87 + t103 * t85;
t1 = [0, 0, 0, 0, 0, 0, -t103, t75, -g(3), -t117, 0, 0, 0, 0, 0, 0, -g(3) * t99 - t103 * t101, -g(3) * t101 + t103 * t99, -t75, -g(1) * (t102 * pkin(1) + t100 * pkin(7)) - g(2) * (t100 * pkin(1) - t102 * pkin(7)) - t117, 0, 0, 0, 0, 0, 0, -t103 * t87 - t116, t72, -t75, -g(1) * t105 - g(2) * t113 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * (t87 * t106 + t111) - g(2) * (t87 * t110 - t107) - t96 * t116, -g(1) * (-t87 * t107 + t110) - g(2) * (-t87 * t111 - t106) + t95 * t116, -t72, -g(1) * (t104 * t102 + t105) - g(2) * (t104 * t100 + t113) - g(3) * (t85 * pkin(3) - t87 * qJ(4) + t114) 0, 0, 0, 0, 0, 0, -g(1) * (t100 * t84 + t86 * t108) - g(2) * (-t102 * t84 + t86 * t112) - t86 * t116, -g(1) * (t100 * t86 - t84 * t108) - g(2) * (-t102 * t86 - t84 * t112) + t84 * t116, -t72, -g(1) * (t81 * t108 - t97 * t109 + t77) - g(2) * (-pkin(4) * t107 + t113) - g(3) * (t85 * t81 + t87 * t97 + t114) + (-g(1) * (-t98 + t115) - g(2) * (t81 * t87 - t85 * t97)) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (t100 * t79 + t80 * t108) - g(2) * (-t102 * t79 + t80 * t112) - t80 * t116, -g(1) * (t100 * t80 - t79 * t108) - g(2) * (-t102 * t80 - t79 * t112) + t79 * t116, -t72, -g(1) * (t73 * t108 - t92 * t109 + t77) - g(2) * (-t102 * t74 + t113) - g(3) * (t85 * t73 + t87 * t92 + t114) + (-g(1) * (t74 - t98) - g(2) * (t73 * t87 - t85 * t92)) * t100;];
U_reg  = t1;

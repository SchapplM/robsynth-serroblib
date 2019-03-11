% Calculate inertial parameters regressor of potential energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:05
% EndTime: 2019-03-09 08:12:05
% DurationCPUTime: 0.15s
% Computational Cost: add. (182->76), mult. (192->91), div. (0->0), fcn. (184->10), ass. (0->45)
t90 = qJ(2) + pkin(9);
t84 = sin(t90);
t107 = qJ(4) * t84;
t86 = cos(t90);
t98 = cos(qJ(1));
t117 = t86 * t98;
t122 = pkin(3) * t117 + t98 * t107;
t121 = g(3) * pkin(6);
t91 = sin(pkin(10));
t120 = pkin(5) * t91;
t119 = g(3) * t86;
t95 = sin(qJ(2));
t118 = t95 * pkin(2) + pkin(6);
t89 = pkin(10) + qJ(6);
t83 = sin(t89);
t96 = sin(qJ(1));
t116 = t96 * t83;
t85 = cos(t89);
t115 = t96 * t85;
t114 = t96 * t91;
t92 = cos(pkin(10));
t113 = t96 * t92;
t112 = t98 * t83;
t111 = t98 * t85;
t110 = t98 * t91;
t109 = t98 * t92;
t97 = cos(qJ(2));
t82 = t97 * pkin(2) + pkin(1);
t94 = -pkin(7) - qJ(3);
t108 = t96 * t82 + t98 * t94;
t106 = qJ(5) * t86;
t105 = t84 * pkin(3) + t118;
t104 = t84 * t110;
t78 = t98 * t82;
t103 = t78 + t122;
t102 = -t96 * t94 + t78;
t101 = t108 + (pkin(3) * t86 + t107) * t96;
t100 = g(1) * t98 + g(2) * t96;
t99 = -t86 * qJ(4) + t105;
t93 = -pkin(8) - qJ(5);
t80 = t92 * pkin(5) + pkin(4);
t74 = g(1) * t96 - g(2) * t98;
t71 = g(3) * t84 + t100 * t86;
t70 = t100 * t84 - t119;
t1 = [0, 0, 0, 0, 0, 0, -t100, t74, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t95 - t100 * t97, -g(3) * t97 + t100 * t95, -t74, -g(1) * (t98 * pkin(1) + t96 * pkin(7)) - g(2) * (t96 * pkin(1) - t98 * pkin(7)) - t121, 0, 0, 0, 0, 0, 0, -t71, t70, -t74, -g(1) * t102 - g(2) * t108 - g(3) * t118, 0, 0, 0, 0, 0, 0, -t74, t71, -t70, -g(1) * (t102 + t122) - g(2) * t101 - g(3) * t99, 0, 0, 0, 0, 0, 0, -g(1) * (t104 + t113) - g(2) * (t84 * t114 - t109) + t91 * t119, -g(1) * (t84 * t109 - t114) - g(2) * (t84 * t113 + t110) + t92 * t119, -t71, -g(1) * (t98 * t106 + (pkin(4) - t94) * t96 + t103) - g(2) * (-t98 * pkin(4) + t96 * t106 + t101) - g(3) * (t84 * qJ(5) + t99) 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t112 + t115) - g(2) * (t84 * t116 - t111) + t83 * t119, -g(1) * (t84 * t111 - t116) - g(2) * (t84 * t115 + t112) + t85 * t119, -t71, -g(1) * (pkin(5) * t104 - t93 * t117 + t103) - g(2) * (-t98 * t80 + t101) - g(3) * (-t84 * t93 + (-qJ(4) - t120) * t86 + t105) + (-g(1) * (t80 - t94) - g(2) * (t84 * t120 - t86 * t93)) * t96;];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:02
% EndTime: 2019-03-09 11:15:02
% DurationCPUTime: 0.21s
% Computational Cost: add. (162->90), mult. (225->107), div. (0->0), fcn. (221->10), ass. (0->50)
t94 = sin(qJ(2));
t105 = qJ(3) * t94;
t95 = sin(qJ(1));
t97 = cos(qJ(2));
t114 = t95 * t97;
t126 = pkin(2) * t114 + t95 * t105;
t125 = g(3) * pkin(6);
t124 = g(3) * t97;
t93 = sin(qJ(4));
t123 = t93 * pkin(4);
t122 = t94 * pkin(2) + pkin(6);
t96 = cos(qJ(4));
t80 = t96 * pkin(4) + pkin(3);
t91 = qJ(4) + pkin(10);
t81 = sin(t91);
t72 = pkin(5) * t81 + t123;
t121 = t72 * t94;
t83 = qJ(6) + t91;
t78 = sin(t83);
t120 = t95 * t78;
t79 = cos(t83);
t119 = t95 * t79;
t118 = t95 * t81;
t82 = cos(t91);
t117 = t95 * t82;
t116 = t95 * t93;
t115 = t95 * t96;
t98 = cos(qJ(1));
t113 = t97 * t98;
t112 = t98 * t78;
t111 = t98 * t79;
t110 = t98 * t81;
t109 = t98 * t82;
t108 = t98 * t93;
t107 = t98 * t96;
t92 = -qJ(5) - pkin(8);
t106 = t98 * pkin(1) + t95 * pkin(7);
t104 = t94 * t116;
t86 = t95 * pkin(1);
t103 = t86 + t126;
t102 = -t98 * pkin(7) + t86;
t101 = pkin(2) * t113 + t98 * t105 + t106;
t100 = -t97 * qJ(3) + t122;
t99 = g(1) * t98 + g(2) * t95;
t90 = -pkin(9) + t92;
t73 = g(1) * t95 - g(2) * t98;
t71 = pkin(5) * t82 + t80;
t70 = g(3) * t94 + t99 * t97;
t69 = t99 * t94 - t124;
t1 = [0, 0, 0, 0, 0, 0, -t99, t73, -g(3), -t125, 0, 0, 0, 0, 0, 0, -t70, t69, -t73, -g(1) * t106 - g(2) * t102 - t125, 0, 0, 0, 0, 0, 0, -t73, t70, -t69, -g(1) * t101 - g(2) * (t102 + t126) - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t108 + t115) - g(2) * (t104 - t107) + t93 * t124, -g(1) * (t94 * t107 - t116) - g(2) * (t94 * t115 + t108) + t96 * t124, -t70, -g(1) * (t95 * pkin(3) + pkin(8) * t113 + t101) - g(2) * (pkin(8) * t114 + (-pkin(3) - pkin(7)) * t98 + t103) - g(3) * (t94 * pkin(8) + t100) 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t110 + t117) - g(2) * (t94 * t118 - t109) + t81 * t124, -g(1) * (t94 * t109 - t118) - g(2) * (t94 * t117 + t110) + t82 * t124, -t70, -g(1) * (t95 * t80 + t101) - g(2) * (pkin(4) * t104 - t92 * t114 + t103) - g(3) * (-t94 * t92 + (-qJ(3) - t123) * t97 + t122) + (-g(1) * (t94 * t123 - t92 * t97) - g(2) * (-pkin(7) - t80)) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t112 + t119) - g(2) * (t94 * t120 - t111) + t78 * t124, -g(1) * (t94 * t111 - t120) - g(2) * (t94 * t119 + t112) + t79 * t124, -t70, -g(1) * (t95 * t71 + t101) - g(2) * (-t90 * t114 + t95 * t121 + t103) - g(3) * (-t94 * t90 + (-qJ(3) - t72) * t97 + t122) + (-g(1) * (-t90 * t97 + t121) - g(2) * (-pkin(7) - t71)) * t98;];
U_reg  = t1;

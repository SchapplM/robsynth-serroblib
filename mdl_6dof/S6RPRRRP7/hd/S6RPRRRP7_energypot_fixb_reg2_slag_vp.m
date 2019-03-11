% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:20:59
% EndTime: 2019-03-09 06:20:59
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->70), mult. (215->87), div. (0->0), fcn. (215->10), ass. (0->44)
t96 = cos(qJ(4));
t82 = t96 * pkin(4) + pkin(3);
t89 = pkin(10) + qJ(3);
t83 = sin(t89);
t84 = cos(t89);
t98 = -pkin(9) - pkin(8);
t120 = t82 * t84 - t83 * t98;
t119 = g(3) * pkin(6);
t118 = g(3) * t83;
t91 = sin(pkin(10));
t117 = t91 * pkin(2) + pkin(6);
t90 = qJ(4) + qJ(5);
t85 = sin(t90);
t95 = sin(qJ(1));
t114 = t95 * t85;
t86 = cos(t90);
t113 = t95 * t86;
t94 = sin(qJ(4));
t112 = t95 * t94;
t111 = t95 * t96;
t97 = cos(qJ(1));
t110 = t97 * t85;
t109 = t97 * t86;
t108 = t97 * t94;
t107 = t97 * t96;
t92 = cos(pkin(10));
t80 = t92 * pkin(2) + pkin(1);
t93 = -pkin(7) - qJ(2);
t106 = t95 * t80 + t97 * t93;
t105 = t97 * t80 - t95 * t93;
t104 = t83 * t82 + t84 * t98 + t117;
t103 = pkin(3) * t84 + pkin(8) * t83;
t102 = g(1) * t97 + g(2) * t95;
t66 = t84 * t114 + t109;
t68 = t84 * t110 - t113;
t101 = g(1) * t68 + g(2) * t66 + t85 * t118;
t100 = pkin(4) * t112 + t120 * t97 + t105;
t99 = -pkin(4) * t108 + t120 * t95 + t106;
t75 = g(1) * t95 - g(2) * t97;
t69 = t84 * t109 + t114;
t67 = t84 * t113 - t110;
t65 = -g(3) * t84 + t102 * t83;
t64 = -g(1) * t69 - g(2) * t67 - t86 * t118;
t1 = [0, 0, 0, 0, 0, 0, -t102, t75, -g(3), -t119, 0, 0, 0, 0, 0, 0, -g(3) * t91 - t102 * t92, -g(3) * t92 + t102 * t91, -t75, -g(1) * (t97 * pkin(1) + t95 * qJ(2)) - g(2) * (t95 * pkin(1) - t97 * qJ(2)) - t119, 0, 0, 0, 0, 0, 0, -t102 * t84 - t118, t65, -t75, -g(1) * t105 - g(2) * t106 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t107 + t112) - g(2) * (t84 * t111 - t108) - t96 * t118, -g(1) * (-t84 * t108 + t111) - g(2) * (-t84 * t112 - t107) + t94 * t118, -t65, -g(1) * (t103 * t97 + t105) - g(2) * (t103 * t95 + t106) - g(3) * (t83 * pkin(3) - t84 * pkin(8) + t117) 0, 0, 0, 0, 0, 0, t64, t101, -t65, -g(1) * t100 - g(2) * t99 - g(3) * t104, 0, 0, 0, 0, 0, 0, t64, -t65, -t101, -g(1) * (t69 * pkin(5) + t68 * qJ(6) + t100) - g(2) * (t67 * pkin(5) + t66 * qJ(6) + t99) - g(3) * ((pkin(5) * t86 + qJ(6) * t85) * t83 + t104);];
U_reg  = t1;

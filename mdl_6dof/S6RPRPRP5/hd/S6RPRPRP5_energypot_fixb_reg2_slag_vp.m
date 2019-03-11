% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:18
% EndTime: 2019-03-09 03:16:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->70), mult. (215->87), div. (0->0), fcn. (215->10), ass. (0->44)
t93 = cos(pkin(10));
t80 = t93 * pkin(4) + pkin(3);
t90 = pkin(9) + qJ(3);
t84 = sin(t90);
t86 = cos(t90);
t95 = -pkin(8) - qJ(4);
t120 = t80 * t86 - t84 * t95;
t119 = g(3) * pkin(6);
t118 = g(3) * t84;
t92 = sin(pkin(9));
t117 = t92 * pkin(2) + pkin(6);
t89 = pkin(10) + qJ(5);
t83 = sin(t89);
t97 = sin(qJ(1));
t114 = t97 * t83;
t85 = cos(t89);
t113 = t97 * t85;
t91 = sin(pkin(10));
t112 = t97 * t91;
t111 = t97 * t93;
t98 = cos(qJ(1));
t110 = t98 * t83;
t109 = t98 * t85;
t108 = t98 * t91;
t107 = t98 * t93;
t94 = cos(pkin(9));
t81 = t94 * pkin(2) + pkin(1);
t96 = -pkin(7) - qJ(2);
t106 = t97 * t81 + t98 * t96;
t105 = t98 * t81 - t97 * t96;
t104 = t84 * t80 + t86 * t95 + t117;
t103 = g(1) * t98 + g(2) * t97;
t102 = pkin(3) * t86 + qJ(4) * t84;
t66 = t86 * t114 + t109;
t68 = t86 * t110 - t113;
t101 = g(1) * t68 + g(2) * t66 + t83 * t118;
t100 = pkin(4) * t112 + t120 * t98 + t105;
t99 = -pkin(4) * t108 + t120 * t97 + t106;
t75 = g(1) * t97 - g(2) * t98;
t69 = t86 * t109 + t114;
t67 = t86 * t113 - t110;
t65 = -g(3) * t86 + t103 * t84;
t64 = -g(1) * t69 - g(2) * t67 - t85 * t118;
t1 = [0, 0, 0, 0, 0, 0, -t103, t75, -g(3), -t119, 0, 0, 0, 0, 0, 0, -g(3) * t92 - t103 * t94, -g(3) * t94 + t103 * t92, -t75, -g(1) * (t98 * pkin(1) + t97 * qJ(2)) - g(2) * (t97 * pkin(1) - t98 * qJ(2)) - t119, 0, 0, 0, 0, 0, 0, -t103 * t86 - t118, t65, -t75, -g(1) * t105 - g(2) * t106 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t107 + t112) - g(2) * (t86 * t111 - t108) - t93 * t118, -g(1) * (-t86 * t108 + t111) - g(2) * (-t86 * t112 - t107) + t91 * t118, -t65, -g(1) * (t102 * t98 + t105) - g(2) * (t102 * t97 + t106) - g(3) * (t84 * pkin(3) - t86 * qJ(4) + t117) 0, 0, 0, 0, 0, 0, t64, t101, -t65, -g(1) * t100 - g(2) * t99 - g(3) * t104, 0, 0, 0, 0, 0, 0, t64, -t65, -t101, -g(1) * (t69 * pkin(5) + t68 * qJ(6) + t100) - g(2) * (t67 * pkin(5) + t66 * qJ(6) + t99) - g(3) * ((pkin(5) * t85 + qJ(6) * t83) * t84 + t104);];
U_reg  = t1;

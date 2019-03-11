% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:45
% EndTime: 2019-03-09 04:40:45
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->70), mult. (215->87), div. (0->0), fcn. (215->10), ass. (0->44)
t96 = cos(qJ(4));
t81 = t96 * pkin(4) + pkin(3);
t88 = pkin(9) + qJ(3);
t82 = sin(t88);
t84 = cos(t88);
t92 = -qJ(5) - pkin(8);
t119 = t81 * t84 - t82 * t92;
t118 = g(3) * pkin(6);
t117 = g(3) * t82;
t90 = sin(pkin(9));
t116 = t90 * pkin(2) + pkin(6);
t89 = qJ(4) + pkin(10);
t83 = sin(t89);
t95 = sin(qJ(1));
t113 = t95 * t83;
t85 = cos(t89);
t112 = t95 * t85;
t94 = sin(qJ(4));
t111 = t95 * t94;
t110 = t95 * t96;
t97 = cos(qJ(1));
t109 = t97 * t83;
t108 = t97 * t85;
t107 = t97 * t94;
t106 = t97 * t96;
t91 = cos(pkin(9));
t79 = t91 * pkin(2) + pkin(1);
t93 = -pkin(7) - qJ(2);
t105 = t95 * t79 + t97 * t93;
t104 = t97 * t79 - t95 * t93;
t103 = t82 * t81 + t84 * t92 + t116;
t102 = pkin(3) * t84 + pkin(8) * t82;
t101 = g(1) * t97 + g(2) * t95;
t65 = t84 * t113 + t108;
t67 = t84 * t109 - t112;
t100 = g(1) * t67 + g(2) * t65 + t83 * t117;
t99 = pkin(4) * t111 + t119 * t97 + t104;
t98 = -pkin(4) * t107 + t119 * t95 + t105;
t74 = g(1) * t95 - g(2) * t97;
t68 = t84 * t108 + t113;
t66 = t84 * t112 - t109;
t64 = -g(3) * t84 + t101 * t82;
t63 = -g(1) * t68 - g(2) * t66 - t85 * t117;
t1 = [0, 0, 0, 0, 0, 0, -t101, t74, -g(3), -t118, 0, 0, 0, 0, 0, 0, -g(3) * t90 - t101 * t91, -g(3) * t91 + t101 * t90, -t74, -g(1) * (t97 * pkin(1) + t95 * qJ(2)) - g(2) * (t95 * pkin(1) - t97 * qJ(2)) - t118, 0, 0, 0, 0, 0, 0, -t101 * t84 - t117, t64, -t74, -g(1) * t104 - g(2) * t105 - g(3) * t116, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t106 + t111) - g(2) * (t84 * t110 - t107) - t96 * t117, -g(1) * (-t84 * t107 + t110) - g(2) * (-t84 * t111 - t106) + t94 * t117, -t64, -g(1) * (t102 * t97 + t104) - g(2) * (t102 * t95 + t105) - g(3) * (t82 * pkin(3) - t84 * pkin(8) + t116) 0, 0, 0, 0, 0, 0, t63, t100, -t64, -g(1) * t99 - g(2) * t98 - g(3) * t103, 0, 0, 0, 0, 0, 0, t63, -t64, -t100, -g(1) * (t68 * pkin(5) + t67 * qJ(6) + t99) - g(2) * (t66 * pkin(5) + t65 * qJ(6) + t98) - g(3) * ((pkin(5) * t85 + qJ(6) * t83) * t82 + t103);];
U_reg  = t1;

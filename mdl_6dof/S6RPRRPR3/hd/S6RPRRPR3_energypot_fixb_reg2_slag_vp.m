% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:14
% EndTime: 2019-03-09 05:07:14
% DurationCPUTime: 0.17s
% Computational Cost: add. (227->67), mult. (248->89), div. (0->0), fcn. (264->10), ass. (0->41)
t92 = sin(qJ(3));
t95 = cos(qJ(4));
t110 = t92 * t95;
t91 = sin(qJ(4));
t112 = t91 * t92;
t118 = pkin(4) * t110 + qJ(5) * t112;
t96 = cos(qJ(3));
t117 = pkin(3) * t96;
t89 = qJ(2) + pkin(6);
t116 = g(3) * t89;
t115 = g(3) * t92;
t88 = qJ(1) + pkin(10);
t82 = sin(t88);
t114 = t82 * t92;
t83 = cos(t88);
t113 = t83 * t92;
t111 = t91 * t96;
t109 = t95 * t96;
t108 = t92 * pkin(3) + t89;
t97 = cos(qJ(1));
t107 = t97 * pkin(1) + t83 * pkin(2) + t82 * pkin(7);
t93 = sin(qJ(1));
t106 = t93 * pkin(1) + t82 * pkin(2) - t83 * pkin(7);
t105 = pkin(8) * t113 + t83 * t117 + t107;
t104 = g(1) * t83 + g(2) * t82;
t103 = -g(1) * t97 - g(2) * t93;
t102 = -t96 * pkin(8) + t108;
t101 = pkin(8) * t114 + t82 * t117 + t106;
t67 = t82 * t111 + t83 * t95;
t69 = t83 * t111 - t82 * t95;
t100 = g(1) * t69 + g(2) * t67 + g(3) * t112;
t70 = t83 * t109 + t82 * t91;
t99 = t70 * pkin(4) + t69 * qJ(5) + t105;
t68 = t82 * t109 - t83 * t91;
t98 = t68 * pkin(4) + t67 * qJ(5) + t101;
t94 = cos(qJ(6));
t90 = sin(qJ(6));
t71 = g(1) * t82 - g(2) * t83;
t64 = -g(3) * t96 + t104 * t92;
t63 = -g(1) * t70 - g(2) * t68 - g(3) * t110;
t1 = [0, 0, 0, 0, 0, 0, t103, g(1) * t93 - g(2) * t97, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t104, t71, -g(3), t103 * pkin(1) - t116, 0, 0, 0, 0, 0, 0, -t104 * t96 - t115, t64, -t71, -g(1) * t107 - g(2) * t106 - t116, 0, 0, 0, 0, 0, 0, t63, t100, -t64, -g(1) * t105 - g(2) * t101 - g(3) * t102, 0, 0, 0, 0, 0, 0, t63, -t64, -t100, -g(1) * t99 - g(2) * t98 - g(3) * (t102 + t118) 0, 0, 0, 0, 0, 0, -g(1) * (t69 * t90 + t70 * t94) - g(2) * (t67 * t90 + t68 * t94) - (t90 * t91 + t94 * t95) * t115, -g(1) * (t69 * t94 - t70 * t90) - g(2) * (t67 * t94 - t68 * t90) - (-t90 * t95 + t91 * t94) * t115, t64, -g(1) * (t70 * pkin(5) - pkin(9) * t113 + t99) - g(2) * (t68 * pkin(5) - pkin(9) * t114 + t98) - g(3) * (pkin(5) * t110 + (-pkin(8) + pkin(9)) * t96 + t108 + t118);];
U_reg  = t1;

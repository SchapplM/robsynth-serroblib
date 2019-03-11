% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP6
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:43
% EndTime: 2019-03-09 06:16:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (202->77), mult. (200->96), div. (0->0), fcn. (196->10), ass. (0->44)
t113 = g(3) * pkin(6);
t94 = -pkin(9) - pkin(8);
t90 = sin(qJ(4));
t112 = pkin(4) * t90;
t85 = pkin(10) + qJ(3);
t77 = sin(t85);
t111 = g(3) * t77;
t87 = sin(pkin(10));
t110 = t87 * pkin(2) + pkin(6);
t92 = cos(qJ(4));
t76 = t92 * pkin(4) + pkin(3);
t84 = -qJ(6) + t94;
t109 = t77 * t84;
t108 = t77 * t94;
t78 = cos(t85);
t93 = cos(qJ(1));
t107 = t78 * t93;
t86 = qJ(4) + qJ(5);
t79 = sin(t86);
t91 = sin(qJ(1));
t106 = t79 * t91;
t105 = t79 * t93;
t80 = cos(t86);
t104 = t80 * t91;
t103 = t90 * t91;
t102 = t90 * t93;
t101 = t91 * t92;
t100 = t92 * t93;
t99 = t93 * t80;
t88 = cos(pkin(10));
t74 = pkin(2) * t88 + pkin(1);
t89 = -pkin(7) - qJ(2);
t98 = t91 * t74 + t93 * t89;
t71 = t93 * t74;
t97 = -t91 * t89 + t71;
t96 = pkin(3) * t78 + pkin(8) * t77;
t95 = g(1) * t93 + g(2) * t91;
t72 = g(1) * t91 - g(2) * t93;
t69 = pkin(5) * t79 + t112;
t68 = pkin(5) * t80 + t76;
t67 = -g(3) * t78 + t95 * t77;
t66 = -g(1) * (t78 * t99 + t106) - g(2) * (t78 * t104 - t105) - t80 * t111;
t65 = -g(1) * (-t78 * t105 + t104) - g(2) * (-t78 * t106 - t99) + t79 * t111;
t1 = [0, 0, 0, 0, 0, 0, -t95, t72, -g(3), -t113, 0, 0, 0, 0, 0, 0, -g(3) * t87 - t95 * t88, -g(3) * t88 + t95 * t87, -t72, -g(1) * (pkin(1) * t93 + qJ(2) * t91) - g(2) * (pkin(1) * t91 - qJ(2) * t93) - t113, 0, 0, 0, 0, 0, 0, -t95 * t78 - t111, t67, -t72, -g(1) * t97 - g(2) * t98 - g(3) * t110, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t100 + t103) - g(2) * (t78 * t101 - t102) - t92 * t111, -g(1) * (-t78 * t102 + t101) - g(2) * (-t78 * t103 - t100) + t90 * t111, -t67, -g(1) * (t96 * t93 + t97) - g(2) * (t96 * t91 + t98) - g(3) * (pkin(3) * t77 - pkin(8) * t78 + t110) 0, 0, 0, 0, 0, 0, t66, t65, -t67, -g(1) * (t76 * t107 - t93 * t108 + t71) - g(2) * (-pkin(4) * t102 + t98) - g(3) * (t77 * t76 + t78 * t94 + t110) + (-g(1) * (-t89 + t112) - g(2) * (t76 * t78 - t108)) * t91, 0, 0, 0, 0, 0, 0, t66, t65, -t67, -g(1) * (t68 * t107 - t93 * t109 + t71) - g(2) * (-t93 * t69 + t98) - g(3) * (t68 * t77 + t78 * t84 + t110) + (-g(1) * (t69 - t89) - g(2) * (t68 * t78 - t109)) * t91;];
U_reg  = t1;

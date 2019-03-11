% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:44
% EndTime: 2019-03-09 03:19:44
% DurationCPUTime: 0.14s
% Computational Cost: add. (172->67), mult. (192->77), div. (0->0), fcn. (184->8), ass. (0->40)
t84 = pkin(9) + qJ(3);
t80 = sin(t84);
t100 = qJ(4) * t80;
t81 = cos(t84);
t92 = cos(qJ(1));
t106 = t81 * t92;
t112 = pkin(3) * t106 + t92 * t100;
t111 = g(3) * pkin(6);
t89 = sin(qJ(5));
t110 = pkin(5) * t89;
t109 = g(3) * t81;
t85 = sin(pkin(9));
t108 = t85 * pkin(2) + pkin(6);
t90 = sin(qJ(1));
t107 = t81 * t90;
t105 = t90 * t89;
t91 = cos(qJ(5));
t104 = t90 * t91;
t103 = t91 * t92;
t102 = t92 * t89;
t86 = cos(pkin(9));
t77 = pkin(2) * t86 + pkin(1);
t88 = -pkin(7) - qJ(2);
t101 = t90 * t77 + t92 * t88;
t99 = t80 * pkin(3) + t108;
t98 = t80 * t102;
t72 = t92 * t77;
t97 = t72 + t112;
t96 = -t90 * t88 + t72;
t95 = pkin(3) * t107 + t90 * t100 + t101;
t94 = g(1) * t92 + g(2) * t90;
t93 = -t81 * qJ(4) + t99;
t87 = -qJ(6) - pkin(8);
t79 = pkin(5) * t91 + pkin(4);
t73 = g(1) * t90 - g(2) * t92;
t68 = g(3) * t80 + t94 * t81;
t67 = t94 * t80 - t109;
t66 = -g(1) * (t98 + t104) - g(2) * (t80 * t105 - t103) + t89 * t109;
t65 = -g(1) * (t80 * t103 - t105) - g(2) * (t80 * t104 + t102) + t91 * t109;
t1 = [0, 0, 0, 0, 0, 0, -t94, t73, -g(3), -t111, 0, 0, 0, 0, 0, 0, -g(3) * t85 - t94 * t86, -g(3) * t86 + t94 * t85, -t73, -g(1) * (pkin(1) * t92 + t90 * qJ(2)) - g(2) * (t90 * pkin(1) - qJ(2) * t92) - t111, 0, 0, 0, 0, 0, 0, -t68, t67, -t73, -g(1) * t96 - g(2) * t101 - g(3) * t108, 0, 0, 0, 0, 0, 0, -t73, t68, -t67, -g(1) * (t96 + t112) - g(2) * t95 - g(3) * t93, 0, 0, 0, 0, 0, 0, t66, t65, -t68, -g(1) * (pkin(8) * t106 + (pkin(4) - t88) * t90 + t97) - g(2) * (-pkin(4) * t92 + pkin(8) * t107 + t95) - g(3) * (pkin(8) * t80 + t93) 0, 0, 0, 0, 0, 0, t66, t65, -t68, -g(1) * (pkin(5) * t98 - t87 * t106 + t97) - g(2) * (-t92 * t79 + t95) - g(3) * (-t80 * t87 + (-qJ(4) - t110) * t81 + t99) + (-g(1) * (t79 - t88) - g(2) * (t80 * t110 - t81 * t87)) * t90;];
U_reg  = t1;

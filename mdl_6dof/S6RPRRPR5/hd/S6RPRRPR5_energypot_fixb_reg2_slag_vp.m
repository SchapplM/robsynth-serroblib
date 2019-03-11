% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR5
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:33
% EndTime: 2019-03-09 05:13:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (193->63), mult. (161->75), div. (0->0), fcn. (149->10), ass. (0->38)
t89 = pkin(10) + qJ(3);
t84 = qJ(4) + t89;
t79 = sin(t84);
t102 = qJ(5) * t79;
t80 = cos(t84);
t96 = cos(qJ(1));
t108 = t80 * t96;
t113 = pkin(4) * t108 + t96 * t102;
t112 = g(3) * pkin(6);
t111 = g(3) * t80;
t90 = sin(pkin(10));
t110 = t90 * pkin(2) + pkin(6);
t91 = cos(pkin(10));
t81 = t91 * pkin(2) + pkin(1);
t94 = sin(qJ(1));
t109 = t80 * t94;
t93 = sin(qJ(6));
t107 = t94 * t93;
t95 = cos(qJ(6));
t106 = t94 * t95;
t105 = t96 * t93;
t104 = t96 * t95;
t92 = -pkin(7) - qJ(2);
t83 = cos(t89);
t70 = pkin(3) * t83 + t81;
t88 = -pkin(8) + t92;
t103 = t94 * t70 + t96 * t88;
t82 = sin(t89);
t101 = pkin(3) * t82 + t110;
t69 = t96 * t70;
t100 = -t94 * t88 + t69;
t99 = pkin(4) * t109 + t94 * t102 + t103;
t98 = g(1) * t96 + g(2) * t94;
t97 = t79 * pkin(4) - t80 * qJ(5) + t101;
t75 = g(1) * t94 - g(2) * t96;
t67 = g(3) * t79 + t98 * t80;
t66 = t98 * t79 - t111;
t1 = [0, 0, 0, 0, 0, 0, -t98, t75, -g(3), -t112, 0, 0, 0, 0, 0, 0, -g(3) * t90 - t98 * t91, -g(3) * t91 + t98 * t90, -t75, -g(1) * (t96 * pkin(1) + t94 * qJ(2)) - g(2) * (t94 * pkin(1) - t96 * qJ(2)) - t112, 0, 0, 0, 0, 0, 0, -g(3) * t82 - t98 * t83, -g(3) * t83 + t98 * t82, -t75, -g(1) * (t96 * t81 - t94 * t92) - g(2) * (t94 * t81 + t96 * t92) - g(3) * t110, 0, 0, 0, 0, 0, 0, -t67, t66, -t75, -g(1) * t100 - g(2) * t103 - g(3) * t101, 0, 0, 0, 0, 0, 0, -t75, t67, -t66, -g(1) * (t100 + t113) - g(2) * t99 - g(3) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t105 + t106) - g(2) * (t79 * t107 - t104) + t93 * t111, -g(1) * (t79 * t104 - t107) - g(2) * (t79 * t106 + t105) + t95 * t111, -t67, -g(1) * (pkin(9) * t108 + t69 + (pkin(5) - t88) * t94 + t113) - g(2) * (-t96 * pkin(5) + pkin(9) * t109 + t99) - g(3) * (t79 * pkin(9) + t97);];
U_reg  = t1;

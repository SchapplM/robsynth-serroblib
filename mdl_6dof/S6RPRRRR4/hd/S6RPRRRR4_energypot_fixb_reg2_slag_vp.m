% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:53
% EndTime: 2019-03-09 07:05:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (207->61), mult. (150->77), div. (0->0), fcn. (138->12), ass. (0->37)
t111 = g(3) * pkin(6);
t91 = pkin(11) + qJ(3);
t85 = qJ(4) + t91;
t82 = qJ(5) + t85;
t75 = sin(t82);
t110 = g(3) * t75;
t92 = sin(pkin(11));
t109 = t92 * pkin(2) + pkin(6);
t93 = cos(pkin(11));
t81 = t93 * pkin(2) + pkin(1);
t95 = sin(qJ(6));
t96 = sin(qJ(1));
t108 = t96 * t95;
t97 = cos(qJ(6));
t107 = t96 * t97;
t98 = cos(qJ(1));
t106 = t98 * t95;
t105 = t98 * t97;
t94 = -pkin(7) - qJ(2);
t84 = cos(t91);
t71 = pkin(3) * t84 + t81;
t80 = cos(t85);
t70 = pkin(4) * t80 + t71;
t90 = -pkin(8) + t94;
t86 = -pkin(9) + t90;
t104 = t96 * t70 + t98 * t86;
t83 = sin(t91);
t103 = pkin(3) * t83 + t109;
t102 = t98 * t70 - t96 * t86;
t79 = sin(t85);
t101 = pkin(4) * t79 + t103;
t76 = cos(t82);
t100 = pkin(5) * t76 + pkin(10) * t75;
t99 = g(1) * t98 + g(2) * t96;
t72 = g(1) * t96 - g(2) * t98;
t67 = -g(3) * t76 + t99 * t75;
t1 = [0, 0, 0, 0, 0, 0, -t99, t72, -g(3), -t111, 0, 0, 0, 0, 0, 0, -g(3) * t92 - t99 * t93, -g(3) * t93 + t99 * t92, -t72, -g(1) * (t98 * pkin(1) + t96 * qJ(2)) - g(2) * (t96 * pkin(1) - t98 * qJ(2)) - t111, 0, 0, 0, 0, 0, 0, -g(3) * t83 - t99 * t84, -g(3) * t84 + t99 * t83, -t72, -g(1) * (t98 * t81 - t96 * t94) - g(2) * (t96 * t81 + t98 * t94) - g(3) * t109, 0, 0, 0, 0, 0, 0, -g(3) * t79 - t99 * t80, -g(3) * t80 + t99 * t79, -t72, -g(1) * (t98 * t71 - t96 * t90) - g(2) * (t96 * t71 + t98 * t90) - g(3) * t103, 0, 0, 0, 0, 0, 0, -t99 * t76 - t110, t67, -t72, -g(1) * t102 - g(2) * t104 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (t76 * t105 + t108) - g(2) * (t76 * t107 - t106) - t97 * t110, -g(1) * (-t76 * t106 + t107) - g(2) * (-t76 * t108 - t105) + t95 * t110, -t67, -g(1) * (t100 * t98 + t102) - g(2) * (t100 * t96 + t104) - g(3) * (t75 * pkin(5) - t76 * pkin(10) + t101);];
U_reg  = t1;

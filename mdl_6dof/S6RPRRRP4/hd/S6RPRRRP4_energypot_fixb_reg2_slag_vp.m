% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:44
% EndTime: 2019-03-09 06:08:44
% DurationCPUTime: 0.17s
% Computational Cost: add. (201->63), mult. (174->78), div. (0->0), fcn. (166->10), ass. (0->37)
t86 = pkin(10) + qJ(3);
t81 = qJ(4) + t86;
t75 = sin(t81);
t76 = cos(t81);
t93 = cos(qJ(5));
t78 = pkin(5) * t93 + pkin(4);
t89 = -qJ(6) - pkin(9);
t109 = -t75 * t89 + t76 * t78;
t108 = g(3) * pkin(6);
t107 = g(3) * t75;
t87 = sin(pkin(10));
t106 = t87 * pkin(2) + pkin(6);
t88 = cos(pkin(10));
t77 = t88 * pkin(2) + pkin(1);
t91 = sin(qJ(5));
t94 = cos(qJ(1));
t103 = t91 * t94;
t92 = sin(qJ(1));
t102 = t92 * t91;
t101 = t92 * t93;
t100 = t94 * t93;
t90 = -pkin(7) - qJ(2);
t80 = cos(t86);
t70 = pkin(3) * t80 + t77;
t85 = -pkin(8) + t90;
t99 = t92 * t70 + t94 * t85;
t79 = sin(t86);
t98 = pkin(3) * t79 + t106;
t69 = t94 * t70;
t97 = -t92 * t85 + t69;
t96 = pkin(4) * t76 + pkin(9) * t75;
t95 = g(1) * t94 + g(2) * t92;
t71 = g(1) * t92 - g(2) * t94;
t67 = -g(3) * t76 + t95 * t75;
t66 = -g(1) * (t76 * t100 + t102) - g(2) * (t76 * t101 - t103) - t93 * t107;
t65 = -g(1) * (-t76 * t103 + t101) - g(2) * (-t76 * t102 - t100) + t91 * t107;
t1 = [0, 0, 0, 0, 0, 0, -t95, t71, -g(3), -t108, 0, 0, 0, 0, 0, 0, -g(3) * t87 - t95 * t88, -g(3) * t88 + t95 * t87, -t71, -g(1) * (pkin(1) * t94 + t92 * qJ(2)) - g(2) * (t92 * pkin(1) - qJ(2) * t94) - t108, 0, 0, 0, 0, 0, 0, -g(3) * t79 - t95 * t80, -g(3) * t80 + t95 * t79, -t71, -g(1) * (t77 * t94 - t92 * t90) - g(2) * (t92 * t77 + t90 * t94) - g(3) * t106, 0, 0, 0, 0, 0, 0, -t95 * t76 - t107, t67, -t71, -g(1) * t97 - g(2) * t99 - g(3) * t98, 0, 0, 0, 0, 0, 0, t66, t65, -t67, -g(1) * (t96 * t94 + t97) - g(2) * (t96 * t92 + t99) - g(3) * (pkin(4) * t75 - pkin(9) * t76 + t98) 0, 0, 0, 0, 0, 0, t66, t65, -t67, -g(1) * (t109 * t94 + t69) - g(2) * (-pkin(5) * t103 + t99) - g(3) * (t75 * t78 + t76 * t89 + t98) + (-g(1) * (pkin(5) * t91 - t85) - g(2) * t109) * t92;];
U_reg  = t1;

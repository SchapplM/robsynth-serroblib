% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:36
% EndTime: 2019-03-09 02:51:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (182->76), mult. (192->91), div. (0->0), fcn. (184->10), ass. (0->45)
t87 = pkin(9) + qJ(3);
t81 = sin(t87);
t104 = qJ(4) * t81;
t83 = cos(t87);
t95 = cos(qJ(1));
t114 = t83 * t95;
t119 = pkin(3) * t114 + t95 * t104;
t118 = g(3) * pkin(6);
t88 = sin(pkin(10));
t117 = pkin(5) * t88;
t116 = g(3) * t83;
t89 = sin(pkin(9));
t115 = t89 * pkin(2) + pkin(6);
t86 = pkin(10) + qJ(6);
t80 = sin(t86);
t94 = sin(qJ(1));
t113 = t94 * t80;
t82 = cos(t86);
t112 = t94 * t82;
t111 = t94 * t88;
t90 = cos(pkin(10));
t110 = t94 * t90;
t109 = t95 * t80;
t108 = t95 * t82;
t107 = t95 * t88;
t106 = t95 * t90;
t91 = cos(pkin(9));
t78 = t91 * pkin(2) + pkin(1);
t93 = -pkin(7) - qJ(2);
t105 = t94 * t78 + t95 * t93;
t103 = qJ(5) * t83;
t102 = t81 * pkin(3) + t115;
t101 = t81 * t107;
t72 = t95 * t78;
t100 = t72 + t119;
t99 = -t94 * t93 + t72;
t98 = t105 + (pkin(3) * t83 + t104) * t94;
t97 = g(1) * t95 + g(2) * t94;
t96 = -t83 * qJ(4) + t102;
t92 = -pkin(8) - qJ(5);
t77 = t90 * pkin(5) + pkin(4);
t73 = g(1) * t94 - g(2) * t95;
t68 = g(3) * t81 + t97 * t83;
t67 = t97 * t81 - t116;
t1 = [0, 0, 0, 0, 0, 0, -t97, t73, -g(3), -t118, 0, 0, 0, 0, 0, 0, -g(3) * t89 - t97 * t91, -g(3) * t91 + t97 * t89, -t73, -g(1) * (t95 * pkin(1) + t94 * qJ(2)) - g(2) * (t94 * pkin(1) - t95 * qJ(2)) - t118, 0, 0, 0, 0, 0, 0, -t68, t67, -t73, -g(1) * t99 - g(2) * t105 - g(3) * t115, 0, 0, 0, 0, 0, 0, -t73, t68, -t67, -g(1) * (t99 + t119) - g(2) * t98 - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t101 + t110) - g(2) * (t81 * t111 - t106) + t88 * t116, -g(1) * (t81 * t106 - t111) - g(2) * (t81 * t110 + t107) + t90 * t116, -t68, -g(1) * (t95 * t103 + (pkin(4) - t93) * t94 + t100) - g(2) * (-t95 * pkin(4) + t94 * t103 + t98) - g(3) * (t81 * qJ(5) + t96) 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t109 + t112) - g(2) * (t81 * t113 - t108) + t80 * t116, -g(1) * (t81 * t108 - t113) - g(2) * (t81 * t112 + t109) + t82 * t116, -t68, -g(1) * (pkin(5) * t101 - t92 * t114 + t100) - g(2) * (-t95 * t77 + t98) - g(3) * (-t81 * t92 + (-qJ(4) - t117) * t83 + t102) + (-g(1) * (t77 - t93) - g(2) * (t81 * t117 - t83 * t92)) * t94;];
U_reg  = t1;

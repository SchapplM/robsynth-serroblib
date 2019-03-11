% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP5
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:24
% EndTime: 2019-03-09 06:12:24
% DurationCPUTime: 0.14s
% Computational Cost: add. (215->60), mult. (187->74), div. (0->0), fcn. (183->10), ass. (0->40)
t93 = pkin(10) + qJ(3);
t88 = qJ(4) + t93;
t82 = sin(t88);
t83 = cos(t88);
t118 = pkin(4) * t83 + pkin(9) * t82;
t117 = g(3) * pkin(6);
t114 = g(3) * t82;
t94 = sin(pkin(10));
t113 = t94 * pkin(2) + pkin(6);
t95 = cos(pkin(10));
t84 = t95 * pkin(2) + pkin(1);
t97 = sin(qJ(5));
t98 = sin(qJ(1));
t112 = t98 * t97;
t99 = cos(qJ(5));
t111 = t98 * t99;
t96 = -pkin(7) - qJ(2);
t100 = cos(qJ(1));
t87 = cos(t93);
t72 = pkin(3) * t87 + t84;
t92 = -pkin(8) + t96;
t110 = t100 * t92 + t98 * t72;
t109 = t100 * t97;
t108 = t100 * t99;
t86 = sin(t93);
t107 = pkin(3) * t86 + t113;
t106 = t100 * t72 - t98 * t92;
t105 = t118 * t98 + t110;
t104 = g(1) * t100 + g(2) * t98;
t103 = t118 * t100 + t106;
t102 = t82 * pkin(4) - t83 * pkin(9) + t107;
t66 = t83 * t112 + t108;
t68 = t83 * t109 - t111;
t101 = g(1) * t68 + g(2) * t66 + t97 * t114;
t77 = g(1) * t98 - g(2) * t100;
t69 = t83 * t108 + t112;
t67 = t83 * t111 - t109;
t65 = -g(3) * t83 + t104 * t82;
t64 = -g(1) * t69 - g(2) * t67 - t99 * t114;
t1 = [0, 0, 0, 0, 0, 0, -t104, t77, -g(3), -t117, 0, 0, 0, 0, 0, 0, -g(3) * t94 - t104 * t95, -g(3) * t95 + t104 * t94, -t77, -g(1) * (t100 * pkin(1) + t98 * qJ(2)) - g(2) * (t98 * pkin(1) - t100 * qJ(2)) - t117, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t104 * t87, -g(3) * t87 + t104 * t86, -t77, -g(1) * (t100 * t84 - t98 * t96) - g(2) * (t100 * t96 + t98 * t84) - g(3) * t113, 0, 0, 0, 0, 0, 0, -t104 * t83 - t114, t65, -t77, -g(1) * t106 - g(2) * t110 - g(3) * t107, 0, 0, 0, 0, 0, 0, t64, t101, -t65, -g(1) * t103 - g(2) * t105 - g(3) * t102, 0, 0, 0, 0, 0, 0, t64, -t65, -t101, -g(1) * (t69 * pkin(5) + t68 * qJ(6) + t103) - g(2) * (pkin(5) * t67 + qJ(6) * t66 + t105) - g(3) * ((pkin(5) * t99 + qJ(6) * t97) * t82 + t102);];
U_reg  = t1;

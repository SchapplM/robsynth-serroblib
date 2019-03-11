% Calculate inertial parameters regressor of potential energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:13
% EndTime: 2019-03-08 19:16:14
% DurationCPUTime: 0.28s
% Computational Cost: add. (359->99), mult. (712->154), div. (0->0), fcn. (880->14), ass. (0->55)
t100 = sin(pkin(11));
t104 = cos(pkin(11));
t109 = sin(qJ(2));
t111 = cos(qJ(2));
t83 = -t109 * t100 + t111 * t104;
t101 = sin(pkin(10));
t105 = cos(pkin(10));
t102 = sin(pkin(6));
t106 = cos(pkin(6));
t126 = t106 * t109;
t81 = pkin(2) * t126 + (-pkin(7) - qJ(3)) * t102;
t93 = t111 * pkin(2) + pkin(1);
t132 = t101 * t93 + t105 * t81;
t99 = sin(pkin(12));
t131 = t106 * t99;
t130 = t106 * pkin(7) + qJ(1);
t129 = t101 * t102;
t128 = t102 * t105;
t127 = t102 * t109;
t125 = t106 * t111;
t122 = t99 * t129;
t121 = t99 * t128;
t120 = -t101 * t81 + t105 * t93;
t119 = pkin(2) * t127 + t106 * qJ(3) + t130;
t118 = g(1) * t101 - g(2) * t105;
t117 = t111 * t100 + t109 * t104;
t107 = -pkin(8) - qJ(4);
t113 = t83 * t106;
t69 = -t101 * t113 - t105 * t117;
t80 = t117 * t106;
t70 = -t101 * t80 + t105 * t83;
t103 = cos(pkin(12));
t92 = t103 * pkin(4) + pkin(3);
t116 = pkin(4) * t122 + t69 * t107 + t70 * t92 + t120;
t78 = t83 * t102;
t79 = t117 * t102;
t115 = pkin(4) * t131 + t78 * t107 + t79 * t92 + t119;
t68 = t101 * t83 + t105 * t80;
t98 = pkin(12) + qJ(5);
t94 = sin(t98);
t95 = cos(t98);
t59 = t95 * t128 + t68 * t94;
t61 = -t95 * t129 + t70 * t94;
t71 = -t106 * t95 + t79 * t94;
t114 = g(1) * t61 + g(2) * t59 + g(3) * t71;
t67 = -t101 * t117 + t105 * t113;
t58 = g(1) * t69 + g(2) * t67 + g(3) * t78;
t112 = -pkin(4) * t121 + t67 * t107 + t68 * t92 + t132;
t110 = cos(qJ(6));
t108 = sin(qJ(6));
t77 = -g(3) * t106 - t118 * t102;
t72 = t106 * t94 + t79 * t95;
t62 = t94 * t129 + t70 * t95;
t60 = -t94 * t128 + t68 * t95;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t101, t118, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t101 * t126 + t105 * t111) - g(2) * (t101 * t111 + t105 * t126) - g(3) * t127, -g(1) * (-t101 * t125 - t105 * t109) - g(2) * (-t101 * t109 + t105 * t125) - g(3) * t102 * t111, t77, -g(1) * (t105 * pkin(1) + pkin(7) * t129) - g(2) * (t101 * pkin(1) - pkin(7) * t128) - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t79, -t58, t77, -g(1) * t120 - g(2) * t132 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t103 + t122) - g(2) * (t68 * t103 - t121) - g(3) * (t79 * t103 + t131) -g(1) * (t103 * t129 - t70 * t99) - g(2) * (-t103 * t128 - t68 * t99) - g(3) * (t106 * t103 - t79 * t99) t58, -g(1) * (t70 * pkin(3) - t69 * qJ(4) + t120) - g(2) * (t68 * pkin(3) - t67 * qJ(4) + t132) - g(3) * (t79 * pkin(3) - t78 * qJ(4) + t119) 0, 0, 0, 0, 0, 0, -g(1) * t62 - g(2) * t60 - g(3) * t72, t114, t58, -g(1) * t116 - g(2) * t112 - g(3) * t115, 0, 0, 0, 0, 0, 0, -g(1) * (-t69 * t108 + t62 * t110) - g(2) * (-t67 * t108 + t60 * t110) - g(3) * (-t78 * t108 + t72 * t110) -g(1) * (-t62 * t108 - t69 * t110) - g(2) * (-t60 * t108 - t67 * t110) - g(3) * (-t72 * t108 - t78 * t110) -t114, -g(1) * (t62 * pkin(5) + t61 * pkin(9) + t116) - g(2) * (pkin(5) * t60 + pkin(9) * t59 + t112) - g(3) * (pkin(5) * t72 + pkin(9) * t71 + t115);];
U_reg  = t1;

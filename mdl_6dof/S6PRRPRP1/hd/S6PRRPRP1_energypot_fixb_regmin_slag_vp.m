% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:47
% EndTime: 2021-01-16 02:38:47
% DurationCPUTime: 0.28s
% Computational Cost: add. (225->73), mult. (379->113), div. (0->0), fcn. (463->12), ass. (0->48)
t111 = sin(qJ(5));
t105 = sin(pkin(10));
t107 = cos(pkin(10));
t113 = sin(qJ(2));
t108 = cos(pkin(6));
t116 = cos(qJ(2));
t124 = t108 * t116;
t87 = t105 * t113 - t107 * t124;
t134 = t87 * t111;
t89 = t105 * t124 + t107 * t113;
t133 = t89 * t111;
t106 = sin(pkin(6));
t132 = t105 * t106;
t131 = t106 * t107;
t112 = sin(qJ(3));
t130 = t106 * t112;
t129 = t106 * t113;
t115 = cos(qJ(3));
t128 = t106 * t115;
t127 = t106 * t116;
t126 = t108 * t112;
t125 = t108 * t113;
t123 = t105 * t130;
t98 = pkin(3) * t115 + pkin(2);
t122 = pkin(3) * t126 + t108 * pkin(7) + t98 * t129 + qJ(1);
t110 = qJ(4) + pkin(8);
t88 = t105 * t116 + t107 * t125;
t121 = t105 * pkin(1) + t87 * t110 + t88 * t98;
t86 = t105 * t125 - t107 * t116;
t120 = t107 * pkin(1) + pkin(3) * t123 + pkin(7) * t132 + t89 * t110 - t86 * t98;
t119 = (-pkin(3) * t112 - pkin(7)) * t107;
t104 = qJ(3) + pkin(11);
t100 = cos(t104);
t99 = sin(t104);
t77 = t100 * t132 + t86 * t99;
t78 = t100 * t131 + t88 * t99;
t83 = -t108 * t100 + t129 * t99;
t118 = g(1) * t77 - g(2) * t78 - g(3) * t83;
t117 = -g(1) * t89 - g(2) * t87 + g(3) * t127;
t114 = cos(qJ(5));
t109 = -qJ(6) - pkin(9);
t97 = pkin(5) * t114 + pkin(4);
t84 = t100 * t129 + t108 * t99;
t80 = -t100 * t86 + t132 * t99;
t79 = t100 * t88 - t131 * t99;
t76 = -g(1) * (t114 * t80 + t133) - g(2) * (t114 * t79 + t134) - g(3) * (-t111 * t127 + t114 * t84);
t75 = -g(1) * (-t111 * t80 + t114 * t89) - g(2) * (-t111 * t79 + t114 * t87) - g(3) * (-t111 * t84 - t114 * t127);
t1 = [-g(3) * qJ(1), 0, g(1) * t86 - g(2) * t88 - g(3) * t129, -t117, 0, 0, 0, 0, 0, -g(1) * (-t115 * t86 + t123) - g(2) * (-t107 * t130 + t115 * t88) - g(3) * (t113 * t128 + t126), -g(1) * (t105 * t128 + t112 * t86) - g(2) * (-t107 * t128 - t112 * t88) - g(3) * (t108 * t115 - t112 * t129), -g(1) * t80 - g(2) * t79 - g(3) * t84, -t118, t117, -g(1) * t120 - g(2) * (t106 * t119 + t121) - g(3) * (-t110 * t127 + t122), 0, 0, 0, 0, 0, t76, t75, t76, t75, t118, -g(1) * (pkin(5) * t133 + t109 * t77 + t80 * t97 + t120) - g(2) * (pkin(5) * t134 - t78 * t109 + t79 * t97 + t121) - g(3) * (-t83 * t109 + t84 * t97 + t122) + (-g(3) * (-pkin(5) * t111 - t110) * t116 - g(2) * t119) * t106;];
U_reg = t1;

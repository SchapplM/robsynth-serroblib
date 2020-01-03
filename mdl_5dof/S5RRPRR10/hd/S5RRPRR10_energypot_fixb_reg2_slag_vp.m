% Calculate inertial parameters regressor of potential energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:38
% EndTime: 2019-12-31 20:26:38
% DurationCPUTime: 0.22s
% Computational Cost: add. (233->76), mult. (538->123), div. (0->0), fcn. (668->12), ass. (0->50)
t96 = cos(pkin(5));
t124 = t96 * pkin(7) + pkin(6);
t94 = sin(pkin(5));
t99 = sin(qJ(2));
t123 = t94 * t99;
t100 = sin(qJ(1));
t104 = cos(qJ(1));
t81 = pkin(2) * t96 * t99 + (-pkin(7) - qJ(3)) * t94;
t103 = cos(qJ(2));
t90 = pkin(2) * t103 + pkin(1);
t122 = t100 * t90 + t104 * t81;
t121 = t100 * t94;
t120 = t100 * t99;
t95 = cos(pkin(10));
t119 = t103 * t95;
t118 = t104 * t94;
t117 = t104 * t99;
t116 = t100 * t103;
t115 = t103 * t104;
t114 = pkin(2) * t123 + t96 * qJ(3) + t124;
t113 = -t100 * t81 + t104 * t90;
t112 = g(1) * t100 - g(2) * t104;
t93 = sin(pkin(10));
t83 = -t93 * t99 + t119;
t111 = t103 * t93 + t95 * t99;
t109 = t83 * t96;
t68 = -t100 * t111 + t104 * t109;
t80 = t111 * t96;
t69 = t100 * t83 + t104 * t80;
t110 = t69 * pkin(3) - pkin(8) * t68 + t122;
t70 = -t100 * t109 - t104 * t111;
t71 = -t100 * t80 + t104 * t83;
t108 = t71 * pkin(3) - pkin(8) * t70 + t113;
t78 = -t94 * t119 + t93 * t123;
t79 = t111 * t94;
t107 = t79 * pkin(3) + pkin(8) * t78 + t114;
t102 = cos(qJ(4));
t98 = sin(qJ(4));
t62 = t102 * t118 + t69 * t98;
t64 = -t102 * t121 + t71 * t98;
t72 = -t96 * t102 + t79 * t98;
t106 = g(1) * t64 + g(2) * t62 + g(3) * t72;
t105 = g(1) * t70 + g(2) * t68 - g(3) * t78;
t101 = cos(qJ(5));
t97 = sin(qJ(5));
t77 = -g(3) * t96 - t112 * t94;
t73 = t102 * t79 + t96 * t98;
t65 = t102 * t71 + t98 * t121;
t63 = t69 * t102 - t98 * t118;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t100, t112, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * (-t96 * t120 + t115) - g(2) * (t96 * t117 + t116) - g(3) * t123, -g(1) * (-t96 * t116 - t117) - g(2) * (t96 * t115 - t120) - g(3) * t94 * t103, t77, -g(1) * (pkin(1) * t104 + pkin(7) * t121) - g(2) * (t100 * pkin(1) - pkin(7) * t118) - g(3) * t124, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t79, -t105, t77, -g(1) * t113 - g(2) * t122 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t63 - g(3) * t73, t106, t105, -g(1) * t108 - g(2) * t110 - g(3) * t107, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t65 - t70 * t97) - g(2) * (t101 * t63 - t68 * t97) - g(3) * (t101 * t73 + t78 * t97), -g(1) * (-t101 * t70 - t65 * t97) - g(2) * (-t101 * t68 - t63 * t97) - g(3) * (t101 * t78 - t73 * t97), -t106, -g(1) * (pkin(4) * t65 + pkin(9) * t64 + t108) - g(2) * (pkin(4) * t63 + pkin(9) * t62 + t110) - g(3) * (pkin(4) * t73 + t72 * pkin(9) + t107);];
U_reg = t1;

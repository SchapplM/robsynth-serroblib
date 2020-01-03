% Calculate inertial parameters regressor of potential energy for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:30
% EndTime: 2019-12-31 22:18:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (202->66), mult. (470->97), div. (0->0), fcn. (575->10), ass. (0->50)
t94 = cos(pkin(5));
t124 = t94 * pkin(7) + pkin(6);
t93 = sin(pkin(5));
t97 = sin(qJ(2));
t123 = t93 * t97;
t98 = sin(qJ(1));
t122 = t93 * t98;
t121 = t98 * t97;
t102 = cos(qJ(1));
t120 = t102 * pkin(1) + pkin(7) * t122;
t100 = cos(qJ(3));
t119 = t100 * t93;
t101 = cos(qJ(2));
t118 = t101 * t93;
t117 = t102 * t93;
t116 = t102 * t97;
t115 = t98 * t101;
t114 = t102 * t101;
t113 = t98 * pkin(1) - pkin(7) * t117;
t112 = g(1) * t98 - g(2) * t102;
t83 = t94 * t115 + t116;
t84 = -t94 * t121 + t114;
t111 = t84 * pkin(2) + t83 * pkin(8) + t120;
t110 = pkin(2) * t123 - pkin(8) * t118 + t124;
t82 = t94 * t116 + t115;
t96 = sin(qJ(3));
t71 = t82 * t100 - t96 * t117;
t81 = -t94 * t114 + t121;
t95 = sin(qJ(4));
t99 = cos(qJ(4));
t62 = t71 * t95 - t81 * t99;
t73 = t84 * t100 + t96 * t122;
t64 = t73 * t95 - t83 * t99;
t80 = t97 * t119 + t94 * t96;
t68 = t99 * t118 + t80 * t95;
t109 = g(1) * t64 + g(2) * t62 + g(3) * t68;
t70 = t100 * t117 + t82 * t96;
t72 = -t98 * t119 + t84 * t96;
t79 = -t94 * t100 + t96 * t123;
t108 = g(1) * t72 + g(2) * t70 + g(3) * t79;
t107 = t82 * pkin(2) + t81 * pkin(8) + t113;
t106 = -g(1) * t83 - g(2) * t81 + g(3) * t118;
t105 = t73 * pkin(3) + t72 * pkin(9) + t111;
t104 = t80 * pkin(3) + t79 * pkin(9) + t110;
t103 = t71 * pkin(3) + t70 * pkin(9) + t107;
t69 = -t95 * t118 + t80 * t99;
t65 = t73 * t99 + t83 * t95;
t63 = t71 * t99 + t81 * t95;
t60 = -g(1) * t65 - g(2) * t63 - g(3) * t69;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t102 - g(2) * t98, t112, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82 - g(3) * t123, -t106, -g(3) * t94 - t112 * t93, -g(1) * t120 - g(2) * t113 - g(3) * t124, 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - g(3) * t80, t108, t106, -g(1) * t111 - g(2) * t107 - g(3) * t110, 0, 0, 0, 0, 0, 0, t60, t109, -t108, -g(1) * t105 - g(2) * t103 - g(3) * t104, 0, 0, 0, 0, 0, 0, t60, -t108, -t109, -g(1) * (t65 * pkin(4) + t64 * qJ(5) + t105) - g(2) * (t63 * pkin(4) + t62 * qJ(5) + t103) - g(3) * (t69 * pkin(4) + t68 * qJ(5) + t104);];
U_reg = t1;

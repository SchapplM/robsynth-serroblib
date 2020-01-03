% Calculate inertial parameters regressor of potential energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:50
% EndTime: 2019-12-31 21:29:50
% DurationCPUTime: 0.18s
% Computational Cost: add. (207->81), mult. (360->121), div. (0->0), fcn. (422->12), ass. (0->47)
t94 = cos(pkin(5));
t122 = t94 * pkin(7) + pkin(6);
t93 = sin(pkin(5));
t98 = sin(qJ(2));
t121 = t93 * t98;
t99 = sin(qJ(1));
t120 = t93 * t99;
t97 = sin(qJ(3));
t119 = t94 * t97;
t118 = t98 * t99;
t103 = cos(qJ(1));
t117 = t103 * pkin(1) + pkin(7) * t120;
t101 = cos(qJ(3));
t116 = t101 * t93;
t102 = cos(qJ(2));
t115 = t102 * t93;
t114 = t103 * t93;
t113 = t103 * t98;
t112 = t99 * t102;
t111 = t102 * t103;
t110 = t97 * t120;
t90 = t99 * pkin(1);
t109 = -pkin(7) * t114 + t90;
t86 = pkin(3) * t101 + pkin(2);
t95 = -qJ(4) - pkin(8);
t108 = pkin(3) * t119 + t95 * t115 + t86 * t121 + t122;
t76 = t94 * t112 + t113;
t77 = -t94 * t118 + t111;
t107 = pkin(3) * t110 - t76 * t95 + t77 * t86 + t117;
t106 = g(1) * t99 - g(2) * t103;
t75 = t94 * t113 + t112;
t92 = qJ(3) + pkin(10);
t87 = sin(t92);
t88 = cos(t92);
t64 = t88 * t114 + t75 * t87;
t66 = -t88 * t120 + t77 * t87;
t70 = t87 * t121 - t94 * t88;
t105 = g(1) * t66 + g(2) * t64 + g(3) * t70;
t74 = -t94 * t111 + t118;
t63 = -g(1) * t76 - g(2) * t74 + g(3) * t115;
t104 = t75 * t86 - t74 * t95 + t90 + (-pkin(3) * t97 - pkin(7)) * t114;
t100 = cos(qJ(5));
t96 = sin(qJ(5));
t71 = t88 * t121 + t87 * t94;
t67 = t87 * t120 + t77 * t88;
t65 = -t87 * t114 + t75 * t88;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t99, t106, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t77 - g(2) * t75 - g(3) * t121, -t63, -g(3) * t94 - t106 * t93, -g(1) * t117 - g(2) * t109 - g(3) * t122, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t77 + t110) - g(2) * (t75 * t101 - t97 * t114) - g(3) * (t98 * t116 + t119), -g(1) * (t99 * t116 - t77 * t97) - g(2) * (-t101 * t114 - t75 * t97) - g(3) * (t101 * t94 - t97 * t121), t63, -g(1) * (pkin(2) * t77 + pkin(8) * t76 + t117) - g(2) * (t75 * pkin(2) + t74 * pkin(8) + t109) - g(3) * ((pkin(2) * t98 - pkin(8) * t102) * t93 + t122), 0, 0, 0, 0, 0, 0, -g(1) * t67 - g(2) * t65 - g(3) * t71, t105, t63, -g(1) * t107 - g(2) * t104 - g(3) * t108, 0, 0, 0, 0, 0, 0, -g(1) * (t100 * t67 + t76 * t96) - g(2) * (t100 * t65 + t74 * t96) - g(3) * (t100 * t71 - t96 * t115), -g(1) * (t100 * t76 - t67 * t96) - g(2) * (t100 * t74 - t65 * t96) - g(3) * (-t100 * t115 - t71 * t96), -t105, -g(1) * (pkin(4) * t67 + pkin(9) * t66 + t107) - g(2) * (t65 * pkin(4) + t64 * pkin(9) + t104) - g(3) * (pkin(4) * t71 + pkin(9) * t70 + t108);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:57
% EndTime: 2019-12-31 22:35:57
% DurationCPUTime: 0.18s
% Computational Cost: add. (207->81), mult. (360->121), div. (0->0), fcn. (422->12), ass. (0->47)
t95 = cos(pkin(5));
t123 = t95 * pkin(7) + pkin(6);
t94 = sin(pkin(5));
t98 = sin(qJ(2));
t122 = t94 * t98;
t99 = sin(qJ(1));
t121 = t94 * t99;
t97 = sin(qJ(3));
t120 = t95 * t97;
t119 = t98 * t99;
t103 = cos(qJ(1));
t118 = t103 * pkin(1) + pkin(7) * t121;
t101 = cos(qJ(3));
t117 = t101 * t94;
t102 = cos(qJ(2));
t116 = t102 * t94;
t115 = t102 * t99;
t114 = t103 * t94;
t113 = t103 * t98;
t112 = t102 * t103;
t111 = t97 * t121;
t91 = t99 * pkin(1);
t110 = -pkin(7) * t114 + t91;
t104 = -pkin(9) - pkin(8);
t87 = pkin(3) * t101 + pkin(2);
t109 = pkin(3) * t120 + t104 * t116 + t87 * t122 + t123;
t77 = t95 * t115 + t113;
t78 = -t95 * t119 + t112;
t108 = pkin(3) * t111 - t77 * t104 + t78 * t87 + t118;
t107 = g(1) * t99 - g(2) * t103;
t76 = t95 * t113 + t115;
t93 = qJ(3) + qJ(4);
t88 = sin(t93);
t89 = cos(t93);
t65 = t89 * t114 + t76 * t88;
t67 = -t89 * t121 + t78 * t88;
t71 = t88 * t122 - t95 * t89;
t106 = g(1) * t67 + g(2) * t65 + g(3) * t71;
t75 = -t95 * t112 + t119;
t64 = -g(1) * t77 - g(2) * t75 + g(3) * t116;
t105 = t76 * t87 - t75 * t104 + t91 + (-pkin(3) * t97 - pkin(7)) * t114;
t100 = cos(qJ(5));
t96 = sin(qJ(5));
t72 = t89 * t122 + t88 * t95;
t68 = t88 * t121 + t78 * t89;
t66 = -t88 * t114 + t76 * t89;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t103 - g(2) * t99, t107, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t122, -t64, -g(3) * t95 - t107 * t94, -g(1) * t118 - g(2) * t110 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t78 + t111) - g(2) * (t101 * t76 - t97 * t114) - g(3) * (t98 * t117 + t120), -g(1) * (t99 * t117 - t78 * t97) - g(2) * (-t101 * t114 - t76 * t97) - g(3) * (t101 * t95 - t97 * t122), t64, -g(1) * (pkin(2) * t78 + pkin(8) * t77 + t118) - g(2) * (pkin(2) * t76 + pkin(8) * t75 + t110) - g(3) * ((pkin(2) * t98 - pkin(8) * t102) * t94 + t123), 0, 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t66 - g(3) * t72, t106, t64, -g(1) * t108 - g(2) * t105 - g(3) * t109, 0, 0, 0, 0, 0, 0, -g(1) * (t100 * t68 + t77 * t96) - g(2) * (t100 * t66 + t75 * t96) - g(3) * (t100 * t72 - t96 * t116), -g(1) * (t100 * t77 - t68 * t96) - g(2) * (t100 * t75 - t66 * t96) - g(3) * (-t100 * t116 - t72 * t96), -t106, -g(1) * (pkin(4) * t68 + pkin(10) * t67 + t108) - g(2) * (pkin(4) * t66 + pkin(10) * t65 + t105) - g(3) * (pkin(4) * t72 + pkin(10) * t71 + t109);];
U_reg = t1;

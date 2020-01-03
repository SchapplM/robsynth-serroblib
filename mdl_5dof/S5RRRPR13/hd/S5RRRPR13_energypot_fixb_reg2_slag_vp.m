% Calculate inertial parameters regressor of potential energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:32
% EndTime: 2019-12-31 21:46:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (176->70), mult. (402->97), div. (0->0), fcn. (481->10), ass. (0->47)
t122 = pkin(4) + pkin(8);
t96 = cos(qJ(2));
t97 = cos(qJ(1));
t110 = t97 * t96;
t92 = sin(qJ(2));
t93 = sin(qJ(1));
t113 = t93 * t92;
t89 = cos(pkin(5));
t75 = -t89 * t110 + t113;
t121 = t75 * pkin(8);
t111 = t97 * t92;
t112 = t93 * t96;
t77 = t89 * t112 + t111;
t120 = t77 * pkin(8);
t119 = t89 * pkin(7) + pkin(6);
t88 = sin(pkin(5));
t118 = t88 * t92;
t117 = t88 * t93;
t95 = cos(qJ(3));
t116 = t88 * t95;
t115 = t88 * t96;
t114 = t88 * t97;
t109 = t97 * pkin(1) + pkin(7) * t117;
t108 = pkin(8) * t115;
t107 = pkin(2) * t118 + t119;
t78 = -t89 * t113 + t110;
t106 = t78 * pkin(2) + t109;
t105 = t93 * pkin(1) - pkin(7) * t114;
t104 = g(1) * t93 - g(2) * t97;
t76 = t89 * t111 + t112;
t103 = t76 * pkin(2) + t105;
t91 = sin(qJ(3));
t73 = t91 * t118 - t89 * t95;
t74 = t92 * t116 + t89 * t91;
t102 = t74 * pkin(3) + t73 * qJ(4) + t107;
t68 = -t93 * t116 + t78 * t91;
t69 = t91 * t117 + t78 * t95;
t101 = t69 * pkin(3) + t68 * qJ(4) + t106;
t66 = t95 * t114 + t76 * t91;
t100 = g(1) * t68 + g(2) * t66 + g(3) * t73;
t67 = -t91 * t114 + t76 * t95;
t99 = g(1) * t69 + g(2) * t67 + g(3) * t74;
t63 = -g(1) * t77 - g(2) * t75 + g(3) * t115;
t98 = t67 * pkin(3) + t66 * qJ(4) + t103;
t94 = cos(qJ(5));
t90 = sin(qJ(5));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t97 - g(2) * t93, t104, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t118, -t63, -g(3) * t89 - t104 * t88, -g(1) * t109 - g(2) * t105 - g(3) * t119, 0, 0, 0, 0, 0, 0, -t99, t100, t63, -g(1) * (t106 + t120) - g(2) * (t103 + t121) - g(3) * (t107 - t108), 0, 0, 0, 0, 0, 0, t63, t99, -t100, -g(1) * (t101 + t120) - g(2) * (t98 + t121) - g(3) * (t102 - t108), 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t90 + t77 * t94) - g(2) * (t66 * t90 + t75 * t94) - g(3) * (-t94 * t115 + t73 * t90), -g(1) * (t68 * t94 - t77 * t90) - g(2) * (t66 * t94 - t75 * t90) - g(3) * (t90 * t115 + t73 * t94), -t99, -g(1) * (t69 * pkin(9) + t122 * t77 + t101) - g(2) * (t67 * pkin(9) + t122 * t75 + t98) - g(3) * (t74 * pkin(9) - t122 * t115 + t102);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RRRPR12
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:34
% EndTime: 2019-12-31 21:40:34
% DurationCPUTime: 0.24s
% Computational Cost: add. (199->82), mult. (424->120), div. (0->0), fcn. (511->12), ass. (0->48)
t92 = cos(pkin(5));
t121 = t92 * pkin(7) + pkin(6);
t90 = sin(pkin(5));
t95 = sin(qJ(2));
t120 = t90 * t95;
t96 = sin(qJ(1));
t119 = t90 * t96;
t97 = cos(qJ(3));
t118 = t90 * t97;
t98 = cos(qJ(2));
t117 = t90 * t98;
t99 = cos(qJ(1));
t116 = t90 * t99;
t115 = t95 * t96;
t114 = t95 * t99;
t113 = t96 * t98;
t112 = t98 * t99;
t111 = t99 * pkin(1) + pkin(7) * t119;
t110 = pkin(2) * t120 + t121;
t76 = -t92 * t115 + t112;
t109 = t76 * pkin(2) + t111;
t89 = sin(pkin(10));
t108 = pkin(4) * t89 + pkin(8);
t107 = t96 * pkin(1) - pkin(7) * t116;
t106 = g(1) * t96 - g(2) * t99;
t74 = t92 * t114 + t113;
t105 = t74 * pkin(2) + t107;
t75 = t92 * t113 + t114;
t104 = pkin(8) * t75 + t109;
t103 = -pkin(8) * t117 + t110;
t94 = sin(qJ(3));
t65 = t97 * t116 + t74 * t94;
t67 = -t96 * t118 + t76 * t94;
t71 = t94 * t120 - t92 * t97;
t102 = g(1) * t67 + g(2) * t65 + g(3) * t71;
t73 = -t92 * t112 + t115;
t101 = t73 * pkin(8) + t105;
t100 = -g(1) * t75 - g(2) * t73 + g(3) * t117;
t93 = -pkin(9) - qJ(4);
t91 = cos(pkin(10));
t88 = pkin(10) + qJ(5);
t84 = cos(t88);
t83 = sin(t88);
t82 = pkin(4) * t91 + pkin(3);
t72 = t95 * t118 + t92 * t94;
t68 = t94 * t119 + t76 * t97;
t66 = -t94 * t116 + t74 * t97;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t99 - g(2) * t96, t106, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t74 - g(3) * t120, -t100, -g(3) * t92 - t106 * t90, -g(1) * t111 - g(2) * t107 - g(3) * t121, 0, 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t66 - g(3) * t72, t102, t100, -g(1) * t104 - g(2) * t101 - g(3) * t103, 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t91 + t75 * t89) - g(2) * (t66 * t91 + t73 * t89) - g(3) * (-t89 * t117 + t72 * t91), -g(1) * (-t68 * t89 + t75 * t91) - g(2) * (-t66 * t89 + t73 * t91) - g(3) * (-t91 * t117 - t72 * t89), -t102, -g(1) * (pkin(3) * t68 + qJ(4) * t67 + t104) - g(2) * (t66 * pkin(3) + t65 * qJ(4) + t101) - g(3) * (pkin(3) * t72 + qJ(4) * t71 + t103), 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t84 + t75 * t83) - g(2) * (t66 * t84 + t73 * t83) - g(3) * (-t83 * t117 + t72 * t84), -g(1) * (-t68 * t83 + t75 * t84) - g(2) * (-t66 * t83 + t73 * t84) - g(3) * (-t84 * t117 - t72 * t83), -t102, -g(1) * (t108 * t75 - t67 * t93 + t68 * t82 + t109) - g(2) * (t108 * t73 - t65 * t93 + t66 * t82 + t105) - g(3) * (-t108 * t117 - t71 * t93 + t72 * t82 + t110);];
U_reg = t1;

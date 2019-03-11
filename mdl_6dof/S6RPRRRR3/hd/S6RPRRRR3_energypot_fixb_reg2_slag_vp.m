% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:28
% EndTime: 2019-03-09 07:02:29
% DurationCPUTime: 0.17s
% Computational Cost: add. (222->83), mult. (186->110), div. (0->0), fcn. (182->12), ass. (0->42)
t94 = -pkin(9) - pkin(8);
t87 = qJ(2) + pkin(6);
t112 = g(3) * t87;
t89 = sin(qJ(3));
t111 = g(3) * t89;
t88 = sin(qJ(4));
t110 = t88 * pkin(4);
t91 = cos(qJ(4));
t74 = t91 * pkin(4) + pkin(3);
t84 = qJ(1) + pkin(11);
t75 = sin(t84);
t109 = t75 * t88;
t92 = cos(qJ(3));
t108 = t75 * t92;
t76 = cos(t84);
t107 = t76 * t92;
t86 = qJ(4) + qJ(5);
t77 = sin(t86);
t106 = t77 * t92;
t78 = cos(t86);
t105 = t78 * t92;
t85 = -pkin(10) + t94;
t104 = t85 * t89;
t103 = t88 * t92;
t102 = t89 * t94;
t101 = t91 * t92;
t90 = sin(qJ(1));
t100 = t90 * pkin(1) + t75 * pkin(2);
t93 = cos(qJ(1));
t99 = t93 * pkin(1) + t76 * pkin(2) + t75 * pkin(7);
t98 = -t76 * pkin(7) + t100;
t97 = pkin(3) * t92 + pkin(8) * t89;
t96 = g(1) * t76 + g(2) * t75;
t95 = -g(1) * t93 - g(2) * t90;
t79 = qJ(6) + t86;
t73 = cos(t79);
t72 = sin(t79);
t68 = pkin(5) * t77 + t110;
t67 = pkin(5) * t78 + t74;
t66 = g(1) * t75 - g(2) * t76;
t65 = -g(3) * t92 + t96 * t89;
t1 = [0, 0, 0, 0, 0, 0, t95, g(1) * t90 - g(2) * t93, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t96, t66, -g(3), t95 * pkin(1) - t112, 0, 0, 0, 0, 0, 0, -t96 * t92 - t111, t65, -t66, -g(1) * t99 - g(2) * t98 - t112, 0, 0, 0, 0, 0, 0, -g(1) * (t76 * t101 + t109) - g(2) * (t75 * t101 - t76 * t88) - t91 * t111, -g(1) * (-t76 * t103 + t75 * t91) - g(2) * (-t75 * t103 - t76 * t91) + t88 * t111, -t65, -g(1) * (t97 * t76 + t99) - g(2) * (t97 * t75 + t98) - g(3) * (t89 * pkin(3) - t92 * pkin(8) + t87) 0, 0, 0, 0, 0, 0, -g(1) * (t76 * t105 + t75 * t77) - g(2) * (t75 * t105 - t76 * t77) - t78 * t111, -g(1) * (-t76 * t106 + t75 * t78) - g(2) * (-t75 * t106 - t76 * t78) + t77 * t111, -t65, -g(1) * (pkin(4) * t109 + t99) - g(2) * (-t75 * t102 + t74 * t108 + t100) - g(3) * (t89 * t74 + t92 * t94 + t87) + (-g(1) * (t74 * t92 - t102) - g(2) * (-pkin(7) - t110)) * t76, 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t107 + t75 * t72) - g(2) * (t73 * t108 - t76 * t72) - t73 * t111, -g(1) * (-t72 * t107 + t75 * t73) - g(2) * (-t72 * t108 - t76 * t73) + t72 * t111, -t65, -g(1) * (t75 * t68 + t99) - g(2) * (-t75 * t104 + t67 * t108 + t100) - g(3) * (t89 * t67 + t92 * t85 + t87) + (-g(1) * (t67 * t92 - t104) - g(2) * (-pkin(7) - t68)) * t76;];
U_reg  = t1;

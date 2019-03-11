% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:25
% EndTime: 2019-03-09 03:42:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (222->83), mult. (186->108), div. (0->0), fcn. (182->12), ass. (0->40)
t88 = qJ(2) + pkin(6);
t109 = g(3) * t88;
t90 = sin(qJ(3));
t108 = g(3) * t90;
t86 = sin(pkin(11));
t107 = t86 * pkin(4);
t87 = cos(pkin(11));
t73 = t87 * pkin(4) + pkin(3);
t85 = qJ(1) + pkin(10);
t75 = sin(t85);
t106 = t75 * t86;
t92 = cos(qJ(3));
t105 = t75 * t92;
t77 = cos(t85);
t104 = t77 * t92;
t89 = -pkin(8) - qJ(4);
t83 = -pkin(9) + t89;
t103 = t83 * t90;
t102 = t86 * t92;
t101 = t87 * t92;
t100 = t89 * t90;
t91 = sin(qJ(1));
t99 = t91 * pkin(1) + t75 * pkin(2);
t84 = pkin(11) + qJ(5);
t93 = cos(qJ(1));
t98 = t93 * pkin(1) + t77 * pkin(2) + t75 * pkin(7);
t97 = -t77 * pkin(7) + t99;
t96 = g(1) * t77 + g(2) * t75;
t95 = -g(1) * t93 - g(2) * t91;
t94 = pkin(3) * t92 + qJ(4) * t90;
t78 = qJ(6) + t84;
t76 = cos(t84);
t74 = sin(t84);
t72 = cos(t78);
t71 = sin(t78);
t67 = pkin(5) * t74 + t107;
t66 = pkin(5) * t76 + t73;
t65 = g(1) * t75 - g(2) * t77;
t64 = -g(3) * t92 + t96 * t90;
t1 = [0, 0, 0, 0, 0, 0, t95, g(1) * t91 - g(2) * t93, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t96, t65, -g(3), t95 * pkin(1) - t109, 0, 0, 0, 0, 0, 0, -t96 * t92 - t108, t64, -t65, -g(1) * t98 - g(2) * t97 - t109, 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t101 + t106) - g(2) * (t75 * t101 - t77 * t86) - t87 * t108, -g(1) * (-t77 * t102 + t75 * t87) - g(2) * (-t75 * t102 - t77 * t87) + t86 * t108, -t64, -g(1) * (t94 * t77 + t98) - g(2) * (t94 * t75 + t97) - g(3) * (t90 * pkin(3) - t92 * qJ(4) + t88) 0, 0, 0, 0, 0, 0, -g(1) * (t76 * t104 + t75 * t74) - g(2) * (t76 * t105 - t77 * t74) - t76 * t108, -g(1) * (-t74 * t104 + t75 * t76) - g(2) * (-t74 * t105 - t77 * t76) + t74 * t108, -t64, -g(1) * (pkin(4) * t106 + t98) - g(2) * (-t75 * t100 + t73 * t105 + t99) - g(3) * (t90 * t73 + t92 * t89 + t88) + (-g(1) * (t73 * t92 - t100) - g(2) * (-pkin(7) - t107)) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t104 + t75 * t71) - g(2) * (t72 * t105 - t77 * t71) - t72 * t108, -g(1) * (-t71 * t104 + t75 * t72) - g(2) * (-t71 * t105 - t77 * t72) + t71 * t108, -t64, -g(1) * (t75 * t67 + t98) - g(2) * (-t75 * t103 + t66 * t105 + t99) - g(3) * (t90 * t66 + t92 * t83 + t88) + (-g(1) * (t66 * t92 - t103) - g(2) * (-pkin(7) - t67)) * t77;];
U_reg  = t1;

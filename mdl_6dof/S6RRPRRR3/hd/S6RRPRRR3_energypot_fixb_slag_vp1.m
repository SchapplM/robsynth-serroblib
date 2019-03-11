% Calculate potential energy for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:40
% EndTime: 2019-03-09 13:21:40
% DurationCPUTime: 0.52s
% Computational Cost: add. (236->104), mult. (212->128), div. (0->0), fcn. (196->12), ass. (0->43)
t106 = rSges(5,3) + pkin(8);
t81 = -pkin(9) - pkin(8);
t105 = rSges(6,3) - t81;
t104 = rSges(7,3) + pkin(10) - t81;
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t103 = g(1) * t80 + g(2) * t77;
t100 = rSges(3,3) + pkin(7);
t76 = sin(qJ(2));
t98 = t76 * pkin(2) + pkin(6);
t78 = cos(qJ(4));
t62 = t78 * pkin(4) + pkin(3);
t71 = qJ(2) + pkin(11);
t64 = sin(t71);
t97 = rSges(4,2) * t64;
t65 = cos(t71);
t96 = t77 * t65;
t73 = qJ(4) + qJ(5);
t66 = sin(t73);
t95 = t77 * t66;
t67 = cos(t73);
t94 = t77 * t67;
t75 = sin(qJ(4));
t93 = t77 * t75;
t92 = t77 * t78;
t91 = t80 * t65;
t90 = t80 * t66;
t89 = t80 * t67;
t88 = t80 * t75;
t87 = t80 * t78;
t79 = cos(qJ(2));
t63 = t79 * pkin(2) + pkin(1);
t74 = -qJ(3) - pkin(7);
t84 = t77 * t63 + t80 * t74;
t58 = t80 * t63;
t83 = -t77 * t74 + t58;
t82 = rSges(3,1) * t79 - rSges(3,2) * t76 + pkin(1);
t68 = qJ(6) + t73;
t61 = cos(t68);
t60 = sin(t68);
t56 = t75 * pkin(4) + pkin(5) * t66;
t55 = pkin(5) * t67 + t62;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t80 * rSges(2,1) - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + t80 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t76 * rSges(3,1) + t79 * rSges(3,2) + pkin(6)) + (g(1) * t82 - g(2) * t100) * t80 + (g(1) * t100 + g(2) * t82) * t77) - m(4) * (g(1) * (rSges(4,1) * t91 - t80 * t97 + t58) + g(2) * (-t80 * rSges(4,3) + t84) + g(3) * (t64 * rSges(4,1) + t65 * rSges(4,2) + t98) + (g(1) * (rSges(4,3) - t74) + g(2) * (rSges(4,1) * t65 - t97)) * t77) - m(5) * (g(1) * (pkin(3) * t91 + (t65 * t87 + t93) * rSges(5,1) + (-t65 * t88 + t92) * rSges(5,2) + t83) + g(2) * (pkin(3) * t96 + (t65 * t92 - t88) * rSges(5,1) + (-t65 * t93 - t87) * rSges(5,2) + t84) + g(3) * (-t106 * t65 + t98) + (g(3) * (rSges(5,1) * t78 - rSges(5,2) * t75 + pkin(3)) + t103 * t106) * t64) - m(6) * (g(1) * (t62 * t91 + pkin(4) * t93 + (t65 * t89 + t95) * rSges(6,1) + (-t65 * t90 + t94) * rSges(6,2) + t83) + g(2) * (t62 * t96 - pkin(4) * t88 + (t65 * t94 - t90) * rSges(6,1) + (-t65 * t95 - t89) * rSges(6,2) + t84) + g(3) * (-t105 * t65 + t98) + (g(3) * (rSges(6,1) * t67 - rSges(6,2) * t66 + t62) + t103 * t105) * t64) - m(7) * (g(1) * (t55 * t91 + t77 * t56 + (t77 * t60 + t61 * t91) * rSges(7,1) + (-t60 * t91 + t77 * t61) * rSges(7,2) + t83) + g(2) * (t55 * t96 - t80 * t56 + (-t80 * t60 + t61 * t96) * rSges(7,1) + (-t60 * t96 - t80 * t61) * rSges(7,2) + t84) + g(3) * (-t104 * t65 + t98) + (g(3) * (rSges(7,1) * t61 - rSges(7,2) * t60 + t55) + t103 * t104) * t64);
U  = t1;

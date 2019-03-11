% Calculate potential energy for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:15
% EndTime: 2019-03-09 05:15:15
% DurationCPUTime: 0.51s
% Computational Cost: add. (236->104), mult. (212->128), div. (0->0), fcn. (196->12), ass. (0->40)
t102 = rSges(5,3) + pkin(8);
t75 = -qJ(5) - pkin(8);
t101 = rSges(6,3) - t75;
t100 = rSges(7,3) + pkin(9) - t75;
t78 = sin(qJ(1));
t80 = cos(qJ(1));
t99 = g(1) * t80 + g(2) * t78;
t73 = sin(pkin(10));
t95 = t73 * pkin(2) + pkin(6);
t79 = cos(qJ(4));
t62 = t79 * pkin(4) + pkin(3);
t71 = pkin(10) + qJ(3);
t63 = sin(t71);
t94 = rSges(4,2) * t63;
t65 = cos(t71);
t93 = t78 * t65;
t72 = qJ(4) + pkin(11);
t66 = cos(t72);
t92 = t78 * t66;
t77 = sin(qJ(4));
t91 = t78 * t77;
t90 = t78 * t79;
t89 = t80 * t65;
t88 = t80 * t66;
t87 = t80 * t77;
t74 = cos(pkin(10));
t60 = pkin(2) * t74 + pkin(1);
t76 = -pkin(7) - qJ(2);
t84 = t78 * t60 + t80 * t76;
t83 = rSges(3,3) + qJ(2);
t57 = t80 * t60;
t82 = -t78 * t76 + t57;
t81 = rSges(3,1) * t74 - rSges(3,2) * t73 + pkin(1);
t67 = qJ(6) + t72;
t64 = sin(t72);
t59 = cos(t67);
t58 = sin(t67);
t55 = t77 * pkin(4) + pkin(5) * t64;
t54 = pkin(5) * t66 + t62;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t80 - t78 * rSges(2,2)) + g(2) * (t78 * rSges(2,1) + rSges(2,2) * t80) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t73 + rSges(3,2) * t74 + pkin(6)) + (g(1) * t81 - g(2) * t83) * t80 + (g(1) * t83 + g(2) * t81) * t78) - m(4) * (g(1) * (rSges(4,1) * t89 - t80 * t94 + t57) + g(2) * (-rSges(4,3) * t80 + t84) + g(3) * (rSges(4,1) * t63 + rSges(4,2) * t65 + t95) + (g(1) * (rSges(4,3) - t76) + g(2) * (rSges(4,1) * t65 - t94)) * t78) - m(5) * (g(1) * (pkin(3) * t89 + (t79 * t89 + t91) * rSges(5,1) + (-t65 * t87 + t90) * rSges(5,2) + t82) + g(2) * (pkin(3) * t93 + (t65 * t90 - t87) * rSges(5,1) + (-t65 * t91 - t79 * t80) * rSges(5,2) + t84) + g(3) * (-t102 * t65 + t95) + (g(3) * (rSges(5,1) * t79 - rSges(5,2) * t77 + pkin(3)) + t99 * t102) * t63) - m(6) * (g(1) * (t62 * t89 + pkin(4) * t91 + (t78 * t64 + t65 * t88) * rSges(6,1) + (-t64 * t89 + t92) * rSges(6,2) + t82) + g(2) * (t62 * t93 - pkin(4) * t87 + (-t80 * t64 + t65 * t92) * rSges(6,1) + (-t64 * t93 - t88) * rSges(6,2) + t84) + g(3) * (-t101 * t65 + t95) + (g(3) * (rSges(6,1) * t66 - rSges(6,2) * t64 + t62) + t99 * t101) * t63) - m(7) * (g(1) * (t54 * t89 + t78 * t55 + (t78 * t58 + t59 * t89) * rSges(7,1) + (-t58 * t89 + t78 * t59) * rSges(7,2) + t82) + g(2) * (t54 * t93 - t80 * t55 + (-t80 * t58 + t59 * t93) * rSges(7,1) + (-t58 * t93 - t80 * t59) * rSges(7,2) + t84) + g(3) * (-t100 * t65 + t95) + (g(3) * (rSges(7,1) * t59 - rSges(7,2) * t58 + t54) + t99 * t100) * t63);
U  = t1;

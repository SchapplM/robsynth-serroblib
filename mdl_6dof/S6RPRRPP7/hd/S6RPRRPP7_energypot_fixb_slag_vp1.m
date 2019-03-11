% Calculate potential energy for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:01
% EndTime: 2019-03-09 04:50:02
% DurationCPUTime: 0.36s
% Computational Cost: add. (143->90), mult. (236->109), div. (0->0), fcn. (230->6), ass. (0->32)
t90 = rSges(7,1) + pkin(5);
t78 = -rSges(7,3) - qJ(6);
t89 = pkin(2) + pkin(6);
t67 = sin(qJ(1));
t88 = g(1) * t67;
t70 = cos(qJ(1));
t87 = g(2) * t70;
t69 = cos(qJ(3));
t86 = rSges(4,2) * t69;
t66 = sin(qJ(3));
t85 = t66 * t67;
t84 = t66 * t70;
t65 = sin(qJ(4));
t83 = t67 * t65;
t68 = cos(qJ(4));
t82 = t67 * t68;
t81 = t70 * t68;
t80 = t70 * pkin(1) + t67 * qJ(2);
t61 = t67 * pkin(1);
t79 = t67 * pkin(7) + t61;
t77 = t70 * pkin(7) + t80;
t76 = t69 * pkin(3) + t66 * pkin(8) + t89;
t75 = pkin(3) * t85 + t77;
t74 = t76 + (pkin(4) * t68 + qJ(5) * t65) * t69;
t49 = t66 * t83 - t81;
t50 = t65 * t70 + t66 * t82;
t73 = t50 * pkin(4) + qJ(5) * t49 + t75;
t72 = -pkin(3) * t84 + t79 + (pkin(8) * t69 - qJ(2)) * t70;
t51 = t65 * t84 + t82;
t52 = -t66 * t81 + t83;
t71 = t52 * pkin(4) - t51 * qJ(5) + t72;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t70 - t67 * rSges(2,2)) + g(2) * (t67 * rSges(2,1) + rSges(2,2) * t70) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t70 + t67 * rSges(3,3) + t80) + g(2) * (-t67 * rSges(3,2) + t61 + (-rSges(3,3) - qJ(2)) * t70) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t85 + t67 * t86 + t77) + g(2) * (t67 * rSges(4,3) + t79) + g(3) * (rSges(4,1) * t69 - rSges(4,2) * t66 + t89) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t66 - qJ(2) - t86)) * t70) - m(5) * (g(1) * (rSges(5,1) * t50 - rSges(5,2) * t49 + t75) + g(2) * (t52 * rSges(5,1) + t51 * rSges(5,2) + t72) + g(3) * (rSges(5,3) * t66 + t76) + (rSges(5,3) * t87 + g(3) * (rSges(5,1) * t68 - rSges(5,2) * t65) + (-rSges(5,3) - pkin(8)) * t88) * t69) - m(6) * (g(1) * (rSges(6,1) * t50 + rSges(6,3) * t49 + t73) + g(2) * (t52 * rSges(6,1) - t51 * rSges(6,3) + t71) + g(3) * (rSges(6,2) * t66 + t74) + (rSges(6,2) * t87 + g(3) * (rSges(6,1) * t68 + rSges(6,3) * t65) + (-rSges(6,2) - pkin(8)) * t88) * t69) - m(7) * (g(1) * (rSges(7,2) * t49 + t90 * t50 + t73) + g(2) * (-t51 * rSges(7,2) + t90 * t52 + t71) + g(3) * (t78 * t66 + t74) + (g(3) * (rSges(7,2) * t65 + t90 * t68) + t78 * t87 + (-pkin(8) - t78) * t88) * t69);
U  = t1;

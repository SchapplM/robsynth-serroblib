% Calculate potential energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:18
% EndTime: 2019-03-09 15:28:19
% DurationCPUTime: 0.40s
% Computational Cost: add. (189->93), mult. (191->112), div. (0->0), fcn. (167->8), ass. (0->33)
t92 = rSges(7,3) + pkin(9);
t71 = -pkin(8) - pkin(7);
t91 = -qJ(5) - t71;
t90 = rSges(3,3) + pkin(7);
t66 = sin(qJ(2));
t88 = t66 * pkin(2) + pkin(6);
t64 = qJ(2) + qJ(3);
t61 = cos(t64);
t67 = sin(qJ(1));
t87 = t67 * t61;
t65 = sin(qJ(6));
t86 = t67 * t65;
t68 = cos(qJ(6));
t85 = t67 * t68;
t60 = sin(t64);
t70 = cos(qJ(1));
t84 = t70 * t60;
t83 = t70 * t61;
t82 = t70 * t65;
t81 = t70 * t68;
t69 = cos(qJ(2));
t58 = t69 * pkin(2) + pkin(1);
t80 = t67 * t58 + t70 * t71;
t79 = qJ(4) * t60;
t78 = t60 * pkin(3) + t88;
t51 = t70 * t58;
t77 = pkin(3) * t83 + t70 * t79 + t51;
t76 = t60 * pkin(4) + t78;
t75 = pkin(3) * t87 + t67 * t79 + t80;
t74 = pkin(4) * t83 + t77;
t73 = pkin(4) * t87 + t70 * qJ(5) + t75;
t72 = rSges(3,1) * t69 - rSges(3,2) * t66 + pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t70 * rSges(2,1) - t67 * rSges(2,2)) + g(2) * (t67 * rSges(2,1) + t70 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t66 * rSges(3,1) + t69 * rSges(3,2) + pkin(6)) + (g(1) * t72 - g(2) * t90) * t70 + (g(1) * t90 + g(2) * t72) * t67) - m(4) * (g(1) * (rSges(4,1) * t83 - rSges(4,2) * t84 + t51) + g(2) * (-t70 * rSges(4,3) + t80) + g(3) * (t60 * rSges(4,1) + t61 * rSges(4,2) + t88) + (g(1) * (rSges(4,3) - t71) + g(2) * (rSges(4,1) * t61 - rSges(4,2) * t60)) * t67) - m(5) * (g(1) * (rSges(5,1) * t83 + rSges(5,3) * t84 + t77) + g(2) * (-t70 * rSges(5,2) + t75) + g(3) * (t60 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t61 + t78) + (g(1) * (rSges(5,2) - t71) + g(2) * (rSges(5,1) * t61 + rSges(5,3) * t60)) * t67) - m(6) * (g(1) * (rSges(6,1) * t84 - rSges(6,2) * t83 + t74) + g(2) * (t70 * rSges(6,3) + t73) + g(3) * (-t60 * rSges(6,2) + (-rSges(6,1) - qJ(4)) * t61 + t76) + (g(1) * (-rSges(6,3) + t91) + g(2) * (rSges(6,1) * t60 - rSges(6,2) * t61)) * t67) - m(7) * (g(1) * (pkin(5) * t84 + (t60 * t81 - t86) * rSges(7,1) + (-t60 * t82 - t85) * rSges(7,2) + t74 + t91 * t67) + g(2) * (t67 * t60 * pkin(5) + (t60 * t85 + t82) * rSges(7,1) + (-t60 * t86 + t81) * rSges(7,2) + t73) + g(3) * (t92 * t60 + t76) + (g(3) * (-rSges(7,1) * t68 + rSges(7,2) * t65 - pkin(5) - qJ(4)) + (g(1) * t70 + g(2) * t67) * t92) * t61);
U  = t1;

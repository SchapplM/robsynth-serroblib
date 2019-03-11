% Calculate potential energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:01
% EndTime: 2019-03-09 02:44:02
% DurationCPUTime: 0.37s
% Computational Cost: add. (192->91), mult. (177->111), div. (0->0), fcn. (153->8), ass. (0->29)
t85 = rSges(7,3) + pkin(8);
t61 = qJ(1) + pkin(9);
t55 = sin(t61);
t63 = sin(qJ(3));
t83 = t55 * t63;
t66 = cos(qJ(3));
t82 = t55 * t66;
t56 = cos(t61);
t81 = t56 * t66;
t62 = sin(qJ(6));
t80 = t62 * t63;
t65 = cos(qJ(6));
t79 = t63 * t65;
t78 = pkin(6) + qJ(2);
t64 = sin(qJ(1));
t59 = t64 * pkin(1);
t77 = t55 * pkin(2) + t59;
t76 = qJ(4) * t63;
t75 = t63 * pkin(3) + t78;
t67 = cos(qJ(1));
t60 = t67 * pkin(1);
t74 = t56 * pkin(2) + t55 * pkin(7) + t60;
t73 = t63 * pkin(4) + t75;
t72 = pkin(3) * t82 + t55 * t76 + t77;
t71 = pkin(3) * t81 + t56 * t76 + t74;
t70 = rSges(6,1) * t63 - rSges(6,2) * t66;
t69 = pkin(4) * t82 + t56 * qJ(5) + t72;
t68 = pkin(4) * t81 + t71;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t67 * rSges(2,1) - t64 * rSges(2,2)) + g(2) * (t64 * rSges(2,1) + t67 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t56 * rSges(3,1) - t55 * rSges(3,2) + t60) + g(2) * (t55 * rSges(3,1) + t56 * rSges(3,2) + t59) + g(3) * (rSges(3,3) + t78)) - m(4) * (g(1) * (t55 * rSges(4,3) + t74) + g(2) * (rSges(4,1) * t82 - rSges(4,2) * t83 + t77) + g(3) * (t63 * rSges(4,1) + t66 * rSges(4,2) + t78) + (g(1) * (rSges(4,1) * t66 - rSges(4,2) * t63) + g(2) * (-rSges(4,3) - pkin(7))) * t56) - m(5) * (g(1) * (t55 * rSges(5,2) + t71) + g(2) * (rSges(5,1) * t82 + rSges(5,3) * t83 + t72) + g(3) * (t63 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t66 + t75) + (g(1) * (rSges(5,1) * t66 + rSges(5,3) * t63) + g(2) * (-rSges(5,2) - pkin(7))) * t56) - m(6) * (g(1) * t68 + g(2) * t69 + g(3) * (-t63 * rSges(6,2) + (-rSges(6,1) - qJ(4)) * t66 + t73) + (g(1) * t70 + g(2) * (rSges(6,3) - pkin(7))) * t56 + (g(1) * (-rSges(6,3) - qJ(5)) + g(2) * t70) * t55) - m(7) * (g(1) * (t56 * t63 * pkin(5) - t55 * qJ(5) + (-t55 * t62 + t56 * t79) * rSges(7,1) + (-t55 * t65 - t56 * t80) * rSges(7,2) + t68) + g(2) * (pkin(5) * t83 - t56 * pkin(7) + (t55 * t79 + t56 * t62) * rSges(7,1) + (-t55 * t80 + t56 * t65) * rSges(7,2) + t69) + g(3) * (t85 * t63 + t73) + (g(3) * (-rSges(7,1) * t65 + rSges(7,2) * t62 - pkin(5) - qJ(4)) + (g(1) * t56 + g(2) * t55) * t85) * t66);
U  = t1;

% Calculate potential energy for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:32
% EndTime: 2019-03-09 04:07:33
% DurationCPUTime: 0.46s
% Computational Cost: add. (174->102), mult. (200->124), div. (0->0), fcn. (184->10), ass. (0->32)
t63 = -pkin(8) - qJ(4);
t86 = rSges(6,3) - t63;
t85 = rSges(7,3) + pkin(9) - t63;
t65 = sin(qJ(1));
t67 = cos(qJ(1));
t84 = -g(1) * t65 + g(2) * t67;
t83 = rSges(5,3) + qJ(4);
t82 = pkin(2) + pkin(6);
t62 = cos(pkin(10));
t49 = t62 * pkin(4) + pkin(3);
t66 = cos(qJ(3));
t79 = rSges(4,2) * t66;
t61 = sin(pkin(10));
t78 = t61 * t67;
t77 = t65 * t61;
t64 = sin(qJ(3));
t76 = t65 * t64;
t75 = t67 * t64;
t72 = t67 * pkin(1) + t65 * qJ(2);
t56 = t65 * pkin(1);
t71 = t65 * pkin(7) + t56;
t60 = pkin(10) + qJ(5);
t69 = t67 * pkin(7) + t72;
t68 = -t67 * qJ(2) + t71;
t52 = qJ(6) + t60;
t51 = cos(t60);
t50 = sin(t60);
t48 = cos(t52);
t47 = sin(t52);
t46 = t61 * pkin(4) + pkin(5) * t50;
t45 = pkin(5) * t51 + t49;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t67 - t65 * rSges(2,2)) + g(2) * (t65 * rSges(2,1) + rSges(2,2) * t67) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t67 + t65 * rSges(3,3) + t72) + g(2) * (-t65 * rSges(3,2) + t56 + (-rSges(3,3) - qJ(2)) * t67) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t76 + t65 * t79 + t69) + g(2) * (t65 * rSges(4,3) + t71) + g(3) * (rSges(4,1) * t66 - rSges(4,2) * t64 + t82) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t64 - qJ(2) - t79)) * t67) - m(5) * (g(1) * (pkin(3) * t76 + (t62 * t76 + t78) * rSges(5,1) + (-t61 * t76 + t62 * t67) * rSges(5,2) + t69) + g(2) * (-pkin(3) * t75 + (-t62 * t75 + t77) * rSges(5,1) + (t61 * t75 + t65 * t62) * rSges(5,2) + t68) + g(3) * (t83 * t64 + t82) + (g(3) * (rSges(5,1) * t62 - rSges(5,2) * t61 + pkin(3)) + t84 * t83) * t66) - m(6) * (g(1) * (t49 * t76 + pkin(4) * t78 + (t50 * t67 + t51 * t76) * rSges(6,1) + (-t50 * t76 + t51 * t67) * rSges(6,2) + t69) + g(2) * (-t49 * t75 + pkin(4) * t77 + (t65 * t50 - t51 * t75) * rSges(6,1) + (t50 * t75 + t65 * t51) * rSges(6,2) + t68) + g(3) * (t64 * t86 + t82) + (g(3) * (rSges(6,1) * t51 - rSges(6,2) * t50 + t49) + t84 * t86) * t66) - m(7) * (g(1) * (t45 * t76 + t67 * t46 + (t47 * t67 + t48 * t76) * rSges(7,1) + (-t47 * t76 + t48 * t67) * rSges(7,2) + t69) + g(2) * (-t45 * t75 + t65 * t46 + (t65 * t47 - t48 * t75) * rSges(7,1) + (t47 * t75 + t65 * t48) * rSges(7,2) + t68) + g(3) * (t85 * t64 + t82) + (g(3) * (rSges(7,1) * t48 - rSges(7,2) * t47 + t45) + t84 * t85) * t66);
U  = t1;

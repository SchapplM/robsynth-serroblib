% Calculate potential energy for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:11:58
% EndTime: 2019-03-09 10:11:58
% DurationCPUTime: 0.39s
% Computational Cost: add. (217->88), mult. (173->105), div. (0->0), fcn. (149->10), ass. (0->36)
t95 = rSges(7,3) + pkin(9);
t94 = rSges(3,3) + pkin(7);
t72 = sin(qJ(2));
t92 = t72 * pkin(2) + pkin(6);
t75 = cos(qJ(2));
t62 = t75 * pkin(2) + pkin(1);
t69 = qJ(2) + pkin(10);
t65 = qJ(4) + t69;
t60 = sin(t65);
t76 = cos(qJ(1));
t91 = t60 * t76;
t71 = sin(qJ(6));
t73 = sin(qJ(1));
t90 = t73 * t71;
t74 = cos(qJ(6));
t89 = t73 * t74;
t61 = cos(t65);
t88 = t76 * t61;
t87 = t76 * t71;
t86 = t76 * t74;
t70 = -qJ(3) - pkin(7);
t85 = rSges(4,3) - t70;
t64 = cos(t69);
t52 = pkin(3) * t64 + t62;
t68 = -pkin(8) + t70;
t84 = t73 * t52 + t76 * t68;
t83 = qJ(5) * t60;
t63 = sin(t69);
t82 = pkin(3) * t63 + t92;
t51 = t76 * t52;
t81 = pkin(4) * t88 + t76 * t83 + t51;
t80 = t60 * pkin(4) + t82;
t79 = t84 + (pkin(4) * t61 + t83) * t73;
t78 = rSges(3,1) * t75 - rSges(3,2) * t72 + pkin(1);
t77 = rSges(4,1) * t64 - rSges(4,2) * t63 + t62;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t76 * rSges(2,1) - t73 * rSges(2,2)) + g(2) * (t73 * rSges(2,1) + t76 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t72 * rSges(3,1) + t75 * rSges(3,2) + pkin(6)) + (g(1) * t78 - g(2) * t94) * t76 + (g(1) * t94 + g(2) * t78) * t73) - m(4) * (g(3) * (t63 * rSges(4,1) + t64 * rSges(4,2) + t92) + (g(1) * t77 - g(2) * t85) * t76 + (g(1) * t85 + g(2) * t77) * t73) - m(5) * (g(1) * (rSges(5,1) * t88 - rSges(5,2) * t91 + t51) + g(2) * (-t76 * rSges(5,3) + t84) + g(3) * (t60 * rSges(5,1) + t61 * rSges(5,2) + t82) + (g(1) * (rSges(5,3) - t68) + g(2) * (rSges(5,1) * t61 - rSges(5,2) * t60)) * t73) - m(6) * (g(1) * (-rSges(6,2) * t88 + rSges(6,3) * t91 + t81) + g(2) * (-t76 * rSges(6,1) + t79) + g(3) * (-t60 * rSges(6,2) + (-rSges(6,3) - qJ(5)) * t61 + t80) + (g(1) * (rSges(6,1) - t68) + g(2) * (-rSges(6,2) * t61 + rSges(6,3) * t60)) * t73) - m(7) * (g(1) * ((t60 * t87 + t89) * rSges(7,1) + (t60 * t86 - t90) * rSges(7,2) + t81 + (pkin(5) - t68) * t73) + g(2) * (-t76 * pkin(5) + (t60 * t90 - t86) * rSges(7,1) + (t60 * t89 + t87) * rSges(7,2) + t79) + g(3) * (t95 * t60 + t80) + (g(3) * (-rSges(7,1) * t71 - rSges(7,2) * t74 - qJ(5)) + (g(1) * t76 + g(2) * t73) * t95) * t61);
U  = t1;

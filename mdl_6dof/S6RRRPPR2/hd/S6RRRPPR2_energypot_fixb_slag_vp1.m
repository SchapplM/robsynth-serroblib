% Calculate potential energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:51
% EndTime: 2019-03-09 15:24:51
% DurationCPUTime: 0.38s
% Computational Cost: add. (217->88), mult. (173->105), div. (0->0), fcn. (149->10), ass. (0->36)
t95 = rSges(7,3) + pkin(9);
t76 = -pkin(8) - pkin(7);
t94 = rSges(3,3) + pkin(7);
t71 = sin(qJ(2));
t92 = t71 * pkin(2) + pkin(6);
t74 = cos(qJ(2));
t62 = t74 * pkin(2) + pkin(1);
t69 = qJ(2) + qJ(3);
t63 = pkin(10) + t69;
t59 = sin(t63);
t75 = cos(qJ(1));
t91 = t59 * t75;
t70 = sin(qJ(6));
t72 = sin(qJ(1));
t90 = t72 * t70;
t73 = cos(qJ(6));
t89 = t72 * t73;
t60 = cos(t63);
t88 = t75 * t60;
t87 = t75 * t70;
t86 = t75 * t73;
t85 = rSges(4,3) - t76;
t65 = cos(t69);
t52 = pkin(3) * t65 + t62;
t68 = -qJ(4) + t76;
t84 = t72 * t52 + t75 * t68;
t83 = qJ(5) * t59;
t64 = sin(t69);
t82 = pkin(3) * t64 + t92;
t51 = t75 * t52;
t81 = pkin(4) * t88 + t75 * t83 + t51;
t80 = t59 * pkin(4) + t82;
t79 = t84 + (pkin(4) * t60 + t83) * t72;
t78 = rSges(3,1) * t74 - rSges(3,2) * t71 + pkin(1);
t77 = rSges(4,1) * t65 - rSges(4,2) * t64 + t62;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t75 * rSges(2,1) - t72 * rSges(2,2)) + g(2) * (t72 * rSges(2,1) + t75 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t71 * rSges(3,1) + t74 * rSges(3,2) + pkin(6)) + (g(1) * t78 - g(2) * t94) * t75 + (g(1) * t94 + g(2) * t78) * t72) - m(4) * (g(3) * (t64 * rSges(4,1) + t65 * rSges(4,2) + t92) + (g(1) * t77 - g(2) * t85) * t75 + (g(1) * t85 + g(2) * t77) * t72) - m(5) * (g(1) * (rSges(5,1) * t88 - rSges(5,2) * t91 + t51) + g(2) * (-t75 * rSges(5,3) + t84) + g(3) * (t59 * rSges(5,1) + t60 * rSges(5,2) + t82) + (g(1) * (rSges(5,3) - t68) + g(2) * (rSges(5,1) * t60 - rSges(5,2) * t59)) * t72) - m(6) * (g(1) * (-rSges(6,2) * t88 + rSges(6,3) * t91 + t81) + g(2) * (-t75 * rSges(6,1) + t79) + g(3) * (-t59 * rSges(6,2) + (-rSges(6,3) - qJ(5)) * t60 + t80) + (g(1) * (rSges(6,1) - t68) + g(2) * (-rSges(6,2) * t60 + rSges(6,3) * t59)) * t72) - m(7) * (g(1) * ((t59 * t87 + t89) * rSges(7,1) + (t59 * t86 - t90) * rSges(7,2) + t81 + (pkin(5) - t68) * t72) + g(2) * (-t75 * pkin(5) + (t59 * t90 - t86) * rSges(7,1) + (t59 * t89 + t87) * rSges(7,2) + t79) + g(3) * (t95 * t59 + t80) + (g(3) * (-rSges(7,1) * t70 - rSges(7,2) * t73 - qJ(5)) + (g(1) * t75 + g(2) * t72) * t95) * t60);
U  = t1;

% Calculate potential energy for
% S6RPRRPR2
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:08
% EndTime: 2019-03-09 05:00:09
% DurationCPUTime: 0.50s
% Computational Cost: add. (245->102), mult. (198->128), div. (0->0), fcn. (182->12), ass. (0->39)
t97 = rSges(5,3) + pkin(8);
t69 = -qJ(5) - pkin(8);
t96 = rSges(6,3) - t69;
t95 = rSges(7,3) + pkin(9) - t69;
t68 = qJ(1) + pkin(10);
t59 = sin(t68);
t61 = cos(t68);
t94 = g(1) * t61 + g(2) * t59;
t73 = cos(qJ(4));
t57 = t73 * pkin(4) + pkin(3);
t71 = sin(qJ(3));
t90 = rSges(4,2) * t71;
t70 = sin(qJ(4));
t89 = t59 * t70;
t74 = cos(qJ(3));
t88 = t59 * t74;
t87 = t61 * t70;
t86 = t61 * t74;
t85 = t70 * t74;
t84 = t73 * t74;
t67 = qJ(4) + pkin(11);
t60 = cos(t67);
t50 = pkin(5) * t60 + t57;
t83 = t74 * t50;
t82 = t74 * t57;
t79 = pkin(6) + qJ(2);
t72 = sin(qJ(1));
t63 = t72 * pkin(1);
t78 = t59 * pkin(2) + t63;
t75 = cos(qJ(1));
t65 = t75 * pkin(1);
t77 = t61 * pkin(2) + t59 * pkin(7) + t65;
t76 = -t61 * pkin(7) + t78;
t62 = qJ(6) + t67;
t58 = sin(t67);
t56 = cos(t62);
t55 = sin(t62);
t51 = t70 * pkin(4) + pkin(5) * t58;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t75 - t72 * rSges(2,2)) + g(2) * (t72 * rSges(2,1) + rSges(2,2) * t75) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t61 - rSges(3,2) * t59 + t65) + g(2) * (rSges(3,1) * t59 + rSges(3,2) * t61 + t63) + g(3) * (rSges(3,3) + t79)) - m(4) * (g(1) * (rSges(4,3) * t59 + t77) + g(2) * (rSges(4,1) * t88 - t59 * t90 + t78) + g(3) * (rSges(4,1) * t71 + rSges(4,2) * t74 + t79) + (g(1) * (rSges(4,1) * t74 - t90) + g(2) * (-rSges(4,3) - pkin(7))) * t61) - m(5) * (g(1) * (pkin(3) * t86 + (t61 * t84 + t89) * rSges(5,1) + (t59 * t73 - t61 * t85) * rSges(5,2) + t77) + g(2) * (pkin(3) * t88 + (t59 * t84 - t87) * rSges(5,1) + (-t59 * t85 - t61 * t73) * rSges(5,2) + t76) + g(3) * (-t97 * t74 + t79) + (g(3) * (rSges(5,1) * t73 - rSges(5,2) * t70 + pkin(3)) + t94 * t97) * t71) - m(6) * (g(1) * (t61 * t82 + pkin(4) * t89 + (t58 * t59 + t60 * t86) * rSges(6,1) + (-t58 * t86 + t59 * t60) * rSges(6,2) + t77) + g(2) * (t59 * t82 - pkin(4) * t87 + (-t58 * t61 + t60 * t88) * rSges(6,1) + (-t58 * t88 - t60 * t61) * rSges(6,2) + t76) + g(3) * (-t96 * t74 + t79) + (g(3) * (rSges(6,1) * t60 - rSges(6,2) * t58 + t57) + t94 * t96) * t71) - m(7) * (g(1) * (t61 * t83 + t59 * t51 + (t55 * t59 + t56 * t86) * rSges(7,1) + (-t55 * t86 + t56 * t59) * rSges(7,2) + t77) + g(2) * (t59 * t83 - t61 * t51 + (-t55 * t61 + t56 * t88) * rSges(7,1) + (-t55 * t88 - t61 * t56) * rSges(7,2) + t76) + g(3) * (-t95 * t74 + t79) + (g(3) * (rSges(7,1) * t56 - rSges(7,2) * t55 + t50) + t94 * t95) * t71);
U  = t1;

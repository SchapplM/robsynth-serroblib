% Calculate potential energy for
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:27
% EndTime: 2019-03-09 03:40:28
% DurationCPUTime: 0.48s
% Computational Cost: add. (245->102), mult. (198->128), div. (0->0), fcn. (182->12), ass. (0->39)
t71 = -pkin(8) - qJ(4);
t97 = rSges(6,3) - t71;
t96 = rSges(7,3) + pkin(9) - t71;
t68 = qJ(1) + pkin(10);
t59 = sin(t68);
t61 = cos(t68);
t95 = g(1) * t61 + g(2) * t59;
t94 = rSges(5,3) + qJ(4);
t70 = cos(pkin(11));
t57 = t70 * pkin(4) + pkin(3);
t72 = sin(qJ(3));
t91 = rSges(4,2) * t72;
t69 = sin(pkin(11));
t90 = t59 * t69;
t74 = cos(qJ(3));
t89 = t59 * t74;
t88 = t61 * t69;
t87 = t61 * t74;
t86 = t69 * t74;
t85 = t70 * t74;
t67 = pkin(11) + qJ(5);
t60 = cos(t67);
t50 = pkin(5) * t60 + t57;
t84 = t74 * t50;
t83 = t74 * t57;
t80 = pkin(6) + qJ(2);
t73 = sin(qJ(1));
t64 = t73 * pkin(1);
t79 = t59 * pkin(2) + t64;
t75 = cos(qJ(1));
t65 = t75 * pkin(1);
t77 = t61 * pkin(2) + t59 * pkin(7) + t65;
t76 = -t61 * pkin(7) + t79;
t62 = qJ(6) + t67;
t58 = sin(t67);
t56 = cos(t62);
t55 = sin(t62);
t51 = t69 * pkin(4) + pkin(5) * t58;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t75 * rSges(2,1) - t73 * rSges(2,2)) + g(2) * (t73 * rSges(2,1) + t75 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t61 * rSges(3,1) - t59 * rSges(3,2) + t65) + g(2) * (t59 * rSges(3,1) + t61 * rSges(3,2) + t64) + g(3) * (rSges(3,3) + t80)) - m(4) * (g(1) * (t59 * rSges(4,3) + t77) + g(2) * (rSges(4,1) * t89 - t59 * t91 + t79) + g(3) * (t72 * rSges(4,1) + t74 * rSges(4,2) + t80) + (g(1) * (rSges(4,1) * t74 - t91) + g(2) * (-rSges(4,3) - pkin(7))) * t61) - m(5) * (g(1) * (pkin(3) * t87 + (t61 * t85 + t90) * rSges(5,1) + (t59 * t70 - t61 * t86) * rSges(5,2) + t77) + g(2) * (pkin(3) * t89 + (t59 * t85 - t88) * rSges(5,1) + (-t59 * t86 - t61 * t70) * rSges(5,2) + t76) + g(3) * (-t94 * t74 + t80) + (g(3) * (rSges(5,1) * t70 - rSges(5,2) * t69 + pkin(3)) + t95 * t94) * t72) - m(6) * (g(1) * (t61 * t83 + pkin(4) * t90 + (t59 * t58 + t60 * t87) * rSges(6,1) + (-t58 * t87 + t59 * t60) * rSges(6,2) + t77) + g(2) * (t59 * t83 - pkin(4) * t88 + (-t61 * t58 + t60 * t89) * rSges(6,1) + (-t58 * t89 - t61 * t60) * rSges(6,2) + t76) + g(3) * (-t97 * t74 + t80) + (g(3) * (rSges(6,1) * t60 - rSges(6,2) * t58 + t57) + t95 * t97) * t72) - m(7) * (g(1) * (t61 * t84 + t59 * t51 + (t59 * t55 + t56 * t87) * rSges(7,1) + (-t55 * t87 + t59 * t56) * rSges(7,2) + t77) + g(2) * (t59 * t84 - t61 * t51 + (-t61 * t55 + t56 * t89) * rSges(7,1) + (-t55 * t89 - t61 * t56) * rSges(7,2) + t76) + g(3) * (-t96 * t74 + t80) + (g(3) * (rSges(7,1) * t56 - rSges(7,2) * t55 + t50) + t95 * t96) * t72);
U  = t1;

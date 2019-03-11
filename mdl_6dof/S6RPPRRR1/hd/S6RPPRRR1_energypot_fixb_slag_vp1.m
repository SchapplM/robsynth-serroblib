% Calculate potential energy for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:52
% EndTime: 2019-03-09 02:17:53
% DurationCPUTime: 0.33s
% Computational Cost: add. (223->85), mult. (148->98), div. (0->0), fcn. (124->12), ass. (0->39)
t93 = rSges(7,3) + pkin(9);
t70 = cos(pkin(11));
t56 = t70 * pkin(3) + pkin(2);
t67 = pkin(11) + qJ(4);
t61 = qJ(5) + t67;
t54 = sin(t61);
t91 = rSges(6,2) * t54;
t68 = qJ(1) + pkin(10);
t58 = sin(t68);
t72 = sin(qJ(6));
t90 = t58 * t72;
t74 = cos(qJ(6));
t89 = t58 * t74;
t55 = cos(t61);
t60 = cos(t68);
t88 = t60 * t55;
t87 = t60 * t72;
t86 = t60 * t74;
t71 = -pkin(7) - qJ(3);
t85 = rSges(5,3) - t71;
t84 = pkin(6) + qJ(2);
t59 = cos(t67);
t51 = pkin(4) * t59 + t56;
t75 = cos(qJ(1));
t65 = t75 * pkin(1);
t83 = t60 * t51 + t65;
t82 = rSges(4,3) + qJ(3);
t69 = sin(pkin(11));
t81 = t69 * pkin(3) + t84;
t73 = sin(qJ(1));
t64 = t73 * pkin(1);
t66 = -pkin(8) + t71;
t80 = t58 * t51 + t60 * t66 + t64;
t57 = sin(t67);
t79 = pkin(4) * t57 + t81;
t78 = g(1) * t65 + g(2) * t64;
t77 = rSges(4,1) * t70 - rSges(4,2) * t69 + pkin(2);
t76 = rSges(5,1) * t59 - rSges(5,2) * t57 + t56;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t75 * rSges(2,1) - t73 * rSges(2,2)) + g(2) * (t73 * rSges(2,1) + t75 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t60 * rSges(3,1) - t58 * rSges(3,2) + t65) + g(2) * (t58 * rSges(3,1) + t60 * rSges(3,2) + t64) + g(3) * (rSges(3,3) + t84)) - m(4) * (g(3) * (t69 * rSges(4,1) + t70 * rSges(4,2) + t84) + (g(1) * t77 - g(2) * t82) * t60 + (g(1) * t82 + g(2) * t77) * t58 + t78) - m(5) * (g(3) * (t57 * rSges(5,1) + t59 * rSges(5,2) + t81) + (g(1) * t76 - g(2) * t85) * t60 + (g(1) * t85 + g(2) * t76) * t58 + t78) - m(6) * (g(1) * (rSges(6,1) * t88 - t60 * t91 + t83) + g(2) * (-t60 * rSges(6,3) + t80) + g(3) * (t54 * rSges(6,1) + t55 * rSges(6,2) + t79) + (g(1) * (rSges(6,3) - t66) + g(2) * (rSges(6,1) * t55 - t91)) * t58) - m(7) * (g(1) * (pkin(5) * t88 - t58 * t66 + (t55 * t86 + t90) * rSges(7,1) + (-t55 * t87 + t89) * rSges(7,2) + t83) + g(2) * (t58 * t55 * pkin(5) + (t55 * t89 - t87) * rSges(7,1) + (-t55 * t90 - t86) * rSges(7,2) + t80) + g(3) * (-t93 * t55 + t79) + (g(3) * (rSges(7,1) * t74 - rSges(7,2) * t72 + pkin(5)) + (g(1) * t60 + g(2) * t58) * t93) * t54);
U  = t1;

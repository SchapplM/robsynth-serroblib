% Calculate potential energy for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:06
% EndTime: 2019-03-09 06:54:06
% DurationCPUTime: 0.34s
% Computational Cost: add. (223->85), mult. (148->98), div. (0->0), fcn. (124->12), ass. (0->39)
t93 = rSges(7,3) + pkin(10);
t75 = -pkin(8) - pkin(7);
t92 = rSges(4,3) + pkin(7);
t73 = cos(qJ(3));
t56 = t73 * pkin(3) + pkin(2);
t68 = qJ(3) + qJ(4);
t61 = qJ(5) + t68;
t54 = sin(t61);
t90 = rSges(6,2) * t54;
t66 = qJ(1) + pkin(11);
t57 = sin(t66);
t69 = sin(qJ(6));
t89 = t57 * t69;
t72 = cos(qJ(6));
t88 = t57 * t72;
t55 = cos(t61);
t58 = cos(t66);
t87 = t58 * t55;
t86 = t58 * t69;
t85 = t58 * t72;
t84 = rSges(5,3) - t75;
t83 = pkin(6) + qJ(2);
t60 = cos(t68);
t51 = pkin(4) * t60 + t56;
t74 = cos(qJ(1));
t65 = t74 * pkin(1);
t82 = t58 * t51 + t65;
t70 = sin(qJ(3));
t81 = t70 * pkin(3) + t83;
t71 = sin(qJ(1));
t63 = t71 * pkin(1);
t67 = -pkin(9) + t75;
t80 = t57 * t51 + t58 * t67 + t63;
t59 = sin(t68);
t79 = pkin(4) * t59 + t81;
t78 = g(1) * t65 + g(2) * t63;
t77 = rSges(4,1) * t73 - rSges(4,2) * t70 + pkin(2);
t76 = rSges(5,1) * t60 - rSges(5,2) * t59 + t56;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t74 * rSges(2,1) - t71 * rSges(2,2)) + g(2) * (t71 * rSges(2,1) + t74 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t58 * rSges(3,1) - t57 * rSges(3,2) + t65) + g(2) * (t57 * rSges(3,1) + t58 * rSges(3,2) + t63) + g(3) * (rSges(3,3) + t83)) - m(4) * (g(3) * (t70 * rSges(4,1) + t73 * rSges(4,2) + t83) + (g(1) * t77 - g(2) * t92) * t58 + (g(1) * t92 + g(2) * t77) * t57 + t78) - m(5) * (g(3) * (t59 * rSges(5,1) + t60 * rSges(5,2) + t81) + (g(1) * t76 - g(2) * t84) * t58 + (g(1) * t84 + g(2) * t76) * t57 + t78) - m(6) * (g(1) * (rSges(6,1) * t87 - t58 * t90 + t82) + g(2) * (-t58 * rSges(6,3) + t80) + g(3) * (t54 * rSges(6,1) + t55 * rSges(6,2) + t79) + (g(1) * (rSges(6,3) - t67) + g(2) * (rSges(6,1) * t55 - t90)) * t57) - m(7) * (g(1) * (pkin(5) * t87 - t57 * t67 + (t55 * t85 + t89) * rSges(7,1) + (-t55 * t86 + t88) * rSges(7,2) + t82) + g(2) * (t57 * t55 * pkin(5) + (t55 * t88 - t86) * rSges(7,1) + (-t55 * t89 - t85) * rSges(7,2) + t80) + g(3) * (-t93 * t55 + t79) + (g(3) * (rSges(7,1) * t72 - rSges(7,2) * t69 + pkin(5)) + (g(1) * t58 + g(2) * t57) * t93) * t54);
U  = t1;

% Calculate potential energy for
% S6RPRPRR1
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:48
% EndTime: 2019-03-09 03:33:48
% DurationCPUTime: 0.35s
% Computational Cost: add. (223->85), mult. (148->98), div. (0->0), fcn. (124->12), ass. (0->39)
t93 = rSges(7,3) + pkin(9);
t92 = rSges(4,3) + pkin(7);
t74 = cos(qJ(3));
t56 = t74 * pkin(3) + pkin(2);
t67 = qJ(3) + pkin(11);
t61 = qJ(5) + t67;
t54 = sin(t61);
t90 = rSges(6,2) * t54;
t68 = qJ(1) + pkin(10);
t58 = sin(t68);
t70 = sin(qJ(6));
t89 = t58 * t70;
t73 = cos(qJ(6));
t88 = t58 * t73;
t55 = cos(t61);
t60 = cos(t68);
t87 = t60 * t55;
t86 = t60 * t70;
t85 = t60 * t73;
t69 = -qJ(4) - pkin(7);
t84 = rSges(5,3) - t69;
t83 = pkin(6) + qJ(2);
t59 = cos(t67);
t51 = pkin(4) * t59 + t56;
t75 = cos(qJ(1));
t65 = t75 * pkin(1);
t82 = t60 * t51 + t65;
t71 = sin(qJ(3));
t81 = t71 * pkin(3) + t83;
t72 = sin(qJ(1));
t63 = t72 * pkin(1);
t66 = -pkin(8) + t69;
t80 = t58 * t51 + t60 * t66 + t63;
t57 = sin(t67);
t79 = pkin(4) * t57 + t81;
t78 = g(1) * t65 + g(2) * t63;
t77 = rSges(4,1) * t74 - rSges(4,2) * t71 + pkin(2);
t76 = rSges(5,1) * t59 - rSges(5,2) * t57 + t56;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t75 - t72 * rSges(2,2)) + g(2) * (t72 * rSges(2,1) + rSges(2,2) * t75) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t60 - rSges(3,2) * t58 + t65) + g(2) * (rSges(3,1) * t58 + rSges(3,2) * t60 + t63) + g(3) * (rSges(3,3) + t83)) - m(4) * (g(3) * (rSges(4,1) * t71 + rSges(4,2) * t74 + t83) + (g(1) * t77 - g(2) * t92) * t60 + (g(1) * t92 + g(2) * t77) * t58 + t78) - m(5) * (g(3) * (rSges(5,1) * t57 + rSges(5,2) * t59 + t81) + (g(1) * t76 - g(2) * t84) * t60 + (g(1) * t84 + g(2) * t76) * t58 + t78) - m(6) * (g(1) * (rSges(6,1) * t87 - t60 * t90 + t82) + g(2) * (-rSges(6,3) * t60 + t80) + g(3) * (rSges(6,1) * t54 + rSges(6,2) * t55 + t79) + (g(1) * (rSges(6,3) - t66) + g(2) * (rSges(6,1) * t55 - t90)) * t58) - m(7) * (g(1) * (pkin(5) * t87 - t58 * t66 + (t55 * t85 + t89) * rSges(7,1) + (-t55 * t86 + t88) * rSges(7,2) + t82) + g(2) * (t58 * t55 * pkin(5) + (t55 * t88 - t86) * rSges(7,1) + (-t55 * t89 - t85) * rSges(7,2) + t80) + g(3) * (-t93 * t55 + t79) + (g(3) * (rSges(7,1) * t73 - rSges(7,2) * t70 + pkin(5)) + (g(1) * t60 + g(2) * t58) * t93) * t54);
U  = t1;

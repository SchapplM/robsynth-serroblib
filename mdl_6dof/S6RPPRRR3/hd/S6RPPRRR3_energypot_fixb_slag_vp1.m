% Calculate potential energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:44
% EndTime: 2019-03-09 02:22:44
% DurationCPUTime: 0.42s
% Computational Cost: add. (193->89), mult. (161->110), div. (0->0), fcn. (141->10), ass. (0->35)
t92 = rSges(6,3) + pkin(8);
t91 = rSges(7,3) + pkin(9) + pkin(8);
t61 = qJ(1) + pkin(10);
t55 = sin(t61);
t56 = cos(t61);
t90 = -g(1) * t55 + g(2) * t56;
t67 = cos(qJ(4));
t86 = rSges(5,2) * t67;
t63 = sin(qJ(5));
t85 = t55 * t63;
t64 = sin(qJ(4));
t84 = t55 * t64;
t83 = t56 * t63;
t62 = qJ(5) + qJ(6);
t57 = sin(t62);
t82 = t57 * t64;
t58 = cos(t62);
t81 = t58 * t64;
t80 = t63 * t64;
t66 = cos(qJ(5));
t54 = t66 * pkin(5) + pkin(4);
t79 = t64 * t54;
t78 = t64 * t66;
t76 = pkin(6) + qJ(2);
t65 = sin(qJ(1));
t59 = t65 * pkin(1);
t75 = t55 * pkin(2) + t59;
t74 = pkin(3) + t76;
t68 = cos(qJ(1));
t60 = t68 * pkin(1);
t73 = t56 * pkin(2) + t55 * qJ(3) + t60;
t72 = t55 * pkin(7) + t75;
t71 = t56 * pkin(7) + t73;
t70 = -t56 * qJ(3) + t72;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t68 * rSges(2,1) - t65 * rSges(2,2)) + g(2) * (t65 * rSges(2,1) + t68 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t56 * rSges(3,1) - t55 * rSges(3,2) + t60) + g(2) * (t55 * rSges(3,1) + t56 * rSges(3,2) + t59) + g(3) * (rSges(3,3) + t76)) - m(4) * (g(1) * (-t56 * rSges(4,2) + t55 * rSges(4,3) + t73) + g(2) * (-t55 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t56 + t75) + g(3) * (rSges(4,1) + t76)) - m(5) * (g(1) * (rSges(5,1) * t84 + t55 * t86 + t71) + g(2) * (t55 * rSges(5,3) + t72) + g(3) * (t67 * rSges(5,1) - t64 * rSges(5,2) + t74) + (g(1) * rSges(5,3) + g(2) * (-rSges(5,1) * t64 - qJ(3) - t86)) * t56) - m(6) * (g(1) * (pkin(4) * t84 + (t55 * t78 + t83) * rSges(6,1) + (-t55 * t80 + t56 * t66) * rSges(6,2) + t71) + g(2) * (-t56 * t64 * pkin(4) + (-t56 * t78 + t85) * rSges(6,1) + (t55 * t66 + t56 * t80) * rSges(6,2) + t70) + g(3) * (t92 * t64 + t74) + (g(3) * (rSges(6,1) * t66 - rSges(6,2) * t63 + pkin(4)) + t90 * t92) * t67) - m(7) * (g(1) * (t55 * t79 + pkin(5) * t83 + (t55 * t81 + t56 * t57) * rSges(7,1) + (-t55 * t82 + t56 * t58) * rSges(7,2) + t71) + g(2) * (-t56 * t79 + pkin(5) * t85 + (t55 * t57 - t56 * t81) * rSges(7,1) + (t55 * t58 + t56 * t82) * rSges(7,2) + t70) + g(3) * (t91 * t64 + t74) + (g(3) * (rSges(7,1) * t58 - rSges(7,2) * t57 + t54) + t90 * t91) * t67);
U  = t1;

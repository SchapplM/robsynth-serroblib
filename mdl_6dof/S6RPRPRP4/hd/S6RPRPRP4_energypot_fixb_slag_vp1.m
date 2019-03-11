% Calculate potential energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:07
% EndTime: 2019-03-09 03:11:07
% DurationCPUTime: 0.34s
% Computational Cost: add. (209->88), mult. (200->107), div. (0->0), fcn. (184->8), ass. (0->34)
t95 = rSges(7,1) + pkin(5);
t94 = rSges(7,3) + qJ(6);
t71 = qJ(1) + pkin(9);
t65 = sin(t71);
t73 = sin(qJ(3));
t93 = t65 * t73;
t76 = cos(qJ(3));
t92 = t65 * t76;
t66 = cos(t71);
t91 = t66 * t76;
t72 = sin(qJ(5));
t90 = t72 * t73;
t75 = cos(qJ(5));
t89 = t73 * t75;
t88 = pkin(6) + qJ(2);
t74 = sin(qJ(1));
t69 = t74 * pkin(1);
t87 = t65 * pkin(2) + t69;
t86 = qJ(4) * t73;
t85 = t73 * pkin(3) + t88;
t77 = cos(qJ(1));
t70 = t77 * pkin(1);
t84 = t66 * pkin(2) + t65 * pkin(7) + t70;
t83 = t73 * pkin(8) + t85;
t82 = pkin(3) * t92 + t65 * t86 + t87;
t81 = pkin(3) * t91 + t66 * t86 + t84;
t80 = g(1) * t66 + g(2) * t65;
t79 = t65 * pkin(4) + pkin(8) * t91 + t81;
t78 = pkin(8) * t92 + t82 + (-pkin(4) - pkin(7)) * t66;
t53 = t65 * t90 - t66 * t75;
t52 = t65 * t89 + t66 * t72;
t51 = t65 * t75 + t66 * t90;
t50 = t65 * t72 - t66 * t89;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t77 - t74 * rSges(2,2)) + g(2) * (t74 * rSges(2,1) + rSges(2,2) * t77) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t66 - rSges(3,2) * t65 + t70) + g(2) * (rSges(3,1) * t65 + rSges(3,2) * t66 + t69) + g(3) * (rSges(3,3) + t88)) - m(4) * (g(1) * (rSges(4,3) * t65 + t84) + g(2) * (rSges(4,1) * t92 - rSges(4,2) * t93 + t87) + g(3) * (rSges(4,1) * t73 + rSges(4,2) * t76 + t88) + (g(1) * (rSges(4,1) * t76 - rSges(4,2) * t73) + g(2) * (-rSges(4,3) - pkin(7))) * t66) - m(5) * (g(1) * (rSges(5,1) * t65 + t81) + g(2) * (-rSges(5,2) * t92 + rSges(5,3) * t93 + t82) + g(3) * (-rSges(5,2) * t73 + (-rSges(5,3) - qJ(4)) * t76 + t85) + (g(1) * (-rSges(5,2) * t76 + rSges(5,3) * t73) + g(2) * (-rSges(5,1) - pkin(7))) * t66) - m(6) * (g(1) * (rSges(6,1) * t51 - rSges(6,2) * t50 + t79) + g(2) * (rSges(6,1) * t53 + rSges(6,2) * t52 + t78) + g(3) * (rSges(6,3) * t73 + t83) + (g(3) * (-rSges(6,1) * t72 - rSges(6,2) * t75 - qJ(4)) + t80 * rSges(6,3)) * t76) - m(7) * (g(1) * (t94 * t50 + t95 * t51 + t79) + g(2) * (-t94 * t52 + t95 * t53 + t78) + g(3) * (rSges(7,2) * t73 + t83) + (g(3) * (-t95 * t72 + t94 * t75 - qJ(4)) + t80 * rSges(7,2)) * t76);
U  = t1;

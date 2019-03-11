% Calculate potential energy for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:45
% EndTime: 2019-03-10 01:05:46
% DurationCPUTime: 0.44s
% Computational Cost: add. (226->99), mult. (212->120), div. (0->0), fcn. (196->10), ass. (0->40)
t100 = rSges(5,3) + pkin(9);
t78 = -pkin(10) - pkin(9);
t99 = rSges(6,3) - t78;
t98 = rSges(7,3) + qJ(6) - t78;
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t97 = g(1) * t77 + g(2) * t74;
t94 = rSges(3,3) + pkin(7);
t73 = sin(qJ(2));
t92 = t73 * pkin(2) + pkin(6);
t75 = cos(qJ(4));
t60 = t75 * pkin(4) + pkin(3);
t71 = qJ(2) + qJ(3);
t64 = sin(t71);
t91 = rSges(4,2) * t64;
t66 = cos(t71);
t90 = t66 * t74;
t89 = t66 * t77;
t72 = sin(qJ(4));
t88 = t74 * t72;
t87 = t74 * t75;
t86 = t77 * t72;
t85 = t77 * t75;
t76 = cos(qJ(2));
t61 = pkin(2) * t76 + pkin(1);
t79 = -pkin(8) - pkin(7);
t82 = t74 * t61 + t77 * t79;
t59 = t77 * t61;
t81 = -t74 * t79 + t59;
t80 = rSges(3,1) * t76 - rSges(3,2) * t73 + pkin(1);
t70 = qJ(4) + qJ(5);
t65 = cos(t70);
t63 = sin(t70);
t57 = pkin(4) * t72 + pkin(5) * t63;
t56 = pkin(5) * t65 + t60;
t55 = t63 * t74 + t65 * t89;
t54 = -t63 * t89 + t65 * t74;
t53 = -t63 * t77 + t65 * t90;
t52 = -t63 * t90 - t65 * t77;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t77 * rSges(2,1) - t74 * rSges(2,2)) + g(2) * (t74 * rSges(2,1) + t77 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t73 * rSges(3,1) + t76 * rSges(3,2) + pkin(6)) + (g(1) * t80 - g(2) * t94) * t77 + (g(1) * t94 + g(2) * t80) * t74) - m(4) * (g(1) * (rSges(4,1) * t89 - t77 * t91 + t59) + g(2) * (-t77 * rSges(4,3) + t82) + g(3) * (t64 * rSges(4,1) + t66 * rSges(4,2) + t92) + (g(1) * (rSges(4,3) - t79) + g(2) * (rSges(4,1) * t66 - t91)) * t74) - m(5) * (g(1) * (pkin(3) * t89 + (t66 * t85 + t88) * rSges(5,1) + (-t66 * t86 + t87) * rSges(5,2) + t81) + g(2) * (pkin(3) * t90 + (t66 * t87 - t86) * rSges(5,1) + (-t66 * t88 - t85) * rSges(5,2) + t82) + g(3) * (-t100 * t66 + t92) + (g(3) * (rSges(5,1) * t75 - rSges(5,2) * t72 + pkin(3)) + t97 * t100) * t64) - m(6) * (g(1) * (t55 * rSges(6,1) + t54 * rSges(6,2) + pkin(4) * t88 + t60 * t89 + t81) + g(2) * (t53 * rSges(6,1) + t52 * rSges(6,2) - pkin(4) * t86 + t60 * t90 + t82) + g(3) * (-t99 * t66 + t92) + (g(3) * (rSges(6,1) * t65 - rSges(6,2) * t63 + t60) + t97 * t99) * t64) - m(7) * (g(1) * (t55 * rSges(7,1) + t54 * rSges(7,2) + t56 * t89 + t74 * t57 + t81) + g(2) * (t53 * rSges(7,1) + t52 * rSges(7,2) + t56 * t90 - t77 * t57 + t82) + g(3) * (-t98 * t66 + t92) + (g(3) * (rSges(7,1) * t65 - rSges(7,2) * t63 + t56) + t97 * t98) * t64);
U  = t1;

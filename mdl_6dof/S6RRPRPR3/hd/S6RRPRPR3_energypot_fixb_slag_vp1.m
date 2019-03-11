% Calculate potential energy for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:44
% EndTime: 2019-03-09 10:15:45
% DurationCPUTime: 0.49s
% Computational Cost: add. (236->104), mult. (212->128), div. (0->0), fcn. (196->12), ass. (0->38)
t99 = rSges(5,3) + pkin(8);
t72 = -qJ(5) - pkin(8);
t98 = rSges(6,3) - t72;
t97 = rSges(7,3) + pkin(9) - t72;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t96 = g(1) * t79 + g(2) * t76;
t93 = rSges(3,3) + pkin(7);
t75 = sin(qJ(2));
t91 = t75 * pkin(2) + pkin(6);
t77 = cos(qJ(4));
t60 = t77 * pkin(4) + pkin(3);
t71 = qJ(2) + pkin(10);
t63 = sin(t71);
t90 = rSges(4,2) * t63;
t74 = sin(qJ(4));
t89 = t74 * t79;
t65 = cos(t71);
t88 = t76 * t65;
t87 = t76 * t74;
t86 = t76 * t77;
t85 = t79 * t65;
t78 = cos(qJ(2));
t61 = pkin(2) * t78 + pkin(1);
t73 = -qJ(3) - pkin(7);
t82 = t76 * t61 + t79 * t73;
t70 = qJ(4) + pkin(11);
t56 = t79 * t61;
t81 = -t76 * t73 + t56;
t80 = rSges(3,1) * t78 - rSges(3,2) * t75 + pkin(1);
t66 = qJ(6) + t70;
t64 = cos(t70);
t62 = sin(t70);
t58 = cos(t66);
t57 = sin(t66);
t54 = t74 * pkin(4) + pkin(5) * t62;
t53 = pkin(5) * t64 + t60;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - t76 * rSges(2,2)) + g(2) * (t76 * rSges(2,1) + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t75 + rSges(3,2) * t78 + pkin(6)) + (g(1) * t80 - g(2) * t93) * t79 + (g(1) * t93 + g(2) * t80) * t76) - m(4) * (g(1) * (rSges(4,1) * t85 - t79 * t90 + t56) + g(2) * (-rSges(4,3) * t79 + t82) + g(3) * (rSges(4,1) * t63 + rSges(4,2) * t65 + t91) + (g(1) * (rSges(4,3) - t73) + g(2) * (rSges(4,1) * t65 - t90)) * t76) - m(5) * (g(1) * (pkin(3) * t85 + (t77 * t85 + t87) * rSges(5,1) + (-t74 * t85 + t86) * rSges(5,2) + t81) + g(2) * (pkin(3) * t88 + (t65 * t86 - t89) * rSges(5,1) + (-t65 * t87 - t77 * t79) * rSges(5,2) + t82) + g(3) * (-t99 * t65 + t91) + (g(3) * (rSges(5,1) * t77 - rSges(5,2) * t74 + pkin(3)) + t96 * t99) * t63) - m(6) * (g(1) * (t60 * t85 + pkin(4) * t87 + (t76 * t62 + t64 * t85) * rSges(6,1) + (-t62 * t85 + t76 * t64) * rSges(6,2) + t81) + g(2) * (t60 * t88 - pkin(4) * t89 + (-t62 * t79 + t64 * t88) * rSges(6,1) + (-t62 * t88 - t64 * t79) * rSges(6,2) + t82) + g(3) * (-t98 * t65 + t91) + (g(3) * (rSges(6,1) * t64 - rSges(6,2) * t62 + t60) + t96 * t98) * t63) - m(7) * (g(1) * (t53 * t85 + t76 * t54 + (t76 * t57 + t58 * t85) * rSges(7,1) + (-t57 * t85 + t76 * t58) * rSges(7,2) + t81) + g(2) * (t53 * t88 - t79 * t54 + (-t57 * t79 + t58 * t88) * rSges(7,1) + (-t57 * t88 - t58 * t79) * rSges(7,2) + t82) + g(3) * (-t97 * t65 + t91) + (g(3) * (rSges(7,1) * t58 - rSges(7,2) * t57 + t53) + t96 * t97) * t63);
U  = t1;

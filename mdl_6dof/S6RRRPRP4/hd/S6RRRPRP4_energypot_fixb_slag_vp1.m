% Calculate potential energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:29
% EndTime: 2019-03-09 16:43:29
% DurationCPUTime: 0.38s
% Computational Cost: add. (201->90), mult. (214->109), div. (0->0), fcn. (198->8), ass. (0->38)
t102 = rSges(7,1) + pkin(5);
t101 = rSges(7,3) + qJ(6);
t100 = rSges(3,3) + pkin(7);
t76 = sin(qJ(2));
t99 = t76 * pkin(2) + pkin(6);
t74 = qJ(2) + qJ(3);
t70 = sin(t74);
t80 = cos(qJ(1));
t98 = t70 * t80;
t71 = cos(t74);
t77 = sin(qJ(1));
t97 = t71 * t77;
t96 = t71 * t80;
t75 = sin(qJ(5));
t95 = t75 * t77;
t94 = t75 * t80;
t78 = cos(qJ(5));
t93 = t77 * t78;
t92 = t78 * t80;
t79 = cos(qJ(2));
t68 = pkin(2) * t79 + pkin(1);
t81 = -pkin(8) - pkin(7);
t91 = t77 * t68 + t80 * t81;
t90 = qJ(4) * t70;
t89 = t70 * pkin(3) + t99;
t61 = t80 * t68;
t88 = pkin(3) * t96 + t80 * t90 + t61;
t87 = t70 * pkin(9) + t89;
t86 = pkin(3) * t97 + t77 * t90 + t91;
t85 = g(1) * t80 + g(2) * t77;
t84 = rSges(3,1) * t79 - rSges(3,2) * t76 + pkin(1);
t83 = -pkin(4) * t80 + pkin(9) * t97 + t86;
t82 = pkin(9) * t96 + t88 + (pkin(4) - t81) * t77;
t56 = t70 * t95 - t92;
t55 = t70 * t93 + t94;
t54 = t70 * t94 + t93;
t53 = -t70 * t92 + t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t80 - rSges(2,2) * t77) + g(2) * (rSges(2,1) * t77 + rSges(2,2) * t80) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t76 + rSges(3,2) * t79 + pkin(6)) + (g(1) * t84 - g(2) * t100) * t80 + (g(1) * t100 + g(2) * t84) * t77) - m(4) * (g(1) * (rSges(4,1) * t96 - rSges(4,2) * t98 + t61) + g(2) * (-rSges(4,3) * t80 + t91) + g(3) * (rSges(4,1) * t70 + rSges(4,2) * t71 + t99) + (g(1) * (rSges(4,3) - t81) + g(2) * (rSges(4,1) * t71 - rSges(4,2) * t70)) * t77) - m(5) * (g(1) * (-rSges(5,2) * t96 + rSges(5,3) * t98 + t88) + g(2) * (-rSges(5,1) * t80 + t86) + g(3) * (-rSges(5,2) * t70 + (-rSges(5,3) - qJ(4)) * t71 + t89) + (g(1) * (rSges(5,1) - t81) + g(2) * (-rSges(5,2) * t71 + rSges(5,3) * t70)) * t77) - m(6) * (g(1) * (t54 * rSges(6,1) - t53 * rSges(6,2) + t82) + g(2) * (rSges(6,1) * t56 + rSges(6,2) * t55 + t83) + g(3) * (rSges(6,3) * t70 + t87) + (g(3) * (-rSges(6,1) * t75 - rSges(6,2) * t78 - qJ(4)) + t85 * rSges(6,3)) * t71) - m(7) * (g(1) * (t101 * t53 + t102 * t54 + t82) + g(2) * (-t101 * t55 + t102 * t56 + t83) + g(3) * (rSges(7,2) * t70 + t87) + (g(3) * (t101 * t78 - t102 * t75 - qJ(4)) + t85 * rSges(7,2)) * t71);
U  = t1;

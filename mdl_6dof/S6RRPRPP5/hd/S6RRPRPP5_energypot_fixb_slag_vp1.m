% Calculate potential energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:24
% EndTime: 2019-03-09 10:03:25
% DurationCPUTime: 0.40s
% Computational Cost: add. (155->95), mult. (267->115), div. (0->0), fcn. (261->6), ass. (0->32)
t101 = rSges(7,1) + pkin(5);
t75 = sin(qJ(1));
t78 = cos(qJ(1));
t83 = g(1) * t78 + g(2) * t75;
t100 = -rSges(7,3) - qJ(6);
t74 = sin(qJ(2));
t97 = pkin(2) * t74 + pkin(6);
t73 = sin(qJ(4));
t96 = t73 * t78;
t95 = t74 * t75;
t76 = cos(qJ(4));
t94 = t75 * t76;
t77 = cos(qJ(2));
t93 = t75 * t77;
t92 = t76 * t78;
t91 = t77 * t78;
t90 = pkin(1) * t78 + pkin(7) * t75;
t89 = qJ(3) * t74;
t87 = pkin(8) * t74 + t97;
t71 = t75 * pkin(1);
t86 = pkin(2) * t93 + t75 * t89 + t71;
t85 = qJ(5) * t76 * t77 + t87;
t84 = pkin(2) * t91 + t78 * t89 + t90;
t82 = pkin(3) * t75 + pkin(8) * t91 + t84;
t81 = pkin(8) * t93 + t86 + (-pkin(3) - pkin(7)) * t78;
t55 = t73 * t75 - t74 * t92;
t56 = t74 * t96 + t94;
t80 = pkin(4) * t56 + t55 * qJ(5) + t82;
t57 = t74 * t94 + t96;
t58 = t73 * t95 - t92;
t79 = pkin(4) * t58 - t57 * qJ(5) + t81;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t78 - rSges(2,2) * t75) + g(2) * (rSges(2,1) * t75 + rSges(2,2) * t78) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t75 * rSges(3,3) + t90) + g(2) * (rSges(3,1) * t93 - rSges(3,2) * t95 + t71) + g(3) * (rSges(3,1) * t74 + rSges(3,2) * t77 + pkin(6)) + (g(1) * (rSges(3,1) * t77 - rSges(3,2) * t74) + g(2) * (-rSges(3,3) - pkin(7))) * t78) - m(4) * (g(1) * (t75 * rSges(4,1) + t84) + g(2) * (-rSges(4,2) * t93 + rSges(4,3) * t95 + t86) + g(3) * (-t74 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t77 + t97) + (g(1) * (-rSges(4,2) * t77 + rSges(4,3) * t74) + g(2) * (-rSges(4,1) - pkin(7))) * t78) - m(5) * (g(1) * (rSges(5,1) * t56 - rSges(5,2) * t55 + t82) + g(2) * (t58 * rSges(5,1) + t57 * rSges(5,2) + t81) + g(3) * (t74 * rSges(5,3) + t87) + (g(3) * (-rSges(5,1) * t73 - rSges(5,2) * t76 - qJ(3)) + t83 * rSges(5,3)) * t77) - m(6) * (g(1) * (t56 * rSges(6,1) + rSges(6,3) * t55 + t80) + g(2) * (t58 * rSges(6,1) - t57 * rSges(6,3) + t79) + g(3) * (t74 * rSges(6,2) + t85) + (g(3) * (rSges(6,3) * t76 - qJ(3) + (-rSges(6,1) - pkin(4)) * t73) + t83 * rSges(6,2)) * t77) - m(7) * (g(1) * (t55 * rSges(7,2) + t101 * t56 + t80) + g(2) * (-t57 * rSges(7,2) + t101 * t58 + t79) + t83 * t77 * t100 + (t85 + (rSges(7,2) * t76 - qJ(3) + (-pkin(4) - t101) * t73) * t77 + t100 * t74) * g(3));
U  = t1;

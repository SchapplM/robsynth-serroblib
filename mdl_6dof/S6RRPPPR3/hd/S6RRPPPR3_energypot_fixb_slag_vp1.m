% Calculate potential energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:46
% EndTime: 2019-03-09 08:13:46
% DurationCPUTime: 0.49s
% Computational Cost: add. (152->103), mult. (226->126), div. (0->0), fcn. (206->8), ass. (0->35)
t100 = rSges(7,3) + pkin(8) + qJ(5);
t71 = sin(qJ(1));
t73 = cos(qJ(1));
t99 = g(1) * t73 + g(2) * t71;
t98 = rSges(6,3) + qJ(5);
t70 = sin(qJ(2));
t95 = t70 * pkin(2) + pkin(6);
t67 = sin(pkin(9));
t94 = t67 * t73;
t93 = t70 * t71;
t92 = t70 * t73;
t66 = pkin(9) + qJ(6);
t58 = sin(t66);
t91 = t71 * t58;
t59 = cos(t66);
t90 = t71 * t59;
t89 = t71 * t67;
t68 = cos(pkin(9));
t88 = t71 * t68;
t72 = cos(qJ(2));
t87 = t71 * t72;
t86 = t72 * t73;
t84 = t73 * pkin(1) + t71 * pkin(7);
t83 = qJ(3) * t70;
t81 = t70 * pkin(3) + t95;
t64 = t71 * pkin(1);
t80 = pkin(2) * t87 + t71 * t83 + t64;
t79 = pkin(2) * t86 + t73 * t83 + t84;
t78 = pkin(3) * t87 + t73 * qJ(4) + t80;
t77 = pkin(3) * t86 + t79;
t76 = rSges(5,1) * t70 - rSges(5,2) * t72;
t75 = -t73 * pkin(7) + t78;
t74 = -t71 * qJ(4) + t77;
t57 = pkin(5) * t68 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t73 - t71 * rSges(2,2)) + g(2) * (t71 * rSges(2,1) + rSges(2,2) * t73) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t71 * rSges(3,3) + t84) + g(2) * (rSges(3,1) * t87 - rSges(3,2) * t93 + t64) + g(3) * (rSges(3,1) * t70 + rSges(3,2) * t72 + pkin(6)) + (g(1) * (rSges(3,1) * t72 - rSges(3,2) * t70) + g(2) * (-rSges(3,3) - pkin(7))) * t73) - m(4) * (g(1) * (t71 * rSges(4,2) + t79) + g(2) * (rSges(4,1) * t87 + rSges(4,3) * t93 + t80) + g(3) * (t70 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t72 + t95) + (g(1) * (rSges(4,1) * t72 + rSges(4,3) * t70) + g(2) * (-rSges(4,2) - pkin(7))) * t73) - m(5) * (g(1) * t77 + g(2) * t78 + g(3) * (-t70 * rSges(5,2) + (-rSges(5,1) - qJ(3)) * t72 + t81) + (g(1) * t76 + g(2) * (rSges(5,3) - pkin(7))) * t73 + (g(1) * (-rSges(5,3) - qJ(4)) + g(2) * t76) * t71) - m(6) * (g(1) * (pkin(4) * t92 + (t68 * t92 - t89) * rSges(6,1) + (-t67 * t92 - t88) * rSges(6,2) + t74) + g(2) * (pkin(4) * t93 + (t70 * t88 + t94) * rSges(6,1) + (t68 * t73 - t70 * t89) * rSges(6,2) + t75) + g(3) * (t98 * t70 + t81) + (g(3) * (-rSges(6,1) * t68 + rSges(6,2) * t67 - pkin(4) - qJ(3)) + t99 * t98) * t72) - m(7) * (g(1) * (t57 * t92 - pkin(5) * t89 + (t59 * t92 - t91) * rSges(7,1) + (-t58 * t92 - t90) * rSges(7,2) + t74) + g(2) * (t57 * t93 + pkin(5) * t94 + (t58 * t73 + t70 * t90) * rSges(7,1) + (t59 * t73 - t70 * t91) * rSges(7,2) + t75) + g(3) * (t100 * t70 + t81) + (g(3) * (-rSges(7,1) * t59 + rSges(7,2) * t58 - qJ(3) - t57) + t99 * t100) * t72);
U  = t1;

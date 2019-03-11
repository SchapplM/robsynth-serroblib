% Calculate potential energy for
% S6RPRPRR2
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:53
% EndTime: 2019-03-09 03:36:53
% DurationCPUTime: 0.43s
% Computational Cost: add. (234->93), mult. (172->111), div. (0->0), fcn. (152->12), ass. (0->42)
t105 = rSges(6,3) + pkin(8);
t104 = rSges(7,3) + pkin(9) + pkin(8);
t71 = qJ(1) + pkin(10);
t62 = sin(t71);
t64 = cos(t71);
t103 = g(1) * t64 + g(2) * t62;
t100 = rSges(4,3) + pkin(7);
t70 = qJ(3) + pkin(11);
t61 = sin(t70);
t98 = rSges(5,2) * t61;
t63 = cos(t70);
t97 = t62 * t63;
t72 = qJ(5) + qJ(6);
t65 = sin(t72);
t96 = t62 * t65;
t66 = cos(t72);
t95 = t62 * t66;
t74 = sin(qJ(5));
t94 = t62 * t74;
t77 = cos(qJ(5));
t93 = t62 * t77;
t92 = t64 * t63;
t91 = t64 * t65;
t90 = t64 * t66;
t89 = t64 * t74;
t88 = t64 * t77;
t86 = pkin(6) + qJ(2);
t78 = cos(qJ(3));
t60 = pkin(3) * t78 + pkin(2);
t79 = cos(qJ(1));
t69 = t79 * pkin(1);
t85 = t60 * t64 + t69;
t75 = sin(qJ(3));
t84 = pkin(3) * t75 + t86;
t76 = sin(qJ(1));
t68 = t76 * pkin(1);
t73 = -qJ(4) - pkin(7);
t83 = t60 * t62 + t64 * t73 + t68;
t82 = -t62 * t73 + t85;
t81 = rSges(4,1) * t78 - rSges(4,2) * t75 + pkin(2);
t59 = pkin(5) * t77 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - rSges(2,2) * t76) + g(2) * (rSges(2,1) * t76 + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t64 - rSges(3,2) * t62 + t69) + g(2) * (rSges(3,1) * t62 + rSges(3,2) * t64 + t68) + g(3) * (rSges(3,3) + t86)) - m(4) * (g(1) * t69 + g(2) * t68 + g(3) * (rSges(4,1) * t75 + rSges(4,2) * t78 + t86) + (g(1) * t81 - g(2) * t100) * t64 + (g(1) * t100 + g(2) * t81) * t62) - m(5) * (g(1) * (rSges(5,1) * t92 - t64 * t98 + t85) + g(2) * (-t64 * rSges(5,3) + t83) + g(3) * (rSges(5,1) * t61 + rSges(5,2) * t63 + t84) + (g(1) * (rSges(5,3) - t73) + g(2) * (rSges(5,1) * t63 - t98)) * t62) - m(6) * (g(1) * (pkin(4) * t92 + (t63 * t88 + t94) * rSges(6,1) + (-t63 * t89 + t93) * rSges(6,2) + t82) + g(2) * (pkin(4) * t97 + (t63 * t93 - t89) * rSges(6,1) + (-t63 * t94 - t88) * rSges(6,2) + t83) + g(3) * (-t105 * t63 + t84) + (g(3) * (rSges(6,1) * t77 - rSges(6,2) * t74 + pkin(4)) + t103 * t105) * t61) - m(7) * (g(1) * (t59 * t92 + pkin(5) * t94 + (t63 * t90 + t96) * rSges(7,1) + (-t63 * t91 + t95) * rSges(7,2) + t82) + g(2) * (t59 * t97 - pkin(5) * t89 + (t63 * t95 - t91) * rSges(7,1) + (-t63 * t96 - t90) * rSges(7,2) + t83) + g(3) * (-t104 * t63 + t84) + (g(3) * (rSges(7,1) * t66 - rSges(7,2) * t65 + t59) + t103 * t104) * t61);
U  = t1;

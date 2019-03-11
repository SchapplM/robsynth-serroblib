% Calculate potential energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:20
% EndTime: 2019-03-09 03:44:20
% DurationCPUTime: 0.47s
% Computational Cost: add. (211->96), mult. (190->120), div. (0->0), fcn. (170->10), ass. (0->35)
t103 = rSges(6,3) + pkin(8);
t102 = rSges(7,3) + pkin(9) + pkin(8);
t72 = qJ(1) + pkin(10);
t65 = sin(t72);
t66 = cos(t72);
t101 = g(1) * t66 + g(2) * t65;
t75 = sin(qJ(3));
t97 = t65 * t75;
t78 = cos(qJ(3));
t96 = t65 * t78;
t73 = qJ(5) + qJ(6);
t67 = sin(t73);
t95 = t67 * t75;
t68 = cos(t73);
t94 = t68 * t75;
t74 = sin(qJ(5));
t93 = t74 * t75;
t77 = cos(qJ(5));
t92 = t75 * t77;
t90 = pkin(6) + qJ(2);
t76 = sin(qJ(1));
t70 = t76 * pkin(1);
t89 = t65 * pkin(2) + t70;
t88 = qJ(4) * t75;
t87 = t65 * t93;
t86 = t66 * t93;
t85 = t75 * pkin(3) + t90;
t79 = cos(qJ(1));
t71 = t79 * pkin(1);
t84 = t66 * pkin(2) + t65 * pkin(7) + t71;
t83 = pkin(3) * t96 + t65 * t88 + t89;
t82 = t84 + (pkin(3) * t78 + t88) * t66;
t81 = -t66 * pkin(7) + t83;
t64 = pkin(5) * t77 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - rSges(2,2) * t76) + g(2) * (rSges(2,1) * t76 + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t66 - rSges(3,2) * t65 + t71) + g(2) * (rSges(3,1) * t65 + rSges(3,2) * t66 + t70) + g(3) * (rSges(3,3) + t90)) - m(4) * (g(1) * (rSges(4,3) * t65 + t84) + g(2) * (rSges(4,1) * t96 - rSges(4,2) * t97 + t89) + g(3) * (rSges(4,1) * t75 + rSges(4,2) * t78 + t90) + (g(1) * (rSges(4,1) * t78 - rSges(4,2) * t75) + g(2) * (-rSges(4,3) - pkin(7))) * t66) - m(5) * (g(1) * (rSges(5,1) * t65 + t82) + g(2) * (-rSges(5,2) * t96 + rSges(5,3) * t97 + t83) + g(3) * (-rSges(5,2) * t75 + (-rSges(5,3) - qJ(4)) * t78 + t85) + (g(1) * (-rSges(5,2) * t78 + rSges(5,3) * t75) + g(2) * (-rSges(5,1) - pkin(7))) * t66) - m(6) * (g(1) * (t65 * pkin(4) + (t65 * t77 + t86) * rSges(6,1) + (-t65 * t74 + t66 * t92) * rSges(6,2) + t82) + g(2) * (-t66 * pkin(4) + (-t66 * t77 + t87) * rSges(6,1) + (t65 * t92 + t66 * t74) * rSges(6,2) + t81) + g(3) * (t103 * t75 + t85) + (g(3) * (-rSges(6,1) * t74 - rSges(6,2) * t77 - qJ(4)) + t101 * t103) * t78) - m(7) * (g(1) * (t65 * t64 + pkin(5) * t86 + (t65 * t68 + t66 * t95) * rSges(7,1) + (-t65 * t67 + t66 * t94) * rSges(7,2) + t82) + g(2) * (-t66 * t64 + pkin(5) * t87 + (t65 * t95 - t66 * t68) * rSges(7,1) + (t65 * t94 + t66 * t67) * rSges(7,2) + t81) + g(3) * (t102 * t75 + t85) + (g(3) * (-rSges(7,1) * t67 - rSges(7,2) * t68 - pkin(5) * t74 - qJ(4)) + t101 * t102) * t78);
U  = t1;

% Calculate potential energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:48
% EndTime: 2019-03-09 04:53:49
% DurationCPUTime: 0.36s
% Computational Cost: add. (143->90), mult. (236->109), div. (0->0), fcn. (230->6), ass. (0->32)
t88 = rSges(7,1) + pkin(5);
t92 = rSges(7,3) + qJ(6);
t91 = pkin(2) + pkin(6);
t69 = sin(qJ(1));
t90 = g(1) * t69;
t72 = cos(qJ(1));
t89 = g(2) * t72;
t71 = cos(qJ(3));
t87 = rSges(4,2) * t71;
t68 = sin(qJ(3));
t86 = t68 * t69;
t85 = t68 * t72;
t67 = sin(qJ(4));
t84 = t69 * t67;
t70 = cos(qJ(4));
t83 = t69 * t70;
t82 = t72 * t70;
t81 = t72 * pkin(1) + t69 * qJ(2);
t63 = t69 * pkin(1);
t80 = t69 * pkin(7) + t63;
t79 = t72 * pkin(7) + t81;
t78 = t71 * pkin(3) + t68 * pkin(8) + t91;
t77 = pkin(3) * t86 + t79;
t76 = t78 + (pkin(4) * t70 + qJ(5) * t67) * t71;
t50 = t68 * t84 - t82;
t51 = t67 * t72 + t68 * t83;
t75 = t51 * pkin(4) + t50 * qJ(5) + t77;
t74 = -pkin(3) * t85 + t80 + (pkin(8) * t71 - qJ(2)) * t72;
t52 = t67 * t85 + t83;
t53 = t68 * t82 - t84;
t73 = -t53 * pkin(4) - t52 * qJ(5) + t74;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t72 - t69 * rSges(2,2)) + g(2) * (t69 * rSges(2,1) + rSges(2,2) * t72) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t72 + t69 * rSges(3,3) + t81) + g(2) * (-t69 * rSges(3,2) + t63 + (-rSges(3,3) - qJ(2)) * t72) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t86 + t69 * t87 + t79) + g(2) * (t69 * rSges(4,3) + t80) + g(3) * (rSges(4,1) * t71 - rSges(4,2) * t68 + t91) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t68 - qJ(2) - t87)) * t72) - m(5) * (g(1) * (rSges(5,1) * t51 - rSges(5,2) * t50 + t77) + g(2) * (-t53 * rSges(5,1) + t52 * rSges(5,2) + t74) + g(3) * (rSges(5,3) * t68 + t78) + (rSges(5,3) * t89 + g(3) * (rSges(5,1) * t70 - rSges(5,2) * t67) + (-rSges(5,3) - pkin(8)) * t90) * t71) - m(6) * (g(1) * (-rSges(6,2) * t51 + rSges(6,3) * t50 + t75) + g(2) * (t53 * rSges(6,2) - t52 * rSges(6,3) + t73) + g(3) * (rSges(6,1) * t68 + t76) + (rSges(6,1) * t89 + g(3) * (-rSges(6,2) * t70 + rSges(6,3) * t67) + (-rSges(6,1) - pkin(8)) * t90) * t71) - m(7) * (g(1) * (t50 * rSges(7,2) + t92 * t51 + t75) + g(2) * (-t52 * rSges(7,2) - t92 * t53 + t73) + g(3) * (t88 * t68 + t76) + (g(3) * (rSges(7,2) * t67 + t92 * t70) + t88 * t89 + (-pkin(8) - t88) * t90) * t71);
U  = t1;

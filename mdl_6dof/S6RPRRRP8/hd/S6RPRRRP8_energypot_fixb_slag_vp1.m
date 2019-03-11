% Calculate potential energy for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:46
% EndTime: 2019-03-09 06:22:46
% DurationCPUTime: 0.38s
% Computational Cost: add. (171->90), mult. (193->105), div. (0->0), fcn. (177->8), ass. (0->40)
t98 = rSges(7,1) + pkin(5);
t97 = rSges(7,3) + qJ(6);
t96 = pkin(2) + pkin(6);
t69 = sin(qJ(3));
t95 = pkin(3) * t69;
t67 = qJ(3) + qJ(4);
t61 = sin(t67);
t94 = pkin(4) * t61;
t70 = sin(qJ(1));
t93 = g(1) * t70;
t74 = -pkin(8) - pkin(7);
t92 = g(1) * t74;
t64 = t70 * pkin(1);
t91 = g(2) * t64;
t90 = rSges(4,3) + pkin(7);
t68 = sin(qJ(5));
t89 = t70 * t68;
t71 = cos(qJ(5));
t88 = t70 * t71;
t73 = cos(qJ(1));
t87 = t73 * t68;
t86 = t73 * t71;
t85 = rSges(5,3) - t74;
t84 = t73 * pkin(1) + t70 * qJ(2);
t72 = cos(qJ(3));
t83 = t72 * pkin(3) + t96;
t82 = t70 * t95 + t84;
t81 = -qJ(2) - t95;
t80 = t70 * t94 + t82;
t62 = cos(t67);
t79 = t62 * pkin(4) + t61 * pkin(9) + t83;
t78 = t73 * t62 * pkin(9) - t70 * t74 + t64;
t77 = rSges(4,1) * t69 + rSges(4,2) * t72;
t76 = rSges(5,1) * t61 + rSges(5,2) * t62;
t75 = t81 - t94;
t54 = -t61 * t86 + t89;
t53 = t61 * t87 + t88;
t52 = t61 * t88 + t87;
t51 = t61 * t89 - t86;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t73 * rSges(2,1) - t70 * rSges(2,2)) + g(2) * (t70 * rSges(2,1) + t73 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t73 * rSges(3,2) + t70 * rSges(3,3) + t84) + g(2) * (-t70 * rSges(3,2) + t64 + (-rSges(3,3) - qJ(2)) * t73) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * t84 + t91 + g(3) * (t72 * rSges(4,1) - t69 * rSges(4,2) + t96) + (g(1) * t77 + g(2) * t90) * t70 + (g(1) * t90 + g(2) * (-qJ(2) - t77)) * t73) - m(5) * (g(1) * t82 + t91 + g(3) * (t62 * rSges(5,1) - t61 * rSges(5,2) + t83) + (g(1) * t76 + g(2) * t85) * t70 + (g(1) * t85 + g(2) * (-t76 + t81)) * t73) - m(6) * (g(1) * (t52 * rSges(6,1) - t51 * rSges(6,2) + t80) + g(2) * (t54 * rSges(6,1) + t53 * rSges(6,2) + t78) + g(3) * (t61 * rSges(6,3) + t79) + (g(3) * (rSges(6,1) * t71 - rSges(6,2) * t68) + (-rSges(6,3) - pkin(9)) * t93) * t62 + (-t92 + g(2) * (rSges(6,3) * t62 + t75)) * t73) - m(7) * (g(1) * (t97 * t51 + t98 * t52 + t80) + g(2) * (-t97 * t53 + t98 * t54 + t78) + g(3) * (t61 * rSges(7,2) + t79) + (g(2) * t75 - t92) * t73 + (g(2) * rSges(7,2) * t73 + g(3) * (t97 * t68 + t98 * t71) + (-rSges(7,2) - pkin(9)) * t93) * t62);
U  = t1;

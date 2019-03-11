% Calculate potential energy for
% S6RRRRRP6
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:27
% EndTime: 2019-03-10 01:26:28
% DurationCPUTime: 0.52s
% Computational Cost: add. (247->106), mult. (255->130), div. (0->0), fcn. (247->10), ass. (0->38)
t106 = rSges(7,1) + pkin(5);
t105 = rSges(4,3) + pkin(8);
t85 = -pkin(9) - pkin(8);
t104 = rSges(5,3) - t85;
t103 = rSges(7,3) + qJ(6);
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t102 = g(1) * t84 + g(2) * t81;
t82 = cos(qJ(3));
t69 = t82 * pkin(3) + pkin(2);
t80 = sin(qJ(2));
t98 = rSges(3,2) * t80;
t79 = sin(qJ(3));
t97 = t81 * t79;
t83 = cos(qJ(2));
t96 = t81 * t83;
t95 = t84 * t79;
t94 = t84 * t83;
t90 = t84 * pkin(1) + t81 * pkin(7);
t78 = qJ(3) + qJ(4);
t71 = cos(t78);
t63 = pkin(4) * t71 + t69;
t77 = -pkin(10) + t85;
t89 = t80 * t63 + t83 * t77 + pkin(6);
t74 = t81 * pkin(1);
t88 = -t84 * pkin(7) + t74;
t70 = sin(t78);
t64 = t79 * pkin(3) + pkin(4) * t70;
t87 = t63 * t94 + t81 * t64 + t90;
t86 = t63 * t96 - t84 * t64 + t88;
t72 = qJ(5) + t78;
t68 = cos(t72);
t67 = sin(t72);
t58 = t81 * t67 + t68 * t94;
t57 = t67 * t94 - t81 * t68;
t56 = -t84 * t67 + t68 * t96;
t55 = t67 * t96 + t84 * t68;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t84 * rSges(2,1) - t81 * rSges(2,2)) + g(2) * (t81 * rSges(2,1) + t84 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t81 * rSges(3,3) + t90) + g(2) * (rSges(3,1) * t96 - t81 * t98 + t74) + g(3) * (t80 * rSges(3,1) + t83 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t83 - t98) + g(2) * (-rSges(3,3) - pkin(7))) * t84) - m(4) * (g(1) * (pkin(2) * t94 + (t82 * t94 + t97) * rSges(4,1) + (-t79 * t94 + t81 * t82) * rSges(4,2) + t90) + g(2) * (pkin(2) * t96 + (t82 * t96 - t95) * rSges(4,1) + (-t79 * t96 - t84 * t82) * rSges(4,2) + t88) + g(3) * (-t105 * t83 + pkin(6)) + (g(3) * (rSges(4,1) * t82 - rSges(4,2) * t79 + pkin(2)) + t102 * t105) * t80) - m(5) * (g(1) * (t69 * t94 + pkin(3) * t97 + (t81 * t70 + t71 * t94) * rSges(5,1) + (-t70 * t94 + t81 * t71) * rSges(5,2) + t90) + g(2) * (t69 * t96 - pkin(3) * t95 + (-t84 * t70 + t71 * t96) * rSges(5,1) + (-t70 * t96 - t84 * t71) * rSges(5,2) + t88) + g(3) * (-t104 * t83 + pkin(6)) + (g(3) * (rSges(5,1) * t71 - rSges(5,2) * t70 + t69) + t102 * t104) * t80) - m(6) * (g(1) * (t58 * rSges(6,1) - t57 * rSges(6,2) + t87) + g(2) * (t56 * rSges(6,1) - t55 * rSges(6,2) + t86) + g(3) * (-t83 * rSges(6,3) + t89) + (g(3) * (rSges(6,1) * t68 - rSges(6,2) * t67) + t102 * (rSges(6,3) - t77)) * t80) - m(7) * (g(1) * (t103 * t57 + t106 * t58 + t87) + g(2) * (t103 * t55 + t106 * t56 + t86) + g(3) * (-t83 * rSges(7,2) + t89) + (g(3) * (t103 * t67 + t106 * t68) + t102 * (rSges(7,2) - t77)) * t80);
U  = t1;

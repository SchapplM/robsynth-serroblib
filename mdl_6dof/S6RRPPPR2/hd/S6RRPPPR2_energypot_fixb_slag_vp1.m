% Calculate potential energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:50
% EndTime: 2019-03-09 08:09:50
% DurationCPUTime: 0.46s
% Computational Cost: add. (206->98), mult. (204->120), div. (0->0), fcn. (184->10), ass. (0->40)
t107 = rSges(7,3) + pkin(8) + qJ(5);
t106 = rSges(6,3) + qJ(5);
t78 = sin(qJ(1));
t80 = cos(qJ(1));
t105 = g(1) * t80 + g(2) * t78;
t102 = rSges(3,3) + pkin(7);
t77 = sin(qJ(2));
t101 = t77 * pkin(2) + pkin(6);
t72 = qJ(2) + pkin(9);
t67 = sin(t72);
t100 = t67 * t80;
t71 = pkin(10) + qJ(6);
t68 = cos(t71);
t99 = t68 * t80;
t69 = cos(t72);
t98 = t69 * t80;
t73 = sin(pkin(10));
t97 = t73 * t80;
t74 = cos(pkin(10));
t96 = t74 * t80;
t66 = sin(t71);
t95 = t78 * t66;
t94 = t78 * t68;
t93 = t78 * t73;
t92 = t78 * t74;
t79 = cos(qJ(2));
t65 = pkin(2) * t79 + pkin(1);
t75 = -qJ(3) - pkin(7);
t90 = t78 * t65 + t80 * t75;
t89 = qJ(4) * t67;
t87 = t67 * pkin(3) + t101;
t86 = t67 * t97;
t85 = t67 * t93;
t61 = t80 * t65;
t84 = pkin(3) * t98 + t80 * t89 + t61;
t83 = t90 + (pkin(3) * t69 + t89) * t78;
t82 = -t78 * t75 + t84;
t81 = rSges(3,1) * t79 - rSges(3,2) * t77 + pkin(1);
t63 = pkin(5) * t74 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t80 - t78 * rSges(2,2)) + g(2) * (t78 * rSges(2,1) + rSges(2,2) * t80) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t77 + rSges(3,2) * t79 + pkin(6)) + (g(1) * t81 - g(2) * t102) * t80 + (g(1) * t102 + g(2) * t81) * t78) - m(4) * (g(1) * (rSges(4,1) * t98 - rSges(4,2) * t100 + t61) + g(2) * (-t80 * rSges(4,3) + t90) + g(3) * (rSges(4,1) * t67 + rSges(4,2) * t69 + t101) + (g(1) * (rSges(4,3) - t75) + g(2) * (rSges(4,1) * t69 - rSges(4,2) * t67)) * t78) - m(5) * (g(1) * (-rSges(5,2) * t98 + rSges(5,3) * t100 + t84) + g(2) * (-rSges(5,1) * t80 + t83) + g(3) * (-rSges(5,2) * t67 + (-rSges(5,3) - qJ(4)) * t69 + t87) + (g(1) * (rSges(5,1) - t75) + g(2) * (-rSges(5,2) * t69 + rSges(5,3) * t67)) * t78) - m(6) * (g(1) * (t78 * pkin(4) + (t86 + t92) * rSges(6,1) + (t67 * t96 - t93) * rSges(6,2) + t82) + g(2) * (-t80 * pkin(4) + (t85 - t96) * rSges(6,1) + (t67 * t92 + t97) * rSges(6,2) + t83) + g(3) * (t106 * t67 + t87) + (g(3) * (-rSges(6,1) * t73 - rSges(6,2) * t74 - qJ(4)) + t105 * t106) * t69) - m(7) * (g(1) * (t78 * t63 + pkin(5) * t86 + (t66 * t100 + t94) * rSges(7,1) + (t67 * t99 - t95) * rSges(7,2) + t82) + g(2) * (-t80 * t63 + pkin(5) * t85 + (t67 * t95 - t99) * rSges(7,1) + (t66 * t80 + t67 * t94) * rSges(7,2) + t83) + g(3) * (t107 * t67 + t87) + (g(3) * (-rSges(7,1) * t66 - rSges(7,2) * t68 - pkin(5) * t73 - qJ(4)) + t105 * t107) * t69);
U  = t1;

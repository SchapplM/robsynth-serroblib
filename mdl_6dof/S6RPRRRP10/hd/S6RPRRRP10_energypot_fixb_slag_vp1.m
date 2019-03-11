% Calculate potential energy for
% S6RPRRRP10
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:12
% EndTime: 2019-03-09 06:30:13
% DurationCPUTime: 0.46s
% Computational Cost: add. (171->98), mult. (215->115), div. (0->0), fcn. (203->8), ass. (0->40)
t99 = rSges(7,1) + pkin(5);
t98 = rSges(5,3) + pkin(8);
t76 = -pkin(9) - pkin(8);
t97 = rSges(7,2) - t76;
t96 = rSges(7,3) + qJ(6);
t95 = pkin(2) + pkin(6);
t72 = sin(qJ(1));
t94 = g(1) * t72;
t75 = cos(qJ(1));
t93 = g(2) * t75;
t70 = sin(qJ(4));
t91 = t70 * t72;
t90 = t70 * t75;
t71 = sin(qJ(3));
t89 = t71 * t72;
t88 = t71 * t75;
t73 = cos(qJ(4));
t87 = t72 * t73;
t74 = cos(qJ(3));
t86 = t72 * t74;
t85 = t73 * t75;
t69 = qJ(4) + qJ(5);
t63 = cos(t69);
t84 = t75 * t63;
t83 = rSges(6,3) - t76;
t82 = t75 * pkin(1) + t72 * qJ(2);
t66 = t72 * pkin(1);
t81 = t72 * pkin(7) + t66;
t61 = pkin(4) * t73 + pkin(3);
t80 = t74 * t61 + t95;
t79 = t75 * pkin(7) + t82;
t78 = -t75 * qJ(2) + t81;
t77 = pkin(4) * t90 + t61 * t89 + t76 * t86 + t79;
t62 = sin(t69);
t59 = pkin(4) * t91;
t54 = t62 * t72 - t71 * t84;
t53 = t62 * t88 + t63 * t72;
t52 = t62 * t75 + t63 * t89;
t51 = t62 * t89 - t84;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t75 - rSges(2,2) * t72) + g(2) * (rSges(2,1) * t72 + rSges(2,2) * t75) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t75 + rSges(3,3) * t72 + t82) + g(2) * (-rSges(3,2) * t72 + t66 + (-rSges(3,3) - qJ(2)) * t75) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t89 + rSges(4,2) * t86 + t79) + g(2) * (rSges(4,3) * t72 + t81) + g(3) * (rSges(4,1) * t74 - rSges(4,2) * t71 + t95) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t71 - rSges(4,2) * t74 - qJ(2))) * t75) - m(5) * (g(1) * (pkin(3) * t89 + (t71 * t87 + t90) * rSges(5,1) + (-t70 * t89 + t85) * rSges(5,2) + t79) + g(2) * (-pkin(3) * t88 + (-t71 * t85 + t91) * rSges(5,1) + (t70 * t88 + t87) * rSges(5,2) + t78) + g(3) * (t98 * t71 + t95) + (g(3) * (rSges(5,1) * t73 - rSges(5,2) * t70 + pkin(3)) + (-t94 + t93) * t98) * t74) - m(6) * (g(1) * (rSges(6,1) * t52 - rSges(6,2) * t51 - rSges(6,3) * t86 + t77) + g(2) * (t54 * rSges(6,1) + t53 * rSges(6,2) + t59 + t81) + g(3) * ((rSges(6,1) * t63 - rSges(6,2) * t62) * t74 + t83 * t71 + t80) + (-t61 * t71 + t83 * t74 - qJ(2)) * t93) - m(7) * (g(1) * (t96 * t51 + t99 * t52 + t77) + g(2) * (-t96 * t53 + t99 * t54 - t61 * t88 + t59 + t78) + g(3) * (t97 * t71 + t80) + (-rSges(7,2) * t94 + g(3) * (t96 * t62 + t99 * t63) + t97 * t93) * t74);
U  = t1;

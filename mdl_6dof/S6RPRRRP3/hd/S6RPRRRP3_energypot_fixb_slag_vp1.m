% Calculate potential energy for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:02:49
% EndTime: 2019-03-09 06:02:49
% DurationCPUTime: 0.42s
% Computational Cost: add. (250->93), mult. (213->115), div. (0->0), fcn. (201->10), ass. (0->41)
t111 = rSges(7,1) + pkin(5);
t110 = rSges(5,3) + pkin(8);
t109 = rSges(7,3) + qJ(6);
t78 = qJ(1) + pkin(10);
t72 = sin(t78);
t73 = cos(t78);
t108 = g(1) * t73 + g(2) * t72;
t81 = sin(qJ(3));
t104 = rSges(4,2) * t81;
t80 = sin(qJ(4));
t103 = t72 * t80;
t84 = cos(qJ(3));
t102 = t72 * t84;
t101 = t73 * t80;
t100 = t73 * t84;
t79 = qJ(4) + qJ(5);
t74 = sin(t79);
t99 = t74 * t84;
t75 = cos(t79);
t98 = t75 * t84;
t97 = t80 * t84;
t83 = cos(qJ(4));
t96 = t83 * t84;
t93 = pkin(6) + qJ(2);
t82 = sin(qJ(1));
t76 = t82 * pkin(1);
t92 = t72 * pkin(2) + t76;
t85 = cos(qJ(1));
t77 = t85 * pkin(1);
t91 = t73 * pkin(2) + t72 * pkin(7) + t77;
t70 = pkin(4) * t83 + pkin(3);
t86 = -pkin(9) - pkin(8);
t90 = t81 * t70 + t84 * t86 + t93;
t89 = -t73 * pkin(7) + t92;
t88 = pkin(4) * t103 + t70 * t100 + t91;
t87 = -pkin(4) * t101 + t70 * t102 + t89;
t61 = t72 * t74 + t73 * t98;
t60 = -t72 * t75 + t73 * t99;
t59 = t72 * t98 - t73 * t74;
t58 = t72 * t99 + t73 * t75;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t85 - rSges(2,2) * t82) + g(2) * (rSges(2,1) * t82 + rSges(2,2) * t85) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t73 - rSges(3,2) * t72 + t77) + g(2) * (rSges(3,1) * t72 + rSges(3,2) * t73 + t76) + g(3) * (rSges(3,3) + t93)) - m(4) * (g(1) * (rSges(4,3) * t72 + t91) + g(2) * (rSges(4,1) * t102 - t72 * t104 + t92) + g(3) * (rSges(4,1) * t81 + rSges(4,2) * t84 + t93) + (g(1) * (rSges(4,1) * t84 - t104) + g(2) * (-rSges(4,3) - pkin(7))) * t73) - m(5) * (g(1) * (pkin(3) * t100 + (t73 * t96 + t103) * rSges(5,1) + (t72 * t83 - t73 * t97) * rSges(5,2) + t91) + g(2) * (pkin(3) * t102 + (t72 * t96 - t101) * rSges(5,1) + (-t72 * t97 - t73 * t83) * rSges(5,2) + t89) + g(3) * (-t110 * t84 + t93) + (g(3) * (rSges(5,1) * t83 - rSges(5,2) * t80 + pkin(3)) + t108 * t110) * t81) - m(6) * (g(1) * (t61 * rSges(6,1) - t60 * rSges(6,2) + t88) + g(2) * (t59 * rSges(6,1) - t58 * rSges(6,2) + t87) + g(3) * (-rSges(6,3) * t84 + t90) + (g(3) * (rSges(6,1) * t75 - rSges(6,2) * t74) + t108 * (rSges(6,3) - t86)) * t81) - m(7) * (g(1) * (t109 * t60 + t111 * t61 + t88) + g(2) * (t109 * t58 + t111 * t59 + t87) + g(3) * (-rSges(7,2) * t84 + t90) + (g(3) * (t109 * t74 + t111 * t75) + t108 * (rSges(7,2) - t86)) * t81);
U  = t1;

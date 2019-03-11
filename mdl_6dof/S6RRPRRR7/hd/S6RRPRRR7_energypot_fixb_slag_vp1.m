% Calculate potential energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:35
% EndTime: 2019-03-09 13:55:36
% DurationCPUTime: 0.39s
% Computational Cost: add. (188->103), mult. (322->124), div. (0->0), fcn. (338->10), ass. (0->36)
t102 = pkin(9) + rSges(6,3);
t96 = pkin(10) + pkin(9) + rSges(7,3);
t81 = sin(qJ(2));
t84 = cos(qJ(4));
t100 = t81 * t84;
t80 = sin(qJ(4));
t85 = cos(qJ(2));
t60 = -t80 * t85 + t100;
t104 = g(3) * t60;
t103 = t81 * pkin(2) + pkin(6);
t82 = sin(qJ(1));
t101 = t81 * t82;
t99 = t82 * t85;
t86 = cos(qJ(1));
t98 = t85 * t86;
t97 = t86 * pkin(1) + t82 * pkin(7);
t95 = qJ(3) * t81;
t75 = t82 * pkin(1);
t94 = pkin(2) * t99 + t82 * t95 + t75;
t93 = pkin(2) * t98 + t86 * t95 + t97;
t92 = pkin(3) * t99 + t86 * pkin(8) + t94;
t91 = pkin(3) * t98 + t93;
t59 = t80 * t81 + t84 * t85;
t90 = t81 * pkin(3) - t85 * qJ(3) + t103;
t78 = qJ(5) + qJ(6);
t70 = sin(t78);
t71 = cos(t78);
t83 = cos(qJ(5));
t89 = t71 * rSges(7,1) - t70 * rSges(7,2) + pkin(5) * t83 + pkin(4);
t79 = sin(qJ(5));
t88 = t70 * rSges(7,1) + t71 * rSges(7,2) + t79 * pkin(5);
t58 = t59 * t86;
t57 = -t86 * t100 + t80 * t98;
t56 = t59 * t82;
t55 = -t82 * t100 + t80 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t86 - rSges(2,2) * t82) + g(2) * (rSges(2,1) * t82 + rSges(2,2) * t86) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t82 + t97) + g(2) * (rSges(3,1) * t99 - rSges(3,2) * t101 + t75) + g(3) * (rSges(3,1) * t81 + rSges(3,2) * t85 + pkin(6)) + (g(1) * (rSges(3,1) * t85 - rSges(3,2) * t81) + g(2) * (-rSges(3,3) - pkin(7))) * t86) - m(4) * (g(1) * (rSges(4,2) * t82 + t93) + g(2) * (rSges(4,1) * t99 + rSges(4,3) * t101 + t94) + g(3) * (rSges(4,1) * t81 + (-rSges(4,3) - qJ(3)) * t85 + t103) + (g(1) * (rSges(4,1) * t85 + rSges(4,3) * t81) + g(2) * (-rSges(4,2) - pkin(7))) * t86) - m(5) * (g(1) * (rSges(5,1) * t58 - rSges(5,2) * t57 + (-rSges(5,3) - pkin(8)) * t82 + t91) + g(2) * (rSges(5,1) * t56 - rSges(5,2) * t55 + (rSges(5,3) - pkin(7)) * t86 + t92) + g(3) * (rSges(5,1) * t60 - rSges(5,2) * t59 + t90)) - m(6) * (g(1) * (t58 * pkin(4) - t82 * pkin(8) + (t58 * t83 - t79 * t82) * rSges(6,1) + (-t58 * t79 - t82 * t83) * rSges(6,2) + t102 * t57 + t91) + g(2) * (t56 * pkin(4) - t86 * pkin(7) + (t56 * t83 + t79 * t86) * rSges(6,1) + (-t56 * t79 + t83 * t86) * rSges(6,2) + t102 * t55 + t92) + g(3) * (t102 * t59 + t90) + (t83 * rSges(6,1) - t79 * rSges(6,2) + pkin(4)) * t104) - m(7) * (g(1) * (t89 * t58 + (-pkin(8) - t88) * t82 + t96 * t57 + t91) + g(2) * (t89 * t56 + (-pkin(7) + t88) * t86 + t96 * t55 + t92) + g(3) * (t96 * t59 + t90) + t89 * t104);
U  = t1;

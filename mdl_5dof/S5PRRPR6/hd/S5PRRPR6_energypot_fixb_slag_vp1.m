% Calculate potential energy for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:00
% EndTime: 2019-12-05 16:30:00
% DurationCPUTime: 0.31s
% Computational Cost: add. (219->95), mult. (435->127), div. (0->0), fcn. (511->12), ass. (0->42)
t81 = sin(pkin(5));
t106 = pkin(6) * t81;
t105 = rSges(4,3) + pkin(7);
t86 = sin(qJ(3));
t104 = t81 * t86;
t87 = sin(qJ(2));
t103 = t81 * t87;
t88 = cos(qJ(3));
t102 = t81 * t88;
t89 = cos(qJ(2));
t101 = t81 * t89;
t84 = cos(pkin(5));
t100 = t84 * t87;
t99 = t84 * t89;
t80 = sin(pkin(9));
t83 = cos(pkin(9));
t98 = t83 * pkin(1) + t80 * t106;
t97 = pkin(8) + qJ(4) + rSges(6,3);
t96 = t84 * pkin(6) + qJ(1);
t95 = qJ(4) + rSges(5,3);
t64 = -t80 * t100 + t83 * t89;
t94 = t64 * pkin(2) + t98;
t93 = pkin(2) * t103 + t96;
t62 = t83 * t100 + t80 * t89;
t75 = t80 * pkin(1);
t92 = t62 * pkin(2) - t83 * t106 + t75;
t78 = pkin(10) + qJ(5);
t73 = sin(t78);
t74 = cos(t78);
t82 = cos(pkin(10));
t91 = t74 * rSges(6,1) - t73 * rSges(6,2) + t82 * pkin(4) + pkin(3);
t79 = sin(pkin(10));
t90 = t73 * rSges(6,1) + t74 * rSges(6,2) + t79 * pkin(4) + pkin(7);
t66 = t87 * t102 + t84 * t86;
t65 = t86 * t103 - t84 * t88;
t63 = t80 * t99 + t83 * t87;
t61 = t80 * t87 - t83 * t99;
t58 = t80 * t104 + t64 * t88;
t57 = -t80 * t102 + t64 * t86;
t56 = -t83 * t104 + t62 * t88;
t55 = t83 * t102 + t62 * t86;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t83 * rSges(2,1) - t80 * rSges(2,2)) + g(2) * (t80 * rSges(2,1) + t83 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t64 * rSges(3,1) - t63 * rSges(3,2) + t98) + g(2) * (t62 * rSges(3,1) - t61 * rSges(3,2) + t75) + g(3) * (t84 * rSges(3,3) + t96) + (g(1) * rSges(3,3) * t80 + g(3) * (rSges(3,1) * t87 + rSges(3,2) * t89) + g(2) * (-rSges(3,3) - pkin(6)) * t83) * t81) - m(4) * (g(1) * (t58 * rSges(4,1) - t57 * rSges(4,2) + t105 * t63 + t94) + g(2) * (t56 * rSges(4,1) - t55 * rSges(4,2) + t105 * t61 + t92) + g(3) * (t66 * rSges(4,1) - t65 * rSges(4,2) - t105 * t101 + t93)) - m(5) * (g(1) * (t58 * pkin(3) + t63 * pkin(7) + (t58 * t82 + t63 * t79) * rSges(5,1) + (-t58 * t79 + t63 * t82) * rSges(5,2) + t95 * t57 + t94) + g(2) * (t56 * pkin(3) + t61 * pkin(7) + (t56 * t82 + t61 * t79) * rSges(5,1) + (-t56 * t79 + t61 * t82) * rSges(5,2) + t95 * t55 + t92) + g(3) * (t66 * pkin(3) - pkin(7) * t101 + (-t79 * t101 + t66 * t82) * rSges(5,1) + (-t82 * t101 - t66 * t79) * rSges(5,2) + t95 * t65 + t93)) - m(6) * (g(1) * (t97 * t57 + t91 * t58 + t90 * t63 + t94) + g(2) * (t97 * t55 + t91 * t56 + t90 * t61 + t92) + g(3) * (-t90 * t101 + t97 * t65 + t91 * t66 + t93));
U = t1;

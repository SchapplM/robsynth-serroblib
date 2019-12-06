% Calculate potential energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:11:55
% EndTime: 2019-12-05 16:11:55
% DurationCPUTime: 0.23s
% Computational Cost: add. (145->83), mult. (279->108), div. (0->0), fcn. (301->8), ass. (0->38)
t104 = rSges(6,1) + pkin(4);
t80 = sin(pkin(7));
t84 = sin(qJ(2));
t103 = t80 * t84;
t86 = cos(qJ(2));
t102 = t80 * t86;
t82 = cos(pkin(7));
t101 = t82 * t84;
t83 = sin(qJ(3));
t100 = t83 * t84;
t99 = t83 * t86;
t85 = cos(qJ(3));
t98 = t84 * t85;
t97 = t85 * t86;
t96 = t82 * pkin(1) + t80 * pkin(5);
t95 = rSges(6,2) + qJ(4);
t94 = rSges(5,3) + qJ(4);
t93 = rSges(6,3) + qJ(5);
t92 = t84 * pkin(2) + qJ(1);
t91 = t82 * t86 * pkin(2) + pkin(6) * t101 + t96;
t65 = t80 * t83 + t82 * t97;
t90 = t65 * pkin(3) + t91;
t76 = t80 * pkin(1);
t89 = pkin(2) * t102 - t82 * pkin(5) + pkin(6) * t103 + t76;
t61 = t80 * t97 - t82 * t83;
t88 = t61 * pkin(3) + t89;
t87 = pkin(3) * t98 - t86 * pkin(6) + qJ(4) * t100 + t92;
t81 = cos(pkin(8));
t79 = sin(pkin(8));
t64 = -t80 * t85 + t82 * t99;
t63 = -t86 * t79 + t81 * t98;
t62 = t79 * t98 + t86 * t81;
t60 = t80 * t99 + t82 * t85;
t57 = t79 * t101 + t65 * t81;
t56 = -t81 * t101 + t65 * t79;
t55 = t79 * t103 + t61 * t81;
t54 = -t81 * t103 + t61 * t79;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t82 * rSges(2,1) - t80 * rSges(2,2)) + g(2) * (t80 * rSges(2,1) + t82 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t80 * rSges(3,3) + t96) + g(2) * (rSges(3,1) * t102 - rSges(3,2) * t103 + t76) + g(3) * (t84 * rSges(3,1) + t86 * rSges(3,2) + qJ(1)) + (g(1) * (rSges(3,1) * t86 - rSges(3,2) * t84) + g(2) * (-rSges(3,3) - pkin(5))) * t82) - m(4) * (g(1) * (t65 * rSges(4,1) - t64 * rSges(4,2) + rSges(4,3) * t101 + t91) + g(2) * (t61 * rSges(4,1) - t60 * rSges(4,2) + rSges(4,3) * t103 + t89) + g(3) * ((-rSges(4,3) - pkin(6)) * t86 + (rSges(4,1) * t85 - rSges(4,2) * t83) * t84 + t92)) - m(5) * (g(1) * (t57 * rSges(5,1) - t56 * rSges(5,2) + t94 * t64 + t90) + g(2) * (t55 * rSges(5,1) - t54 * rSges(5,2) + t94 * t60 + t88) + g(3) * (t63 * rSges(5,1) - t62 * rSges(5,2) + rSges(5,3) * t100 + t87)) - m(6) * (g(1) * (t104 * t57 + t93 * t56 + t95 * t64 + t90) + g(2) * (t104 * t55 + t93 * t54 + t95 * t60 + t88) + g(3) * (rSges(6,2) * t100 + t104 * t63 + t93 * t62 + t87));
U = t1;

% Calculate potential energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:20
% EndTime: 2019-12-05 16:02:20
% DurationCPUTime: 0.39s
% Computational Cost: add. (175->97), mult. (358->133), div. (0->0), fcn. (405->10), ass. (0->42)
t79 = sin(pkin(9));
t108 = g(1) * t79;
t81 = cos(pkin(9));
t107 = g(2) * t81;
t106 = pkin(8) + rSges(6,3);
t80 = sin(pkin(5));
t105 = t79 * t80;
t84 = sin(qJ(4));
t104 = t80 * t84;
t85 = sin(qJ(2));
t103 = t80 * t85;
t87 = cos(qJ(4));
t102 = t80 * t87;
t88 = cos(qJ(2));
t101 = t80 * t88;
t82 = cos(pkin(5));
t100 = t82 * t85;
t99 = t82 * t88;
t98 = t81 * pkin(1) + pkin(6) * t105;
t97 = t88 * qJ(3);
t96 = t82 * pkin(6) + qJ(1);
t95 = pkin(2) * t103 + t96;
t94 = (-pkin(3) - pkin(6)) * t81;
t63 = t79 * t85 - t81 * t99;
t64 = t81 * t100 + t79 * t88;
t75 = t79 * pkin(1);
t93 = t64 * pkin(2) + t63 * qJ(3) + t75;
t92 = t82 * pkin(3) + pkin(7) * t103 + t95;
t65 = t79 * t99 + t81 * t85;
t66 = -t79 * t100 + t81 * t88;
t91 = t66 * pkin(2) + t65 * qJ(3) + t98;
t90 = pkin(3) * t105 + t91;
t89 = t64 * pkin(7) + t93;
t86 = cos(qJ(5));
t83 = sin(qJ(5));
t68 = -t84 * t101 + t82 * t87;
t67 = t87 * t101 + t82 * t84;
t59 = -t81 * t102 + t63 * t84;
t58 = t81 * t104 + t63 * t87;
t57 = t79 * t102 + t65 * t84;
t56 = t79 * t104 - t65 * t87;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t81 - rSges(2,2) * t79) + g(2) * (rSges(2,1) * t79 + rSges(2,2) * t81) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t66 - rSges(3,2) * t65 + t98) + g(2) * (rSges(3,1) * t64 - rSges(3,2) * t63 + t75) + g(3) * (t82 * rSges(3,3) + t96) + (rSges(3,3) * t108 + g(3) * (rSges(3,1) * t85 + rSges(3,2) * t88) + (-rSges(3,3) - pkin(6)) * t107) * t80) - m(4) * (g(1) * (-rSges(4,2) * t66 + rSges(4,3) * t65 + t91) + g(2) * (-rSges(4,2) * t64 + rSges(4,3) * t63 + t93) + g(3) * (t82 * rSges(4,1) + t95) + (rSges(4,1) * t108 + g(3) * (-rSges(4,2) * t85 - rSges(4,3) * t88 - t97) + (-rSges(4,1) - pkin(6)) * t107) * t80) - m(5) * (g(1) * (rSges(5,1) * t57 - rSges(5,2) * t56 + (rSges(5,3) + pkin(7)) * t66 + t90) + g(2) * (rSges(5,1) * t59 + rSges(5,2) * t58 + rSges(5,3) * t64 + t89) + g(3) * (t68 * rSges(5,1) - t67 * rSges(5,2) + t92) + (g(3) * (rSges(5,3) * t85 - t97) + g(2) * t94) * t80) - m(6) * (g(1) * (t57 * pkin(4) + t66 * pkin(7) + (t57 * t86 + t66 * t83) * rSges(6,1) + (-t57 * t83 + t66 * t86) * rSges(6,2) + t106 * t56 + t90) + g(2) * (t59 * pkin(4) + (t59 * t86 + t64 * t83) * rSges(6,1) + (-t59 * t83 + t64 * t86) * rSges(6,2) + t80 * t94 - t106 * t58 + t89) + g(3) * (t68 * pkin(4) - t80 * t97 + (t83 * t103 + t68 * t86) * rSges(6,1) + (t86 * t103 - t68 * t83) * rSges(6,2) + t106 * t67 + t92));
U = t1;

% Calculate potential energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR12_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:36:58
% EndTime: 2019-12-31 21:36:59
% DurationCPUTime: 0.34s
% Computational Cost: add. (219->95), mult. (435->125), div. (0->0), fcn. (511->12), ass. (0->44)
t110 = rSges(4,3) + pkin(8);
t84 = cos(pkin(5));
t109 = t84 * pkin(7) + pkin(6);
t82 = sin(pkin(5));
t87 = sin(qJ(2));
t108 = t82 * t87;
t88 = sin(qJ(1));
t107 = t82 * t88;
t89 = cos(qJ(3));
t106 = t82 * t89;
t90 = cos(qJ(2));
t105 = t82 * t90;
t91 = cos(qJ(1));
t104 = t82 * t91;
t103 = t88 * t87;
t102 = t88 * t90;
t101 = t91 * t87;
t100 = t91 * t90;
t99 = t91 * pkin(1) + pkin(7) * t107;
t98 = pkin(9) + qJ(4) + rSges(6,3);
t97 = qJ(4) + rSges(5,3);
t96 = pkin(2) * t108 + t109;
t68 = -t84 * t103 + t100;
t95 = t68 * pkin(2) + t99;
t66 = t84 * t101 + t102;
t78 = t88 * pkin(1);
t94 = t66 * pkin(2) - pkin(7) * t104 + t78;
t80 = pkin(10) + qJ(5);
t75 = sin(t80);
t76 = cos(t80);
t83 = cos(pkin(10));
t93 = t76 * rSges(6,1) - t75 * rSges(6,2) + t83 * pkin(4) + pkin(3);
t81 = sin(pkin(10));
t92 = t75 * rSges(6,1) + t76 * rSges(6,2) + t81 * pkin(4) + pkin(8);
t86 = sin(qJ(3));
t67 = t84 * t102 + t101;
t65 = -t84 * t100 + t103;
t64 = t87 * t106 + t84 * t86;
t63 = t86 * t108 - t84 * t89;
t60 = t86 * t107 + t68 * t89;
t59 = -t88 * t106 + t68 * t86;
t58 = -t86 * t104 + t66 * t89;
t57 = t89 * t104 + t66 * t86;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t91 * rSges(2,1) - t88 * rSges(2,2)) + g(2) * (t88 * rSges(2,1) + t91 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t68 * rSges(3,1) - t67 * rSges(3,2) + t99) + g(2) * (t66 * rSges(3,1) - t65 * rSges(3,2) + t78) + g(3) * (t84 * rSges(3,3) + t109) + (g(1) * rSges(3,3) * t88 + g(3) * (rSges(3,1) * t87 + rSges(3,2) * t90) + g(2) * (-rSges(3,3) - pkin(7)) * t91) * t82) - m(4) * (g(1) * (t60 * rSges(4,1) - t59 * rSges(4,2) + t110 * t67 + t95) + g(2) * (t58 * rSges(4,1) - t57 * rSges(4,2) + t110 * t65 + t94) + g(3) * (t64 * rSges(4,1) - t63 * rSges(4,2) - t110 * t105 + t96)) - m(5) * (g(1) * (t60 * pkin(3) + t67 * pkin(8) + (t60 * t83 + t67 * t81) * rSges(5,1) + (-t60 * t81 + t67 * t83) * rSges(5,2) + t97 * t59 + t95) + g(2) * (t58 * pkin(3) + t65 * pkin(8) + (t58 * t83 + t65 * t81) * rSges(5,1) + (-t58 * t81 + t65 * t83) * rSges(5,2) + t97 * t57 + t94) + g(3) * (t64 * pkin(3) - pkin(8) * t105 + (-t81 * t105 + t64 * t83) * rSges(5,1) + (-t83 * t105 - t64 * t81) * rSges(5,2) + t97 * t63 + t96)) - m(6) * (g(1) * (t98 * t59 + t93 * t60 + t92 * t67 + t95) + g(2) * (t98 * t57 + t93 * t58 + t92 * t65 + t94) + g(3) * (-t92 * t105 + t98 * t63 + t93 * t64 + t96));
U = t1;

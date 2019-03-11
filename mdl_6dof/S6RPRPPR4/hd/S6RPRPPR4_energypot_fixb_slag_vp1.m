% Calculate potential energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:41
% EndTime: 2019-03-09 02:46:41
% DurationCPUTime: 0.49s
% Computational Cost: add. (239->101), mult. (274->127), div. (0->0), fcn. (278->10), ass. (0->39)
t104 = -rSges(7,3) - pkin(8);
t77 = sin(pkin(9));
t103 = t77 * pkin(2) + pkin(6);
t75 = pkin(9) + qJ(3);
t73 = cos(t75);
t84 = cos(qJ(1));
t102 = t73 * t84;
t72 = sin(t75);
t82 = sin(qJ(1));
t101 = t82 * t72;
t76 = sin(pkin(10));
t100 = t82 * t76;
t78 = cos(pkin(10));
t99 = t82 * t78;
t98 = t84 * t72;
t97 = t84 * t76;
t96 = t84 * t78;
t79 = cos(pkin(9));
t70 = t79 * pkin(2) + pkin(1);
t80 = -pkin(7) - qJ(2);
t95 = t82 * t70 + t84 * t80;
t94 = qJ(4) * t72;
t93 = rSges(3,3) + qJ(2);
t92 = rSges(6,3) + qJ(5);
t91 = t72 * pkin(3) + t103;
t90 = t95 + (pkin(3) * t73 + t94) * t82;
t89 = t91 + (pkin(4) * t78 + qJ(5) * t76) * t72;
t57 = t73 * t99 - t97;
t88 = t57 * pkin(4) + t90;
t65 = t84 * t70;
t87 = pkin(3) * t102 - t82 * t80 + t84 * t94 + t65;
t86 = rSges(3,1) * t79 - rSges(3,2) * t77 + pkin(1);
t59 = t73 * t96 + t100;
t85 = t59 * pkin(4) + t87;
t83 = cos(qJ(6));
t81 = sin(qJ(6));
t58 = t73 * t97 - t99;
t56 = t73 * t100 + t96;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t84 * rSges(2,1) - t82 * rSges(2,2)) + g(2) * (t82 * rSges(2,1) + t84 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t77 * rSges(3,1) + t79 * rSges(3,2) + pkin(6)) + (g(1) * t86 - g(2) * t93) * t84 + (g(1) * t93 + g(2) * t86) * t82) - m(4) * (g(1) * (rSges(4,1) * t102 - rSges(4,2) * t98 + t65) + g(2) * (-t84 * rSges(4,3) + t95) + g(3) * (t72 * rSges(4,1) + t73 * rSges(4,2) + t103) + (g(1) * (rSges(4,3) - t80) + g(2) * (rSges(4,1) * t73 - rSges(4,2) * t72)) * t82) - m(5) * (g(1) * (t59 * rSges(5,1) - t58 * rSges(5,2) + rSges(5,3) * t98 + t87) + g(2) * (t57 * rSges(5,1) - t56 * rSges(5,2) + rSges(5,3) * t101 + t90) + g(3) * ((-rSges(5,3) - qJ(4)) * t73 + (rSges(5,1) * t78 - rSges(5,2) * t76) * t72 + t91)) - m(6) * (g(1) * (t59 * rSges(6,1) + rSges(6,2) * t98 + t92 * t58 + t85) + g(2) * (t57 * rSges(6,1) + rSges(6,2) * t101 + t92 * t56 + t88) + g(3) * ((-rSges(6,2) - qJ(4)) * t73 + (rSges(6,1) * t78 + rSges(6,3) * t76) * t72 + t89)) - m(7) * (g(1) * (t59 * pkin(5) + t58 * qJ(5) + (t58 * t81 + t59 * t83) * rSges(7,1) + (t58 * t83 - t59 * t81) * rSges(7,2) + t85) + g(2) * (t57 * pkin(5) + t56 * qJ(5) + (t56 * t81 + t57 * t83) * rSges(7,1) + (t56 * t83 - t57 * t81) * rSges(7,2) + t88) + (g(1) * t84 + g(2) * t82) * t72 * t104 + (t89 + (-qJ(4) - t104) * t73 + (t78 * pkin(5) + (t76 * t81 + t78 * t83) * rSges(7,1) + (t76 * t83 - t78 * t81) * rSges(7,2)) * t72) * g(3));
U  = t1;

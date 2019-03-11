% Calculate potential energy for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:11
% EndTime: 2019-03-09 21:05:11
% DurationCPUTime: 0.47s
% Computational Cost: add. (224->99), mult. (280->122), div. (0->0), fcn. (278->8), ass. (0->34)
t104 = rSges(7,1) + pkin(5);
t103 = rSges(4,3) + pkin(8);
t102 = rSges(7,3) + qJ(6);
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t101 = g(1) * t80 + g(2) * t77;
t76 = sin(qJ(2));
t97 = rSges(3,2) * t76;
t75 = sin(qJ(3));
t96 = t77 * t75;
t79 = cos(qJ(2));
t95 = t77 * t79;
t94 = t80 * t75;
t93 = t80 * t79;
t90 = t80 * pkin(1) + t77 * pkin(7);
t78 = cos(qJ(3));
t67 = t78 * pkin(3) + pkin(2);
t81 = -pkin(9) - pkin(8);
t89 = t76 * t67 + t79 * t81 + pkin(6);
t72 = t77 * pkin(1);
t87 = -t80 * pkin(7) + t72;
t86 = pkin(3) * t96 + t67 * t93 + t90;
t74 = qJ(3) + qJ(4);
t69 = sin(t74);
t70 = cos(t74);
t85 = t89 + (pkin(4) * t70 + qJ(5) * t69) * t76;
t58 = t69 * t93 - t77 * t70;
t59 = t77 * t69 + t70 * t93;
t84 = t59 * pkin(4) + t58 * qJ(5) + t86;
t83 = -pkin(3) * t94 + t67 * t95 + t87;
t56 = t69 * t95 + t80 * t70;
t57 = -t80 * t69 + t70 * t95;
t82 = t57 * pkin(4) + t56 * qJ(5) + t83;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t80 * rSges(2,1) - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + t80 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t77 * rSges(3,3) + t90) + g(2) * (rSges(3,1) * t95 - t77 * t97 + t72) + g(3) * (t76 * rSges(3,1) + t79 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t79 - t97) + g(2) * (-rSges(3,3) - pkin(7))) * t80) - m(4) * (g(1) * (pkin(2) * t93 + (t78 * t93 + t96) * rSges(4,1) + (-t75 * t93 + t77 * t78) * rSges(4,2) + t90) + g(2) * (pkin(2) * t95 + (t78 * t95 - t94) * rSges(4,1) + (-t75 * t95 - t80 * t78) * rSges(4,2) + t87) + g(3) * (-t103 * t79 + pkin(6)) + (g(3) * (rSges(4,1) * t78 - rSges(4,2) * t75 + pkin(2)) + t101 * t103) * t76) - m(5) * (g(1) * (t59 * rSges(5,1) - t58 * rSges(5,2) + t86) + g(2) * (t57 * rSges(5,1) - t56 * rSges(5,2) + t83) + g(3) * (-t79 * rSges(5,3) + t89) + (g(3) * (rSges(5,1) * t70 - rSges(5,2) * t69) + t101 * (rSges(5,3) - t81)) * t76) - m(6) * (g(1) * (t59 * rSges(6,1) + t58 * rSges(6,3) + t84) + g(2) * (t57 * rSges(6,1) + t56 * rSges(6,3) + t82) + g(3) * (-t79 * rSges(6,2) + t85) + (g(3) * (rSges(6,1) * t70 + rSges(6,3) * t69) + t101 * (rSges(6,2) - t81)) * t76) - m(7) * (g(1) * (t58 * rSges(7,2) + t104 * t59 + t84) + g(2) * (t56 * rSges(7,2) + t104 * t57 + t82) + g(3) * (t102 * t79 + t85) + (g(3) * (rSges(7,2) * t69 + t104 * t70) + t101 * (-t81 - t102)) * t76);
U  = t1;

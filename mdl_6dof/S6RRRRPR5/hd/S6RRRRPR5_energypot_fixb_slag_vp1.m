% Calculate potential energy for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:11:56
% EndTime: 2019-03-09 22:11:56
% DurationCPUTime: 0.49s
% Computational Cost: add. (239->101), mult. (274->127), div. (0->0), fcn. (278->10), ass. (0->38)
t102 = rSges(3,3) + pkin(7);
t101 = -rSges(7,3) - pkin(10);
t77 = sin(qJ(2));
t100 = t77 * pkin(2) + pkin(6);
t74 = qJ(2) + qJ(3);
t72 = cos(t74);
t82 = cos(qJ(1));
t99 = t72 * t82;
t71 = sin(t74);
t78 = sin(qJ(1));
t98 = t78 * t71;
t76 = sin(qJ(4));
t97 = t78 * t76;
t80 = cos(qJ(4));
t96 = t78 * t80;
t95 = t82 * t71;
t94 = t82 * t76;
t93 = t82 * t80;
t81 = cos(qJ(2));
t69 = t81 * pkin(2) + pkin(1);
t83 = -pkin(8) - pkin(7);
t92 = t78 * t69 + t82 * t83;
t91 = rSges(6,3) + qJ(5);
t90 = t71 * pkin(3) + t100;
t89 = t78 * t72 * pkin(3) + pkin(9) * t98 + t92;
t88 = t90 + (pkin(4) * t80 + qJ(5) * t76) * t71;
t56 = t72 * t96 - t94;
t87 = t56 * pkin(4) + t89;
t61 = t82 * t69;
t86 = pkin(3) * t99 + pkin(9) * t95 - t78 * t83 + t61;
t85 = rSges(3,1) * t81 - rSges(3,2) * t77 + pkin(1);
t58 = t72 * t93 + t97;
t84 = t58 * pkin(4) + t86;
t79 = cos(qJ(6));
t75 = sin(qJ(6));
t57 = t72 * t94 - t96;
t55 = t72 * t97 + t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t82 * rSges(2,1) - t78 * rSges(2,2)) + g(2) * (t78 * rSges(2,1) + t82 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t77 * rSges(3,1) + t81 * rSges(3,2) + pkin(6)) + (g(1) * t85 - g(2) * t102) * t82 + (g(1) * t102 + g(2) * t85) * t78) - m(4) * (g(1) * (rSges(4,1) * t99 - rSges(4,2) * t95 + t61) + g(2) * (-t82 * rSges(4,3) + t92) + g(3) * (t71 * rSges(4,1) + t72 * rSges(4,2) + t100) + (g(1) * (rSges(4,3) - t83) + g(2) * (rSges(4,1) * t72 - rSges(4,2) * t71)) * t78) - m(5) * (g(1) * (t58 * rSges(5,1) - t57 * rSges(5,2) + rSges(5,3) * t95 + t86) + g(2) * (t56 * rSges(5,1) - t55 * rSges(5,2) + rSges(5,3) * t98 + t89) + g(3) * ((-rSges(5,3) - pkin(9)) * t72 + (rSges(5,1) * t80 - rSges(5,2) * t76) * t71 + t90)) - m(6) * (g(1) * (t58 * rSges(6,1) + rSges(6,2) * t95 + t91 * t57 + t84) + g(2) * (t56 * rSges(6,1) + rSges(6,2) * t98 + t91 * t55 + t87) + g(3) * ((-rSges(6,2) - pkin(9)) * t72 + (rSges(6,1) * t80 + rSges(6,3) * t76) * t71 + t88)) - m(7) * (g(1) * (t58 * pkin(5) + t57 * qJ(5) + (t57 * t75 + t58 * t79) * rSges(7,1) + (t57 * t79 - t58 * t75) * rSges(7,2) + t84) + g(2) * (t56 * pkin(5) + t55 * qJ(5) + (t55 * t75 + t56 * t79) * rSges(7,1) + (t55 * t79 - t56 * t75) * rSges(7,2) + t87) + (g(1) * t82 + g(2) * t78) * t71 * t101 + (t88 + (-pkin(9) - t101) * t72 + (t80 * pkin(5) + (t75 * t76 + t79 * t80) * rSges(7,1) + (-t75 * t80 + t76 * t79) * rSges(7,2)) * t71) * g(3));
U  = t1;

% Calculate potential energy for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:47
% EndTime: 2019-03-10 03:27:47
% DurationCPUTime: 0.35s
% Computational Cost: add. (231->84), mult. (162->100), div. (0->0), fcn. (138->12), ass. (0->40)
t97 = rSges(7,3) + pkin(11);
t79 = -pkin(8) - pkin(7);
t96 = rSges(3,3) + pkin(7);
t74 = sin(qJ(2));
t94 = t74 * pkin(2) + pkin(6);
t77 = cos(qJ(2));
t63 = t77 * pkin(2) + pkin(1);
t72 = qJ(2) + qJ(3);
t67 = qJ(4) + t72;
t64 = qJ(5) + t67;
t57 = sin(t64);
t93 = rSges(6,2) * t57;
t73 = sin(qJ(6));
t75 = sin(qJ(1));
t92 = t75 * t73;
t76 = cos(qJ(6));
t91 = t75 * t76;
t58 = cos(t64);
t78 = cos(qJ(1));
t90 = t78 * t58;
t89 = t78 * t73;
t88 = t78 * t76;
t87 = rSges(4,3) - t79;
t71 = -pkin(9) + t79;
t86 = rSges(5,3) - t71;
t66 = cos(t72);
t54 = pkin(3) * t66 + t63;
t62 = cos(t67);
t53 = pkin(4) * t62 + t54;
t68 = -pkin(10) + t71;
t85 = t75 * t53 + t78 * t68;
t65 = sin(t72);
t84 = pkin(3) * t65 + t94;
t61 = sin(t67);
t83 = pkin(4) * t61 + t84;
t82 = rSges(3,1) * t77 - rSges(3,2) * t74 + pkin(1);
t81 = rSges(4,1) * t66 - rSges(4,2) * t65 + t63;
t80 = rSges(5,1) * t62 - rSges(5,2) * t61 + t54;
t52 = t78 * t53;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t78 * rSges(2,1) - t75 * rSges(2,2)) + g(2) * (t75 * rSges(2,1) + t78 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t74 * rSges(3,1) + t77 * rSges(3,2) + pkin(6)) + (g(1) * t82 - g(2) * t96) * t78 + (g(1) * t96 + g(2) * t82) * t75) - m(4) * (g(3) * (t65 * rSges(4,1) + t66 * rSges(4,2) + t94) + (g(1) * t81 - g(2) * t87) * t78 + (g(1) * t87 + g(2) * t81) * t75) - m(5) * (g(3) * (t61 * rSges(5,1) + t62 * rSges(5,2) + t84) + (g(1) * t80 - g(2) * t86) * t78 + (g(1) * t86 + g(2) * t80) * t75) - m(6) * (g(1) * (rSges(6,1) * t90 - t78 * t93 + t52) + g(2) * (-t78 * rSges(6,3) + t85) + g(3) * (t57 * rSges(6,1) + t58 * rSges(6,2) + t83) + (g(1) * (rSges(6,3) - t68) + g(2) * (rSges(6,1) * t58 - t93)) * t75) - m(7) * (g(1) * (pkin(5) * t90 + t52 - t75 * t68 + (t58 * t88 + t92) * rSges(7,1) + (-t58 * t89 + t91) * rSges(7,2)) + g(2) * (t75 * t58 * pkin(5) + (t58 * t91 - t89) * rSges(7,1) + (-t58 * t92 - t88) * rSges(7,2) + t85) + g(3) * (-t97 * t58 + t83) + (g(3) * (rSges(7,1) * t76 - rSges(7,2) * t73 + pkin(5)) + (g(1) * t78 + g(2) * t75) * t97) * t57);
U  = t1;

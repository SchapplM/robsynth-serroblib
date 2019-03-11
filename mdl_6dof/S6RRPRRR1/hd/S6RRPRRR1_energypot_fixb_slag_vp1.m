% Calculate potential energy for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:22
% EndTime: 2019-03-09 13:12:22
% DurationCPUTime: 0.39s
% Computational Cost: add. (231->84), mult. (162->100), div. (0->0), fcn. (138->12), ass. (0->38)
t95 = rSges(7,3) + pkin(10);
t94 = rSges(3,3) + pkin(7);
t75 = sin(qJ(2));
t92 = t75 * pkin(2) + pkin(6);
t78 = cos(qJ(2));
t64 = t78 * pkin(2) + pkin(1);
t72 = qJ(2) + pkin(11);
t67 = qJ(4) + t72;
t63 = qJ(5) + t67;
t57 = sin(t63);
t91 = rSges(6,2) * t57;
t74 = sin(qJ(6));
t76 = sin(qJ(1));
t90 = t76 * t74;
t77 = cos(qJ(6));
t89 = t76 * t77;
t58 = cos(t63);
t79 = cos(qJ(1));
t88 = t79 * t58;
t73 = -qJ(3) - pkin(7);
t87 = rSges(4,3) - t73;
t71 = -pkin(8) + t73;
t86 = rSges(5,3) - t71;
t66 = cos(t72);
t54 = pkin(3) * t66 + t64;
t62 = cos(t67);
t53 = pkin(4) * t62 + t54;
t68 = -pkin(9) + t71;
t85 = t76 * t53 + t79 * t68;
t65 = sin(t72);
t84 = pkin(3) * t65 + t92;
t61 = sin(t67);
t83 = pkin(4) * t61 + t84;
t82 = rSges(3,1) * t78 - rSges(3,2) * t75 + pkin(1);
t81 = rSges(4,1) * t66 - rSges(4,2) * t65 + t64;
t80 = rSges(5,1) * t62 - rSges(5,2) * t61 + t54;
t52 = t79 * t53;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - t76 * rSges(2,2)) + g(2) * (t76 * rSges(2,1) + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t75 + rSges(3,2) * t78 + pkin(6)) + (g(1) * t82 - g(2) * t94) * t79 + (g(1) * t94 + g(2) * t82) * t76) - m(4) * (g(3) * (rSges(4,1) * t65 + rSges(4,2) * t66 + t92) + (g(1) * t81 - g(2) * t87) * t79 + (g(1) * t87 + g(2) * t81) * t76) - m(5) * (g(3) * (rSges(5,1) * t61 + rSges(5,2) * t62 + t84) + (g(1) * t80 - g(2) * t86) * t79 + (g(1) * t86 + g(2) * t80) * t76) - m(6) * (g(1) * (rSges(6,1) * t88 - t79 * t91 + t52) + g(2) * (-rSges(6,3) * t79 + t85) + g(3) * (rSges(6,1) * t57 + rSges(6,2) * t58 + t83) + (g(1) * (rSges(6,3) - t68) + g(2) * (rSges(6,1) * t58 - t91)) * t76) - m(7) * (g(1) * (pkin(5) * t88 + t52 - t76 * t68 + (t77 * t88 + t90) * rSges(7,1) + (-t74 * t88 + t89) * rSges(7,2)) + g(2) * (t76 * t58 * pkin(5) + (t58 * t89 - t74 * t79) * rSges(7,1) + (-t58 * t90 - t77 * t79) * rSges(7,2) + t85) + g(3) * (-t95 * t58 + t83) + (g(3) * (rSges(7,1) * t77 - rSges(7,2) * t74 + pkin(5)) + (g(1) * t79 + g(2) * t76) * t95) * t57);
U  = t1;

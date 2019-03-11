% Calculate potential energy for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:43
% EndTime: 2019-03-09 07:04:43
% DurationCPUTime: 0.35s
% Computational Cost: add. (231->84), mult. (162->100), div. (0->0), fcn. (138->12), ass. (0->40)
t97 = rSges(7,3) + pkin(10);
t73 = sin(pkin(11));
t95 = t73 * pkin(2) + pkin(6);
t74 = cos(pkin(11));
t63 = t74 * pkin(2) + pkin(1);
t72 = pkin(11) + qJ(3);
t67 = qJ(4) + t72;
t64 = qJ(5) + t67;
t57 = sin(t64);
t94 = rSges(6,2) * t57;
t76 = sin(qJ(6));
t77 = sin(qJ(1));
t93 = t77 * t76;
t78 = cos(qJ(6));
t92 = t77 * t78;
t58 = cos(t64);
t79 = cos(qJ(1));
t91 = t79 * t58;
t90 = t79 * t76;
t89 = t79 * t78;
t75 = -pkin(7) - qJ(2);
t88 = rSges(4,3) - t75;
t71 = -pkin(8) + t75;
t87 = rSges(5,3) - t71;
t66 = cos(t72);
t54 = pkin(3) * t66 + t63;
t62 = cos(t67);
t53 = pkin(4) * t62 + t54;
t68 = -pkin(9) + t71;
t86 = t77 * t53 + t79 * t68;
t85 = rSges(3,3) + qJ(2);
t65 = sin(t72);
t84 = pkin(3) * t65 + t95;
t61 = sin(t67);
t83 = pkin(4) * t61 + t84;
t82 = rSges(3,1) * t74 - rSges(3,2) * t73 + pkin(1);
t81 = rSges(4,1) * t66 - rSges(4,2) * t65 + t63;
t80 = rSges(5,1) * t62 - rSges(5,2) * t61 + t54;
t52 = t79 * t53;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t79 * rSges(2,1) - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + t79 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t73 * rSges(3,1) + t74 * rSges(3,2) + pkin(6)) + (g(1) * t82 - g(2) * t85) * t79 + (g(1) * t85 + g(2) * t82) * t77) - m(4) * (g(3) * (t65 * rSges(4,1) + t66 * rSges(4,2) + t95) + (g(1) * t81 - g(2) * t88) * t79 + (g(1) * t88 + g(2) * t81) * t77) - m(5) * (g(3) * (t61 * rSges(5,1) + t62 * rSges(5,2) + t84) + (g(1) * t80 - g(2) * t87) * t79 + (g(1) * t87 + g(2) * t80) * t77) - m(6) * (g(1) * (rSges(6,1) * t91 - t79 * t94 + t52) + g(2) * (-t79 * rSges(6,3) + t86) + g(3) * (t57 * rSges(6,1) + t58 * rSges(6,2) + t83) + (g(1) * (rSges(6,3) - t68) + g(2) * (rSges(6,1) * t58 - t94)) * t77) - m(7) * (g(1) * (pkin(5) * t91 + t52 - t77 * t68 + (t58 * t89 + t93) * rSges(7,1) + (-t58 * t90 + t92) * rSges(7,2)) + g(2) * (t77 * t58 * pkin(5) + (t58 * t92 - t90) * rSges(7,1) + (-t58 * t93 - t89) * rSges(7,2) + t86) + g(3) * (-t97 * t58 + t83) + (g(3) * (rSges(7,1) * t78 - rSges(7,2) * t76 + pkin(5)) + (g(1) * t79 + g(2) * t77) * t97) * t57);
U  = t1;

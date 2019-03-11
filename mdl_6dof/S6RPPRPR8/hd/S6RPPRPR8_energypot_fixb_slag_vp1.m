% Calculate potential energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:13
% EndTime: 2019-03-09 01:55:14
% DurationCPUTime: 0.35s
% Computational Cost: add. (155->84), mult. (167->91), div. (0->0), fcn. (143->8), ass. (0->32)
t59 = sin(qJ(6));
t61 = cos(qJ(6));
t66 = rSges(7,1) * t59 + rSges(7,2) * t61;
t84 = pkin(2) + pkin(6);
t60 = sin(qJ(1));
t53 = t60 * pkin(1);
t83 = g(2) * t53;
t55 = pkin(9) + qJ(4);
t49 = sin(t55);
t82 = t49 * pkin(4);
t56 = sin(pkin(9));
t81 = t56 * pkin(3);
t80 = rSges(7,3) + pkin(8);
t58 = -pkin(7) - qJ(3);
t77 = rSges(6,1) - t58;
t76 = rSges(5,3) - t58;
t62 = cos(qJ(1));
t75 = t62 * pkin(1) + t60 * qJ(2);
t50 = cos(t55);
t74 = t50 * qJ(5);
t73 = rSges(4,3) + qJ(3);
t57 = cos(pkin(9));
t72 = t57 * pkin(3) + t84;
t71 = t60 * t81 + t75;
t70 = -qJ(2) - t81;
t69 = t50 * pkin(4) + t49 * qJ(5) + t72;
t68 = rSges(4,1) * t56 + rSges(4,2) * t57;
t67 = rSges(5,1) * t49 + rSges(5,2) * t50;
t65 = -rSges(6,2) * t49 - rSges(6,3) * t50;
t64 = t61 * rSges(7,1) - t59 * rSges(7,2) + pkin(5) - t58;
t63 = g(1) * (t60 * t82 + t71) + g(2) * (t62 * t74 + t53);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t62 - t60 * rSges(2,2)) + g(2) * (t60 * rSges(2,1) + rSges(2,2) * t62) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t62 + t60 * rSges(3,3) + t75) + g(2) * (-t60 * rSges(3,2) + t53 + (-rSges(3,3) - qJ(2)) * t62) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * t75 + t83 + g(3) * (rSges(4,1) * t57 - rSges(4,2) * t56 + t84) + (g(1) * t68 + g(2) * t73) * t60 + (g(1) * t73 + g(2) * (-qJ(2) - t68)) * t62) - m(5) * (g(1) * t71 + t83 + g(3) * (rSges(5,1) * t50 - rSges(5,2) * t49 + t72) + (g(1) * t67 + g(2) * t76) * t60 + (g(1) * t76 + g(2) * (-t67 + t70)) * t62) - m(6) * (g(3) * (-rSges(6,2) * t50 + rSges(6,3) * t49 + t69) + (g(1) * (t65 - t74) + g(2) * t77) * t60 + (g(1) * t77 + g(2) * (-t65 + t70 - t82)) * t62 + t63) - m(7) * (g(3) * (t49 * t66 + t50 * t80 + t69) + (g(2) * t64 + (t80 * t49 + (-qJ(5) - t66) * t50) * g(1)) * t60 + (g(1) * t64 + (t70 + (-pkin(4) - t80) * t49 + t66 * t50) * g(2)) * t62 + t63);
U  = t1;

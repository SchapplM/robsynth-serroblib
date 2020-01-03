% Calculate potential energy for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR15_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:04
% EndTime: 2019-12-31 20:41:05
% DurationCPUTime: 0.41s
% Computational Cost: add. (116->84), mult. (172->107), div. (0->0), fcn. (156->8), ass. (0->28)
t82 = rSges(5,3) + pkin(7);
t81 = rSges(6,3) + pkin(8) + pkin(7);
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t80 = g(1) * t61 + g(2) * t58;
t57 = sin(qJ(2));
t76 = t57 * pkin(2) + pkin(5);
t75 = t57 * t58;
t74 = t57 * t61;
t59 = cos(qJ(4));
t73 = t58 * t59;
t60 = cos(qJ(2));
t72 = t58 * t60;
t71 = t59 * t61;
t69 = t61 * pkin(1) + t58 * pkin(6);
t68 = qJ(3) * t57;
t56 = sin(qJ(4));
t67 = t56 * t75;
t66 = t56 * t74;
t53 = t58 * pkin(1);
t65 = pkin(2) * t72 + t58 * t68 + t53;
t64 = t69 + (pkin(2) * t60 + t68) * t61;
t63 = -t61 * pkin(6) + t65;
t55 = qJ(4) + qJ(5);
t50 = cos(t55);
t49 = sin(t55);
t48 = pkin(4) * t59 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t61 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t61) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t58 + t69) + g(2) * (rSges(3,1) * t72 - rSges(3,2) * t75 + t53) + g(3) * (rSges(3,1) * t57 + rSges(3,2) * t60 + pkin(5)) + (g(1) * (rSges(3,1) * t60 - rSges(3,2) * t57) + g(2) * (-rSges(3,3) - pkin(6))) * t61) - m(4) * (g(1) * (rSges(4,1) * t58 + t64) + g(2) * (-rSges(4,2) * t72 + rSges(4,3) * t75 + t65) + g(3) * (-rSges(4,2) * t57 + (-rSges(4,3) - qJ(3)) * t60 + t76) + (g(1) * (-rSges(4,2) * t60 + rSges(4,3) * t57) + g(2) * (-rSges(4,1) - pkin(6))) * t61) - m(5) * (g(1) * (t58 * pkin(3) + (t66 + t73) * rSges(5,1) + (-t56 * t58 + t57 * t71) * rSges(5,2) + t64) + g(2) * (-t61 * pkin(3) + (t67 - t71) * rSges(5,1) + (t56 * t61 + t57 * t73) * rSges(5,2) + t63) + g(3) * (t82 * t57 + t76) + (g(3) * (-rSges(5,1) * t56 - rSges(5,2) * t59 - qJ(3)) + t80 * t82) * t60) - m(6) * (g(1) * (t58 * t48 + pkin(4) * t66 + (t49 * t74 + t50 * t58) * rSges(6,1) + (-t49 * t58 + t50 * t74) * rSges(6,2) + t64) + g(2) * (-t61 * t48 + pkin(4) * t67 + (t49 * t75 - t50 * t61) * rSges(6,1) + (t49 * t61 + t50 * t75) * rSges(6,2) + t63) + g(3) * (t81 * t57 + t76) + (g(3) * (-rSges(6,1) * t49 - rSges(6,2) * t50 - pkin(4) * t56 - qJ(3)) + t80 * t81) * t60);
U = t1;

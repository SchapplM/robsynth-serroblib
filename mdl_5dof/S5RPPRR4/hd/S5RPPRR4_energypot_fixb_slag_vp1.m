% Calculate potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:20
% EndTime: 2020-01-03 11:30:21
% DurationCPUTime: 0.44s
% Computational Cost: add. (152->78), mult. (180->96), div. (0->0), fcn. (168->10), ass. (0->34)
t62 = -pkin(6) - qJ(3);
t88 = rSges(5,3) - t62;
t87 = rSges(6,3) + pkin(7) - t62;
t63 = sin(qJ(1));
t64 = cos(qJ(1));
t86 = g(2) * t63 - g(3) * t64;
t85 = rSges(4,3) + qJ(3);
t58 = sin(pkin(9));
t84 = -pkin(3) * t58 - qJ(2);
t83 = g(3) * pkin(1);
t55 = t63 * pkin(1);
t81 = g(2) * t55;
t61 = cos(pkin(8));
t80 = g(2) * t61;
t78 = g(3) * t61;
t60 = cos(pkin(9));
t50 = t60 * pkin(3) + pkin(2);
t76 = t61 * t64;
t75 = t63 * t58;
t74 = t63 * t60;
t71 = -rSges(3,3) - qJ(2);
t57 = pkin(9) + qJ(4);
t59 = sin(pkin(8));
t69 = rSges(3,1) * t61 - rSges(3,2) * t59;
t51 = sin(t57);
t52 = cos(t57);
t68 = rSges(5,1) * t52 - rSges(5,2) * t51 + t50;
t53 = qJ(5) + t57;
t48 = sin(t53);
t49 = cos(t53);
t67 = rSges(6,1) * t49 - rSges(6,2) * t48 + pkin(4) * t52 + t50;
t66 = -rSges(6,1) * t48 - rSges(6,2) * t49 - pkin(4) * t51 + t84;
t65 = -rSges(5,1) * t51 - rSges(5,2) * t52 + t84;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (t63 * rSges(2,1) + rSges(2,2) * t64) + g(3) * (-rSges(2,1) * t64 + t63 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,1) * t59 + rSges(3,2) * t61 + pkin(5)) + t81 + (g(2) * t69 + g(3) * t71) * t63 + (g(2) * t71 + g(3) * (-pkin(1) - t69)) * t64) - m(4) * (g(1) * (-t85 * t61 + pkin(5)) + g(2) * (t63 * t61 * pkin(2) + t55 - t64 * qJ(2) + (-t58 * t64 + t61 * t74) * rSges(4,1) + (-t60 * t64 - t61 * t75) * rSges(4,2)) + g(3) * (-pkin(2) * t76 - t64 * pkin(1) - t63 * qJ(2) + (-t60 * t76 - t75) * rSges(4,1) + (t58 * t76 - t74) * rSges(4,2)) + (g(1) * (rSges(4,1) * t60 - rSges(4,2) * t58 + pkin(2)) + t86 * t85) * t59) - m(5) * (g(1) * (-t88 * t61 + pkin(5)) + t81 + (g(3) * t65 + t68 * t80) * t63 + (g(2) * t65 - t68 * t78 - t83) * t64 + (g(1) * t68 + t86 * t88) * t59) - m(6) * (g(1) * (-t87 * t61 + pkin(5)) + t81 + (g(3) * t66 + t67 * t80) * t63 + (g(2) * t66 - t67 * t78 - t83) * t64 + (g(1) * t67 + t86 * t87) * t59);
U = t1;

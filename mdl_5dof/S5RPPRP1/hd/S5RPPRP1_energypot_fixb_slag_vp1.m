% Calculate potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:08
% EndTime: 2022-01-23 09:12:09
% DurationCPUTime: 0.31s
% Computational Cost: add. (147->72), mult. (141->89), div. (0->0), fcn. (125->8), ass. (0->31)
t81 = rSges(5,3) + pkin(6);
t80 = rSges(6,3) + qJ(5) + pkin(6);
t56 = qJ(1) + pkin(7);
t52 = sin(t56);
t53 = cos(t56);
t79 = g(1) * t53 + g(2) * t52;
t57 = sin(pkin(8));
t75 = rSges(4,2) * t57;
t58 = cos(pkin(8));
t74 = t52 * t58;
t60 = sin(qJ(4));
t73 = t52 * t60;
t72 = t53 * t58;
t71 = t53 * t60;
t70 = t58 * t60;
t62 = cos(qJ(4));
t69 = t58 * t62;
t67 = pkin(5) + qJ(2);
t61 = sin(qJ(1));
t54 = t61 * pkin(1);
t66 = t52 * pkin(2) + t54;
t63 = cos(qJ(1));
t55 = t63 * pkin(1);
t65 = t53 * pkin(2) + t52 * qJ(3) + t55;
t64 = -t53 * qJ(3) + t66;
t51 = t62 * pkin(4) + pkin(3);
t47 = t53 * t69 + t73;
t46 = t52 * t62 - t53 * t70;
t45 = t52 * t69 - t71;
t44 = -t52 * t70 - t53 * t62;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t63 * rSges(2,1) - t61 * rSges(2,2)) + g(2) * (t61 * rSges(2,1) + t63 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t53 * rSges(3,1) - t52 * rSges(3,2) + t55) + g(2) * (t52 * rSges(3,1) + t53 * rSges(3,2) + t54) + g(3) * (rSges(3,3) + t67)) - m(4) * (g(1) * (t52 * rSges(4,3) + t65) + g(2) * (rSges(4,1) * t74 - t52 * t75 + t66) + g(3) * (t57 * rSges(4,1) + t58 * rSges(4,2) + t67) + (g(1) * (rSges(4,1) * t58 - t75) + g(2) * (-rSges(4,3) - qJ(3))) * t53) - m(5) * (g(1) * (t47 * rSges(5,1) + t46 * rSges(5,2) + pkin(3) * t72 + t65) + g(2) * (t45 * rSges(5,1) + t44 * rSges(5,2) + pkin(3) * t74 + t64) + g(3) * (-t81 * t58 + t67) + (g(3) * (rSges(5,1) * t62 - rSges(5,2) * t60 + pkin(3)) + t79 * t81) * t57) - m(6) * (g(1) * (t47 * rSges(6,1) + t46 * rSges(6,2) + pkin(4) * t73 + t51 * t72 + t65) + g(2) * (t45 * rSges(6,1) + t44 * rSges(6,2) - pkin(4) * t71 + t51 * t74 + t64) + g(3) * (-t80 * t58 + t67) + (g(3) * (rSges(6,1) * t62 - rSges(6,2) * t60 + t51) + t79 * t80) * t57);
U = t1;

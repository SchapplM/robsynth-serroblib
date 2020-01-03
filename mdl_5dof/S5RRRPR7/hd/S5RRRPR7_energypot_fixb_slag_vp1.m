% Calculate potential energy for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:02
% EndTime: 2019-12-31 21:16:02
% DurationCPUTime: 0.37s
% Computational Cost: add. (151->79), mult. (154->99), div. (0->0), fcn. (138->10), ass. (0->30)
t78 = rSges(6,3) + pkin(8) + qJ(4);
t58 = sin(qJ(1));
t60 = cos(qJ(1));
t77 = g(1) * t60 + g(2) * t58;
t76 = rSges(5,3) + qJ(4);
t73 = rSges(3,3) + pkin(6);
t57 = sin(qJ(2));
t72 = t57 * pkin(2) + pkin(5);
t53 = qJ(2) + qJ(3);
t49 = sin(t53);
t71 = rSges(4,2) * t49;
t54 = sin(pkin(9));
t70 = t54 * t58;
t69 = t54 * t60;
t50 = cos(t53);
t68 = t58 * t50;
t67 = t60 * t50;
t59 = cos(qJ(2));
t45 = pkin(2) * t59 + pkin(1);
t61 = -pkin(7) - pkin(6);
t65 = t58 * t45 + t60 * t61;
t43 = t60 * t45;
t63 = -t58 * t61 + t43;
t62 = rSges(3,1) * t59 - rSges(3,2) * t57 + pkin(1);
t55 = cos(pkin(9));
t52 = pkin(9) + qJ(5);
t48 = cos(t52);
t47 = sin(t52);
t44 = pkin(4) * t55 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t60 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t60) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t57 + rSges(3,2) * t59 + pkin(5)) + (g(1) * t62 - g(2) * t73) * t60 + (g(1) * t73 + g(2) * t62) * t58) - m(4) * (g(1) * (rSges(4,1) * t67 - t60 * t71 + t43) + g(2) * (-rSges(4,3) * t60 + t65) + g(3) * (rSges(4,1) * t49 + rSges(4,2) * t50 + t72) + (g(1) * (rSges(4,3) - t61) + g(2) * (rSges(4,1) * t50 - t71)) * t58) - m(5) * (g(1) * (pkin(3) * t67 + (t55 * t67 + t70) * rSges(5,1) + (-t54 * t67 + t55 * t58) * rSges(5,2) + t63) + g(2) * (pkin(3) * t68 + (t55 * t68 - t69) * rSges(5,1) + (-t54 * t68 - t55 * t60) * rSges(5,2) + t65) + g(3) * (-t76 * t50 + t72) + (g(3) * (rSges(5,1) * t55 - rSges(5,2) * t54 + pkin(3)) + t77 * t76) * t49) - m(6) * (g(1) * (t44 * t67 + pkin(4) * t70 + (t47 * t58 + t48 * t67) * rSges(6,1) + (-t47 * t67 + t48 * t58) * rSges(6,2) + t63) + g(2) * (t44 * t68 - pkin(4) * t69 + (-t47 * t60 + t48 * t68) * rSges(6,1) + (-t47 * t68 - t48 * t60) * rSges(6,2) + t65) + g(3) * (-t78 * t50 + t72) + (g(3) * (rSges(6,1) * t48 - rSges(6,2) * t47 + t44) + t77 * t78) * t49);
U = t1;

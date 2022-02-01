% Calculate potential energy for
% S5RRRPR3
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:23
% EndTime: 2022-01-20 11:42:23
% DurationCPUTime: 0.17s
% Computational Cost: add. (135->59), mult. (97->64), div. (0->0), fcn. (73->10), ass. (0->27)
t63 = pkin(5) + pkin(6);
t62 = rSges(4,3) + pkin(7);
t53 = cos(qJ(3));
t37 = t53 * pkin(3) + pkin(2);
t50 = -qJ(4) - pkin(7);
t61 = rSges(5,3) - t50;
t60 = rSges(6,3) + pkin(8) - t50;
t48 = qJ(3) + pkin(9);
t51 = sin(qJ(3));
t59 = t51 * pkin(3) + t63;
t52 = sin(qJ(1));
t44 = t52 * pkin(1);
t54 = cos(qJ(1));
t46 = t54 * pkin(1);
t58 = g(1) * t46 + g(2) * t44;
t57 = rSges(4,1) * t53 - rSges(4,2) * t51 + pkin(2);
t38 = sin(t48);
t39 = cos(t48);
t56 = rSges(5,1) * t39 - rSges(5,2) * t38 + t37;
t40 = qJ(5) + t48;
t35 = sin(t40);
t36 = cos(t40);
t55 = rSges(6,1) * t36 - rSges(6,2) * t35 + pkin(4) * t39 + t37;
t49 = qJ(1) + qJ(2);
t42 = cos(t49);
t41 = sin(t49);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t54 - t52 * rSges(2,2)) + g(2) * (t52 * rSges(2,1) + rSges(2,2) * t54) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t42 - rSges(3,2) * t41 + t46) + g(2) * (rSges(3,1) * t41 + rSges(3,2) * t42 + t44) + g(3) * (rSges(3,3) + t63)) - m(4) * (g(3) * (rSges(4,1) * t51 + rSges(4,2) * t53 + t63) + (g(1) * t57 - g(2) * t62) * t42 + (g(1) * t62 + g(2) * t57) * t41 + t58) - m(5) * (g(3) * (rSges(5,1) * t38 + rSges(5,2) * t39 + t59) + (g(1) * t56 - g(2) * t61) * t42 + (g(1) * t61 + g(2) * t56) * t41 + t58) - m(6) * (g(3) * (rSges(6,1) * t35 + rSges(6,2) * t36 + pkin(4) * t38 + t59) + (g(1) * t55 - g(2) * t60) * t42 + (g(1) * t60 + g(2) * t55) * t41 + t58);
U = t1;

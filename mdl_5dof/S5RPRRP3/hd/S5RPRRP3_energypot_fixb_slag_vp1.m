% Calculate potential energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:29
% EndTime: 2022-01-23 09:29:30
% DurationCPUTime: 0.21s
% Computational Cost: add. (129->57), mult. (97->62), div. (0->0), fcn. (73->8), ass. (0->25)
t59 = rSges(6,1) + pkin(4);
t49 = -pkin(7) - pkin(6);
t58 = rSges(4,3) + pkin(6);
t47 = cos(qJ(3));
t33 = t47 * pkin(3) + pkin(2);
t57 = rSges(5,3) - t49;
t56 = rSges(6,3) + qJ(5) - t49;
t55 = pkin(5) + qJ(2);
t45 = sin(qJ(3));
t54 = t45 * pkin(3) + t55;
t46 = sin(qJ(1));
t39 = t46 * pkin(1);
t48 = cos(qJ(1));
t41 = t48 * pkin(1);
t53 = g(1) * t41 + g(2) * t39;
t52 = rSges(4,1) * t47 - rSges(4,2) * t45 + pkin(2);
t44 = qJ(3) + qJ(4);
t36 = sin(t44);
t37 = cos(t44);
t51 = rSges(5,1) * t37 - rSges(5,2) * t36 + t33;
t50 = -rSges(6,2) * t36 + t59 * t37 + t33;
t43 = qJ(1) + pkin(8);
t35 = cos(t43);
t34 = sin(t43);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t48 * rSges(2,1) - t46 * rSges(2,2)) + g(2) * (t46 * rSges(2,1) + t48 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t35 * rSges(3,1) - t34 * rSges(3,2) + t41) + g(2) * (t34 * rSges(3,1) + t35 * rSges(3,2) + t39) + g(3) * (rSges(3,3) + t55)) - m(4) * (g(3) * (t45 * rSges(4,1) + t47 * rSges(4,2) + t55) + (g(1) * t52 - g(2) * t58) * t35 + (g(1) * t58 + g(2) * t52) * t34 + t53) - m(5) * (g(3) * (t36 * rSges(5,1) + t37 * rSges(5,2) + t54) + (g(1) * t51 - g(2) * t57) * t35 + (g(1) * t57 + g(2) * t51) * t34 + t53) - m(6) * (g(3) * (t37 * rSges(6,2) + t59 * t36 + t54) + (g(1) * t50 - g(2) * t56) * t35 + (g(1) * t56 + g(2) * t50) * t34 + t53);
U = t1;

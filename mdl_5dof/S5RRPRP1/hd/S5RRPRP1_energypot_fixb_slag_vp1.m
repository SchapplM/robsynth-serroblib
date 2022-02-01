% Calculate potential energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:29
% EndTime: 2022-01-20 10:19:29
% DurationCPUTime: 0.21s
% Computational Cost: add. (130->54), mult. (85->58), div. (0->0), fcn. (61->8), ass. (0->23)
t55 = rSges(6,1) + pkin(4);
t54 = pkin(5) + pkin(6);
t53 = rSges(5,3) + pkin(7);
t52 = rSges(6,3) + qJ(5) + pkin(7);
t40 = qJ(1) + qJ(2);
t36 = sin(t40);
t43 = sin(qJ(1));
t38 = t43 * pkin(1);
t51 = pkin(2) * t36 + t38;
t37 = cos(t40);
t45 = cos(qJ(1));
t39 = t45 * pkin(1);
t50 = pkin(2) * t37 + t39;
t49 = qJ(3) + t54;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t48 = rSges(5,1) * t44 - rSges(5,2) * t42 + pkin(3);
t47 = -rSges(6,2) * t42 + t55 * t44 + pkin(3);
t46 = g(1) * t50 + g(2) * t51;
t35 = pkin(8) + t40;
t31 = cos(t35);
t30 = sin(t35);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t45 * rSges(2,1) - t43 * rSges(2,2)) + g(2) * (t43 * rSges(2,1) + t45 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t37 * rSges(3,1) - t36 * rSges(3,2) + t39) + g(2) * (t36 * rSges(3,1) + t37 * rSges(3,2) + t38) + g(3) * (rSges(3,3) + t54)) - m(4) * (g(1) * (t31 * rSges(4,1) - t30 * rSges(4,2) + t50) + g(2) * (t30 * rSges(4,1) + t31 * rSges(4,2) + t51) + g(3) * (rSges(4,3) + t49)) - m(5) * (g(3) * (t42 * rSges(5,1) + t44 * rSges(5,2) + t49) + (g(1) * t48 - g(2) * t53) * t31 + (g(1) * t53 + g(2) * t48) * t30 + t46) - m(6) * (g(3) * (t44 * rSges(6,2) + t55 * t42 + t49) + (g(1) * t47 - g(2) * t52) * t31 + (g(1) * t52 + g(2) * t47) * t30 + t46);
U = t1;

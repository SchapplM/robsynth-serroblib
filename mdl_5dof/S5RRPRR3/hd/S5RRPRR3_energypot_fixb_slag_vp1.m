% Calculate potential energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:55
% EndTime: 2022-01-20 10:33:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (134->54), mult. (74->56), div. (0->0), fcn. (50->10), ass. (0->25)
t59 = pkin(5) + pkin(6);
t58 = rSges(6,3) + pkin(8);
t46 = qJ(1) + qJ(2);
t42 = sin(t46);
t48 = sin(qJ(1));
t44 = t48 * pkin(1);
t57 = pkin(2) * t42 + t44;
t43 = cos(t46);
t50 = cos(qJ(1));
t45 = t50 * pkin(1);
t56 = pkin(2) * t43 + t45;
t55 = qJ(3) + t59;
t41 = pkin(9) + t46;
t36 = sin(t41);
t54 = pkin(3) * t36 + t57;
t37 = cos(t41);
t53 = pkin(3) * t37 + t56;
t52 = pkin(7) + t55;
t47 = sin(qJ(5));
t49 = cos(qJ(5));
t51 = rSges(6,1) * t49 - rSges(6,2) * t47 + pkin(4);
t40 = qJ(4) + t41;
t35 = cos(t40);
t34 = sin(t40);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t50 - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 + rSges(2,2) * t50) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t43 - rSges(3,2) * t42 + t45) + g(2) * (rSges(3,1) * t42 + rSges(3,2) * t43 + t44) + g(3) * (rSges(3,3) + t59)) - m(4) * (g(1) * (rSges(4,1) * t37 - rSges(4,2) * t36 + t56) + g(2) * (rSges(4,1) * t36 + rSges(4,2) * t37 + t57) + g(3) * (rSges(4,3) + t55)) - m(5) * (g(1) * (rSges(5,1) * t35 - rSges(5,2) * t34 + t53) + g(2) * (rSges(5,1) * t34 + rSges(5,2) * t35 + t54) + g(3) * (rSges(5,3) + t52)) - m(6) * (g(1) * t53 + g(2) * t54 + g(3) * (rSges(6,1) * t47 + rSges(6,2) * t49 + t52) + (g(1) * t51 - g(2) * t58) * t35 + (g(1) * t58 + g(2) * t51) * t34);
U = t1;

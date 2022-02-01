% Calculate potential energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:16
% EndTime: 2022-01-20 09:51:16
% DurationCPUTime: 0.18s
% Computational Cost: add. (136->56), mult. (85->60), div. (0->0), fcn. (61->10), ass. (0->25)
t59 = pkin(5) + pkin(6);
t58 = rSges(6,3) + pkin(7) + qJ(4);
t45 = qJ(1) + qJ(2);
t40 = sin(t45);
t49 = sin(qJ(1));
t42 = t49 * pkin(1);
t57 = pkin(2) * t40 + t42;
t41 = cos(t45);
t50 = cos(qJ(1));
t43 = t50 * pkin(1);
t56 = pkin(2) * t41 + t43;
t55 = rSges(5,3) + qJ(4);
t54 = qJ(3) + t59;
t46 = sin(pkin(9));
t47 = cos(pkin(9));
t53 = rSges(5,1) * t47 - rSges(5,2) * t46 + pkin(3);
t44 = pkin(9) + qJ(5);
t37 = sin(t44);
t38 = cos(t44);
t52 = rSges(6,1) * t38 - rSges(6,2) * t37 + pkin(4) * t47 + pkin(3);
t51 = g(1) * t56 + g(2) * t57;
t39 = pkin(8) + t45;
t33 = cos(t39);
t32 = sin(t39);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t50 - t49 * rSges(2,2)) + g(2) * (t49 * rSges(2,1) + rSges(2,2) * t50) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t41 - rSges(3,2) * t40 + t43) + g(2) * (rSges(3,1) * t40 + rSges(3,2) * t41 + t42) + g(3) * (rSges(3,3) + t59)) - m(4) * (g(1) * (rSges(4,1) * t33 - rSges(4,2) * t32 + t56) + g(2) * (rSges(4,1) * t32 + rSges(4,2) * t33 + t57) + g(3) * (rSges(4,3) + t54)) - m(5) * (g(3) * (rSges(5,1) * t46 + rSges(5,2) * t47 + t54) + (g(1) * t53 - g(2) * t55) * t33 + (g(1) * t55 + g(2) * t53) * t32 + t51) - m(6) * (g(3) * (rSges(6,1) * t37 + rSges(6,2) * t38 + pkin(4) * t46 + t54) + (g(1) * t52 - g(2) * t58) * t33 + (g(1) * t58 + g(2) * t52) * t32 + t51);
U = t1;

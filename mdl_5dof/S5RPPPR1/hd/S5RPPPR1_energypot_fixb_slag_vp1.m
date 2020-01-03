% Calculate potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:01
% EndTime: 2020-01-03 11:20:02
% DurationCPUTime: 0.43s
% Computational Cost: add. (157->72), mult. (141->84), div. (0->0), fcn. (125->10), ass. (0->28)
t48 = sin(pkin(8));
t50 = cos(pkin(8));
t45 = pkin(9) + qJ(5);
t40 = sin(t45);
t42 = cos(t45);
t49 = cos(pkin(9));
t56 = rSges(6,1) * t42 - rSges(6,2) * t40 + pkin(4) * t49 + pkin(3);
t67 = rSges(6,3) + pkin(6) + qJ(4);
t69 = t48 * t67 + t56 * t50;
t66 = rSges(5,3) + qJ(4);
t65 = t50 * pkin(3);
t53 = cos(qJ(1));
t64 = t53 * pkin(1);
t47 = sin(pkin(9));
t63 = t47 * t50;
t62 = t49 * t50;
t61 = pkin(5) + qJ(2);
t46 = qJ(1) + pkin(7);
t41 = sin(t46);
t52 = sin(qJ(1));
t44 = t52 * pkin(1);
t60 = t41 * pkin(2) + t44;
t59 = -rSges(4,3) - qJ(3);
t57 = rSges(4,1) * t50 - rSges(4,2) * t48;
t55 = g(2) * t60 - g(3) * t64;
t54 = -t40 * rSges(6,1) - t42 * rSges(6,2) - t47 * pkin(4) - qJ(3);
t43 = cos(t46);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (t52 * rSges(2,1) + rSges(2,2) * t53) + g(3) * (-rSges(2,1) * t53 + t52 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t61) + g(2) * (rSges(3,1) * t41 + rSges(3,2) * t43 + t44) + g(3) * (-t43 * rSges(3,1) + t41 * rSges(3,2) - t64)) - m(4) * (g(1) * (rSges(4,1) * t48 + rSges(4,2) * t50 + t61) + (g(2) * t57 + g(3) * t59) * t41 + (g(2) * t59 + g(3) * (-pkin(2) - t57)) * t43 + t55) - m(5) * (g(1) * (-t66 * t50 + t61) + g(2) * (t41 * t65 - t43 * qJ(3) + (t41 * t62 - t43 * t47) * rSges(5,1) + (-t41 * t63 - t43 * t49) * rSges(5,2) + t60) + g(3) * (-t64 + (-t47 * rSges(5,1) - t49 * rSges(5,2) - qJ(3)) * t41 + (-t62 * rSges(5,1) + t63 * rSges(5,2) - pkin(2) - t65) * t43) + (g(1) * (rSges(5,1) * t49 - rSges(5,2) * t47 + pkin(3)) + (g(2) * t41 - g(3) * t43) * t66) * t48) - m(6) * ((t69 * g(2) + g(3) * t54) * t41 + (g(2) * t54 + (-pkin(2) - t69) * g(3)) * t43 + t55 + (t56 * t48 - t50 * t67 + t61) * g(1));
U = t1;

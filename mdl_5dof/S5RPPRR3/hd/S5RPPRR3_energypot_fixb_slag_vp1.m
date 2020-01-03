% Calculate potential energy for
% S5RPPRR3
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:23
% EndTime: 2020-01-03 11:27:23
% DurationCPUTime: 0.21s
% Computational Cost: add. (135->59), mult. (97->64), div. (0->0), fcn. (73->10), ass. (0->27)
t52 = cos(qJ(1));
t62 = pkin(1) * t52;
t49 = cos(pkin(9));
t36 = t49 * pkin(3) + pkin(2);
t50 = -pkin(6) - qJ(3);
t61 = -rSges(5,3) + t50;
t60 = -rSges(6,3) - pkin(7) + t50;
t59 = pkin(5) + qJ(2);
t58 = -rSges(4,3) - qJ(3);
t46 = pkin(9) + qJ(4);
t48 = sin(pkin(9));
t57 = t48 * pkin(3) + t59;
t51 = sin(qJ(1));
t44 = t51 * pkin(1);
t56 = g(2) * t44 - g(3) * t62;
t55 = rSges(4,1) * t49 - rSges(4,2) * t48 + pkin(2);
t37 = sin(t46);
t39 = cos(t46);
t54 = rSges(5,1) * t39 - rSges(5,2) * t37 + t36;
t41 = qJ(5) + t46;
t34 = sin(t41);
t35 = cos(t41);
t53 = rSges(6,1) * t35 - rSges(6,2) * t34 + pkin(4) * t39 + t36;
t47 = qJ(1) + pkin(8);
t40 = cos(t47);
t38 = sin(t47);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (t51 * rSges(2,1) + rSges(2,2) * t52) + g(3) * (-rSges(2,1) * t52 + t51 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t59) + g(2) * (rSges(3,1) * t38 + rSges(3,2) * t40 + t44) + g(3) * (-t40 * rSges(3,1) + t38 * rSges(3,2) - t62)) - m(4) * (g(1) * (rSges(4,1) * t48 + rSges(4,2) * t49 + t59) + (g(2) * t58 - g(3) * t55) * t40 + (g(2) * t55 + g(3) * t58) * t38 + t56) - m(5) * (g(1) * (rSges(5,1) * t37 + rSges(5,2) * t39 + t57) + (g(2) * t61 - g(3) * t54) * t40 + (g(2) * t54 + g(3) * t61) * t38 + t56) - m(6) * (g(1) * (rSges(6,1) * t34 + rSges(6,2) * t35 + pkin(4) * t37 + t57) + (g(2) * t60 - g(3) * t53) * t40 + (g(2) * t53 + g(3) * t60) * t38 + t56);
U = t1;

% Calculate Gravitation load on the joints for
% S5RRPRR5
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
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:31
% EndTime: 2020-01-03 12:03:32
% DurationCPUTime: 0.29s
% Computational Cost: add. (276->70), mult. (199->89), div. (0->0), fcn. (157->10), ass. (0->41)
t44 = -pkin(7) - qJ(3);
t72 = rSges(5,3) - t44;
t71 = rSges(6,3) + pkin(8) - t44;
t40 = pkin(9) + qJ(4);
t33 = qJ(5) + t40;
t26 = sin(t33);
t27 = cos(t33);
t55 = t27 * rSges(6,1) - rSges(6,2) * t26;
t31 = sin(t40);
t32 = cos(t40);
t64 = rSges(5,1) * t32;
t70 = -rSges(5,2) * t31 + t64;
t69 = rSges(4,3) + qJ(3);
t43 = cos(pkin(9));
t68 = -rSges(4,2) * sin(pkin(9)) + rSges(4,1) * t43 + pkin(2);
t41 = qJ(1) + qJ(2);
t35 = cos(t41);
t59 = t27 * t35;
t60 = t26 * t35;
t67 = rSges(6,1) * t60 + rSges(6,2) * t59;
t34 = sin(t41);
t66 = g(2) * t34;
t28 = t43 * pkin(3) + pkin(2);
t58 = t31 * t35;
t57 = t34 * rSges(3,1) + t35 * rSges(3,2);
t56 = t35 * rSges(3,1) - t34 * rSges(3,2);
t54 = rSges(5,1) * t31 + rSges(5,2) * t32;
t53 = -rSges(6,1) * t26 - rSges(6,2) * t27;
t52 = t69 * t34 + t68 * t35;
t19 = pkin(4) * t32;
t4 = t19 + t28;
t51 = -t71 * t35 + (t4 + t55) * t34;
t50 = rSges(6,1) * t59 - rSges(6,2) * t60 + t71 * t34 + t35 * t4;
t49 = -t72 * t35 + (t28 + t70) * t34;
t48 = -rSges(5,2) * t58 + (t28 + t64) * t35 + t72 * t34;
t47 = t68 * t34 - t69 * t35;
t46 = cos(qJ(1));
t45 = sin(qJ(1));
t38 = t46 * pkin(1);
t37 = t45 * pkin(1);
t1 = [-m(2) * (g(2) * (t46 * rSges(2,1) - t45 * rSges(2,2)) + g(3) * (t45 * rSges(2,1) + t46 * rSges(2,2))) - m(3) * (g(2) * (t38 + t56) + g(3) * (t37 + t57)) - m(4) * (g(2) * (t38 + t52) + g(3) * (t37 + t47)) - m(5) * (g(2) * (t38 + t48) + g(3) * (t37 + t49)) - m(6) * (g(2) * (t38 + t50) + g(3) * (t37 + t51)), -m(3) * (g(2) * t56 + g(3) * t57) - m(4) * (g(2) * t52 + g(3) * t47) - m(5) * (g(2) * t48 + g(3) * t49) - m(6) * (g(2) * t50 + g(3) * t51), (-m(4) - m(5) - m(6)) * (-g(2) * t35 - g(3) * t34), -m(5) * (g(3) * t54 * t35 + g(1) * t70) - m(6) * (g(1) * (t19 + t55) + g(3) * (pkin(4) * t58 + t67)) + (m(5) * t54 - m(6) * (-pkin(4) * t31 + t53)) * t66, -m(6) * (g(1) * t55 + g(3) * t67 + t53 * t66)];
taug = t1(:);

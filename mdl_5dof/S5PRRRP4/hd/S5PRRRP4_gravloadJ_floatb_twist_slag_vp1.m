% Calculate Gravitation load on the joints for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:37
% EndTime: 2019-12-05 16:45:39
% DurationCPUTime: 0.47s
% Computational Cost: add. (240->82), mult. (313->114), div. (0->0), fcn. (285->8), ass. (0->45)
t66 = rSges(6,1) + pkin(4);
t33 = cos(qJ(4));
t65 = -t33 * t66 - pkin(3);
t29 = sin(pkin(8));
t28 = qJ(2) + qJ(3);
t25 = sin(t28);
t31 = sin(qJ(4));
t56 = t25 * t31;
t44 = rSges(5,2) * t56;
t26 = cos(t28);
t55 = t26 * t29;
t64 = rSges(5,3) * t55 + t29 * t44;
t30 = cos(pkin(8));
t54 = t26 * t30;
t63 = rSges(5,3) * t54 + t30 * t44;
t62 = g(1) * t30 + g(2) * t29;
t45 = rSges(6,3) + qJ(5);
t61 = t62 * t25;
t32 = sin(qJ(2));
t60 = pkin(2) * t32;
t53 = t26 * t31;
t52 = t26 * t33;
t51 = t29 * t31;
t50 = t29 * t33;
t49 = t30 * t31;
t48 = t30 * t33;
t47 = t26 * pkin(3) + t25 * pkin(7);
t46 = qJ(5) * t31;
t43 = -rSges(5,1) * t33 - pkin(3);
t11 = pkin(7) * t55;
t42 = -t29 * t60 + t11;
t13 = pkin(7) * t54;
t41 = -t30 * t60 + t13;
t39 = t25 * rSges(6,2) + rSges(6,3) * t53 + t26 * t46 + t66 * t52 + t47;
t38 = rSges(5,1) * t52 - rSges(5,2) * t53 + t25 * rSges(5,3) + t47;
t34 = cos(qJ(2));
t27 = t34 * pkin(2);
t21 = t26 * rSges(4,1);
t10 = rSges(6,2) * t54;
t8 = rSges(6,2) * t55;
t4 = t26 * t48 + t51;
t3 = t26 * t49 - t50;
t2 = t26 * t50 - t49;
t1 = t26 * t51 + t48;
t5 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(3) * (t34 * rSges(3,1) - t32 * rSges(3,2)) + t62 * (-rSges(3,1) * t32 - rSges(3,2) * t34)) - m(4) * (g(3) * (-t25 * rSges(4,2) + t21 + t27) + t62 * (-rSges(4,1) * t25 - rSges(4,2) * t26 - t60)) - m(5) * (g(1) * (t41 + t63) + g(2) * (t42 + t64) + g(3) * (t27 + t38) + t43 * t61) - m(6) * (g(1) * (t10 + t41) + g(2) * (t42 + t8) + g(3) * (t27 + t39) + (-t45 * t31 + t65) * t61), -m(4) * (g(3) * t21 + (-g(1) * t54 - g(2) * t55) * rSges(4,2)) - m(5) * (g(1) * (t13 + t63) + g(2) * (t11 + t64) + g(3) * t38) - m(6) * (g(1) * (t10 + t13) + g(2) * (t11 + t8) + g(3) * t39) + (m(4) * g(3) * rSges(4,2) + t62 * (m(4) * rSges(4,1) - m(5) * t43 - m(6) * (-rSges(6,3) * t31 - t46 + t65))) * t25, -m(5) * (g(1) * (-t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) - t2 * rSges(5,2))) - m(6) * (g(1) * (-t3 * t66 + t45 * t4) + g(2) * (-t1 * t66 + t45 * t2)) + (-m(5) * (-rSges(5,1) * t31 - rSges(5,2) * t33) - m(6) * (-t31 * t66 + t33 * t45)) * g(3) * t25, -m(6) * (g(1) * t3 + g(2) * t1 + g(3) * t56)];
taug = t5(:);

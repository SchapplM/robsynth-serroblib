% Calculate Gravitation load on the joints for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:00
% EndTime: 2019-12-05 16:12:02
% DurationCPUTime: 0.54s
% Computational Cost: add. (176->93), mult. (437->132), div. (0->0), fcn. (456->8), ass. (0->46)
t24 = sin(pkin(8));
t26 = cos(pkin(8));
t42 = rSges(6,3) + qJ(5);
t53 = rSges(6,1) + pkin(4);
t32 = -t42 * t24 - t53 * t26;
t25 = sin(pkin(7));
t27 = cos(pkin(7));
t60 = g(1) * t27 + g(2) * t25;
t29 = sin(qJ(2));
t59 = t60 * t29;
t58 = -m(5) - m(6);
t30 = cos(qJ(3));
t48 = t29 * t30;
t55 = g(3) * qJ(4) * t48;
t54 = g(3) * t29;
t31 = cos(qJ(2));
t52 = t25 * t31;
t51 = t27 * t31;
t28 = sin(qJ(3));
t50 = t28 * t31;
t49 = t29 * t26;
t47 = t30 * t31;
t46 = t31 * t26;
t45 = t31 * pkin(2) + t29 * pkin(6);
t44 = rSges(6,2) + qJ(4);
t43 = rSges(5,3) + qJ(4);
t41 = -pkin(3) * t30 - pkin(2);
t40 = pkin(3) * t47 + qJ(4) * t50 + t45;
t38 = -rSges(5,1) * t26 + rSges(5,2) * t24;
t36 = t24 * t31 - t26 * t48;
t35 = t24 * t48 + t46;
t16 = pkin(6) * t51;
t14 = pkin(6) * t52;
t13 = t29 * t24 + t30 * t46;
t12 = t24 * t47 - t49;
t11 = t25 * t28 + t27 * t47;
t10 = -t25 * t30 + t27 * t50;
t9 = t25 * t47 - t27 * t28;
t8 = t25 * t50 + t27 * t30;
t7 = t10 * pkin(3);
t6 = t8 * pkin(3);
t5 = t36 * t27;
t4 = t35 * t27;
t3 = t36 * t25;
t2 = t35 * t25;
t1 = [(-m(2) - m(3) - m(4) + t58) * g(3), -m(3) * (g(3) * (t31 * rSges(3,1) - t29 * rSges(3,2)) + t60 * (-rSges(3,1) * t29 - rSges(3,2) * t31)) - m(4) * (g(1) * (rSges(4,3) * t51 + t16) + g(2) * (rSges(4,3) * t52 + t14) + g(3) * (rSges(4,1) * t47 - rSges(4,2) * t50 + t45) + (g(3) * rSges(4,3) + t60 * (-rSges(4,1) * t30 + rSges(4,2) * t28 - pkin(2))) * t29) - m(5) * (g(1) * (t5 * rSges(5,1) + t4 * rSges(5,2) + t16) + g(2) * (t3 * rSges(5,1) + t2 * rSges(5,2) + t14) + g(3) * (t13 * rSges(5,1) - t12 * rSges(5,2) + rSges(5,3) * t50 + t40) + (-t43 * t28 + t41) * t59) - m(6) * (g(1) * (-t42 * t4 + t53 * t5 + t16) + g(2) * (-t42 * t2 + t53 * t3 + t14) + g(3) * (rSges(6,2) * t50 + t42 * t12 + t53 * t13 + t40) + (-t44 * t28 + t41) * t59), -m(4) * (g(1) * (-t10 * rSges(4,1) - t11 * rSges(4,2)) + g(2) * (-t8 * rSges(4,1) - t9 * rSges(4,2))) - m(5) * (g(1) * (t38 * t10 + t43 * t11 - t7) + g(2) * (t38 * t8 + t43 * t9 - t6) + t55) - m(6) * (g(1) * (t32 * t10 + t44 * t11 - t7) + g(2) * (t32 * t8 + t44 * t9 - t6) + t55) + ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t30 + (m(4) * rSges(4,1) - m(5) * (-pkin(3) + t38) - m(6) * (-pkin(3) + t32)) * t28) * t54, t58 * (g(1) * t10 + g(2) * t8 + t28 * t54), -m(6) * (g(1) * (t11 * t24 - t27 * t49) + g(2) * (t9 * t24 - t25 * t49) + g(3) * t35)];
taug = t1(:);

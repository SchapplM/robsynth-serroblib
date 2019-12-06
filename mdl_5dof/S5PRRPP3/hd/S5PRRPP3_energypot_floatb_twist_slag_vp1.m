% Calculate potential energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:11:54
% EndTime: 2019-12-05 16:11:55
% DurationCPUTime: 0.28s
% Computational Cost: add. (163->91), mult. (279->108), div. (0->0), fcn. (301->8), ass. (0->39)
t53 = rSges(6,1) + pkin(4);
t27 = sin(pkin(7));
t31 = sin(qJ(2));
t52 = t27 * t31;
t33 = cos(qJ(2));
t51 = t27 * t33;
t29 = cos(pkin(7));
t50 = t29 * t31;
t30 = sin(qJ(3));
t49 = t30 * t31;
t48 = t30 * t33;
t32 = cos(qJ(3));
t47 = t31 * t32;
t46 = t32 * t33;
t45 = rSges(6,2) + qJ(4);
t44 = rSges(5,3) + qJ(4);
t43 = rSges(6,3) + qJ(5);
t42 = t27 * pkin(1) + r_base(2);
t41 = qJ(1) + r_base(3);
t40 = t29 * pkin(1) + t27 * pkin(5) + r_base(1);
t39 = t31 * pkin(2) + t41;
t38 = t29 * t33 * pkin(2) + pkin(6) * t50 + t40;
t12 = t27 * t30 + t29 * t46;
t37 = t12 * pkin(3) + t38;
t36 = pkin(2) * t51 - t29 * pkin(5) + pkin(6) * t52 + t42;
t8 = t27 * t46 - t29 * t30;
t35 = t8 * pkin(3) + t36;
t34 = pkin(3) * t47 - t33 * pkin(6) + qJ(4) * t49 + t39;
t28 = cos(pkin(8));
t26 = sin(pkin(8));
t11 = -t27 * t32 + t29 * t48;
t10 = -t33 * t26 + t28 * t47;
t9 = t26 * t47 + t33 * t28;
t7 = t27 * t48 + t29 * t32;
t4 = t12 * t28 + t26 * t50;
t3 = t12 * t26 - t28 * t50;
t2 = t26 * t52 + t8 * t28;
t1 = t8 * t26 - t28 * t52;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t27 * rSges(2,2) + r_base(1)) + g(2) * (t27 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t27 * rSges(3,3) + t40) + g(2) * (rSges(3,1) * t51 - rSges(3,2) * t52 + t42) + g(3) * (t31 * rSges(3,1) + t33 * rSges(3,2) + t41) + (g(1) * (rSges(3,1) * t33 - rSges(3,2) * t31) + g(2) * (-rSges(3,3) - pkin(5))) * t29) - m(4) * (g(1) * (t12 * rSges(4,1) - t11 * rSges(4,2) + rSges(4,3) * t50 + t38) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + rSges(4,3) * t52 + t36) + g(3) * ((-rSges(4,3) - pkin(6)) * t33 + (rSges(4,1) * t32 - rSges(4,2) * t30) * t31 + t39)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t44 * t11 + t37) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t44 * t7 + t35) + g(3) * (t10 * rSges(5,1) - t9 * rSges(5,2) + rSges(5,3) * t49 + t34)) - m(6) * (g(1) * (t45 * t11 + t43 * t3 + t53 * t4 + t37) + g(2) * (t43 * t1 + t53 * t2 + t45 * t7 + t35) + g(3) * (rSges(6,2) * t49 + t53 * t10 + t43 * t9 + t34));
U = t5;

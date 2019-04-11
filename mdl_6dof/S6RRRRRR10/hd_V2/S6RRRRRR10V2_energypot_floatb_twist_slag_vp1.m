% Calculate potential energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10V2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10V2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:31
% EndTime: 2019-04-11 14:41:31
% DurationCPUTime: 0.30s
% Computational Cost: add. (244->100), mult. (277->125), div. (0->0), fcn. (295->12), ass. (0->45)
t24 = qJ(2) + qJ(3);
t22 = cos(t24);
t54 = pkin(3) * t22;
t53 = pkin(6) + rSges(7,3);
t21 = sin(t24);
t27 = sin(qJ(4));
t52 = t21 * t27;
t29 = sin(qJ(1));
t51 = t21 * t29;
t32 = cos(qJ(4));
t50 = t21 * t32;
t34 = cos(qJ(1));
t49 = t21 * t34;
t48 = t29 * t27;
t47 = t29 * t32;
t46 = t34 * t27;
t45 = t34 * t32;
t44 = pkin(4) + r_base(3);
t33 = cos(qJ(2));
t20 = t33 * pkin(2) + pkin(1);
t43 = t29 * t20 + r_base(2);
t42 = t34 * t20 + r_base(1);
t28 = sin(qJ(2));
t41 = t28 * pkin(2) + t44;
t40 = t21 * pkin(3) + t41;
t39 = pkin(5) * t51 + t29 * t54 + t43;
t38 = pkin(5) * t49 + t34 * t54 + t42;
t37 = rSges(4,1) * t22 - rSges(4,2) * t21;
t36 = rSges(3,1) * t33 - rSges(3,2) * t28 + pkin(1);
t35 = -t22 * pkin(5) + t40;
t31 = cos(qJ(5));
t30 = cos(qJ(6));
t26 = sin(qJ(5));
t25 = sin(qJ(6));
t10 = t22 * t45 + t48;
t9 = t22 * t46 - t47;
t8 = t22 * t47 - t46;
t7 = t22 * t48 + t45;
t6 = -t22 * t26 + t31 * t50;
t5 = t22 * t31 + t26 * t50;
t4 = t10 * t31 + t26 * t49;
t3 = t10 * t26 - t31 * t49;
t2 = t26 * t51 + t8 * t31;
t1 = t8 * t26 - t31 * t51;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t29 * rSges(2,2) + r_base(1)) + g(2) * (t29 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t44)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t28 * rSges(3,1) + t33 * rSges(3,2) + t44) + (-g(2) * rSges(3,3) + g(1) * t36) * t34 + (g(1) * rSges(3,3) + g(2) * t36) * t29) - m(4) * (g(1) * (t29 * rSges(4,3) + t37 * t34 + t42) + g(2) * (-t34 * rSges(4,3) + t37 * t29 + t43) + g(3) * (t21 * rSges(4,1) + t22 * rSges(4,2) + t41)) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + rSges(5,3) * t49 + t38) + g(2) * (t8 * rSges(5,1) - t7 * rSges(5,2) + rSges(5,3) * t51 + t39) + g(3) * ((-rSges(5,3) - pkin(5)) * t22 + (rSges(5,1) * t32 - rSges(5,2) * t27) * t21 + t40)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t9 * rSges(6,3) + t38) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t7 * rSges(6,3) + t39) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + rSges(6,3) * t52 + t35)) - m(7) * (g(1) * ((t9 * t25 + t4 * t30) * rSges(7,1) + (-t4 * t25 + t9 * t30) * rSges(7,2) + t53 * t3 + t38) + g(2) * ((t2 * t30 + t7 * t25) * rSges(7,1) + (-t2 * t25 + t7 * t30) * rSges(7,2) + t53 * t1 + t39) + g(3) * ((t25 * t52 + t6 * t30) * rSges(7,1) + (-t6 * t25 + t30 * t52) * rSges(7,2) + t53 * t5 + t35));
U  = t11;

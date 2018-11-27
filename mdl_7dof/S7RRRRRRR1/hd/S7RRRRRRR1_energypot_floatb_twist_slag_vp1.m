% Calculate potential energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% rSges [8x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S7RRRRRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [7x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp1: m has to be [8x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [8,3]), ...
  'S7RRRRRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [8x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 19:14:15
% EndTime: 2018-11-26 19:14:16
% DurationCPUTime: 0.47s
% Computational Cost: add. (318->121), mult. (689->160), div. (0->0), fcn. (857->14), ass. (0->54)
t56 = rSges(6,3) + pkin(3);
t55 = pkin(4) + rSges(8,3);
t54 = cos(qJ(4));
t30 = sin(qJ(3));
t31 = sin(qJ(2));
t53 = t30 * t31;
t32 = sin(qJ(1));
t52 = t32 * t31;
t37 = cos(qJ(2));
t51 = t32 * t37;
t38 = cos(qJ(1));
t50 = t38 * t30;
t49 = t38 * t31;
t36 = cos(qJ(3));
t48 = t38 * t36;
t47 = pkin(1) + r_base(3);
t46 = t31 * t54;
t45 = t37 * pkin(2) + t47;
t44 = rSges(3,1) * t37 - rSges(3,2) * t31;
t43 = -pkin(2) * t52 + r_base(2);
t42 = -pkin(2) * t49 + r_base(1);
t29 = sin(qJ(4));
t17 = t31 * t36 * t29 + t37 * t54;
t41 = t17 * pkin(3) + t45;
t20 = t36 * t51 + t50;
t13 = t20 * t29 - t32 * t46;
t40 = t13 * pkin(3) + t43;
t22 = -t32 * t30 + t37 * t48;
t15 = t22 * t29 - t38 * t46;
t39 = t15 * pkin(3) + t42;
t35 = cos(qJ(5));
t34 = cos(qJ(6));
t33 = cos(qJ(7));
t28 = sin(qJ(5));
t27 = sin(qJ(6));
t26 = sin(qJ(7));
t21 = -t32 * t36 - t37 * t50;
t19 = -t30 * t51 + t48;
t18 = -t37 * t29 + t36 * t46;
t16 = t22 * t54 + t29 * t49;
t14 = t20 * t54 + t29 * t52;
t12 = t18 * t35 - t28 * t53;
t11 = t18 * t28 + t35 * t53;
t10 = t16 * t35 + t21 * t28;
t9 = t16 * t28 - t21 * t35;
t8 = t14 * t35 + t19 * t28;
t7 = t14 * t28 - t19 * t35;
t6 = t12 * t34 + t17 * t27;
t5 = -t12 * t27 + t17 * t34;
t4 = t10 * t34 + t15 * t27;
t3 = -t10 * t27 + t15 * t34;
t2 = t13 * t27 + t8 * t34;
t1 = t13 * t34 - t8 * t27;
t23 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t38 * rSges(2,1) - t32 * rSges(2,2) + r_base(1)) + g(2) * (t32 * rSges(2,1) + t38 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t47)) - m(3) * (g(1) * (t32 * rSges(3,3) + t44 * t38 + r_base(1)) + g(2) * (-t38 * rSges(3,3) + t44 * t32 + r_base(2)) + g(3) * (t31 * rSges(3,1) + t37 * rSges(3,2) + t47)) - m(4) * (g(1) * (t22 * rSges(4,1) + t21 * rSges(4,2) + r_base(1)) + g(2) * (t20 * rSges(4,1) + t19 * rSges(4,2) + r_base(2)) + g(3) * (t37 * rSges(4,3) + t45) + (g(3) * (rSges(4,1) * t36 - rSges(4,2) * t30) + (g(1) * t38 + g(2) * t32) * (-rSges(4,3) - pkin(2))) * t31) - m(5) * (g(1) * (t16 * rSges(5,1) - t15 * rSges(5,2) + t21 * rSges(5,3) + t42) + g(2) * (t14 * rSges(5,1) - t13 * rSges(5,2) + t19 * rSges(5,3) + t43) + g(3) * (t18 * rSges(5,1) - t17 * rSges(5,2) - rSges(5,3) * t53 + t45)) - m(6) * (g(1) * (t10 * rSges(6,1) - t9 * rSges(6,2) + t56 * t15 + t42) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t56 * t13 + t43) + g(3) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t56 * t17 + t45)) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t9 * rSges(7,3) + t39) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t7 * rSges(7,3) + t40) + g(3) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t11 * rSges(7,3) + t41)) - m(8) * (g(1) * ((-t9 * t26 + t4 * t33) * rSges(8,1) + (-t4 * t26 - t9 * t33) * rSges(8,2) + t55 * t3 + t39) + g(2) * ((t2 * t33 - t7 * t26) * rSges(8,1) + (-t2 * t26 - t7 * t33) * rSges(8,2) + t55 * t1 + t40) + g(3) * ((-t11 * t26 + t6 * t33) * rSges(8,1) + (-t11 * t33 - t6 * t26) * rSges(8,2) + t55 * t5 + t41));
U  = t23;

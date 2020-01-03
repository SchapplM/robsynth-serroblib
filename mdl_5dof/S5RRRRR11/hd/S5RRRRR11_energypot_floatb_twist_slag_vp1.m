% Calculate potential energy for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:42
% EndTime: 2019-12-31 22:38:43
% DurationCPUTime: 0.34s
% Computational Cost: add. (237->103), mult. (435->125), div. (0->0), fcn. (511->12), ass. (0->45)
t56 = rSges(4,3) + pkin(8);
t55 = pkin(9) + rSges(5,3);
t25 = sin(pkin(5));
t29 = sin(qJ(2));
t54 = t25 * t29;
t30 = sin(qJ(1));
t53 = t25 * t30;
t32 = cos(qJ(3));
t52 = t25 * t32;
t33 = cos(qJ(2));
t51 = t25 * t33;
t34 = cos(qJ(1));
t50 = t25 * t34;
t49 = t30 * t29;
t48 = t30 * t33;
t47 = t34 * t29;
t46 = t34 * t33;
t45 = pkin(10) + pkin(9) + rSges(6,3);
t44 = pkin(6) + r_base(3);
t43 = t30 * pkin(1) + r_base(2);
t26 = cos(pkin(5));
t42 = t26 * pkin(7) + t44;
t41 = t34 * pkin(1) + pkin(7) * t53 + r_base(1);
t40 = pkin(2) * t54 + t42;
t12 = -t26 * t49 + t46;
t39 = t12 * pkin(2) + t41;
t24 = qJ(4) + qJ(5);
t19 = sin(t24);
t20 = cos(t24);
t31 = cos(qJ(4));
t38 = t20 * rSges(6,1) - t19 * rSges(6,2) + t31 * pkin(4) + pkin(3);
t10 = t26 * t47 + t48;
t37 = t10 * pkin(2) - pkin(7) * t50 + t43;
t27 = sin(qJ(4));
t36 = t19 * rSges(6,1) + t20 * rSges(6,2) + t27 * pkin(4) + pkin(8);
t28 = sin(qJ(3));
t11 = t26 * t48 + t47;
t9 = -t26 * t46 + t49;
t8 = t26 * t28 + t29 * t52;
t7 = -t26 * t32 + t28 * t54;
t4 = t12 * t32 + t28 * t53;
t3 = t12 * t28 - t30 * t52;
t2 = t10 * t32 - t28 * t50;
t1 = t10 * t28 + t32 * t50;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t34 * rSges(2,1) - t30 * rSges(2,2) + r_base(1)) + g(2) * (t30 * rSges(2,1) + t34 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t44)) - m(3) * (g(1) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t41) + g(2) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t43) + g(3) * (t26 * rSges(3,3) + t42) + (g(1) * rSges(3,3) * t30 + g(3) * (rSges(3,1) * t29 + rSges(3,2) * t33) + g(2) * (-rSges(3,3) - pkin(7)) * t34) * t25) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t56 * t11 + t39) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t56 * t9 + t37) + g(3) * (t8 * rSges(4,1) - t7 * rSges(4,2) - t56 * t51 + t40)) - m(5) * (g(1) * (t4 * pkin(3) + t11 * pkin(8) + (t11 * t27 + t4 * t31) * rSges(5,1) + (t11 * t31 - t4 * t27) * rSges(5,2) + t55 * t3 + t39) + g(2) * (t2 * pkin(3) + t9 * pkin(8) + (t2 * t31 + t9 * t27) * rSges(5,1) + (-t2 * t27 + t9 * t31) * rSges(5,2) + t55 * t1 + t37) + g(3) * (t8 * pkin(3) - pkin(8) * t51 + (-t27 * t51 + t8 * t31) * rSges(5,1) + (-t8 * t27 - t31 * t51) * rSges(5,2) + t55 * t7 + t40)) - m(6) * (g(1) * (t36 * t11 + t45 * t3 + t38 * t4 + t39) + g(2) * (t45 * t1 + t38 * t2 + t36 * t9 + t37) + g(3) * (-t36 * t51 + t38 * t8 + t45 * t7 + t40));
U = t5;

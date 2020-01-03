% Calculate potential energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP10_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:03
% EndTime: 2019-12-31 22:08:04
% DurationCPUTime: 0.32s
% Computational Cost: add. (225->104), mult. (435->127), div. (0->0), fcn. (511->10), ass. (0->48)
t58 = rSges(4,3) + pkin(8);
t57 = rSges(5,3) + pkin(9);
t28 = sin(pkin(5));
t33 = sin(qJ(2));
t56 = t28 * t33;
t34 = sin(qJ(1));
t55 = t28 * t34;
t36 = cos(qJ(3));
t54 = t28 * t36;
t37 = cos(qJ(2));
t53 = t28 * t37;
t38 = cos(qJ(1));
t52 = t28 * t38;
t51 = t34 * t33;
t50 = t34 * t37;
t49 = t38 * t33;
t48 = t38 * t37;
t47 = rSges(6,3) + qJ(5) + pkin(9);
t46 = pkin(6) + r_base(3);
t45 = t34 * pkin(1) + r_base(2);
t31 = sin(qJ(4));
t44 = pkin(4) * t31 + pkin(8);
t29 = cos(pkin(5));
t43 = t29 * pkin(7) + t46;
t42 = t38 * pkin(1) + pkin(7) * t55 + r_base(1);
t41 = pkin(2) * t56 + t43;
t18 = -t29 * t51 + t48;
t40 = t18 * pkin(2) + t42;
t16 = t29 * t49 + t50;
t39 = t16 * pkin(2) - pkin(7) * t52 + t45;
t35 = cos(qJ(4));
t32 = sin(qJ(3));
t24 = t35 * pkin(4) + pkin(3);
t17 = t29 * t50 + t49;
t15 = -t29 * t48 + t51;
t14 = t29 * t32 + t33 * t54;
t13 = -t29 * t36 + t32 * t56;
t10 = t18 * t36 + t32 * t55;
t9 = t18 * t32 - t34 * t54;
t8 = t16 * t36 - t32 * t52;
t7 = t16 * t32 + t36 * t52;
t6 = t14 * t35 - t31 * t53;
t5 = -t14 * t31 - t35 * t53;
t4 = t10 * t35 + t17 * t31;
t3 = -t10 * t31 + t17 * t35;
t2 = t15 * t31 + t8 * t35;
t1 = t15 * t35 - t8 * t31;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t38 * rSges(2,1) - t34 * rSges(2,2) + r_base(1)) + g(2) * (t34 * rSges(2,1) + t38 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t46)) - m(3) * (g(1) * (t18 * rSges(3,1) - t17 * rSges(3,2) + t42) + g(2) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t45) + g(3) * (t29 * rSges(3,3) + t43) + (g(1) * rSges(3,3) * t34 + g(3) * (rSges(3,1) * t33 + rSges(3,2) * t37) + g(2) * (-rSges(3,3) - pkin(7)) * t38) * t28) - m(4) * (g(1) * (t10 * rSges(4,1) - t9 * rSges(4,2) + t58 * t17 + t40) + g(2) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t58 * t15 + t39) + g(3) * (t14 * rSges(4,1) - t13 * rSges(4,2) - t58 * t53 + t41)) - m(5) * (g(1) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t10 * pkin(3) + t17 * pkin(8) + t57 * t9 + t40) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t8 * pkin(3) + t15 * pkin(8) + t57 * t7 + t39) + g(3) * (t6 * rSges(5,1) + t5 * rSges(5,2) + t14 * pkin(3) - pkin(8) * t53 + t57 * t13 + t41)) - m(6) * (g(1) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t10 * t24 + t44 * t17 + t47 * t9 + t40) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t44 * t15 + t8 * t24 + t47 * t7 + t39) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t47 * t13 + t14 * t24 - t44 * t53 + t41));
U = t11;

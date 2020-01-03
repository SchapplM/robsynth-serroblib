% Calculate potential energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:18
% EndTime: 2019-12-31 17:29:18
% DurationCPUTime: 0.27s
% Computational Cost: add. (145->83), mult. (272->106), div. (0->0), fcn. (311->10), ass. (0->39)
t48 = rSges(4,3) + pkin(7);
t47 = pkin(8) + rSges(5,3);
t21 = sin(pkin(4));
t25 = sin(qJ(2));
t46 = t21 * t25;
t26 = sin(qJ(1));
t45 = t21 * t26;
t28 = cos(qJ(3));
t44 = t21 * t28;
t29 = cos(qJ(2));
t43 = t21 * t29;
t30 = cos(qJ(1));
t42 = t21 * t30;
t41 = t26 * t25;
t40 = t26 * t29;
t39 = t30 * t25;
t38 = t30 * t29;
t37 = pkin(5) + r_base(3);
t36 = t26 * pkin(1) + r_base(2);
t22 = cos(pkin(4));
t35 = t22 * pkin(6) + t37;
t34 = t30 * pkin(1) + pkin(6) * t45 + r_base(1);
t33 = pkin(2) * t46 + t35;
t12 = -t22 * t41 + t38;
t32 = t12 * pkin(2) + t34;
t10 = t22 * t39 + t40;
t31 = t10 * pkin(2) - pkin(6) * t42 + t36;
t27 = cos(qJ(4));
t24 = sin(qJ(3));
t23 = sin(qJ(4));
t11 = t22 * t40 + t39;
t9 = -t22 * t38 + t41;
t8 = t22 * t24 + t25 * t44;
t7 = -t22 * t28 + t24 * t46;
t4 = t12 * t28 + t24 * t45;
t3 = t12 * t24 - t26 * t44;
t2 = t10 * t28 - t24 * t42;
t1 = t10 * t24 + t28 * t42;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t30 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t30 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(1) * (t12 * rSges(3,1) - t11 * rSges(3,2) + t34) + g(2) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t36) + g(3) * (t22 * rSges(3,3) + t35) + (g(1) * rSges(3,3) * t26 + g(3) * (rSges(3,1) * t25 + rSges(3,2) * t29) + g(2) * (-rSges(3,3) - pkin(6)) * t30) * t21) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t48 * t11 + t32) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t48 * t9 + t31) + g(3) * (t8 * rSges(4,1) - t7 * rSges(4,2) - t48 * t43 + t33)) - m(5) * (g(1) * (t4 * pkin(3) + t11 * pkin(7) + (t11 * t23 + t4 * t27) * rSges(5,1) + (t11 * t27 - t4 * t23) * rSges(5,2) + t47 * t3 + t32) + g(2) * (t2 * pkin(3) + t9 * pkin(7) + (t2 * t27 + t9 * t23) * rSges(5,1) + (-t2 * t23 + t9 * t27) * rSges(5,2) + t47 * t1 + t31) + g(3) * (t8 * pkin(3) - pkin(7) * t43 + (-t23 * t43 + t8 * t27) * rSges(5,1) + (-t8 * t23 - t27 * t43) * rSges(5,2) + t47 * t7 + t33));
U = t5;

% Calculate potential energy for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:58:59
% EndTime: 2022-01-23 08:58:59
% DurationCPUTime: 0.38s
% Computational Cost: add. (176->116), mult. (263->148), div. (0->0), fcn. (280->10), ass. (0->48)
t24 = sin(pkin(7));
t50 = t24 * qJ(3) + pkin(1);
t23 = sin(pkin(8));
t49 = t23 * qJ(4) + pkin(2);
t28 = sin(qJ(5));
t48 = t23 * t28;
t30 = cos(qJ(5));
t47 = t23 * t30;
t26 = cos(pkin(8));
t46 = t24 * t26;
t29 = sin(qJ(1));
t45 = t24 * t29;
t31 = cos(qJ(1));
t44 = t24 * t31;
t25 = cos(pkin(9));
t27 = cos(pkin(7));
t43 = t27 * t25;
t42 = t29 * t23;
t41 = t29 * t26;
t40 = t31 * t23;
t39 = t31 * t26;
t38 = pkin(5) + r_base(3);
t37 = t29 * qJ(2) + r_base(1);
t36 = qJ(4) * t26 - qJ(2);
t35 = -t27 * qJ(3) + t38;
t34 = rSges(3,1) * t27 - rSges(3,2) * t24 + pkin(1);
t22 = sin(pkin(9));
t33 = t22 * pkin(4) - t25 * pkin(6) + qJ(3);
t18 = t25 * pkin(4) + t22 * pkin(6) + pkin(3);
t32 = -t18 * t23 + t36;
t17 = t26 * pkin(3) + t49;
t16 = pkin(2) * t27 + t50;
t15 = -t23 * pkin(3) + t36;
t14 = t27 * t39 + t42;
t13 = t27 * t40 - t41;
t12 = t27 * t41 - t40;
t11 = t27 * t42 + t39;
t10 = t25 * t48 + t30 * t26;
t9 = t24 * t22 + t26 * t43;
t8 = -t27 * t22 + t25 * t46;
t7 = t22 * t46 + t43;
t6 = t18 * t26 + t49;
t5 = t17 * t27 + t50;
t4 = t14 * t22 - t25 * t44;
t3 = t12 * t22 - t25 * t45;
t2 = t27 * t47 - t9 * t28;
t1 = t24 * t33 + t6 * t27 + pkin(1);
t19 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t31 * rSges(2,1) - t29 * rSges(2,2) + r_base(1)) + g(2) * (t29 * rSges(2,1) + t31 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t38)) - m(3) * (g(1) * t37 + g(2) * r_base(2) + g(3) * (t24 * rSges(3,1) + t27 * rSges(3,2) + t38) + (g(1) * rSges(3,3) + g(2) * t34) * t29 + (g(1) * t34 + g(2) * (-rSges(3,3) - qJ(2))) * t31) - m(4) * (g(1) * (t14 * rSges(4,1) - t13 * rSges(4,2) + t16 * t31 + t37) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) - t31 * qJ(2) + t16 * t29 + r_base(2)) + g(3) * (-t27 * rSges(4,3) + t35) + (g(3) * (rSges(4,1) * t26 - rSges(4,2) * t23 + pkin(2)) + (g(1) * t31 + g(2) * t29) * rSges(4,3)) * t24) - m(5) * (g(1) * (t5 * t31 - t15 * t29 + r_base(1) + (t14 * t25 + t22 * t44) * rSges(5,1) - t4 * rSges(5,2) + t13 * rSges(5,3)) + g(2) * (t5 * t29 + t15 * t31 + r_base(2) + (t12 * t25 + t22 * t45) * rSges(5,1) - t3 * rSges(5,2) + t11 * rSges(5,3)) + g(3) * (t8 * rSges(5,1) - t7 * rSges(5,2) + (rSges(5,3) * t23 + t17) * t24 + t35)) - m(6) * (g(1) * (t1 * t31 - t32 * t29 + r_base(1) + ((t27 * t48 + t9 * t30) * t31 + t29 * (t25 * t47 - t28 * t26)) * rSges(6,1) + (-t29 * t10 + t2 * t31) * rSges(6,2) + t4 * rSges(6,3)) + g(2) * (t1 * t29 + t32 * t31 + r_base(2) + ((-t25 * t40 + t9 * t29) * t30 + t11 * t28) * rSges(6,1) + (t31 * t10 + t2 * t29) * rSges(6,2) + t3 * rSges(6,3)) + g(3) * (t6 * t24 - t33 * t27 + (t24 * t48 + t8 * t30) * rSges(6,1) + (t24 * t47 - t8 * t28) * rSges(6,2) + t7 * rSges(6,3) + t38));
U = t19;

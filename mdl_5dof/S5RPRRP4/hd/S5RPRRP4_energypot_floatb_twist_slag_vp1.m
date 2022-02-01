% Calculate potential energy for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:42
% EndTime: 2022-01-23 09:31:43
% DurationCPUTime: 0.43s
% Computational Cost: add. (160->94), mult. (173->101), div. (0->0), fcn. (161->8), ass. (0->39)
t48 = rSges(4,3) + pkin(6);
t25 = pkin(7) + pkin(6);
t47 = rSges(6,3) + qJ(5) + t25;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t46 = g(1) * t24 + g(2) * t22;
t21 = sin(qJ(3));
t43 = t21 * pkin(3);
t23 = cos(qJ(3));
t42 = t23 * pkin(3);
t20 = cos(pkin(8));
t41 = pkin(2) * t20 + pkin(1);
t19 = sin(pkin(8));
t40 = rSges(3,2) * t19;
t39 = t20 * t22;
t18 = qJ(3) + qJ(4);
t10 = sin(t18);
t38 = t22 * t10;
t11 = cos(t18);
t37 = t22 * t11;
t36 = t24 * t10;
t35 = t24 * t11;
t33 = pkin(5) + r_base(3);
t32 = t22 * qJ(2) + r_base(1);
t31 = t22 * pkin(1) + r_base(2);
t30 = t19 * pkin(2) + t33;
t29 = t24 * pkin(1) + t32;
t28 = rSges(4,1) * t23 - rSges(4,2) * t21;
t27 = t21 * rSges(4,1) + t23 * rSges(4,2);
t26 = t48 * t19 + t28 * t20 + t41;
t9 = qJ(2) + t43;
t7 = pkin(4) * t10 + t43;
t6 = pkin(4) * t11 + pkin(2) + t42;
t5 = t20 * t35 + t38;
t4 = -t20 * t36 + t37;
t3 = t20 * t37 - t36;
t2 = -t20 * t38 - t35;
t1 = t25 * t19 + t20 * t42 + t41;
t8 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t33)) - m(3) * (g(1) * (t22 * rSges(3,3) + t29) + g(2) * (rSges(3,1) * t39 - t22 * t40 + t31) + g(3) * (t19 * rSges(3,1) + t20 * rSges(3,2) + t33) + (g(1) * (rSges(3,1) * t20 - t40) + g(2) * (-rSges(3,3) - qJ(2))) * t24) - m(4) * (g(1) * (t27 * t22 + t26 * t24 + t32) + g(2) * (r_base(2) + (-qJ(2) - t27) * t24 + t26 * t22) + g(3) * (t28 * t19 - t48 * t20 + t30)) - m(5) * (g(1) * (t5 * rSges(5,1) + t4 * rSges(5,2) + t1 * t24 + t9 * t22 + r_base(1)) + g(2) * (t3 * rSges(5,1) + t2 * rSges(5,2) + t1 * t22 - t9 * t24 + r_base(2)) + g(3) * (t30 + (-rSges(5,3) - t25) * t20) + (g(3) * (rSges(5,1) * t11 - rSges(5,2) * t10 + t42) + t46 * rSges(5,3)) * t19) - m(6) * (g(1) * (t24 * t20 * t6 + t5 * rSges(6,1) + t4 * rSges(6,2) + t22 * t7 + t29) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t6 * t39 + t31 + (-qJ(2) - t7) * t24) + g(3) * (-t47 * t20 + t33) + (g(3) * (rSges(6,1) * t11 - rSges(6,2) * t10 + t6) + t46 * t47) * t19);
U = t8;

% Calculate potential energy for
% S5RPRRP3
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:29
% EndTime: 2022-01-23 09:29:29
% DurationCPUTime: 0.24s
% Computational Cost: add. (147->65), mult. (97->62), div. (0->0), fcn. (73->8), ass. (0->26)
t31 = rSges(6,1) + pkin(4);
t18 = -pkin(7) - pkin(6);
t16 = cos(qJ(3));
t2 = t16 * pkin(3) + pkin(2);
t30 = rSges(4,3) + pkin(6);
t29 = rSges(5,3) - t18;
t28 = rSges(6,3) + qJ(5) - t18;
t27 = pkin(5) + r_base(3);
t15 = sin(qJ(1));
t26 = t15 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t25 = t17 * pkin(1) + r_base(1);
t24 = qJ(2) + t27;
t14 = sin(qJ(3));
t23 = t14 * pkin(3) + t24;
t13 = qJ(3) + qJ(4);
t5 = sin(t13);
t6 = cos(t13);
t22 = rSges(5,1) * t6 - rSges(5,2) * t5 + t2;
t21 = -rSges(6,2) * t5 + t31 * t6 + t2;
t20 = rSges(4,1) * t16 - rSges(4,2) * t14 + pkin(2);
t19 = g(1) * t25 + g(2) * t26;
t12 = qJ(1) + pkin(8);
t4 = cos(t12);
t3 = sin(t12);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t17 * rSges(2,1) - t15 * rSges(2,2) + r_base(1)) + g(2) * (t15 * rSges(2,1) + t17 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t27)) - m(3) * (g(1) * (t4 * rSges(3,1) - t3 * rSges(3,2) + t25) + g(2) * (t3 * rSges(3,1) + t4 * rSges(3,2) + t26) + g(3) * (rSges(3,3) + t24)) - m(4) * (g(3) * (t14 * rSges(4,1) + t16 * rSges(4,2) + t24) + (g(1) * t20 - g(2) * t30) * t4 + (g(1) * t30 + g(2) * t20) * t3 + t19) - m(5) * (g(3) * (t5 * rSges(5,1) + t6 * rSges(5,2) + t23) + (g(1) * t22 - g(2) * t29) * t4 + (g(1) * t29 + g(2) * t22) * t3 + t19) - m(6) * (g(3) * (t6 * rSges(6,2) + t31 * t5 + t23) + (g(1) * t21 - g(2) * t28) * t4 + (g(1) * t28 + g(2) * t21) * t3 + t19);
U = t1;

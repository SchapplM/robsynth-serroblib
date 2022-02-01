% Calculate potential energy for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:28
% EndTime: 2022-01-20 10:19:29
% DurationCPUTime: 0.25s
% Computational Cost: add. (148->62), mult. (85->58), div. (0->0), fcn. (61->8), ass. (0->24)
t29 = rSges(6,1) + pkin(4);
t28 = rSges(5,3) + pkin(7);
t27 = rSges(6,3) + qJ(5) + pkin(7);
t11 = qJ(1) + qJ(2);
t26 = pkin(5) + r_base(3);
t14 = sin(qJ(1));
t25 = t14 * pkin(1) + r_base(2);
t16 = cos(qJ(1));
t24 = t16 * pkin(1) + r_base(1);
t23 = pkin(6) + t26;
t7 = sin(t11);
t22 = pkin(2) * t7 + t25;
t8 = cos(t11);
t21 = pkin(2) * t8 + t24;
t20 = qJ(3) + t23;
t13 = sin(qJ(4));
t15 = cos(qJ(4));
t19 = rSges(5,1) * t15 - rSges(5,2) * t13 + pkin(3);
t18 = -rSges(6,2) * t13 + t29 * t15 + pkin(3);
t17 = g(1) * t21 + g(2) * t22;
t6 = pkin(8) + t11;
t2 = cos(t6);
t1 = sin(t6);
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t16 * rSges(2,1) - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + t16 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t24) + g(2) * (t7 * rSges(3,1) + t8 * rSges(3,2) + t25) + g(3) * (rSges(3,3) + t23)) - m(4) * (g(1) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t21) + g(2) * (t1 * rSges(4,1) + t2 * rSges(4,2) + t22) + g(3) * (rSges(4,3) + t20)) - m(5) * (g(3) * (t13 * rSges(5,1) + t15 * rSges(5,2) + t20) + (g(1) * t19 - g(2) * t28) * t2 + (g(1) * t28 + g(2) * t19) * t1 + t17) - m(6) * (g(3) * (t15 * rSges(6,2) + t29 * t13 + t20) + (g(1) * t18 - g(2) * t27) * t2 + (g(1) * t27 + g(2) * t18) * t1 + t17);
U = t3;

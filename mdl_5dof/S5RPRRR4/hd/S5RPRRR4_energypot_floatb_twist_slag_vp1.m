% Calculate potential energy for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:25
% EndTime: 2022-01-23 09:34:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (152->62), mult. (74->56), div. (0->0), fcn. (50->10), ass. (0->26)
t31 = rSges(6,3) + pkin(8);
t30 = pkin(5) + r_base(3);
t15 = qJ(1) + pkin(9);
t17 = sin(qJ(1));
t29 = t17 * pkin(1) + r_base(2);
t19 = cos(qJ(1));
t28 = t19 * pkin(1) + r_base(1);
t10 = sin(t15);
t27 = pkin(2) * t10 + t29;
t11 = cos(t15);
t26 = pkin(2) * t11 + t28;
t25 = qJ(2) + t30;
t12 = qJ(3) + t15;
t7 = sin(t12);
t24 = pkin(3) * t7 + t27;
t8 = cos(t12);
t23 = pkin(3) * t8 + t26;
t22 = pkin(6) + t25;
t21 = pkin(7) + t22;
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t20 = rSges(6,1) * t18 - rSges(6,2) * t16 + pkin(4);
t9 = qJ(4) + t12;
t4 = cos(t9);
t3 = sin(t9);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t28) + g(2) * (t10 * rSges(3,1) + t11 * rSges(3,2) + t29) + g(3) * (rSges(3,3) + t25)) - m(4) * (g(1) * (t8 * rSges(4,1) - t7 * rSges(4,2) + t26) + g(2) * (t7 * rSges(4,1) + t8 * rSges(4,2) + t27) + g(3) * (rSges(4,3) + t22)) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t23) + g(2) * (t3 * rSges(5,1) + t4 * rSges(5,2) + t24) + g(3) * (rSges(5,3) + t21)) - m(6) * (g(1) * t23 + g(2) * t24 + g(3) * (t16 * rSges(6,1) + t18 * rSges(6,2) + t21) + (g(1) * t20 - g(2) * t31) * t4 + (g(1) * t31 + g(2) * t20) * t3);
U = t1;

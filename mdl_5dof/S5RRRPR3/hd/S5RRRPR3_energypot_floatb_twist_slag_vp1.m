% Calculate potential energy for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:23
% EndTime: 2022-01-20 11:42:23
% DurationCPUTime: 0.21s
% Computational Cost: add. (153->67), mult. (97->64), div. (0->0), fcn. (73->10), ass. (0->28)
t33 = rSges(4,3) + pkin(7);
t20 = cos(qJ(3));
t4 = t20 * pkin(3) + pkin(2);
t17 = -qJ(4) - pkin(7);
t32 = rSges(5,3) - t17;
t31 = rSges(6,3) + pkin(8) - t17;
t30 = pkin(5) + r_base(3);
t15 = qJ(3) + pkin(9);
t19 = sin(qJ(1));
t29 = t19 * pkin(1) + r_base(2);
t21 = cos(qJ(1));
t28 = t21 * pkin(1) + r_base(1);
t27 = pkin(6) + t30;
t18 = sin(qJ(3));
t26 = t18 * pkin(3) + t27;
t5 = sin(t15);
t6 = cos(t15);
t25 = rSges(5,1) * t6 - rSges(5,2) * t5 + t4;
t7 = qJ(5) + t15;
t2 = sin(t7);
t3 = cos(t7);
t24 = rSges(6,1) * t3 - rSges(6,2) * t2 + pkin(4) * t6 + t4;
t23 = rSges(4,1) * t20 - rSges(4,2) * t18 + pkin(2);
t22 = g(1) * t28 + g(2) * t29;
t16 = qJ(1) + qJ(2);
t9 = cos(t16);
t8 = sin(t16);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t21 * rSges(2,1) - t19 * rSges(2,2) + r_base(1)) + g(2) * (t19 * rSges(2,1) + t21 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t30)) - m(3) * (g(1) * (t9 * rSges(3,1) - t8 * rSges(3,2) + t28) + g(2) * (t8 * rSges(3,1) + t9 * rSges(3,2) + t29) + g(3) * (rSges(3,3) + t27)) - m(4) * (g(3) * (t18 * rSges(4,1) + t20 * rSges(4,2) + t27) + (g(1) * t23 - g(2) * t33) * t9 + (g(1) * t33 + g(2) * t23) * t8 + t22) - m(5) * (g(3) * (t5 * rSges(5,1) + t6 * rSges(5,2) + t26) + (g(1) * t25 - g(2) * t32) * t9 + (g(1) * t32 + g(2) * t25) * t8 + t22) - m(6) * (g(3) * (t2 * rSges(6,1) + t3 * rSges(6,2) + pkin(4) * t5 + t26) + (g(1) * t24 - g(2) * t31) * t9 + (g(1) * t31 + g(2) * t24) * t8 + t22);
U = t1;

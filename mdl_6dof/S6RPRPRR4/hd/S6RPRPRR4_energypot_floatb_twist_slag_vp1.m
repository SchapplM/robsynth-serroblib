% Calculate potential energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:19
% EndTime: 2019-03-09 03:44:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (232->101), mult. (190->115), div. (0->0), fcn. (170->10), ass. (0->36)
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t56 = pkin(3) * t22 + qJ(4) * t19;
t55 = rSges(6,3) + pkin(8);
t54 = rSges(7,3) + pkin(9) + pkin(8);
t53 = -rSges(5,2) * t22 + rSges(5,3) * t19;
t52 = rSges(4,1) * t22 - rSges(4,2) * t19;
t16 = qJ(1) + pkin(10);
t10 = cos(t16);
t9 = sin(t16);
t51 = g(1) * t10 + g(2) * t9;
t17 = qJ(5) + qJ(6);
t11 = sin(t17);
t42 = t11 * t19;
t12 = cos(t17);
t41 = t12 * t19;
t18 = sin(qJ(5));
t40 = t18 * t19;
t21 = cos(qJ(5));
t39 = t19 * t21;
t36 = pkin(6) + r_base(3);
t35 = t9 * t40;
t20 = sin(qJ(1));
t34 = t20 * pkin(1) + r_base(2);
t23 = cos(qJ(1));
t33 = t23 * pkin(1) + r_base(1);
t32 = t10 * t40;
t31 = t9 * pkin(2) + t34;
t30 = qJ(2) + t36;
t29 = t10 * pkin(2) + t9 * pkin(7) + t33;
t28 = t19 * pkin(3) + t30;
t27 = t56 * t9 + t31;
t26 = t56 * t10 + t29;
t25 = -t10 * pkin(7) + t27;
t8 = t21 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t33) + g(2) * (t9 * rSges(3,1) + t10 * rSges(3,2) + t34) + g(3) * (rSges(3,3) + t30)) - m(4) * (g(1) * (t9 * rSges(4,3) + t29) + g(2) * (t52 * t9 + t31) + g(3) * (t19 * rSges(4,1) + t22 * rSges(4,2) + t30) + (g(1) * t52 + g(2) * (-rSges(4,3) - pkin(7))) * t10) - m(5) * (g(1) * (t9 * rSges(5,1) + t26) + g(2) * (t53 * t9 + t27) + g(3) * (-t19 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t22 + t28) + (g(1) * t53 + g(2) * (-rSges(5,1) - pkin(7))) * t10) - m(6) * (g(1) * (t9 * pkin(4) + (t9 * t21 + t32) * rSges(6,1) + (t10 * t39 - t9 * t18) * rSges(6,2) + t26) + g(2) * (-t10 * pkin(4) + (-t10 * t21 + t35) * rSges(6,1) + (t10 * t18 + t9 * t39) * rSges(6,2) + t25) + g(3) * (t55 * t19 + t28) + (g(3) * (-rSges(6,1) * t18 - rSges(6,2) * t21 - qJ(4)) + t51 * t55) * t22) - m(7) * (g(1) * (t9 * t8 + pkin(5) * t32 + (t10 * t42 + t9 * t12) * rSges(7,1) + (t10 * t41 - t9 * t11) * rSges(7,2) + t26) + g(2) * (-t10 * t8 + pkin(5) * t35 + (-t10 * t12 + t9 * t42) * rSges(7,1) + (t10 * t11 + t9 * t41) * rSges(7,2) + t25) + g(3) * (t54 * t19 + t28) + (g(3) * (-rSges(7,1) * t11 - rSges(7,2) * t12 - pkin(5) * t18 - qJ(4)) + t51 * t54) * t22);
U  = t1;

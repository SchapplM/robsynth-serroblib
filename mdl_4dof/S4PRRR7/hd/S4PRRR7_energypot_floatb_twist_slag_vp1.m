% Calculate potential energy for
% S4PRRR7
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
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR7_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:00
% EndTime: 2019-12-31 16:36:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (145->83), mult. (272->108), div. (0->0), fcn. (311->10), ass. (0->37)
t22 = sin(pkin(4));
t46 = pkin(5) * t22;
t45 = rSges(4,3) + pkin(6);
t44 = pkin(7) + rSges(5,3);
t26 = sin(qJ(3));
t43 = t22 * t26;
t27 = sin(qJ(2));
t42 = t22 * t27;
t29 = cos(qJ(3));
t41 = t22 * t29;
t30 = cos(qJ(2));
t40 = t22 * t30;
t24 = cos(pkin(4));
t39 = t24 * t27;
t38 = t24 * t30;
t21 = sin(pkin(8));
t37 = t21 * pkin(1) + r_base(2);
t36 = qJ(1) + r_base(3);
t23 = cos(pkin(8));
t35 = t23 * pkin(1) + t21 * t46 + r_base(1);
t34 = t24 * pkin(5) + t36;
t10 = -t21 * t39 + t23 * t30;
t33 = t10 * pkin(2) + t35;
t32 = pkin(2) * t42 + t34;
t8 = t21 * t30 + t23 * t39;
t31 = t8 * pkin(2) - t23 * t46 + t37;
t28 = cos(qJ(4));
t25 = sin(qJ(4));
t12 = t24 * t26 + t27 * t41;
t11 = -t24 * t29 + t26 * t42;
t9 = t21 * t38 + t23 * t27;
t7 = t21 * t27 - t23 * t38;
t4 = t10 * t29 + t21 * t43;
t3 = t10 * t26 - t21 * t41;
t2 = -t23 * t43 + t8 * t29;
t1 = t23 * t41 + t8 * t26;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t35) + g(2) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t37) + g(3) * (t24 * rSges(3,3) + t34) + (g(1) * rSges(3,3) * t21 + g(3) * (rSges(3,1) * t27 + rSges(3,2) * t30) + g(2) * (-rSges(3,3) - pkin(5)) * t23) * t22) - m(4) * (g(1) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t45 * t9 + t33) + g(2) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t45 * t7 + t31) + g(3) * (t12 * rSges(4,1) - t11 * rSges(4,2) - t45 * t40 + t32)) - m(5) * (g(1) * (t4 * pkin(3) + t9 * pkin(6) + (t9 * t25 + t4 * t28) * rSges(5,1) + (-t4 * t25 + t9 * t28) * rSges(5,2) + t44 * t3 + t33) + g(2) * (t2 * pkin(3) + t7 * pkin(6) + (t2 * t28 + t7 * t25) * rSges(5,1) + (-t2 * t25 + t7 * t28) * rSges(5,2) + t44 * t1 + t31) + g(3) * (t12 * pkin(3) - pkin(6) * t40 + (t12 * t28 - t25 * t40) * rSges(5,1) + (-t12 * t25 - t28 * t40) * rSges(5,2) + t44 * t11 + t32));
U = t5;

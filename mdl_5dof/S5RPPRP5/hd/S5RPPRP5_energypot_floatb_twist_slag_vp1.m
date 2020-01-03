% Calculate potential energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:17
% EndTime: 2019-12-31 17:53:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (130->84), mult. (194->91), div. (0->0), fcn. (188->6), ass. (0->29)
t42 = rSges(6,1) + pkin(4);
t21 = sin(pkin(7));
t24 = sin(qJ(1));
t41 = t21 * t24;
t25 = cos(qJ(4));
t40 = t21 * t25;
t22 = cos(pkin(7));
t39 = t22 * t24;
t26 = cos(qJ(1));
t38 = t22 * t26;
t37 = qJ(3) * t21;
t36 = rSges(6,3) + qJ(5);
t35 = pkin(5) + r_base(3);
t34 = t24 * pkin(1) + r_base(2);
t33 = t21 * pkin(2) + t35;
t32 = t26 * pkin(1) + t24 * qJ(2) + r_base(1);
t31 = pkin(2) * t39 + t24 * t37 + t34;
t23 = sin(qJ(4));
t5 = t21 * t23 + t22 * t25;
t30 = pkin(2) * t38 + t26 * t37 + t32;
t29 = pkin(3) * t39 + t26 * pkin(6) + t31;
t28 = pkin(3) * t38 + t30;
t27 = t21 * pkin(3) - t22 * qJ(3) + t33;
t6 = -t22 * t23 + t40;
t4 = t5 * t26;
t3 = t23 * t38 - t26 * t40;
t2 = t5 * t24;
t1 = t23 * t39 - t24 * t40;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t26 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t26 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t35)) - m(3) * (g(1) * (t24 * rSges(3,3) + t32) + g(2) * (rSges(3,1) * t39 - rSges(3,2) * t41 + t34) + g(3) * (t21 * rSges(3,1) + t22 * rSges(3,2) + t35) + (g(1) * (rSges(3,1) * t22 - rSges(3,2) * t21) + g(2) * (-rSges(3,3) - qJ(2))) * t26) - m(4) * (g(1) * (t24 * rSges(4,2) + t30) + g(2) * (rSges(4,1) * t39 + rSges(4,3) * t41 + t31) + g(3) * (t21 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t22 + t33) + (g(1) * (rSges(4,1) * t22 + rSges(4,3) * t21) + g(2) * (-rSges(4,2) - qJ(2))) * t26) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + (-rSges(5,3) - pkin(6)) * t24 + t28) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + (rSges(5,3) - qJ(2)) * t26 + t29) + g(3) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t27)) - m(6) * (g(1) * (t42 * t4 + t36 * t3 + (-rSges(6,2) - pkin(6)) * t24 + t28) + g(2) * ((rSges(6,2) - qJ(2)) * t26 + t42 * t2 + t36 * t1 + t29) + g(3) * (t36 * t5 + t42 * t6 + t27));
U = t7;

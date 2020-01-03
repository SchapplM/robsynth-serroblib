% Calculate potential energy for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:32
% EndTime: 2019-12-31 17:47:32
% DurationCPUTime: 0.37s
% Computational Cost: add. (134->93), mult. (203->110), div. (0->0), fcn. (199->8), ass. (0->32)
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t33 = t23 * pkin(1) + r_base(2);
t19 = sin(pkin(7));
t36 = qJ(3) * t19;
t21 = cos(pkin(7));
t42 = t21 * t23;
t29 = pkin(2) * t42 + t23 * t36 + t33;
t35 = qJ(4) * t21;
t45 = t23 * t35 + t29 + (-pkin(3) - qJ(2)) * t25;
t44 = pkin(6) + rSges(6,3);
t43 = t19 * t23;
t41 = t21 * t25;
t18 = sin(pkin(8));
t40 = t23 * t18;
t20 = cos(pkin(8));
t39 = t23 * t20;
t38 = t25 * t18;
t37 = t25 * t20;
t34 = pkin(5) + r_base(3);
t32 = t19 * pkin(2) + t34;
t31 = t25 * pkin(1) + t23 * qJ(2) + r_base(1);
t30 = t19 * qJ(4) + t32;
t27 = pkin(2) * t41 + t25 * t36 + t31;
t26 = t23 * pkin(3) + t25 * t35 + t27;
t24 = cos(qJ(5));
t22 = sin(qJ(5));
t4 = t19 * t40 - t37;
t3 = t19 * t39 + t38;
t2 = t19 * t38 + t39;
t1 = -t19 * t37 + t40;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t23 * rSges(2,2) + r_base(1)) + g(2) * (t23 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t23 * rSges(3,3) + t31) + g(2) * (rSges(3,1) * t42 - rSges(3,2) * t43 + t33) + g(3) * (t19 * rSges(3,1) + t21 * rSges(3,2) + t34) + (g(1) * (rSges(3,1) * t21 - rSges(3,2) * t19) + g(2) * (-rSges(3,3) - qJ(2))) * t25) - m(4) * (g(1) * (t23 * rSges(4,1) + t27) + g(2) * (-rSges(4,2) * t42 + rSges(4,3) * t43 + t29) + g(3) * (-t19 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t21 + t32) + (g(1) * (-rSges(4,2) * t21 + rSges(4,3) * t19) + g(2) * (-rSges(4,1) - qJ(2))) * t25) - m(5) * (g(1) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t26) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t45) + g(3) * (t19 * rSges(5,3) + t30) + (g(3) * (-rSges(5,1) * t18 - rSges(5,2) * t20 - qJ(3)) + (g(1) * t25 + g(2) * t23) * rSges(5,3)) * t21) - m(6) * (g(1) * (t2 * pkin(4) + (t2 * t24 + t22 * t41) * rSges(6,1) + (-t2 * t22 + t24 * t41) * rSges(6,2) + t44 * t1 + t26) + g(2) * (t4 * pkin(4) + (t22 * t42 + t4 * t24) * rSges(6,1) + (-t4 * t22 + t24 * t42) * rSges(6,2) - t44 * t3 + t45) + g(3) * ((t22 * rSges(6,1) + t24 * rSges(6,2)) * t19 + (-qJ(3) + t44 * t20 + (-t24 * rSges(6,1) + t22 * rSges(6,2) - pkin(4)) * t18) * t21 + t30));
U = t5;

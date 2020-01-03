% Calculate potential energy for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:28:56
% EndTime: 2019-12-31 20:28:56
% DurationCPUTime: 0.46s
% Computational Cost: add. (138->90), mult. (215->104), div. (0->0), fcn. (216->8), ass. (0->30)
t44 = pkin(8) + rSges(6,3);
t22 = sin(qJ(4));
t27 = cos(qJ(2));
t23 = sin(qJ(2));
t26 = cos(qJ(4));
t42 = t23 * t26;
t45 = t27 * t22 - t42;
t24 = sin(qJ(1));
t43 = t23 * t24;
t41 = t24 * t27;
t28 = cos(qJ(1));
t39 = t27 * t28;
t38 = qJ(3) * t23;
t37 = pkin(5) + r_base(3);
t36 = t24 * pkin(1) + r_base(2);
t35 = t23 * pkin(2) + t37;
t34 = t28 * pkin(1) + t24 * pkin(6) + r_base(1);
t33 = pkin(2) * t41 + t24 * t38 + t36;
t5 = t23 * t22 + t27 * t26;
t32 = pkin(2) * t39 + t28 * t38 + t34;
t31 = pkin(3) * t41 + t28 * pkin(7) + t33;
t30 = pkin(3) * t39 + t32;
t29 = t23 * pkin(3) - t27 * qJ(3) + t35;
t25 = cos(qJ(5));
t21 = sin(qJ(5));
t4 = t5 * t28;
t3 = t22 * t39 - t28 * t42;
t2 = t5 * t24;
t1 = t45 * t24;
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t28 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t28 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(1) * (t24 * rSges(3,3) + t34) + g(2) * (rSges(3,1) * t41 - rSges(3,2) * t43 + t36) + g(3) * (t23 * rSges(3,1) + t27 * rSges(3,2) + t37) + (g(1) * (rSges(3,1) * t27 - rSges(3,2) * t23) + g(2) * (-rSges(3,3) - pkin(6))) * t28) - m(4) * (g(1) * (t24 * rSges(4,2) + t32) + g(2) * (rSges(4,1) * t41 + rSges(4,3) * t43 + t33) + g(3) * (t23 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t27 + t35) + (g(1) * (rSges(4,1) * t27 + rSges(4,3) * t23) + g(2) * (-rSges(4,2) - pkin(6))) * t28) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + (-rSges(5,3) - pkin(7)) * t24 + t30) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + (rSges(5,3) - pkin(6)) * t28 + t31) + g(3) * (-rSges(5,1) * t45 - t5 * rSges(5,2) + t29)) - m(6) * (g(1) * (t4 * pkin(4) - t24 * pkin(7) + (-t24 * t21 + t4 * t25) * rSges(6,1) + (-t4 * t21 - t24 * t25) * rSges(6,2) + t44 * t3 + t30) + g(2) * (t2 * pkin(4) - t28 * pkin(6) + (t2 * t25 + t28 * t21) * rSges(6,1) + (-t2 * t21 + t28 * t25) * rSges(6,2) + t44 * t1 + t31) + (t29 - (t25 * rSges(6,1) - t21 * rSges(6,2) + pkin(4)) * t45 + t44 * t5) * g(3));
U = t6;

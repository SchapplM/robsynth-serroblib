% Calculate potential energy for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:57
% EndTime: 2019-12-31 17:37:58
% DurationCPUTime: 0.45s
% Computational Cost: add. (138->89), mult. (215->103), div. (0->0), fcn. (216->8), ass. (0->29)
t44 = pkin(6) + rSges(6,3);
t21 = sin(pkin(8));
t23 = cos(pkin(8));
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t45 = t28 * t21 - t26 * t23;
t22 = sin(pkin(7));
t43 = t22 * t26;
t42 = t22 * t28;
t24 = cos(pkin(7));
t41 = t24 * t28;
t38 = qJ(3) * t26;
t37 = t22 * pkin(1) + r_base(2);
t36 = qJ(1) + r_base(3);
t35 = t24 * pkin(1) + t22 * pkin(5) + r_base(1);
t34 = t26 * pkin(2) + t36;
t33 = pkin(2) * t42 + t22 * t38 + t37;
t5 = t26 * t21 + t28 * t23;
t32 = pkin(2) * t41 + t24 * t38 + t35;
t31 = pkin(3) * t42 + t24 * qJ(4) + t33;
t30 = pkin(3) * t41 + t32;
t29 = t26 * pkin(3) - t28 * qJ(3) + t34;
t27 = cos(qJ(5));
t25 = sin(qJ(5));
t4 = t5 * t24;
t3 = t45 * t24;
t2 = t5 * t22;
t1 = t45 * t22;
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * (t22 * rSges(3,3) + t35) + g(2) * (rSges(3,1) * t42 - rSges(3,2) * t43 + t37) + g(3) * (t26 * rSges(3,1) + t28 * rSges(3,2) + t36) + (g(1) * (rSges(3,1) * t28 - rSges(3,2) * t26) + g(2) * (-rSges(3,3) - pkin(5))) * t24) - m(4) * (g(1) * (t22 * rSges(4,2) + t32) + g(2) * (rSges(4,1) * t42 + rSges(4,3) * t43 + t33) + g(3) * (t26 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t28 + t34) + (g(1) * (rSges(4,1) * t28 + rSges(4,3) * t26) + g(2) * (-rSges(4,2) - pkin(5))) * t24) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t22 + t30) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + (rSges(5,3) - pkin(5)) * t24 + t31) + g(3) * (-rSges(5,1) * t45 - t5 * rSges(5,2) + t29)) - m(6) * (g(1) * (t4 * pkin(4) - t22 * qJ(4) + (-t22 * t25 + t4 * t27) * rSges(6,1) + (-t22 * t27 - t4 * t25) * rSges(6,2) + t44 * t3 + t30) + g(2) * (t2 * pkin(4) - t24 * pkin(5) + (t2 * t27 + t24 * t25) * rSges(6,1) + (-t2 * t25 + t24 * t27) * rSges(6,2) + t44 * t1 + t31) + (t29 - (t27 * rSges(6,1) - t25 * rSges(6,2) + pkin(4)) * t45 + t44 * t5) * g(3));
U = t6;

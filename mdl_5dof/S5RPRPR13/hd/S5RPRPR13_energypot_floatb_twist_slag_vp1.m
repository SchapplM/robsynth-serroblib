% Calculate potential energy for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR13_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR13_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR13_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:51
% EndTime: 2019-12-31 18:31:52
% DurationCPUTime: 0.39s
% Computational Cost: add. (151->84), mult. (141->93), div. (0->0), fcn. (121->8), ass. (0->29)
t38 = rSges(6,3) + pkin(7);
t13 = pkin(8) + qJ(3);
t10 = sin(t13);
t20 = cos(qJ(1));
t36 = t10 * t20;
t17 = sin(qJ(5));
t18 = sin(qJ(1));
t35 = t18 * t17;
t19 = cos(qJ(5));
t34 = t18 * t19;
t11 = cos(t13);
t33 = t20 * t11;
t32 = t20 * t17;
t31 = t20 * t19;
t30 = qJ(4) * t10;
t29 = rSges(3,3) + qJ(2);
t28 = pkin(5) + r_base(3);
t15 = cos(pkin(8));
t8 = pkin(2) * t15 + pkin(1);
t27 = t20 * t8 + r_base(1);
t16 = -pkin(6) - qJ(2);
t26 = t20 * t16 + t18 * t8 + r_base(2);
t14 = sin(pkin(8));
t25 = t14 * pkin(2) + t28;
t24 = pkin(3) * t33 + t20 * t30 + t27;
t23 = t10 * pkin(3) + t25;
t22 = t26 + (pkin(3) * t11 + t30) * t18;
t21 = rSges(3,1) * t15 - rSges(3,2) * t14 + pkin(1);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t20 - rSges(2,2) * t18 + r_base(1)) + g(2) * (rSges(2,1) * t18 + rSges(2,2) * t20 + r_base(2)) + g(3) * (rSges(2,3) + t28)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t14 + rSges(3,2) * t15 + t28) + (g(1) * t21 - g(2) * t29) * t20 + (g(1) * t29 + g(2) * t21) * t18) - m(4) * (g(1) * (rSges(4,1) * t33 - rSges(4,2) * t36 + t27) + g(2) * (-t20 * rSges(4,3) + t26) + g(3) * (rSges(4,1) * t10 + rSges(4,2) * t11 + t25) + (g(1) * (rSges(4,3) - t16) + g(2) * (rSges(4,1) * t11 - rSges(4,2) * t10)) * t18) - m(5) * (g(1) * (-rSges(5,2) * t33 + rSges(5,3) * t36 + t24) + g(2) * (-t20 * rSges(5,1) + t22) + g(3) * (-t10 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t11 + t23) + (g(1) * (rSges(5,1) - t16) + g(2) * (-rSges(5,2) * t11 + rSges(5,3) * t10)) * t18) - m(6) * (g(1) * ((t10 * t32 + t34) * rSges(6,1) + (t10 * t31 - t35) * rSges(6,2) + t24 + (pkin(4) - t16) * t18) + g(2) * (-t20 * pkin(4) + (t10 * t35 - t31) * rSges(6,1) + (t10 * t34 + t32) * rSges(6,2) + t22) + g(3) * (t38 * t10 + t23) + (g(3) * (-rSges(6,1) * t17 - rSges(6,2) * t19 - qJ(4)) + (g(1) * t20 + g(2) * t18) * t38) * t11);
U = t1;

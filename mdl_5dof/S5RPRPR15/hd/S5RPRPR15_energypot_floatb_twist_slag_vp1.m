% Calculate potential energy for
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR15_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:33
% EndTime: 2019-12-31 18:36:33
% DurationCPUTime: 0.39s
% Computational Cost: add. (124->85), mult. (143->95), div. (0->0), fcn. (127->8), ass. (0->26)
t35 = rSges(6,3) + pkin(7) + qJ(4);
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t34 = -g(1) * t14 + g(2) * t16;
t33 = rSges(5,3) + qJ(4);
t15 = cos(qJ(3));
t30 = rSges(4,2) * t15;
t10 = sin(pkin(8));
t29 = t14 * t10;
t13 = sin(qJ(3));
t28 = t14 * t13;
t27 = t16 * t10;
t26 = t16 * t13;
t23 = pkin(5) + r_base(3);
t22 = t14 * pkin(1) + r_base(2);
t21 = pkin(2) + t23;
t20 = t16 * pkin(1) + t14 * qJ(2) + r_base(1);
t19 = t14 * pkin(6) + t22;
t18 = t16 * pkin(6) + t20;
t17 = -t16 * qJ(2) + t19;
t11 = cos(pkin(8));
t9 = pkin(8) + qJ(5);
t3 = cos(t9);
t2 = sin(t9);
t1 = t11 * pkin(4) + pkin(3);
t4 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t16 * rSges(2,1) - t14 * rSges(2,2) + r_base(1)) + g(2) * (t14 * rSges(2,1) + t16 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t23)) - m(3) * (g(1) * (-t16 * rSges(3,2) + t14 * rSges(3,3) + t20) + g(2) * (-t14 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t16 + t22) + g(3) * (rSges(3,1) + t23)) - m(4) * (g(1) * (rSges(4,1) * t28 + t14 * t30 + t18) + g(2) * (t14 * rSges(4,3) + t19) + g(3) * (t15 * rSges(4,1) - t13 * rSges(4,2) + t21) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t13 - qJ(2) - t30)) * t16) - m(5) * (g(1) * (pkin(3) * t28 + (t11 * t28 + t27) * rSges(5,1) + (-t10 * t28 + t16 * t11) * rSges(5,2) + t18) + g(2) * (-pkin(3) * t26 + (-t11 * t26 + t29) * rSges(5,1) + (t10 * t26 + t14 * t11) * rSges(5,2) + t17) + g(3) * (t33 * t13 + t21) + (g(3) * (rSges(5,1) * t11 - rSges(5,2) * t10 + pkin(3)) + t34 * t33) * t15) - m(6) * (g(1) * (t1 * t28 + pkin(4) * t27 + (t16 * t2 + t3 * t28) * rSges(6,1) + (t16 * t3 - t2 * t28) * rSges(6,2) + t18) + g(2) * (-t1 * t26 + pkin(4) * t29 + (t14 * t2 - t3 * t26) * rSges(6,1) + (t14 * t3 + t2 * t26) * rSges(6,2) + t17) + g(3) * (t35 * t13 + t21) + (g(3) * (rSges(6,1) * t3 - rSges(6,2) * t2 + t1) + t34 * t35) * t15);
U = t4;

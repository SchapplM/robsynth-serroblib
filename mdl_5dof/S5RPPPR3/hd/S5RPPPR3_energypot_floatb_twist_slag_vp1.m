% Calculate potential energy for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:47
% EndTime: 2019-12-31 17:43:48
% DurationCPUTime: 0.34s
% Computational Cost: add. (156->73), mult. (131->77), div. (0->0), fcn. (113->8), ass. (0->26)
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t41 = pkin(3) * t15 + qJ(4) * t14;
t40 = rSges(5,1) * t15 + rSges(5,3) * t14;
t39 = rSges(4,1) * t15 - rSges(4,2) * t14;
t37 = -pkin(6) - rSges(6,3);
t31 = pkin(5) + r_base(3);
t17 = sin(qJ(1));
t30 = t17 * pkin(1) + r_base(2);
t19 = cos(qJ(1));
t29 = t19 * pkin(1) + r_base(1);
t13 = qJ(1) + pkin(7);
t8 = sin(t13);
t28 = t8 * pkin(2) + t30;
t27 = qJ(2) + t31;
t9 = cos(t13);
t26 = t9 * pkin(2) + t8 * qJ(3) + t29;
t25 = t14 * pkin(3) + t27;
t24 = t41 * t8 + t28;
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t23 = t14 * t18 - t15 * t16;
t22 = t14 * t16 + t15 * t18;
t21 = t41 * t9 + t26;
t20 = t22 * rSges(6,1) + t23 * rSges(6,2) + t15 * pkin(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t17 * rSges(2,2) + r_base(1)) + g(2) * (t17 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t31)) - m(3) * (g(1) * (t9 * rSges(3,1) - t8 * rSges(3,2) + t29) + g(2) * (t8 * rSges(3,1) + t9 * rSges(3,2) + t30) + g(3) * (rSges(3,3) + t27)) - m(4) * (g(1) * (t8 * rSges(4,3) + t26) + g(2) * (t39 * t8 + t28) + g(3) * (t14 * rSges(4,1) + t15 * rSges(4,2) + t27) + (g(1) * t39 + g(2) * (-rSges(4,3) - qJ(3))) * t9) - m(5) * (g(1) * (t8 * rSges(5,2) + t21) + g(2) * (t40 * t8 + t24) + g(3) * (t14 * rSges(5,1) + (-rSges(5,3) - qJ(4)) * t15 + t25) + (g(1) * t40 + g(2) * (-rSges(5,2) - qJ(3))) * t9) - m(6) * (g(1) * t21 + g(2) * t24 + g(3) * (t23 * rSges(6,1) - t22 * rSges(6,2) + t14 * pkin(4) - t15 * qJ(4) + t25) + (g(1) * t37 + g(2) * t20) * t8 + (g(1) * t20 + g(2) * (-qJ(3) - t37)) * t9);
U = t1;

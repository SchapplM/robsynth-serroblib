% Calculate potential energy for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:13
% EndTime: 2019-12-31 17:32:14
% DurationCPUTime: 0.23s
% Computational Cost: add. (124->66), mult. (141->66), div. (0->0), fcn. (141->8), ass. (0->23)
t32 = cos(qJ(3));
t31 = sin(qJ(3));
t30 = rSges(6,3) + pkin(6) + qJ(4);
t29 = rSges(5,3) + qJ(4);
t28 = sin(pkin(7));
t27 = t28 * pkin(1) + r_base(2);
t26 = qJ(1) + r_base(3);
t17 = cos(pkin(7));
t25 = t17 * pkin(1) + t28 * qJ(2) + r_base(1);
t24 = -pkin(5) + t26;
t23 = t17 * pkin(2) + t25;
t16 = cos(pkin(8));
t14 = pkin(8) + qJ(5);
t7 = sin(t14);
t8 = cos(t14);
t22 = -rSges(6,1) * t8 + rSges(6,2) * t7 - t16 * pkin(4) - pkin(3);
t15 = sin(pkin(8));
t21 = -rSges(5,1) * t16 + rSges(5,2) * t15 - pkin(3);
t20 = t28 * pkin(2) - t17 * qJ(2) + t27;
t19 = g(1) * t23 + g(2) * t20;
t2 = t17 * t31 - t28 * t32;
t1 = -t17 * t32 - t28 * t31;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t17 * rSges(2,1) - t28 * rSges(2,2) + r_base(1)) + g(2) * (t28 * rSges(2,1) + t17 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * (t17 * rSges(3,1) + t28 * rSges(3,3) + t25) + g(2) * (t28 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t17 + t27) + g(3) * (rSges(3,2) + t26)) - m(4) * (g(1) * (-t1 * rSges(4,1) - t2 * rSges(4,2) + t23) + g(2) * (-t2 * rSges(4,1) + t1 * rSges(4,2) + t20) + g(3) * (-rSges(4,3) + t24)) - m(5) * (g(3) * (-t15 * rSges(5,1) - t16 * rSges(5,2) + t24) + (g(1) * t29 + g(2) * t21) * t2 + (g(1) * t21 - g(2) * t29) * t1 + t19) - m(6) * (g(3) * (-t7 * rSges(6,1) - t8 * rSges(6,2) - t15 * pkin(4) + t24) + (g(1) * t30 + g(2) * t22) * t2 + (g(1) * t22 - g(2) * t30) * t1 + t19);
U = t3;

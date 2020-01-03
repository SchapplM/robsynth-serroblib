% Calculate potential energy for
% S5PPRPR5
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:22
% EndTime: 2019-12-31 17:33:22
% DurationCPUTime: 0.17s
% Computational Cost: add. (113->64), mult. (132->62), div. (0->0), fcn. (132->6), ass. (0->20)
t29 = -rSges(6,3) - pkin(6);
t28 = cos(qJ(3));
t27 = sin(qJ(3));
t26 = rSges(5,3) + qJ(4);
t25 = sin(pkin(7));
t24 = t25 * pkin(1) + r_base(2);
t23 = qJ(1) + r_base(3);
t13 = cos(pkin(7));
t22 = t13 * pkin(1) + t25 * qJ(2) + r_base(1);
t21 = -pkin(5) + t23;
t20 = t13 * pkin(2) + t22;
t3 = -t13 * t28 - t25 * t27;
t19 = -t3 * pkin(3) + t20;
t14 = sin(qJ(5));
t15 = cos(qJ(5));
t18 = rSges(6,1) * t14 + rSges(6,2) * t15 + qJ(4);
t17 = t25 * pkin(2) - t13 * qJ(2) + t24;
t4 = t13 * t27 - t25 * t28;
t16 = -t4 * pkin(3) + t17;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t13 * rSges(2,1) - t25 * rSges(2,2) + r_base(1)) + g(2) * (t25 * rSges(2,1) + t13 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t23)) - m(3) * (g(1) * (t13 * rSges(3,1) + t25 * rSges(3,3) + t22) + g(2) * (t25 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t13 + t24) + g(3) * (rSges(3,2) + t23)) - m(4) * (g(1) * (-t3 * rSges(4,1) - t4 * rSges(4,2) + t20) + g(2) * (-t4 * rSges(4,1) + t3 * rSges(4,2) + t17) + g(3) * (-rSges(4,3) + t21)) - m(5) * (g(1) * (t3 * rSges(5,2) + t26 * t4 + t19) + g(2) * (t4 * rSges(5,2) - t26 * t3 + t16) + g(3) * (-rSges(5,1) + t21)) - m(6) * (g(1) * t19 + g(2) * t16 + g(3) * (-t15 * rSges(6,1) + t14 * rSges(6,2) - pkin(4) + t21) + (g(1) * t18 + g(2) * t29) * t4 + (g(1) * t29 - g(2) * t18) * t3);
U = t1;

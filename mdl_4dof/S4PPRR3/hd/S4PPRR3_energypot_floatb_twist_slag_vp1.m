% Calculate potential energy for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPRR3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:19
% EndTime: 2019-12-31 16:17:19
% DurationCPUTime: 0.14s
% Computational Cost: add. (80->52), mult. (92->52), div. (0->0), fcn. (86->6), ass. (0->17)
t24 = rSges(5,3) + pkin(5);
t23 = cos(qJ(3));
t22 = sin(qJ(3));
t21 = sin(pkin(6));
t20 = t21 * pkin(1) + r_base(2);
t19 = qJ(1) + r_base(3);
t11 = cos(pkin(6));
t18 = t11 * pkin(1) + t21 * qJ(2) + r_base(1);
t17 = -pkin(4) + t19;
t16 = t11 * pkin(2) + t18;
t12 = sin(qJ(4));
t13 = cos(qJ(4));
t15 = -rSges(5,1) * t13 + rSges(5,2) * t12 - pkin(3);
t14 = t21 * pkin(2) - t11 * qJ(2) + t20;
t2 = t11 * t22 - t21 * t23;
t1 = -t11 * t23 - t21 * t22;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t11 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t11 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t19)) - m(3) * (g(1) * (t11 * rSges(3,1) + t21 * rSges(3,3) + t18) + g(2) * (t21 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t11 + t20) + g(3) * (rSges(3,2) + t19)) - m(4) * (g(1) * (-rSges(4,1) * t1 - rSges(4,2) * t2 + t16) + g(2) * (-rSges(4,1) * t2 + rSges(4,2) * t1 + t14) + g(3) * (-rSges(4,3) + t17)) - m(5) * (g(1) * t16 + g(2) * t14 + g(3) * (-t12 * rSges(5,1) - t13 * rSges(5,2) + t17) + (g(1) * t24 + g(2) * t15) * t2 + (g(1) * t15 - g(2) * t24) * t1);
U = t3;

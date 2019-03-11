% Calculate potential energy for
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:13
% EndTime: 2019-03-08 18:26:14
% DurationCPUTime: 0.27s
% Computational Cost: add. (108->70), mult. (175->79), div. (0->0), fcn. (177->6), ass. (0->25)
t30 = rSges(5,1) + pkin(3);
t33 = rSges(5,3) + qJ(4);
t18 = sin(qJ(1));
t32 = g(1) * t18;
t19 = cos(qJ(1));
t31 = g(2) * t19;
t14 = sin(pkin(6));
t29 = t14 * t18;
t17 = cos(pkin(4));
t28 = t17 * t19;
t16 = cos(pkin(6));
t27 = t18 * t16;
t26 = pkin(5) + r_base(3);
t25 = t18 * pkin(1) + r_base(2);
t24 = t17 * qJ(2) + t26;
t15 = sin(pkin(4));
t23 = t18 * t15 * qJ(2) + t19 * pkin(1) + r_base(1);
t22 = t15 * t14 * pkin(2) + t24;
t3 = -t16 * t28 + t29;
t4 = t14 * t28 + t27;
t21 = t4 * pkin(2) + t3 * qJ(3) + t25;
t5 = t14 * t19 + t17 * t27;
t6 = t16 * t19 - t17 * t29;
t20 = t6 * pkin(2) + qJ(3) * t5 + t23;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t18 * rSges(2,2) + r_base(1)) + g(2) * (t18 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t26)) - m(3) * (g(1) * (rSges(3,1) * t6 - rSges(3,2) * t5 + t23) + g(2) * (t4 * rSges(3,1) - t3 * rSges(3,2) + t25) + g(3) * (rSges(3,3) * t17 + t24) + (rSges(3,3) * t32 + g(3) * (rSges(3,1) * t14 + rSges(3,2) * t16) + (-rSges(3,3) - qJ(2)) * t31) * t15) - m(4) * (g(1) * (-rSges(4,2) * t6 + rSges(4,3) * t5 + t20) + g(2) * (-t4 * rSges(4,2) + t3 * rSges(4,3) + t21) + g(3) * (rSges(4,1) * t17 + t22) + (rSges(4,1) * t32 + g(3) * (-rSges(4,2) * t14 + (-rSges(4,3) - qJ(3)) * t16) + (-rSges(4,1) - qJ(2)) * t31) * t15) - m(5) * (g(1) * (rSges(5,2) * t5 + t33 * t6 + t20) + g(2) * (t3 * rSges(5,2) + t33 * t4 + t21) + g(3) * (t30 * t17 + t22) + (g(3) * ((-rSges(5,2) - qJ(3)) * t16 + t33 * t14) + t30 * t32 + (-qJ(2) - t30) * t31) * t15);
U  = t1;

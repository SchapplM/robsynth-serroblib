% Calculate potential energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:50
% EndTime: 2019-07-18 13:24:50
% DurationCPUTime: 0.24s
% Computational Cost: add. (80->69), mult. (123->87), div. (0->0), fcn. (115->8), ass. (0->22)
t12 = cos(qJ(3));
t8 = sin(qJ(3));
t24 = rSges(4,1) * t12 - rSges(4,2) * t8;
t9 = sin(qJ(1));
t23 = t8 * t9;
t7 = sin(qJ(4));
t22 = t9 * t7;
t11 = cos(qJ(4));
t20 = t11 * t8;
t13 = cos(qJ(1));
t19 = t13 * t8;
t18 = t9 * t11;
t16 = t12 * t13;
t15 = t9 * qJ(2) + r_base(1);
t14 = -t13 * qJ(2) + r_base(2);
t10 = cos(qJ(5));
t6 = sin(qJ(5));
t4 = t11 * t16 + t22;
t3 = t7 * t16 - t18;
t2 = t12 * t18 - t13 * t7;
t1 = t11 * t13 + t12 * t22;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t13 - t9 * rSges(2,2) + r_base(1)) + g(2) * (t9 * rSges(2,1) + rSges(2,2) * t13 + r_base(2)) + g(3) * (r_base(3) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t13 + t9 * rSges(3,3) + t15) + g(2) * (t9 * rSges(3,1) + r_base(2) + (-rSges(3,3) - qJ(2)) * t13) + g(3) * (r_base(3) + rSges(3,2))) - m(4) * (g(1) * (t9 * rSges(4,3) + t15) + g(2) * (t24 * t9 + r_base(2)) + g(3) * (rSges(4,1) * t8 + rSges(4,2) * t12 + r_base(3)) + (g(1) * t24 + g(2) * (-rSges(4,3) - qJ(2))) * t13) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t19 + t15) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + rSges(5,3) * t23 + t14) + g(3) * (-rSges(5,3) * t12 + r_base(3) + (rSges(5,1) * t11 - rSges(5,2) * t7) * t8)) - m(6) * (g(1) * ((t4 * t10 + t6 * t19) * rSges(6,1) + (t10 * t19 - t4 * t6) * rSges(6,2) + t3 * rSges(6,3) + t15) + g(2) * ((t10 * t2 + t6 * t23) * rSges(6,1) + (t10 * t23 - t2 * t6) * rSges(6,2) + t1 * rSges(6,3) + t14) + g(3) * (r_base(3) + (t10 * t20 - t12 * t6) * rSges(6,1) + (-t10 * t12 - t6 * t20) * rSges(6,2) + t8 * t7 * rSges(6,3)));
U  = t5;

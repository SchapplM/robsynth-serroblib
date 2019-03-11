% Calculate potential energy for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR7_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:07
% EndTime: 2019-03-09 02:56:07
% DurationCPUTime: 0.34s
% Computational Cost: add. (176->92), mult. (167->91), div. (0->0), fcn. (143->8), ass. (0->33)
t14 = sin(qJ(6));
t17 = cos(qJ(6));
t23 = rSges(7,1) * t14 + rSges(7,2) * t17;
t12 = qJ(3) + pkin(9);
t6 = sin(t12);
t43 = t6 * pkin(4);
t15 = sin(qJ(3));
t42 = t15 * pkin(3);
t41 = rSges(4,3) + pkin(7);
t40 = rSges(7,3) + pkin(8);
t7 = cos(t12);
t37 = t7 * qJ(5);
t13 = -qJ(4) - pkin(7);
t36 = rSges(6,1) - t13;
t35 = rSges(5,3) - t13;
t34 = pkin(6) + r_base(3);
t16 = sin(qJ(1));
t33 = t16 * pkin(1) + r_base(2);
t32 = pkin(2) + t34;
t19 = cos(qJ(1));
t31 = t19 * pkin(1) + t16 * qJ(2) + r_base(1);
t30 = -qJ(2) - t42;
t29 = g(2) * t33;
t18 = cos(qJ(3));
t28 = t18 * pkin(3) + t32;
t27 = t16 * t42 + t31;
t26 = rSges(5,1) * t6 + rSges(5,2) * t7;
t25 = -rSges(6,2) * t6 - rSges(6,3) * t7;
t24 = rSges(4,1) * t15 + rSges(4,2) * t18;
t22 = t7 * pkin(4) + t6 * qJ(5) + t28;
t21 = t17 * rSges(7,1) - t14 * rSges(7,2) + pkin(5) - t13;
t20 = g(1) * (t16 * t43 + t27) + g(2) * (t19 * t37 + t33);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t19 * rSges(2,1) - t16 * rSges(2,2) + r_base(1)) + g(2) * (t16 * rSges(2,1) + t19 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (-t19 * rSges(3,2) + t16 * rSges(3,3) + t31) + g(2) * (-t16 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t19 + t33) + g(3) * (rSges(3,1) + t34)) - m(4) * (g(1) * t31 + t29 + g(3) * (t18 * rSges(4,1) - t15 * rSges(4,2) + t32) + (g(1) * t24 + g(2) * t41) * t16 + (g(1) * t41 + g(2) * (-qJ(2) - t24)) * t19) - m(5) * (g(1) * t27 + t29 + g(3) * (t7 * rSges(5,1) - t6 * rSges(5,2) + t28) + (g(1) * t26 + g(2) * t35) * t16 + (g(1) * t35 + g(2) * (-t26 + t30)) * t19) - m(6) * (g(3) * (-t7 * rSges(6,2) + t6 * rSges(6,3) + t22) + (g(1) * (t25 - t37) + g(2) * t36) * t16 + (g(1) * t36 + g(2) * (-t25 + t30 - t43)) * t19 + t20) - m(7) * (g(3) * (t23 * t6 + t40 * t7 + t22) + (g(2) * t21 + (t40 * t6 + (-qJ(5) - t23) * t7) * g(1)) * t16 + (g(1) * t21 + (t30 + (-pkin(4) - t40) * t6 + t23 * t7) * g(2)) * t19 + t20);
U  = t1;

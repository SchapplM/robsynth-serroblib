% Calculate potential energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:40
% EndTime: 2019-03-09 02:46:41
% DurationCPUTime: 0.56s
% Computational Cost: add. (260->111), mult. (274->129), div. (0->0), fcn. (278->10), ass. (0->40)
t53 = -rSges(7,3) - pkin(8);
t22 = pkin(9) + qJ(3);
t20 = cos(t22);
t31 = cos(qJ(1));
t52 = t20 * t31;
t19 = sin(t22);
t29 = sin(qJ(1));
t51 = t29 * t19;
t23 = sin(pkin(10));
t50 = t29 * t23;
t25 = cos(pkin(10));
t49 = t29 * t25;
t48 = t31 * t19;
t47 = t31 * t23;
t46 = t31 * t25;
t45 = qJ(4) * t19;
t44 = rSges(3,3) + qJ(2);
t43 = rSges(6,3) + qJ(5);
t42 = pkin(6) + r_base(3);
t26 = cos(pkin(9));
t17 = pkin(2) * t26 + pkin(1);
t41 = t31 * t17 + r_base(1);
t24 = sin(pkin(9));
t40 = pkin(2) * t24 + t42;
t27 = -pkin(7) - qJ(2);
t39 = t17 * t29 + t27 * t31 + r_base(2);
t38 = pkin(3) * t19 + t40;
t37 = t39 + (pkin(3) * t20 + t45) * t29;
t36 = rSges(3,1) * t26 - rSges(3,2) * t24 + pkin(1);
t35 = t38 + (pkin(4) * t25 + qJ(5) * t23) * t19;
t4 = t20 * t49 - t47;
t34 = pkin(4) * t4 + t37;
t33 = pkin(3) * t52 - t29 * t27 + t31 * t45 + t41;
t6 = t20 * t46 + t50;
t32 = pkin(4) * t6 + t33;
t30 = cos(qJ(6));
t28 = sin(qJ(6));
t5 = t20 * t47 - t49;
t3 = t20 * t50 + t46;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t31 - rSges(2,2) * t29 + r_base(1)) + g(2) * (rSges(2,1) * t29 + rSges(2,2) * t31 + r_base(2)) + g(3) * (rSges(2,3) + t42)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(3,1) * t24 + rSges(3,2) * t26 + t42) + (g(1) * t36 - g(2) * t44) * t31 + (g(1) * t44 + g(2) * t36) * t29) - m(4) * (g(1) * (rSges(4,1) * t52 - rSges(4,2) * t48 + t41) + g(2) * (-t31 * rSges(4,3) + t39) + g(3) * (rSges(4,1) * t19 + rSges(4,2) * t20 + t40) + (g(1) * (rSges(4,3) - t27) + g(2) * (rSges(4,1) * t20 - rSges(4,2) * t19)) * t29) - m(5) * (g(1) * (rSges(5,1) * t6 - rSges(5,2) * t5 + rSges(5,3) * t48 + t33) + g(2) * (rSges(5,1) * t4 - rSges(5,2) * t3 + rSges(5,3) * t51 + t37) + g(3) * ((-rSges(5,3) - qJ(4)) * t20 + (rSges(5,1) * t25 - rSges(5,2) * t23) * t19 + t38)) - m(6) * (g(1) * (t6 * rSges(6,1) + rSges(6,2) * t48 + t43 * t5 + t32) + g(2) * (t4 * rSges(6,1) + rSges(6,2) * t51 + t3 * t43 + t34) + g(3) * ((-rSges(6,2) - qJ(4)) * t20 + (rSges(6,1) * t25 + rSges(6,3) * t23) * t19 + t35)) - m(7) * (g(1) * (t6 * pkin(5) + t5 * qJ(5) + (t28 * t5 + t30 * t6) * rSges(7,1) + (-t28 * t6 + t30 * t5) * rSges(7,2) + t32) + g(2) * (t4 * pkin(5) + t3 * qJ(5) + (t28 * t3 + t30 * t4) * rSges(7,1) + (-t28 * t4 + t3 * t30) * rSges(7,2) + t34) + (g(1) * t31 + g(2) * t29) * t19 * t53 + (t35 + (-qJ(4) - t53) * t20 + (t25 * pkin(5) + (t23 * t28 + t25 * t30) * rSges(7,1) + (t23 * t30 - t25 * t28) * rSges(7,2)) * t19) * g(3));
U  = t1;

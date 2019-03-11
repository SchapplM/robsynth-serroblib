% Calculate potential energy for
% S6RPRPPR2
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:09
% EndTime: 2019-03-09 02:41:10
% DurationCPUTime: 0.40s
% Computational Cost: add. (235->96), mult. (159->103), div. (0->0), fcn. (135->10), ass. (0->35)
t46 = rSges(7,3) + pkin(8);
t45 = rSges(4,3) + pkin(7);
t17 = qJ(3) + pkin(10);
t10 = sin(t17);
t18 = qJ(1) + pkin(9);
t13 = cos(t18);
t43 = t10 * t13;
t11 = sin(t18);
t20 = sin(qJ(6));
t42 = t11 * t20;
t23 = cos(qJ(6));
t41 = t11 * t23;
t12 = cos(t17);
t40 = t13 * t12;
t39 = t13 * t20;
t38 = t13 * t23;
t37 = qJ(5) * t10;
t36 = pkin(6) + r_base(3);
t22 = sin(qJ(1));
t35 = t22 * pkin(1) + r_base(2);
t25 = cos(qJ(1));
t34 = t25 * pkin(1) + r_base(1);
t24 = cos(qJ(3));
t9 = t24 * pkin(3) + pkin(2);
t33 = t13 * t9 + t34;
t32 = qJ(2) + t36;
t19 = -qJ(4) - pkin(7);
t31 = t11 * t9 + t13 * t19 + t35;
t21 = sin(qJ(3));
t30 = t21 * pkin(3) + t32;
t29 = pkin(4) * t40 + t13 * t37 + t33;
t28 = t10 * pkin(4) + t30;
t27 = t31 + (pkin(4) * t12 + t37) * t11;
t26 = rSges(4,1) * t24 - rSges(4,2) * t21 + pkin(2);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t36)) - m(3) * (g(1) * (t13 * rSges(3,1) - t11 * rSges(3,2) + t34) + g(2) * (t11 * rSges(3,1) + t13 * rSges(3,2) + t35) + g(3) * (rSges(3,3) + t32)) - m(4) * (g(1) * t34 + g(2) * t35 + g(3) * (t21 * rSges(4,1) + t24 * rSges(4,2) + t32) + (g(1) * t26 - g(2) * t45) * t13 + (g(1) * t45 + g(2) * t26) * t11) - m(5) * (g(1) * (rSges(5,1) * t40 - rSges(5,2) * t43 + t33) + g(2) * (-t13 * rSges(5,3) + t31) + g(3) * (t10 * rSges(5,1) + t12 * rSges(5,2) + t30) + (g(1) * (rSges(5,3) - t19) + g(2) * (rSges(5,1) * t12 - rSges(5,2) * t10)) * t11) - m(6) * (g(1) * (-rSges(6,2) * t40 + rSges(6,3) * t43 + t29) + g(2) * (-t13 * rSges(6,1) + t27) + g(3) * (-t10 * rSges(6,2) + (-rSges(6,3) - qJ(5)) * t12 + t28) + (g(1) * (rSges(6,1) - t19) + g(2) * (-rSges(6,2) * t12 + rSges(6,3) * t10)) * t11) - m(7) * (g(1) * ((t10 * t39 + t41) * rSges(7,1) + (t10 * t38 - t42) * rSges(7,2) + t29 + (pkin(5) - t19) * t11) + g(2) * (-t13 * pkin(5) + (t10 * t42 - t38) * rSges(7,1) + (t10 * t41 + t39) * rSges(7,2) + t27) + g(3) * (t46 * t10 + t28) + (g(3) * (-rSges(7,1) * t20 - rSges(7,2) * t23 - qJ(5)) + (g(1) * t13 + g(2) * t11) * t46) * t12);
U  = t1;

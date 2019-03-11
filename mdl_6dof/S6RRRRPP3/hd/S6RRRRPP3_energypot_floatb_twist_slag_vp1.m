% Calculate potential energy for
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:38
% EndTime: 2019-03-09 20:53:38
% DurationCPUTime: 0.45s
% Computational Cost: add. (244->102), mult. (248->112), div. (0->0), fcn. (242->8), ass. (0->39)
t53 = rSges(7,2) + qJ(5);
t52 = rSges(7,3) + qJ(6);
t51 = rSges(7,1) + pkin(5);
t50 = rSges(3,3) + pkin(7);
t23 = qJ(2) + qJ(3);
t20 = sin(t23);
t26 = sin(qJ(1));
t49 = t20 * t26;
t29 = cos(qJ(1));
t48 = t20 * t29;
t21 = cos(t23);
t47 = t21 * t29;
t24 = sin(qJ(4));
t46 = t26 * t24;
t27 = cos(qJ(4));
t45 = t26 * t27;
t44 = t29 * t24;
t43 = t29 * t27;
t42 = rSges(6,3) + qJ(5);
t41 = pkin(6) + r_base(3);
t28 = cos(qJ(2));
t18 = t28 * pkin(2) + pkin(1);
t40 = t29 * t18 + r_base(1);
t25 = sin(qJ(2));
t39 = t25 * pkin(2) + t41;
t30 = -pkin(8) - pkin(7);
t38 = t26 * t18 + t29 * t30 + r_base(2);
t37 = t20 * pkin(3) + t39;
t36 = t26 * t21 * pkin(3) + pkin(9) * t49 + t38;
t35 = rSges(3,1) * t28 - rSges(3,2) * t25 + pkin(1);
t34 = t37 + (pkin(4) * t27 + qJ(5) * t24) * t20;
t4 = t21 * t45 - t44;
t33 = t4 * pkin(4) + t36;
t32 = pkin(3) * t47 + pkin(9) * t48 - t26 * t30 + t40;
t6 = t21 * t43 + t46;
t31 = t6 * pkin(4) + t32;
t5 = t21 * t44 - t45;
t3 = t21 * t46 + t43;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (t25 * rSges(3,1) + t28 * rSges(3,2) + t41) + (g(1) * t35 - g(2) * t50) * t29 + (g(1) * t50 + g(2) * t35) * t26) - m(4) * (g(1) * (rSges(4,1) * t47 - rSges(4,2) * t48 + t40) + g(2) * (-t29 * rSges(4,3) + t38) + g(3) * (t20 * rSges(4,1) + t21 * rSges(4,2) + t39) + (g(1) * (rSges(4,3) - t30) + g(2) * (rSges(4,1) * t21 - rSges(4,2) * t20)) * t26) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + rSges(5,3) * t48 + t32) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t49 + t36) + g(3) * ((-rSges(5,3) - pkin(9)) * t21 + (rSges(5,1) * t27 - rSges(5,2) * t24) * t20 + t37)) - m(6) * (g(1) * (rSges(6,1) * t48 - t6 * rSges(6,2) + t42 * t5 + t31) + g(2) * (rSges(6,1) * t49 - t4 * rSges(6,2) + t3 * t42 + t33) + g(3) * ((-rSges(6,1) - pkin(9)) * t21 + (-rSges(6,2) * t27 + rSges(6,3) * t24) * t20 + t34)) - m(7) * (g(1) * (t53 * t5 + t52 * t6 + t31) + g(2) * (t53 * t3 + t52 * t4 + t33) + (g(1) * t29 + g(2) * t26) * t20 * t51 + (t34 + (-pkin(9) - t51) * t21 + (rSges(7,2) * t24 + t52 * t27) * t20) * g(3));
U  = t1;

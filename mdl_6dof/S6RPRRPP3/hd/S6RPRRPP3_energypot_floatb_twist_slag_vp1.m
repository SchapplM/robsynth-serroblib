% Calculate potential energy for
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP3_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:34:57
% EndTime: 2019-03-09 04:34:57
% DurationCPUTime: 0.42s
% Computational Cost: add. (253->98), mult. (234->108), div. (0->0), fcn. (228->8), ass. (0->35)
t50 = rSges(7,2) + qJ(5);
t49 = rSges(7,3) + qJ(6);
t48 = rSges(7,1) + pkin(5);
t23 = qJ(1) + pkin(9);
t18 = sin(t23);
t25 = sin(qJ(3));
t47 = t18 * t25;
t28 = cos(qJ(3));
t46 = t18 * t28;
t19 = cos(t23);
t45 = t19 * t25;
t24 = sin(qJ(4));
t44 = t24 * t28;
t27 = cos(qJ(4));
t43 = t27 * t28;
t42 = rSges(6,3) + qJ(5);
t41 = pkin(6) + r_base(3);
t26 = sin(qJ(1));
t40 = t26 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t39 = t29 * pkin(1) + r_base(1);
t38 = qJ(2) + t41;
t37 = t18 * pkin(2) + t40;
t36 = t25 * pkin(3) + t38;
t35 = t19 * pkin(2) + t18 * pkin(7) + t39;
t34 = t19 * t28 * pkin(3) + pkin(8) * t45 + t35;
t33 = t36 + (pkin(4) * t27 + qJ(5) * t24) * t25;
t6 = t18 * t24 + t19 * t43;
t32 = t6 * pkin(4) + t34;
t31 = pkin(3) * t46 - t19 * pkin(7) + pkin(8) * t47 + t37;
t4 = t18 * t43 - t19 * t24;
t30 = t4 * pkin(4) + t31;
t5 = -t18 * t27 + t19 * t44;
t3 = t18 * t44 + t19 * t27;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t19 * rSges(3,1) - t18 * rSges(3,2) + t39) + g(2) * (t18 * rSges(3,1) + t19 * rSges(3,2) + t40) + g(3) * (rSges(3,3) + t38)) - m(4) * (g(1) * (t18 * rSges(4,3) + t35) + g(2) * (rSges(4,1) * t46 - rSges(4,2) * t47 + t37) + g(3) * (t25 * rSges(4,1) + t28 * rSges(4,2) + t38) + (g(1) * (rSges(4,1) * t28 - rSges(4,2) * t25) + g(2) * (-rSges(4,3) - pkin(7))) * t19) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + rSges(5,3) * t45 + t34) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t47 + t31) + g(3) * ((-rSges(5,3) - pkin(8)) * t28 + (rSges(5,1) * t27 - rSges(5,2) * t24) * t25 + t36)) - m(6) * (g(1) * (rSges(6,1) * t45 - t6 * rSges(6,2) + t42 * t5 + t32) + g(2) * (rSges(6,1) * t47 - t4 * rSges(6,2) + t42 * t3 + t30) + g(3) * ((-rSges(6,1) - pkin(8)) * t28 + (-rSges(6,2) * t27 + rSges(6,3) * t24) * t25 + t33)) - m(7) * (g(1) * (t49 * t6 + t50 * t5 + t32) + g(2) * (t50 * t3 + t49 * t4 + t30) + (g(1) * t19 + g(2) * t18) * t25 * t48 + (t33 + (-pkin(8) - t48) * t28 + (rSges(7,2) * t24 + t49 * t27) * t25) * g(3));
U  = t1;

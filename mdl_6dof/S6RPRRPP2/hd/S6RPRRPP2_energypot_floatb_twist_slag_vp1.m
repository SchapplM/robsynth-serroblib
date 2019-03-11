% Calculate potential energy for
% S6RPRRPP2
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:23
% EndTime: 2019-03-09 04:31:23
% DurationCPUTime: 0.42s
% Computational Cost: add. (253->98), mult. (234->108), div. (0->0), fcn. (228->8), ass. (0->35)
t49 = rSges(7,1) + pkin(5);
t48 = rSges(7,2) + qJ(5);
t22 = qJ(1) + pkin(9);
t17 = sin(t22);
t24 = sin(qJ(3));
t47 = t17 * t24;
t27 = cos(qJ(3));
t46 = t17 * t27;
t18 = cos(t22);
t45 = t18 * t24;
t23 = sin(qJ(4));
t44 = t23 * t27;
t26 = cos(qJ(4));
t43 = t26 * t27;
t42 = rSges(6,3) + qJ(5);
t41 = -rSges(7,3) - qJ(6);
t40 = pkin(6) + r_base(3);
t25 = sin(qJ(1));
t39 = pkin(1) * t25 + r_base(2);
t28 = cos(qJ(1));
t38 = pkin(1) * t28 + r_base(1);
t37 = qJ(2) + t40;
t36 = pkin(2) * t17 + t39;
t35 = pkin(3) * t24 + t37;
t34 = pkin(2) * t18 + pkin(7) * t17 + t38;
t33 = pkin(3) * t18 * t27 + pkin(8) * t45 + t34;
t32 = t35 + (pkin(4) * t26 + qJ(5) * t23) * t24;
t6 = t17 * t23 + t18 * t43;
t31 = pkin(4) * t6 + t33;
t30 = pkin(3) * t46 - t18 * pkin(7) + pkin(8) * t47 + t36;
t4 = t17 * t43 - t18 * t23;
t29 = pkin(4) * t4 + t30;
t5 = -t17 * t26 + t18 * t44;
t3 = t17 * t44 + t18 * t26;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t28 - rSges(2,2) * t25 + r_base(1)) + g(2) * (rSges(2,1) * t25 + rSges(2,2) * t28 + r_base(2)) + g(3) * (rSges(2,3) + t40)) - m(3) * (g(1) * (rSges(3,1) * t18 - rSges(3,2) * t17 + t38) + g(2) * (rSges(3,1) * t17 + rSges(3,2) * t18 + t39) + g(3) * (rSges(3,3) + t37)) - m(4) * (g(1) * (t17 * rSges(4,3) + t34) + g(2) * (rSges(4,1) * t46 - rSges(4,2) * t47 + t36) + g(3) * (rSges(4,1) * t24 + rSges(4,2) * t27 + t37) + (g(1) * (rSges(4,1) * t27 - rSges(4,2) * t24) + g(2) * (-rSges(4,3) - pkin(7))) * t18) - m(5) * (g(1) * (rSges(5,1) * t6 - rSges(5,2) * t5 + rSges(5,3) * t45 + t33) + g(2) * (rSges(5,1) * t4 - rSges(5,2) * t3 + rSges(5,3) * t47 + t30) + g(3) * ((-rSges(5,3) - pkin(8)) * t27 + (rSges(5,1) * t26 - rSges(5,2) * t23) * t24 + t35)) - m(6) * (g(1) * (t6 * rSges(6,1) + rSges(6,2) * t45 + t42 * t5 + t31) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t47 + t3 * t42 + t29) + g(3) * ((-rSges(6,2) - pkin(8)) * t27 + (rSges(6,1) * t26 + rSges(6,3) * t23) * t24 + t32)) - m(7) * (g(1) * (t48 * t5 + t49 * t6 + t31) + g(2) * (t3 * t48 + t4 * t49 + t29) + (g(1) * t18 + g(2) * t17) * t24 * t41 + (t32 + (-pkin(8) - t41) * t27 + (rSges(7,2) * t23 + t26 * t49) * t24) * g(3));
U  = t1;
